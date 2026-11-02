#!/usr/bin/env python3
"""
Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.

Check and/or apply clang-format and AMD copyright banners on source files.

Run ``python3 scripts/apply_format_changes.py -h`` for flags and examples.
"""

from __future__ import print_function

import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile

# Repo-root prefixes (generated / excluded trees).
SKIP_PREFIXES = ("tmp/", "old/", ".git/")
# Skip these directory names anywhere in the relative path.
SKIP_DIR_COMPONENTS = {".git"}
SKIP_WALK_DIRS = {".git", "tmp", "old"}
FORMAT_EXTS = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx"}
COPYRIGHT_EXTS = FORMAT_EXTS | {".py", ".sh"}
# BUILD.md requires clang-format 10 or newer. Try the plain name first, then
# each versioned name from newest to oldest (clang-format-20 .. clang-format-10).
CLANG_FORMAT_MIN_VERSION = 10
CLANG_FORMAT_MAX_VERSION = 20
CLANG_FORMAT_NAMES = ("clang-format",) + tuple(
    "clang-format-%d" % version
    for version in range(CLANG_FORMAT_MAX_VERSION, CLANG_FORMAT_MIN_VERSION - 1, -1)
)

COPYRIGHT_WORD_RE = re.compile(r"\bcopyright\b", re.IGNORECASE)
COPYRIGHT_PREFIX_RE = re.compile(r"\bcopyright\b(?:\s*\([Cc]\))?", re.IGNORECASE)
YEAR_RE = re.compile(r"\b(?:19|20)\d{2}\b")
AMD_RE = re.compile(r"Advanced Micro Devices|\bAMD,\s*Inc", re.IGNORECASE)
HEADER_LINES = 60
COPYRIGHT_MODE = "amd-if-present"
AMD_BANNER = "Copyright (C) %d, Advanced Micro Devices, Inc. All rights reserved."
AMD_MOD_BANNER = "Modifications Copyright (C) %d, Advanced Micro Devices, Inc. All rights reserved."
MODIFICATIONS_COPYRIGHT_RE = re.compile(
    r"\bmodifications?\s+copyright\b", re.IGNORECASE
)
RANGE_BEFORE_YEAR_RE = re.compile(r"(?:19|20)\d{2}\s*[-–]\s*$")
# PEP 263: optional whitespace, a # comment, then coding[:=] somewhere on that line.
SHEBANG_CODING_RE = re.compile(r"^[ \t\f]*#.*coding[:=]\s*[-\w.]+")
# A module docstring is a plain, r, or u triple-quoted string. f and b
# prefixes are executable literals; ast.get_docstring ignores them.
_DOCSTRING_OPENER_RE = re.compile(r"^(?:[rR]|[uU])?(\"\"\"|''')")


def repo_root():
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            universal_newlines=True,
        )
        return out.strip()
    except (subprocess.CalledProcessError, OSError):
        return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def _nul_paths(argv):
    """Return Git path names from a command that prints NUL-terminated paths."""
    out = subprocess.check_output(argv)
    return [item.decode("utf-8", "surrogateescape") for item in out.split(b"\0") if item]


def git_files(args):
    if args.git_changed:
        names = set(
            _nul_paths(
                ["git", "diff", "--name-only", "-z", "--diff-filter=ACMR", "HEAD"]
            )
        )
        names.update(
            _nul_paths(["git", "ls-files", "-z", "--others", "--exclude-standard"])
        )
    else:
        names = set(
            _nul_paths(
                ["git", "diff", "--name-only", "-z", "--diff-filter=ACMR", "--cached"]
            )
        )
    return sorted(n for n in names if not should_skip(n.replace("\\", "/")))


def staged_rewrite_groups():
    """Return staged source paths the pre-commit hook may rewrite.

    The first list matches the index. The second list is staged and also has
    unstaged edits. Skipped directories and non-source files are omitted, so
    the hook does not stage them again.
    """
    # T is a staged type change. ACMR skips it, so a regular file replaced by
    # a symlink would never reach the hook's symlink check.
    out = subprocess.check_output(
        ["git", "diff", "--name-only", "--diff-filter=ACMRT", "--cached", "-z"],
    )
    full = []
    partial = []
    for name in out.split(b"\0"):
        if not name:
            continue
        rel = name.decode("utf-8", "surrogateescape").replace("\\", "/")
        if should_skip(rel):
            continue
        if not copyright_eligible(rel):
            continue
        proc = subprocess.run(
            ["git", "diff", "--quiet", "--", name],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        text = name.decode("utf-8", "surrogateescape")
        if proc.returncode == 0:
            full.append(text)
        elif proc.returncode == 1:
            partial.append(text)
        else:
            raise subprocess.CalledProcessError(proc.returncode, ["git", "diff", "--quiet", "--", text])
    return full, partial


def rel_to_root(path, root):
    rel = os.path.relpath(path, root)
    return rel.replace("\\", "/")


def should_skip(rel):
    rel = rel.replace("\\", "/")
    if any(rel == p.rstrip("/") or rel.startswith(p) for p in SKIP_PREFIXES):
        return True
    parts = [p for p in rel.split("/") if p and p != "."]
    return any(part in SKIP_DIR_COMPONENTS for part in parts)


def resolve_user_path(raw, root, start_cwd):
    if os.path.isabs(raw):
        return os.path.normpath(raw)
    for base in (start_cwd, root, os.getcwd()):
        candidate = os.path.normpath(os.path.join(base, raw))
        if os.path.exists(candidate):
            return candidate
    return os.path.normpath(os.path.join(root, raw))


def collect_source_files(paths, root, start_cwd=None):
    if start_cwd is None:
        start_cwd = os.getcwd()
    seen = set()
    files = []
    errors = []
    for raw in paths:
        path = resolve_user_path(raw, root, start_cwd)
        if not os.path.exists(path):
            errors.append("%s: path not found" % raw)
            continue
        if os.path.isdir(path):
            if os.path.islink(path):
                errors.append("%s: refusing to walk a symlink directory" % raw)
                continue
            for dirpath, dirnames, filenames in os.walk(path):
                dirnames[:] = [d for d in dirnames if d not in SKIP_WALK_DIRS]
                for name in filenames:
                    full = os.path.join(dirpath, name)
                    rel = rel_to_root(full, root)
                    if should_skip(rel):
                        continue
                    if not copyright_eligible(full):
                        continue
                    if full not in seen:
                        seen.add(full)
                        files.append(full)
        else:
            if path not in seen:
                seen.add(path)
                files.append(path)
    return files, errors


def clang_format_bin():
    """Return the clang-format executable path if it is on PATH and runnable."""
    for name in CLANG_FORMAT_NAMES:
        found = shutil.which(name)
        if not found:
            continue
        try:
            subprocess.check_output([found, "--version"], stderr=subprocess.DEVNULL)
            return found
        except (subprocess.CalledProcessError, OSError):
            continue
    return None


def notify_missing_clang_format(applying):
    names = "clang-format, clang-format-%d .. clang-format-%d" % (
        CLANG_FORMAT_MAX_VERSION,
        CLANG_FORMAT_MIN_VERSION,
    )
    action = "apply" if applying else "check"
    print(
        "error: clang-format executable was not found; cannot %s formatting." % action,
        file=sys.stderr,
    )
    print("  Looked for: %s" % names, file=sys.stderr)
    print(
        "  Install clang-format and add its directory to PATH, then retry.",
        file=sys.stderr,
    )
    print('  Example: export PATH="/usr/bin:$PATH"', file=sys.stderr)
    print(
        "  To skip formatting, pass --skip-clang-format (copyright can still run).",
        file=sys.stderr,
    )


def needs_clang_format(files):
    for path in files:
        if os.path.splitext(path)[1].lower() in FORMAT_EXTS:
            return True
    return False


def symlink_blocks_rewrite(path, root=None):
    """True when rewriting path would follow a symlink.

    The file itself is always refused. When root is given, a symlink
    component between root and path is refused too, as is a resolved path
    that leaves root. A symlink above root, such as a symlinked /tmp, is
    not this path escaping the repository.
    """
    if os.path.islink(path):
        return True
    if root is None:
        return False
    abs_path = os.path.abspath(path)
    root_abs = os.path.abspath(root)
    if abs_path != root_abs and not abs_path.startswith(root_abs + os.sep):
        return False
    current = root_abs
    rel = abs_path[len(root_abs) :].lstrip(os.sep)
    for part in rel.split(os.sep):
        if part in ("", "."):
            continue
        current = os.path.join(current, part)
        if os.path.islink(current):
            return True
    if not os.path.lexists(path):
        return False
    real = os.path.realpath(path)
    root_real = os.path.realpath(root)
    return real != root_real and not real.startswith(root_real + os.sep)


def check_format(path, cf, root, fix):
    ext = os.path.splitext(path)[1].lower()
    if ext not in FORMAT_EXTS:
        return None
    if not cf:
        return "clang-format executable was not found; add it to PATH"
    style = os.path.join(root, ".clang-format")
    cmd = [cf, "-style=file"]
    cwd = root if os.path.isfile(style) else os.path.dirname(path) or "."
    if fix:
        if symlink_blocks_rewrite(path, root):
            return "refusing to rewrite through a symlink"
        cmd.extend(["-i", path])
        mode = file_permissions(path)
        try:
            subprocess.check_call(cmd, cwd=cwd)
        except subprocess.CalledProcessError as exc:
            return "clang-format -i failed (exit %s)" % exc.returncode
        finally:
            if os.path.exists(path):
                restore_permissions(path, mode)
        return None
    cmd.extend(["--dry-run", "--Werror", path])
    proc = subprocess.run(
        cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True
    )
    if proc.returncode != 0:
        snippet = (proc.stdout or "").strip().splitlines()
        extra = snippet[0] if snippet else "needs clang-format"
        return extra
    return None


def years_on_copyright_line(line):
    if not COPYRIGHT_WORD_RE.search(line):
        return []
    return [int(match.group(0)) for match in YEAR_RE.finditer(line)]


def detect_newline(text):
    """Match the first line ending so inserted text keeps the file's format."""
    for index, char in enumerate(text):
        if char == "\r":
            if index + 1 < len(text) and text[index + 1] == "\n":
                return "\r\n"
            return "\r"
        if char == "\n":
            return "\n"
    return "\n"


def is_cmake_lists(path):
    """True for CMakeLists.txt. Other .txt files are not copyright sources."""
    return os.path.basename(path).lower() == "cmakelists.txt"


def is_makefile(path):
    """True for Makefile, which uses shell-style # comments."""
    return os.path.basename(path).lower() == "makefile"


def copyright_eligible(path):
    """True when this path is checked for an AMD copyright banner.

    C and C++ sources, Python, and shell are included by extension.
    Makefile and CMakeLists.txt are also included. Plain-text files such as
    LICENSE.txt and NOTICE.txt are not.
    """
    ext = os.path.splitext(path)[1].lower()
    return ext in COPYRIGHT_EXTS or is_cmake_lists(path) or is_makefile(path)


def hash_comment_ext(path):
    ext = os.path.splitext(path)[1].lower()
    return ext in {".py", ".sh"} or is_cmake_lists(path) or is_makefile(path)


def decode_source(data):
    """Decode source bytes without substituting U+FFFD.

    Valid UTF-8 stays UTF-8. Any other byte sequence is decoded as Latin-1,
    which maps each byte to one character and encodes back to the same bytes.
    """
    try:
        return data.decode("utf-8"), "utf-8"
    except UnicodeDecodeError:
        return data.decode("latin-1"), "latin-1"


def split_header_bytes(data, max_lines):
    """Split off the first max_lines records. The remainder is not decoded."""
    count = 0
    index = 0
    size = len(data)
    while index < size and count < max_lines:
        byte = data[index]
        if byte == 0x0D:
            index += 2 if index + 1 < size and data[index + 1] == 0x0A else 1
            count += 1
        elif byte == 0x0A:
            index += 1
            count += 1
        else:
            index += 1
    return data[:index], data[index:]


def read_source(path):
    with open(path, "rb") as handle:
        data = handle.read()
    header, body = split_header_bytes(data, HEADER_LINES)
    text, encoding = decode_source(header)
    return text, encoding, body


def read_header_lines(path):
    text, _encoding, _body = read_source(path)
    return text.splitlines(keepends=True)


def file_permissions(path):
    return os.stat(path).st_mode & 0o7777


def restore_permissions(path, mode):
    if os.stat(path).st_mode & 0o7777 != mode:
        os.chmod(path, mode)


def write_source(path, text, encoding, body):
    """Replace path only after the full contents are on disk.

    The original file is left untouched if encoding or the temporary write
    fails. A symlink, including a symlinked parent directory, is refused so
    the link target is not rewritten.
    """
    if symlink_blocks_rewrite(path, repo_root()):
        raise OSError("refusing to rewrite through a symlink: %s" % path)
    payload = text.encode(encoding) + body
    mode = file_permissions(path)
    directory = os.path.dirname(path) or "."
    fd, temp_name = tempfile.mkstemp(prefix=".apply-format-", dir=directory)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temp_name, mode)
        os.replace(temp_name, path)
        temp_name = None
    finally:
        if temp_name is not None:
            try:
                os.unlink(temp_name)
            except OSError:
                pass


_CPP_RAW_MARKERS = ("u8R\"", "UR\"", "LR\"", "uR\"", "R\"")


def _cpp_raw_string_open(line, i):
    """Return (delimiter, content index) if a C++ raw string starts at i."""
    if i > 0 and (line[i - 1].isalnum() or line[i - 1] == "_"):
        return None
    marker = None
    for candidate in _CPP_RAW_MARKERS:
        if line.startswith(candidate, i):
            marker = candidate
            break
    if marker is None:
        return None
    j = i + len(marker)
    start = j
    n = len(line)
    while j < n and (j - start) < 16:
        ch = line[j]
        if ch == "(":
            return line[start:j], j + 1
        if ch in " ()\\" or ord(ch) < 32:
            return None
        j += 1
    if j < n and line[j] == "(" and (j - start) <= 16:
        return line[start:j], j + 1
    return None


def _raw_string_close(line, content_at, delim):
    """Index just after a raw string's closing quote, or None if it stays open."""
    end = line.find(")" + delim + '"', content_at)
    if end == -1:
        return None
    return end + len(delim) + 2


def comment_spans(line, hash_style, in_block, raw_delim=None):
    """Return (spans, still_in_block, raw_delim) covering comment text on this line.

    raw_delim is the C++ raw-string delimiter still open after this line, or None.
    Text inside a raw string is not a comment, even when it contains a quote
    followed by // or /*.
    """
    spans = []
    i = 0
    n = len(line)
    if raw_delim is not None and not hash_style:
        end = _raw_string_close(line, 0, raw_delim)
        if end is None:
            return spans, in_block, raw_delim
        i = end
        raw_delim = None
    if in_block and not hash_style and i == 0:
        end = line.find("*/")
        if end == -1:
            return [(0, n)], True, None
        spans.append((0, end + 2))
        i = end + 2
        in_block = False
    while i < n:
        if not hash_style:
            opened = _cpp_raw_string_open(line, i)
            if opened:
                delim, content_at = opened
                end = _raw_string_close(line, content_at, delim)
                if end is None:
                    return spans, in_block, delim
                i = end
                continue
        if hash_style:
            if line[i] in "\"'":
                quote = line[i]
                i += 1
                while i < n and line[i] != quote:
                    i += 2 if line[i] == "\\" else 1
                i += 1
                continue
            if line[i] == "#":
                spans.append((i, n))
                break
            i += 1
            continue
        if line.startswith("//", i):
            spans.append((i, n))
            break
        if line.startswith("/*", i):
            end = line.find("*/", i + 2)
            if end == -1:
                spans.append((i, n))
                return spans, True, None
            spans.append((i, end + 2))
            i = end + 2
            continue
        if line[i] in "\"'":
            quote = line[i]
            i += 1
            while i < n and line[i] != quote:
                i += 2 if line[i] == "\\" else 1
            i += 1
            continue
        i += 1
    return spans, in_block, raw_delim


def _trailing_splice(content):
    """True when an odd trailing backslash joins this line to the next."""
    body = content.rstrip("\r\n")
    count = 0
    i = len(body) - 1
    while i >= 0 and body[i] == "\\":
        count += 1
        i -= 1
    return count % 2 == 1


def _leading_comment_kind(content, kinds):
    """Return the comment introducer at the first non-whitespace, or None."""
    i = 0
    n = len(content)
    while i < n and content[i] in " \t":
        i += 1
    for kind in kinds:
        if content.startswith(kind, i):
            return kind
    return None


def _span_opens_header(content, spans):
    """True when a comment span starts at the first non-whitespace character."""
    i = 0
    n = len(content)
    while i < n and content[i] in " \t":
        i += 1
    return any(start == i for start, _end in spans)


def _backslash_outside_comment(content, spans):
    """True when a trailing backslash is not inside a comment span."""
    if not _trailing_splice(content):
        return False
    pos = len(content.rstrip("\r\n")) - 1
    for start, end in spans:
        if start <= pos < end:
            return False
    return True


def copyright_line_role(line, hash_style, in_block=False, raw_delim=None, at_header=False):
    """Classify Copyright matches as comment, code, ambiguous, or none.

    at_header is true when the comment begins the logical line, or this line
    continues a block comment that began that way. An inline comment such as
    code followed by // Copyright is not a banner.
    """
    matches = list(COPYRIGHT_WORD_RE.finditer(line))
    if not matches:
        return "none"
    spans, _still, _raw = comment_spans(line, hash_style, in_block, raw_delim)
    inside = [
        any(start <= match.start() and match.end() <= end for start, end in spans)
        for match in matches
    ]
    if any(inside) and not all(inside):
        return "ambiguous"
    if all(inside) and at_header and _spanned_copyright_is_notice(line, spans):
        return "comment"
    if all(inside):
        return "none"
    return "code"


def _shell_read_heredoc(line, i):
    """Read a here-document operator at i. Return (end, delim, strip_tabs).

    delim is the delimiter word after quote removal, or None when << has no
    word. The word is not limited to letters: <<END-OF-FILE closes on
    END-OF-FILE. end is the index just after that word.
    """
    n = len(line)
    j = i + 2
    strip_tabs = j < n and line[j] == "-"
    if strip_tabs:
        j += 1
    while j < n and line[j] in " \t":
        j += 1
    chars = []
    got_word = False
    while j < n:
        ch = line[j]
        if ch in "\"'":
            quote = ch
            j += 1
            while j < n and line[j] != quote:
                chars.append(line[j])
                j += 1
            if j < n:
                j += 1
            got_word = True
            continue
        if ch == "\\" and j + 1 < n:
            chars.append(line[j + 1])
            j += 2
            got_word = True
            continue
        if not got_word and ch == "#":
            break
        if ch in " \t\r\n|&;()<>":
            break
        chars.append(ch)
        j += 1
        got_word = True
    if not got_word:
        return j, None, False
    return j, "".join(chars), strip_tabs


def _shell_line_scan(line):
    """Return (comment spans, heredocs) for a shell line outside a heredoc.

    # starts a comment only at the beginning of a word. It is not a comment
    inside a word (foo#Copyright) or a parameter expansion (${value#Copyright}).
    heredocs is a list of (delimiter, strip_tabs) opened before any comment.
    A here-string (<<<) does not open a heredoc.
    """
    spans = []
    heredocs = []
    i = 0
    n = len(line)
    at_word_start = True
    # 'p' is ${...}; 'c' is $(...) nested inside one. # is not a comment in 'p'.
    stack = []
    while i < n:
        ch = line[i]
        if ch in "\"'":
            quote = ch
            i += 1
            while i < n and line[i] != quote:
                if quote == '"' and line[i] == "\\":
                    i += 2
                else:
                    i += 1
            if i < n:
                i += 1
            at_word_start = False
            continue
        if ch == "\\" and i + 1 < n:
            i += 2
            at_word_start = False
            continue
        in_param = bool(stack) and stack[-1] == "p"
        if ch == "#" and not in_param and at_word_start:
            spans.append((i, n))
            break
        if not in_param and line.startswith("<<<", i):
            i += 3
            at_word_start = True
            continue
        if not in_param and line.startswith("<<", i):
            j, delim, strip_tabs = _shell_read_heredoc(line, i)
            if delim is not None:
                heredocs.append((delim, strip_tabs))
            i = j
            at_word_start = True
            continue
        if line.startswith("${", i):
            stack.append("p")
            i += 2
            at_word_start = False
            continue
        if line.startswith("$(", i) and in_param:
            stack.append("c")
            i += 2
            at_word_start = True
            continue
        if ch == "}" and stack and stack[-1] == "p":
            stack.pop()
            i += 1
            at_word_start = False
            continue
        if ch == ")" and stack and stack[-1] == "c":
            stack.pop()
            i += 1
            at_word_start = False
            continue
        if not in_param and ch in " \t":
            at_word_start = True
            i += 1
            continue
        if not in_param and ch in "|&;()<>":
            at_word_start = True
            i += 1
            continue
        if ch in "\r\n":
            break
        at_word_start = False
        i += 1
    return spans, heredocs


def _iter_shell_lines(lines):
    """Yield (idx, content, in_body, depth_after, started, spans) for a shell header.

    in_body is true for heredoc payload, including the closing delimiter line.
    A # line in that payload is data, not a comment.
    """
    queue = []
    for idx, raw in enumerate(lines):
        content = split_leading_bom(raw)[1]
        if queue:
            delim, strip_tabs = queue[0]
            check = content.lstrip("\t") if strip_tabs else content
            if check.rstrip("\r\n") == delim:
                queue.pop(0)
            yield idx, content, True, len(queue), 0, []
            continue
        spans, started = _shell_line_scan(content)
        queue.extend(started)
        yield idx, content, False, len(queue), len(started), spans


def _blank_line(content):
    return not content.strip()


def _shell_copyright_indexes(lines):
    result = []
    code_splice = False
    in_header = True
    for idx, content, in_body, _depth, _started, spans in _iter_shell_lines(lines):
        if in_body:
            code_splice = False
            in_header = False
            continue
        # A # comment counts only when it begins the logical line and that line
        # is still in the leading header. A comment after code, including one
        # inside a function or $(...), is not a banner.
        comment_line = (not code_splice) and _span_opens_header(content, spans)
        if in_header and comment_line and _copyright_inside_spans(content, spans) == "comment":
            result.append(idx)
        if in_header and not code_splice and not _blank_line(content) and not comment_line:
            in_header = False
        code_splice = _backslash_outside_comment(content, spans)
    return result


def _shell_comment_end(lines, copy_idx):
    """Last line to keep before an inserted banner, or None if a heredoc stays open."""
    rows = list(_iter_shell_lines(lines))
    if copy_idx >= len(rows) or rows[copy_idx][2]:
        return copy_idx
    end = copy_idx
    while end + 1 < len(rows) and rows[end][4] == 0:
        _idx, content, in_body, _depth, _started, _spans = rows[end + 1]
        if content.lstrip().startswith("#") and not in_body:
            end += 1
            continue
        break
    if rows[end][4] == 0:
        return end
    depth_before = rows[end][3] - rows[end][4]
    for j in range(end + 1, len(rows)):
        if rows[j][3] == depth_before:
            return j
    return None


def leading_docstring_bounds(lines):
    """Return (start, end) of the module docstring.

    Comments and blank lines before the first statement are skipped, so a
    coding cookie on line 2 does not hide the docstring. end is None when the
    opener is present but the closer is not in lines. (None, None) means there
    is no module docstring.
    """
    i = 0
    n = len(lines)
    while i < n:
        content = split_leading_bom(lines[i])[1].strip()
        if content == "" or content.startswith("#"):
            i += 1
            continue
        break
    if i >= n:
        return None, None
    stripped = split_leading_bom(lines[i])[1].lstrip()
    opener = _DOCSTRING_OPENER_RE.match(stripped)
    if not opener:
        return None, None
    quote = opener.group(1)
    if _unescaped_triple(stripped, quote, opener.end()) is not None:
        return i, i
    j = i + 1
    while j < n:
        if _unescaped_triple(split_leading_bom(lines[j])[1], quote) is not None:
            return i, j
        j += 1
    return i, None


def _unescaped_triple(line, quote, start=0):
    """Index of quote not consumed by a backslash, or None.

    A backslash skips the next character, matching python_comment_spans, so an
    escaped triple quote inside a module docstring is not the closer.
    """
    i = start
    n = len(line)
    while i < n:
        if line.startswith("\\", i):
            i += 2
            continue
        if line.startswith(quote, i):
            return i
        i += 1
    return None


def _python_string_open(line, index):
    """Return (delimiter, length) if a Python string starts at index."""
    if index > 0 and (line[index - 1].isalnum() or line[index - 1] == "_"):
        return None
    rest = line[index:]
    lower = rest.lower()
    for prefix in ("rf", "fr", "rb", "br", "r", "u", "b", "f", ""):
        if prefix and not lower.startswith(prefix):
            continue
        after = rest[len(prefix) :]
        if after.startswith(('"""', "'''")):
            return after[:3], len(prefix) + 3
        if after[:1] in ("'", '"'):
            return after[0], len(prefix) + 1
    return None


def python_comment_spans(line, in_triple):
    """Return (# comment spans, triple-quote state) for one Python line.

    A # inside a single or triple-quoted string is not a comment. Triple-quoted
    strings carry state onto the next line. in_triple is ''' or \"\"\" or None.
    """
    spans = []
    i = 0
    n = len(line)
    if in_triple:
        while i < n:
            if line.startswith("\\", i):
                i += 2
                continue
            if line.startswith(in_triple, i):
                i += 3
                in_triple = None
                break
            i += 1
        if in_triple:
            return spans, in_triple
    while i < n:
        opened = _python_string_open(line, i)
        if opened:
            delim, length = opened
            i += length
            if len(delim) == 3:
                closed = False
                while i < n:
                    if line.startswith("\\", i):
                        i += 2
                        continue
                    if line.startswith(delim, i):
                        i += 3
                        closed = True
                        break
                    i += 1
                if not closed:
                    return spans, delim
            else:
                while i < n and line[i] != delim:
                    i += 2 if line[i] == "\\" else 1
                i += 1
            continue
        if line[i] == "#":
            spans.append((i, n))
            break
        i += 1
    return spans, None


def _copyright_inside_spans(content, spans):
    matches = list(COPYRIGHT_WORD_RE.finditer(content))
    if not matches:
        return "none"
    inside = [
        any(start <= match.start() and match.end() <= end for start, end in spans)
        for match in matches
    ]
    if any(inside) and not all(inside):
        return "ambiguous"
    if all(inside) and _spanned_copyright_is_notice(content, spans):
        return "comment"
    if all(inside):
        return "none"
    return "code"


def _is_copyright_notice(text):
    """True when text is a copyright notice, not prose that only mentions the word.

    A notice contains Copyright followed by (C) or a 19xx/20xx year, including
    later in the line ("This file is Copyright (C) 2019, ..."). A line that
    starts with Copyright is a notice when the rest is a banner tail, such as
    (C), a year, or Advanced Micro Devices. Prose such as "copyright banners"
    is not a notice.
    """
    text = text.strip()
    text = re.sub(r"^(?://|/\*)\s*", "", text)
    text = re.sub(r"\*/\s*$", "", text)
    text = text.strip().lstrip("*#-").strip()
    notice_tail = r"(?:\(\s*[Cc]\s*\)|\b(?:19|20)\d{2}\b)"
    if re.search(
        r"\b(?:modifications?\s+)?copyright\b\s*" + notice_tail,
        text,
        re.IGNORECASE,
    ):
        return True
    match = re.match(r"(?:modifications?\s+)?copyright\b\s*(.*)$", text, re.IGNORECASE)
    if not match:
        return False
    rest = match.group(1).strip()
    if not rest:
        return True
    return bool(
        re.match(
            r"(?:\(\s*[Cc]\s*\)|(?:19|20)\d{2}\b|,|Advanced Micro Devices\b|AMD\b|by\b)",
            rest,
            re.IGNORECASE,
        )
    )


def _spanned_copyright_is_notice(content, spans):
    snippet = _span_text_for_copyright(content, spans)
    return snippet is not None and _is_copyright_notice(snippet)


def _docstring_notice(content):
    """True when a docstring line is a copyright notice, not prose mentioning the word.

    A notice may sit later in the line, as in "This file is Copyright (C) 2019,
    ...". Prose that only uses the word, such as "copyright banners", does not
    count. A one-line module docstring such as \"\"\"Copyright (C) ...\"\"\"
    counts. The opening quotes and an optional r or u prefix are ignored.
    An f or b prefix is not a docstring.
    """
    text = content.strip()
    opener = _DOCSTRING_OPENER_RE.match(text)
    if opener:
        quote = opener.group(1)
        text = text[opener.end() :]
        if text.rstrip().endswith(quote):
            text = text.rstrip()[: -len(quote)]
    return _is_copyright_notice(text)


def header_copyright_indexes(lines, ext, path=""):
    """Header lines whose Copyright text is only inside a leading comment.

    A match in source text is ignored. A line with both a comment match and a
    source match is ambiguous and is not rewritten. A leading Python module
    docstring counts as the header comment. A # inside a later Python string
    does not. A comment after executable code has started is not a header, so
    a later # Copyright, // Copyright, or shell comment is left unchanged.
    """
    if ext == ".sh" or is_cmake_lists(path) or is_makefile(path):
        return _shell_copyright_indexes(lines)
    hash_style = ext == ".py"
    result = []
    in_block = False
    header_block = False
    raw_delim = None
    in_triple = None
    code_splice = False
    slash_header = False
    in_header = True
    doc_start, doc_end = leading_docstring_bounds(lines) if ext == ".py" else (None, None)
    for idx, raw in enumerate(lines):
        content = split_leading_bom(raw)[1]
        if doc_start is not None and doc_start <= idx and (doc_end is None or idx <= doc_end):
            if in_header and _docstring_notice(content):
                result.append(idx)
            continue
        if ext == ".py":
            spans, in_triple = python_comment_spans(content, in_triple)
            comment_line = (not code_splice) and _span_opens_header(content, spans)
            if in_header and comment_line and _copyright_inside_spans(content, spans) == "comment":
                result.append(idx)
            if in_header and not code_splice and not _blank_line(content) and not comment_line:
                in_header = False
            code_splice = _backslash_outside_comment(content, spans) and in_triple is None
            continue
        in_block_before = in_block
        raw_before = raw_delim
        opens = (
            in_header
            and not in_block_before
            and raw_before is None
            and not code_splice
            and _leading_comment_kind(content, ("//", "/*"))
        )
        continues_slash = code_splice and slash_header and not in_block_before
        at_header = bool(opens) or (in_block_before and header_block) or continues_slash
        role = copyright_line_role(content, hash_style, in_block, raw_delim, at_header=at_header)
        _spans, in_block, raw_delim = comment_spans(content, hash_style, in_block, raw_delim)
        if role == "comment" or (continues_slash and in_header and _is_copyright_notice(content)):
            result.append(idx)
        if in_block_before and header_block:
            header_block = in_block
        elif opens == "/*" and in_block:
            header_block = True
        else:
            header_block = False
        slash_header = (opens == "//" or continues_slash) and _trailing_splice(content) and raw_delim is None
        if raw_delim is None and not in_block:
            code_splice = _trailing_splice(content)
        else:
            code_splice = False
        if in_header and not (in_block_before or opens or continues_slash or _blank_line(content)):
            in_header = False
    return result


def analyze_copyright(header, ext, path=""):
    banner_idxs = header_copyright_indexes(header, ext, path)
    amd_years = []
    copy_idxs = []
    amd_idxs = []
    any_year = False
    has_other = False
    amd_line_idxs = []
    for idx in banner_idxs:
        line = header[idx]
        copy_idxs.append(idx)
        years = years_on_copyright_line(line)
        if AMD_RE.search(line):
            amd_line_idxs.append(idx)
            if years:
                any_year = True
                amd_years.append(max(years))
                amd_idxs.append(idx)
        else:
            has_other = True
            if years:
                any_year = True
    return {
        "has_word": bool(copy_idxs),
        "any_year": any_year,
        "has_other": has_other,
        "amd_years": amd_years,
        "amd_line_idxs": amd_line_idxs,
        "copy_idxs": copy_idxs,
        "amd_idxs": amd_idxs,
    }


def check_copyright(path, year, mode):
    ext = os.path.splitext(path)[1].lower()
    if not copyright_eligible(path):
        return None
    if mode == "off":
        return None
    try:
        header = read_header_lines(path)
    except OSError as exc:
        return "cannot read: %s" % exc
    info = analyze_copyright(header, ext, path)
    if not info["has_word"]:
        return "missing Copyright banner in the first %d lines" % HEADER_LINES

    amd_years = info["amd_years"]
    # Every AMD banner needs its own 19xx/20xx year. One dated line does not
    # cover a second AMD line that has no year.
    if any(idx not in info["amd_idxs"] for idx in info["amd_line_idxs"]):
        return "Copyright text found but year could not be parsed"

    if amd_years:
        if info["has_other"]:
            for idx in info["amd_idxs"]:
                if not MODIFICATIONS_COPYRIGHT_RE.search(header[idx]):
                    return (
                        "AMD copyright must use Modifications Copyright when another "
                        "copyright is present"
                    )
        stale = [value for value in amd_years if value < year]
        if stale:
            oldest = min(stale)
            return "AMD copyright last year is %d; touch requires %d (e.g. 2024-%d)" % (
                oldest,
                year,
                year,
            )
        return None

    if info["has_other"]:
        return (
            "no AMD copyright banner (add: Modifications Copyright (C) %d, "
            "Advanced Micro Devices, Inc.)" % year
        )
    return (
        "no AMD copyright banner (add: Copyright (C) %d, Advanced Micro Devices, Inc.)" % year
    )


def bump_last_year(line, year):
    matches = list(YEAR_RE.finditer(line))
    if not matches:
        return line, False
    last = matches[-1]
    last_year = int(last.group(0))
    if last_year >= year:
        return line, False
    prefix = line[: last.start()]
    if RANGE_BEFORE_YEAR_RE.search(prefix):
        new_line = prefix + str(year) + line[last.end() :]
    else:
        new_line = prefix + "%d-%d" % (last_year, year) + line[last.end() :]
    return new_line, new_line != line


def inject_year_after_copyright(line, year):
    if years_on_copyright_line(line):
        return line, False

    def repl(match):
        return "%s %d," % (match.group(0), year)

    new_line, n = COPYRIGHT_PREFIX_RE.subn(repl, line, count=1)
    return new_line, n > 0


UTF8_BOM = "\ufeff"
UTF8_BOM_LATIN1 = "\xef\xbb\xbf"


def split_leading_bom(line):
    """Return (bom, rest). A UTF-8 BOM may be U+FEFF or three Latin-1 bytes."""
    if line.startswith(UTF8_BOM):
        return UTF8_BOM, line[len(UTF8_BOM) :]
    if line.startswith(UTF8_BOM_LATIN1):
        return UTF8_BOM_LATIN1, line[len(UTF8_BOM_LATIN1) :]
    return "", line


def shebang_insert_index(lines):
    """Index after a shebang and a PEP 263 coding cookie in the first two lines."""
    idx = 0
    if lines:
        _bom, content = split_leading_bom(lines[0])
        if content.startswith("#!"):
            idx = 1
    # PEP 263 reads a coding cookie on line 1 or line 2 only.
    for cookie_idx in range(min(len(lines), 2)):
        _bom, content = split_leading_bom(lines[cookie_idx])
        if SHEBANG_CODING_RE.search(content):
            return cookie_idx + 1
    return idx


def new_amd_banner_lines(path, year, newline, modifications=False):
    template = AMD_MOD_BANNER if modifications else AMD_BANNER
    body = template % year
    if hash_comment_ext(path):
        return ["# %s%s" % (body, newline)]
    return ["/*%s" % newline, "    %s%s" % (body, newline), "*/%s" % newline]


def ensure_modifications_copyright(line):
    """Prefix Copyright with Modifications on an AMD line. Leave it when already present."""
    if MODIFICATIONS_COPYRIGHT_RE.search(line):
        return line, False
    new_line, count = COPYRIGHT_WORD_RE.subn("Modifications Copyright", line, count=1)
    return new_line, count > 0


def insert_lines_at(lines, idx, new_lines, newline):
    new_lines = list(new_lines)
    if idx == 0 and lines and new_lines:
        bom, rest = split_leading_bom(lines[0])
        if bom and not split_leading_bom(new_lines[0])[0]:
            lines[0] = rest
            new_lines[0] = bom + new_lines[0]
    if idx > 0 and lines[idx - 1] and not lines[idx - 1].endswith(("\n", "\r")):
        lines[idx - 1] += newline
    for offset, extra in enumerate(new_lines):
        lines.insert(idx + offset, extra)
    return lines


def _span_text_for_copyright(content, spans):
    for match in COPYRIGHT_WORD_RE.finditer(content):
        for start, end in spans:
            if start <= match.start() and match.end() <= end:
                return content[start:end]
    return None


def _line_comment_run_end(lines, idx):
    end = idx
    while end + 1 < len(lines):
        current = split_leading_bom(lines[end])[1]
        nxt = split_leading_bom(lines[end + 1])[1]
        if _trailing_splice(current) or nxt.lstrip().startswith("//"):
            end += 1
            continue
        break
    return end


def _continued_line_comment(lines, idx):
    """True when lines[idx] continues a // comment joined by a backslash."""
    prev = idx - 1
    while prev >= 0:
        content = split_leading_bom(lines[prev])[1]
        if not _trailing_splice(content):
            return False
        if content.lstrip().startswith("//"):
            return True
        prev -= 1
    return False


def _hash_comment_run_end(lines, idx):
    end = idx
    while end + 1 < len(lines):
        content = split_leading_bom(lines[end + 1])[1]
        if content.lstrip().startswith("#"):
            end += 1
            continue
        break
    return end


def _block_comment_end(lines, start):
    """Line that closes a block comment opened on start, or None if it stays open."""
    in_block = False
    for idx in range(start, len(lines)):
        content = split_leading_bom(lines[idx])[1]
        _spans, in_block, _raw = comment_spans(content, False, in_block if idx > start else False)
        if not in_block:
            return idx
    return None


def copyright_comment_end(lines, ext, copy_idx, path=""):
    """Last line of the comment span that contains the copyright at copy_idx."""
    if ext == ".sh" or is_cmake_lists(path) or is_makefile(path):
        return _shell_comment_end(lines, copy_idx)
    if ext == ".py":
        doc_start, doc_end = leading_docstring_bounds(lines)
        if doc_start is not None and doc_start <= copy_idx and (doc_end is None or copy_idx <= doc_end):
            return doc_end
        return _hash_comment_run_end(lines, copy_idx)

    in_block = False
    raw_delim = None
    block_start = None
    for idx, raw in enumerate(lines):
        content = split_leading_bom(raw)[1]
        spans, still, raw_delim = comment_spans(content, False, in_block, raw_delim)
        if in_block:
            if not still and block_start is not None and block_start <= copy_idx <= idx:
                return idx
            in_block = still
            if not still:
                block_start = None
            continue
        if idx == copy_idx:
            snippet = _span_text_for_copyright(content, spans)
            if snippet is not None and snippet.startswith("//"):
                return _line_comment_run_end(lines, idx)
            if still:
                return _block_comment_end(lines, idx)
            if _continued_line_comment(lines, idx):
                return _line_comment_run_end(lines, idx)
            return idx
        if still:
            in_block = True
            block_start = idx
    if in_block and block_start is not None and block_start <= copy_idx:
        return None
    return copy_idx


def _banner_fits_header(insert_at, block):
    """True when the inserted Copyright line stays inside the checked header."""
    for offset, line in enumerate(block):
        if COPYRIGHT_WORD_RE.search(line):
            return insert_at + offset < HEADER_LINES
    return insert_at < HEADER_LINES


def insert_amd_modifications(lines, path, year, newline, copy_idxs):
    """Insert the AMD banner after the comment that holds the detected copyright.

    Returns (lines, error). error is set when the third-party comment is still
    open, or when the AMD Copyright line would land past the first HEADER_LINES
    lines. The file must not be written in either case.
    """
    block = new_amd_banner_lines(path, year, newline, modifications=True)
    ext = os.path.splitext(path)[1].lower()
    if copy_idxs:
        end = copyright_comment_end(lines, ext, min(copy_idxs), path)
        if end is None:
            return lines, (
                "third-party comment continues past the first %d lines; "
                "copyright was not changed" % HEADER_LINES
            )
        insert_at = end + 1
    else:
        insert_at = shebang_insert_index(lines) if hash_comment_ext(path) else 0
    if not _banner_fits_header(insert_at, block):
        return lines, (
            "AMD copyright would fall outside the first %d lines; "
            "copyright was not changed" % HEADER_LINES
        )
    return insert_lines_at(lines, insert_at, block, newline), None


def fix_copyright(path, year, mode):
    ext = os.path.splitext(path)[1].lower()
    if not copyright_eligible(path):
        return None
    if mode == "off":
        return None
    try:
        text, encoding, body = read_source(path)
    except OSError as exc:
        return "cannot read: %s" % exc

    newline = detect_newline(text)
    lines = text.splitlines(keepends=True)
    header = lines[:HEADER_LINES]
    info = analyze_copyright(header, ext, path)
    changed = False

    if not info["has_word"]:
        idx = shebang_insert_index(lines) if hash_comment_ext(path) else 0
        lines = insert_lines_at(lines, idx, new_amd_banner_lines(path, year, newline), newline)
        changed = True
    else:
        # Inject a year only onto an AMD banner that has none. Third-party
        # notices are not rewritten; a missing AMD banner is added later.
        for idx in info["amd_line_idxs"]:
            if years_on_copyright_line(lines[idx]):
                continue
            new_line, did = inject_year_after_copyright(lines[idx], year)
            if did:
                lines[idx] = new_line
                changed = True
        if info["amd_line_idxs"]:
            header = lines[:HEADER_LINES]
            info = analyze_copyright(header, ext, path)
            if any(idx not in info["amd_idxs"] for idx in info["amd_line_idxs"]):
                return "Copyright text found but year could not be parsed"

        header = lines[:HEADER_LINES]
        info = analyze_copyright(header, ext, path)
        if info["has_other"]:
            for idx in info["amd_idxs"]:
                new_line, did = ensure_modifications_copyright(lines[idx])
                if did:
                    lines[idx] = new_line
                    changed = True
            header = lines[:HEADER_LINES]
            info = analyze_copyright(header, ext, path)
        amd_years = info["amd_years"]
        if amd_years:
            for idx in info["amd_idxs"]:
                new_line, did = bump_last_year(lines[idx], year)
                if did:
                    lines[idx] = new_line
                    changed = True
        else:
            # Another copyright is present and there is no AMD line yet.
            lines, err = insert_amd_modifications(
                lines, path, year, newline, info["copy_idxs"]
            )
            if err:
                return err
            changed = True

    if not changed:
        return None
    try:
        write_source(path, "".join(lines), encoding, body)
    except (OSError, UnicodeEncodeError) as exc:
        return "cannot write: %s" % exc
    return None


HELP_DESCRIPTION = """\
Check or apply clang-format and AMD copyright banners.

Pass one or more paths. A path that is a file is checked as that file. A path
that is a directory is walked for source files. Paths you pass are resolved
from the current directory or the repository root. --git-staged and
--git-changed names are resolved from the repository root. The default is
check-only.

clang-format uses the repository .clang-format file (WebKit-based AOCL style)
on .c/.cc/.cpp/.cxx/.h/.hh/.hpp/.hxx. A check runs clang-format --dry-run
--Werror. --apply-clang-format and --apply rewrite those files with
clang-format -i. .py and .sh files are copyright-only and do not need clang-format.
Makefile and CMakeLists.txt are also copyright-only. Other .txt files, including
LICENSE.txt and NOTICE.txt, are left unchanged.

clang-format is required only when a selected file has a C or C++ extension
(.c/.cc/.cpp/.cxx/.h/.hh/.hpp/.hxx). The script looks for clang-format, then
clang-format-20 down to clang-format-10. If one of those files
is selected and no executable is on PATH, the script stops and tells you to
install clang-format and set PATH. --skip-clang-format skips that lookup.

Copyright is checked in the first 60 lines, and only inside a header comment
(a // line, a /* */ block, a leading # comment, or a leading Python module
docstring). Makefile and CMakeLists.txt use leading # comments. Other .txt
files, such as LICENSE.txt and NOTICE.txt, are not rewritten.
The comment must begin the logical line, so int x = 0; // Copyright
is not a banner. A later line of that block comment or module docstring
still counts. A comment after executable code has started is not a header.

Only the whole word Copyright counts, so prose such as copyrighted is not a
banner and is not rewritten. Prose that only mentions the word, such as
"copyright banners" or /* copyright banners */, is not a notice in C, C++,
Python, or shell. A # Copyright line inside a Python string, other than the
module docstring, is not a banner. A # line inside a shell heredoc is not a
banner. A shell # starts a comment only at the beginning of a word, so
foo#Copyright and ${value#Copyright} are not banners. A heredoc delimiter is
the whole shell word, so <<END-OF-FILE closes on END-OF-FILE. A one-line
module docstring that begins with Copyright counts. A plain, r, or u
triple-quoted string can be that docstring. An f-string or bytes literal,
such as f\"\"\"...\"\"\" or b\"\"\"...\"\"\", is not. A docstring notice later
in the line, such as "This file is Copyright (C) 2019, ...", also counts.
A Copyright word in source text is ignored. A C++ raw string is source text,
so a quote inside it does not make a later // Copyright into a banner. A line
that has the word both in a comment and in source text is left unchanged.

A Copyright banner is required. If no copyright is present, the AMD banner is:

    /*
        Copyright (C) YEAR, Advanced Micro Devices, Inc. All rights reserved.
    */

Python and shell use one line, after a #! line and a following encoding cookie:

    # Copyright (C) YEAR, Advanced Micro Devices, Inc. All rights reserved.

If another party's banner is already present (for example UT Austin or Intel),
a separate AMD comment is required and uses Modifications Copyright:

    /*
        Modifications Copyright (C) YEAR, Advanced Micro Devices, Inc. All rights reserved.
    */

If an AMD banner is present, its last 19xx/20xx year must be the current
calendar year or later. --apply-copyright and --apply bump 2024 to
2024-<current year>, insert a missing AMD banner, or add that Modifications
Copyright block after the comment that contains the other party's copyright.
An existing AMD Copyright line is prefixed with Modifications when another
copyright is already present. A different /* block, a // banner, or the
characters /* inside source text do not choose the insertion point. The AMD
line is not inserted inside the other banner. If that comment is still open
after the first 60 lines, or the new AMD Copyright line would land after
those 60 lines, the file is not rewritten. A third-party notice is not given
a year. A year is added only to an AMD line that lacks one. Every AMD line
needs its own 19xx/20xx year. One dated AMD line does not excuse another AMD
line that has no year.

A Python coding cookie must be a # comment in the first or second line.
Source text such as print('coding: utf-8') is not a cookie. The cookie stays
on that line; a new banner is inserted after it. An existing Modifications
Copyright line keeps that wording, and its last year is bumped when that year
is older than the current year. A current Copyright line does not leave a
stale Modifications Copyright year unchanged. --apply-copyright cannot be
combined with --skip-copyright.

Directory walks, --git-staged, and --git-changed skip tmp/, old/, 
.git/ directories. A path given as a file is still
checked. Copyright edits rewrite only the header and copy the rest of the
file unchanged, including bytes that are not UTF-8. New lines use the file's
first line ending. The replacement is written to a temporary file in the same
directory and takes the original name only after that write finishes, so a
failed write does not leave the source truncated. Permission bits are kept.
A symlink is not rewritten, and neither is a file reached through a symlinked
parent directory. A symlink to a directory is not walked, so the files in its
target stay unchanged.
"""

HELP_EPILOG = """\
examples:
  Check one file (format and copyright):
    python3 scripts/apply_format_changes.py src/map/lapack2flamec/FLA_bdsqr.c

  Check a directory:
    python3 scripts/apply_format_changes.py src/map/lapack2flamec

  Check a file and a directory together:
    python3 scripts/apply_format_changes.py src/foo.c src/map/lapack2flamec

  Apply clang-format and copyright:
    python3 scripts/apply_format_changes.py --apply src/map/lapack2flamec/FLA_bdsqr.c

  Apply only copyright:
    python3 scripts/apply_format_changes.py --apply-copyright --skip-clang-format path/to/file.c

  Apply only clang-format:
    python3 scripts/apply_format_changes.py --apply-clang-format --skip-copyright path/to/file.c

  Several files:
    python3 scripts/apply_format_changes.py --apply file1.c include/foo.h

  Unstaged, staged, and untracked files:
    python3 scripts/apply_format_changes.py --git-changed

  Pre-commit (files already added with git add):
    python3 scripts/apply_format_changes.py --git-staged

  The hook scripts/githooks/pre-commit rewrites eligible staged source files
  in the working tree on commit and leaves the index untouched. A partly
  staged file is formatted in full, including its unstaged edits. Deleted
  worktree files are skipped. If any file is rewritten, the commit is
  rejected so the changes can be reviewed and staged manually. A non-source
  file or a skipped path such as tmp/, old/, or CMakeFiles/ is not rewritten.
  A staged symlink is rejected and left as a link. A file under a symlinked
  directory is rejected before anything is rewritten.

  --git-staged, when you run this script yourself, reads the working-tree
  file. The hook does not. git commit --no-verify skips the hook.

  Git does not run scripts/githooks/pre-commit after a clone. From the
  repository root, once per clone:

    ln -sfn ../../scripts/githooks/pre-commit .git/hooks/pre-commit

  That links .git/hooks/pre-commit to the tracked hook. On Windows, if ln
  cannot create the link:

    cp -f scripts/githooks/pre-commit .git/hooks/pre-commit

  core.hooksPath is left unset. The hook runs python3, then python, then
  py -3.
"""


def _file_bytes(path):
    """Return the file contents, or None when it cannot be read."""
    try:
        with open(path, "rb") as handle:
            return handle.read()
    except OSError:
        return None


def main(argv=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=HELP_DESCRIPTION,
        epilog=HELP_EPILOG,
    )
    parser.add_argument(
        "paths",
        nargs="*",
        help="one or more files or directories; a file is checked as itself, a directory is walked",
    )
    parser.add_argument(
        "--git-staged",
        action="store_true",
        help="check files already added with git add (pre-commit). Does not include unstaged files that were not added, or untracked files",
    )
    parser.add_argument(
        "--git-changed",
        action="store_true",
        help="check unstaged changes, staged changes, and untracked files versus HEAD",
    )
    parser.add_argument(
        "--apply-clang-format",
        action="store_true",
        help="rewrite files with clang-format -i (copyright is not changed unless "
        "--apply-copyright or --apply is also set)",
    )
    parser.add_argument(
        "--apply-copyright",
        action="store_true",
        help="rewrite copyright: bump AMD years, insert a missing AMD banner, or add a "
        "separate Modifications Copyright block when another copyright "
        "is already present",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="apply both clang-format and copyright (same as "
        "--apply-clang-format --apply-copyright)",
    )
    parser.add_argument(
        "--skip-clang-format",
        action="store_true",
        help="do not run clang-format checks or --apply-clang-format",
    )
    parser.add_argument(
        "--skip-copyright",
        action="store_true",
        help="do not run copyright checks or --apply-copyright",
    )
    args = parser.parse_args(argv)
    year = datetime.date.today().year
    if args.apply:
        args.apply_clang_format = True
        args.apply_copyright = True
    if args.apply_copyright and args.skip_copyright:
        print("error: --apply-copyright cannot be used with --skip-copyright", file=sys.stderr)
        return 1
    start_cwd = os.getcwd()
    root = repo_root()
    os.chdir(root)

    paths = list(args.paths)
    if args.git_staged or args.git_changed:
        try:
            # Git names are relative to the repository root. Joining them here
            # keeps a same-named path in the current directory from winning.
            paths.extend(os.path.normpath(os.path.join(root, name)) for name in git_files(args))
        except (subprocess.CalledProcessError, OSError) as exc:
            print("git listing failed: %s" % exc, file=sys.stderr)
            return 1
    if not paths:
        print("no files to check")
        print("hint: pass a file or directory, or use --git-staged / --git-changed", file=sys.stderr)
        return 0

    files, path_errors = collect_source_files(paths, root, start_cwd=start_cwd)
    cf = None
    if not args.skip_clang_format and needs_clang_format(files):
        cf = clang_format_bin()
        if cf is None:
            notify_missing_clang_format(args.apply_clang_format)
            return 1
    problems = list(path_errors)
    bad_files = set()
    checked = 0
    changed_files = []
    applying = args.apply_clang_format or args.apply_copyright
    for path in files:
        rel = os.path.relpath(path, root)
        checked += 1
        before = _file_bytes(path) if applying else None
        if (args.apply_copyright or args.apply_clang_format) and symlink_blocks_rewrite(path, root):
            problems.append("%s: refusing to rewrite through a symlink" % rel)
            bad_files.add(rel)
            continue
        if not args.skip_copyright:
            if args.apply_copyright:
                err = fix_copyright(path, year, COPYRIGHT_MODE)
                if err:
                    problems.append("%s: copyright: %s" % (rel, err))
                    bad_files.add(rel)
            err = check_copyright(path, year, COPYRIGHT_MODE)
            if err:
                problems.append("%s: copyright: %s" % (rel, err))
                bad_files.add(rel)

        if cf and not args.skip_clang_format:
            err = check_format(path, cf, root, args.apply_clang_format)
            if err:
                problems.append("%s: format: %s" % (rel, err))
                bad_files.add(rel)

        if applying and _file_bytes(path) != before:
            changed_files.append(rel)

    if not checked and not path_errors:
        print("no matching source files")
        return 0

    if problems:
        nfiles = len(bad_files) + sum(1 for item in path_errors)
        print("%d issue(s) in %d file(s):" % (len(problems), nfiles))
        for item in problems:
            print("  " + item)
        return 1

    if applying:
        for rel in changed_files:
            print("  changed: %s" % rel)
        print("fixed: %d of %d file(s) changed" % (len(changed_files), checked))
    else:
        print("ok: %d file(s)" % checked)
    return 0


if __name__ == "__main__":
    sys.exit(main())
