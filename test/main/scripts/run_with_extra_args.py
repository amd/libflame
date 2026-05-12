#!/usr/bin/env python3
###############################################################################
# Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
"""Runtime CTest wrapper for the main LAPACK test executable.

CTest runs: ``python3 run_with_extra_args.py <exe> <per-test args...>``.
Read ``test/main/extra_args.txt`` (parent of ``scripts/``), skip ``#`` comment
lines, append parsed tokens after the per-test arguments, then ``exec`` the
binary—same contract as a thin shell wrapper: global flags such as
``--interface=`` or ``--test-mode=`` apply to yaml_generated and other main-suite
CTests without CMake reconfigure (see ``test/main/ReadMe.txt``).
"""

import os
import shlex
import sys
from pathlib import Path


def _read_extra_args(path: Path) -> list:
    """Load optional CLI tokens from ``extra_args.txt``.

    If ``path`` is missing or not a regular file, returns an empty list.
    Otherwise reads UTF-8 text, ignores blank lines and lines whose first
    non-whitespace character is ``#``, and extends the result with tokens
    from each remaining line using POSIX ``shlex`` rules (``posix=True``).

    Args:
        path: Filesystem path to ``extra_args.txt`` (typically beside ``scripts/``).

    Returns:
        Flat list of argument strings to append after per-test argv.
    """
    if not path.is_file():
        return []
    out = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            out.extend(shlex.split(line, posix=True))
    return out


def main() -> int:
    if len(sys.argv) < 2:
        sys.stderr.write(
            "usage: run_with_extra_args.py <executable> [test arguments...]\n"
        )
        return 2

    exe = sys.argv[1]
    test_args = sys.argv[2:]

    extra_path = Path(__file__).resolve().parent.parent / "extra_args.txt"
    full_argv = [exe] + test_args + _read_extra_args(extra_path)

    try:
        os.execv(exe, full_argv)
    except OSError as exc:
        sys.stderr.write(f"run_with_extra_args: cannot exec {exe!r}: {exc}\n")
        return 127


if __name__ == "__main__":
    raise SystemExit(main())
