"""
Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
Name: auto_config.py
Purpose: To check CPU ISA support for auto configuration
"""

#Global Imports
import os
import platform
import shutil
import subprocess
import sys
import tempfile


def helper_source():
    """Return the C source used to detect host AVX ISA support."""
    return """
#include <stdio.h>

#define XCR0_XMM       (1ull << 1)
#define XCR0_YMM       (1ull << 2)
#define XCR0_OPMASK    (1ull << 5)
#define XCR0_ZMM_HI256 (1ull << 6)
#define XCR0_HI16_ZMM  (1ull << 7)

#define XCR0_AVX2_MASK   (XCR0_XMM | XCR0_YMM)
#define XCR0_AVX512_MASK (XCR0_XMM | XCR0_YMM | XCR0_OPMASK | \\
                          XCR0_ZMM_HI256 | XCR0_HI16_ZMM)

#if defined(_MSC_VER)
#include <intrin.h>
static void cpuid_count(unsigned int leaf, unsigned int subleaf,
                        unsigned int regs[4])
{
    int cpu_info[4];
    __cpuidex(cpu_info, (int)leaf, (int)subleaf);
    regs[0] = (unsigned int)cpu_info[0];
    regs[1] = (unsigned int)cpu_info[1];
    regs[2] = (unsigned int)cpu_info[2];
    regs[3] = (unsigned int)cpu_info[3];
}

static unsigned long long xgetbv0(void)
{
    return _xgetbv(0);
}
#else
#include <cpuid.h>
static void cpuid_count(unsigned int leaf, unsigned int subleaf,
                        unsigned int regs[4])
{
    __cpuid_count(leaf, subleaf, regs[0], regs[1], regs[2], regs[3]);
}

static unsigned long long xgetbv0(void)
{
    unsigned int eax;
    unsigned int edx;
    __asm__ volatile("xgetbv" : "=a"(eax), "=d"(edx) : "c"(0));
    return ((unsigned long long)edx << 32) | eax;
}
#endif

int main(void)
{
    unsigned int regs[4] = {0, 0, 0, 0};
    unsigned int max_leaf;
    unsigned long long xcr0;

    cpuid_count(0, 0, regs);
    max_leaf = regs[0];
    if(max_leaf < 7)
    {
        puts("none");
        return 0;
    }

    /* CPUID.1:ECX[27] = OSXSAVE. Without it, XCR0 is not valid here. */
    cpuid_count(1, 0, regs);
    if(!(regs[2] & (1u << 27)))
    {
        puts("none");
        return 0;
    }

    /* CPUID reports hardware ISA support; XCR0 reports OS-enabled state. */
    xcr0 = xgetbv0();
    cpuid_count(7, 0, regs);
    /* CPUID.7.0:EBX[16] = AVX512F, require XCR0 bits 1,2,5,6,7. */
    if((regs[1] & (1u << 16)) &&
       ((xcr0 & XCR0_AVX512_MASK) == XCR0_AVX512_MASK))
    {
        puts("avx512");
        return 0;
    }
    /* CPUID.7.0:EBX[5] = AVX2, require XCR0 bits 1,2. */
    if((regs[1] & (1u << 5)) &&
       ((xcr0 & XCR0_AVX2_MASK) == XCR0_AVX2_MASK))
    {
        puts("avx2");
        return 0;
    }

    puts("none");
    return 0;
}
"""


def build_helper(cc, source, output_path):
    """Compile the CPUID helper with the selected C compiler."""
    compiler_name = os.path.basename(cc).lower()

    if compiler_name in ("cl", "cl.exe", "clang-cl", "clang-cl.exe"):
        command = [cc, "/nologo", source, "/Fe:%s" % output_path]
    else:
        command = [cc, source, "-o", output_path]

    return subprocess.run(command, cwd=os.path.dirname(output_path),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True)


def host_is_x86():
    """Return whether the build host CPU is an x86-family processor."""
    machine = platform.machine().lower()
    return machine in ("x86_64", "amd64", "i386", "i486", "i586", "i686",
                       "x86")


def config_check():                                                             #Function to detect the ISA token
    """Detect and return the best supported ISA token for this host."""
    if not host_is_x86():
        return "none"

    cc = os.environ.get("CC") or os.environ.get("CMAKE_C_COMPILER") or \
        shutil.which("gcc") or shutil.which("clang") or shutil.which("cc") \
        or shutil.which("cl") or shutil.which("clang-cl")

    if not cc:
        raise RuntimeError("Unable to auto-detect ISA: no C compiler found")

    with tempfile.TemporaryDirectory(prefix="libflame-cpuid-") as tmpdir:
        source = os.path.join(tmpdir, "detect_isa.c")
        exe_name = "detect_isa.exe" if 'Windows' in platform.system() \
            else "detect_isa"
        output_path = os.path.join(tmpdir, exe_name)

        with open(source, "w", encoding="utf-8") as source_file:
            source_file.write(helper_source())

        build_result = build_helper(cc, source, output_path)
        if build_result.returncode != 0:
            raise RuntimeError("Unable to build CPUID helper with %s\n%s" %
                               (cc, build_result.stderr))

        run_result = subprocess.run([output_path], stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)
        if run_result.returncode != 0:
            raise RuntimeError("CPUID helper failed\n%s" % run_result.stderr)

        return run_result.stdout.strip()


if __name__ == "__main__":
    try:
        print(config_check())                                                   #Function call for ISA token
    except Exception as e:
        print("Exception due to %s" % e, file=sys.stderr)
        sys.exit(1)
