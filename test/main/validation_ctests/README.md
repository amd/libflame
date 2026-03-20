##########################################################################################
Validation CTests - YAML-Based Test Generation
##########################################################################################

# Introduction

This document explains how to use the YAML-based test generation system for creating
CMake CTest definitions for LAPACK API validation tests. This system allows developers
to define test cases in YAML format, which are then automatically converted to CMake
test files during the build process.

The validation_ctests directory has the following contents:
   1. README.md - This file
   2. auto_generate_tests.py - Python script that generates CMake test files from YAML
   3. auto_generate_label_groups.yaml - Global label group definitions and size threshold configuration
   4. <API>.yaml - Test definitions for each API. Each API has its own YAML file.

## Overview

The test generation workflow follows this pattern:

   YAML Files → Python Generator Script → CMake CTest Files → CTest Execution

1. YAML Files (*.yaml): Define test parameters and code paths for each API
2. Generator Script (auto_generate_tests.py): Reads YAML files and generates CMake test files
3. CMake Integration: Tests are generated at configure time in the build directory
4. CTest Execution: Generated tests are executed via CMake's CTest framework

## YAML File Format

Each API has its own YAML file (e.g., getrf.yaml, syev.yaml). The structure follows
this pattern:

### Structure

   api_name:
     description: "Brief description of the API"
     api_group: "API group (e.g., LIN, EIG)"
     labels: ["presubmit"]  # Optional: API-level labels
     precisions: [s, d, c, z]  # Supported precisions
     param_order: [param1, param2, ...]  # Order of parameters for test command

     custom_paths:
       # Code path: src/path/to/file.c (s)
       - name: test_name
         precisions: sdcz  # works with individual or Combined precision syntax
         labels: []  # Optional: Path-specific labels (overrides API-level if non-empty)
         params: {param1: value1, param2: value2, ...}

## Labels

Labels allow filtering and organizing tests. Hierarchy: path-level (takes precedence) >
API-level (if path labels empty) > default labels (API name, API group, precision).

Common labels: presubmit, avx2, avx512, precision_s/d/c/z (auto-generated).

Filter tests:
   # Single label (case-insensitive for API name/group)
   $ ctest -L presubmit
   $ ctest -L avx2
   $ ctest -L GETRF        # Uppercase works
   $ ctest -L getrf        # Lowercase also works (both cases supported)

   # Multiple labels (use multiple -L flags, not semicolon)
   $ ctest -L GETRF -L small    # Tests with both GETRF and small labels
   $ ctest -L avx2 -L precision_d  # AVX2 tests with double precision

### Label Groups

Define label groups in auto_generate_label_groups.yaml to automatically include multiple labels across
all API files. This avoids manually adding multiple labels to each test.

Example in auto_generate_label_groups.yaml:
   label_groups:
     postsubmit: [avx2, avx512]

When any test uses the "postsubmit" label, it automatically includes both "postsubmit"
and its constituent labels (avx2, avx512). This allows managing labels centrally for
multiple API files without editing each file individually. API-specific label_groups can
override global groups if needed.

### Automatic Size-Based Labeling

Tests are automatically assigned size labels (small, medium, large) based on parameter
values defined in auto_generate_label_groups.yaml. Group labels (e.g., "short") are then automatically
assigned to tests matching certain size criteria.

Configuration in auto_generate_label_groups.yaml:
   size_thresholds:
     parameter: max  # Check 'n', 'm', or 'max' (max of m and n)
     small: 16      # n < 16 or max(m,n) < 16
     medium: 100    # 16 <= n < 100 or 16 <= max(m,n) < 100
     # large is automatically assigned if >= medium threshold

   auto_assign_groups:
     short: [small, medium]  # Automatically add "short" label to tests with small or medium size

Benefits:
- Zero manual edits: Tests automatically get size and group labels based on parameters
- Centralized configuration: All thresholds in auto_generate_label_groups.yaml
- Works across all API files automatically
- Flexible: Can define different thresholds for different parameters

Example: A test with m=16, n=16 automatically gets "medium" and "short" labels without
any manual editing in the API YAML file.


## Adding New APIs

To add tests for a new API, follow these steps:

1. Create YAML File
   Create <api_name>.yaml in the validation_ctests folder. See getrf.yaml for reference.

2. Determine Parameter Order
   Check test/main/src/test_<api>.c to determine the correct parameter order.

3. Define Test Paths
   Add test cases covering different code paths: single element, small/large sizes,
   architecture-specific optimizations (AVX2, AVX512), and edge cases.

4. Verify Generation
   $ cd test/main/validation_ctests
   $ python3 auto_generate_tests.py --api <api_name> /tmp/test_output

5. Build and Test
   The CMake build system will automatically generate test files during configuration
   when BUILD_TEST=ON. Tests are generated in <build_dir>/test/main/validation_ctests/

## Adding Test Paths to Existing APIs

To add a new test case to an existing API:

1. Open the YAML file (e.g., getrf.yaml)
2. Add code path comments documenting which source files are exercised
3. Add the test definition in custom_paths section
4. Verify by regenerating tests:
   $ python3 auto_generate_tests.py --api getrf /tmp/test_output

## Running Tests

Tests are automatically generated during CMake configuration when BUILD_TEST=ON.
Generated files are in <build_dir>/test/main/validation_ctests/*_tests.cmake

For manual generation:
   $ cd test/main/validation_ctests
   $ python3 auto_generate_tests.py --api getrf /tmp/output
   $ python3 auto_generate_tests.py --all /tmp/output

Run tests:
   # Filter by label (case-insensitive for API name and group)
   $ ctest -L GETRF        # Matches tests with GETRF label
   $ ctest -L getrf        # Also matches (both cases supported)
   $ ctest -L GETRF -L small  # Multiple labels (use multiple -L flags)
   $ ctest -L avx2         # Architecture-specific tests

   # Filter by test name regex (uppercase API name)
   $ ctest -R GETRF        # Matches all GETRF tests
   $ ctest -R "GETRF.*kernel"  # Matches GETRF kernel tests

   # Combined filtering
   $ ctest -L GETRF -R ".*kernel.*"  # GETRF tests with "kernel" in name

## Example

Generated CMake test files are in
<build_dir>/test/main/validation_ctests/getrf_tests.cmake after configuration.

## Integration with CMake

The test generation is integrated into the CMake build system via test/main/CMakeLists.txt.
Tests are generated at configure time in the build directory and are not committed to
the source tree.
