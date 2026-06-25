# ########################################################################
# Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
# ########################################################################

# CompilerFlags.cmake
# This module handles all compiler flag configuration for libflame
# including ISA settings, optimization, security, and platform-specific flags

# ============================================================================
# ISA Configuration Functions
# ============================================================================

# Function to validate and normalize ISA configuration
function(validate_and_normalize_isa_config)
    # Make LF_ISA_CONFIG case insensitive
    string(TOLOWER "${LF_ISA_CONFIG}" lf_isa_config_lower)
    
    set(VALID_ISA_CONFIGS "none" "avx" "avx2" "avx2-strict" "avx512" "avx512-strict" "auto")
    
    if(NOT lf_isa_config_lower IN_LIST VALID_ISA_CONFIGS)
        message(FATAL_ERROR
            "Machine ISA configuration error, valid values are "
            "none, avx, avx2, avx2-strict, avx512, avx512-strict, auto")
    endif()
    
    message(STATUS "Machine ISA configuration (LF_ISA_CONFIG) set to ${LF_ISA_CONFIG}")
    
    # Export to parent scope
    set(lf_isa_config_lower ${lf_isa_config_lower} PARENT_SCOPE)
endfunction()

# Function to auto-detect ISA configuration
function(auto_detect_isa_config)
    # Unified-build fix: use ${FLAME_ROOT} (libflame's own source root) instead
    # of ${CMAKE_SOURCE_DIR}, which points at the top-level AOCL project when
    # libflame is brought in via add_subdirectory.
    set(AUTO_CONFIG_PY "${FLAME_ROOT}/build/auto_config.py")
    set(AUTO_CONFIG_ENV)

    if(CMAKE_C_COMPILER)
        list(APPEND AUTO_CONFIG_ENV "CC=${CMAKE_C_COMPILER}")
        list(APPEND AUTO_CONFIG_ENV "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
    endif()
    
    # Run python script to detect ISA support using CPUID.
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env ${AUTO_CONFIG_ENV} ${Python3_EXECUTABLE} ${AUTO_CONFIG_PY}
        RESULT_VARIABLE CMD_RESULT
        OUTPUT_VARIABLE CMD_OUTPUT
        ERROR_VARIABLE CMD_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    # Check if the command failed
    if(NOT CMD_RESULT EQUAL 0)
        message(FATAL_ERROR
            "Auto-detecting ISA configuration failed with exit code ${CMD_RESULT}.\n"
            "Command: ${CMAKE_COMMAND} -E env ${AUTO_CONFIG_ENV} ${Python3_EXECUTABLE} ${AUTO_CONFIG_PY}\n"
            "stdout: ${CMD_OUTPUT}\n"
            "stderr: ${CMD_ERROR}"
        )
    endif()
    
    message(STATUS "Auto configuring ISA using CPUID helper: ${CMD_OUTPUT}")
    
    if("${CMD_OUTPUT}" STREQUAL "avx512")
        set(lf_isa_config_lower "avx512" PARENT_SCOPE)
    elseif("${CMD_OUTPUT}" STREQUAL "avx2")
        set(lf_isa_config_lower "avx2" PARENT_SCOPE)
    elseif("${CMD_OUTPUT}" STREQUAL "none")
        set(lf_isa_config_lower "none" PARENT_SCOPE)
    else()
        message(FATAL_ERROR
            "Auto-detecting ISA configuration returned an unknown value: ${CMD_OUTPUT}")
    endif()
endfunction()

# Function to handle strict ISA variants
function(handle_strict_isa_variants ISA_CONFIG_VAR)
    set(isa_config ${${ISA_CONFIG_VAR}})
    
    # Accept fixed/forced variants and map them to a compile-time definition
    if("${isa_config}" STREQUAL "avx2-strict")
        message(STATUS "LF_ISA_CONFIG requests a compile-time forced AVX2 ISA")
        # Define compile-time forced ISA macro for C code (FLA_Context.c expects FLA_STRICT_ARCH)
        add_definitions(-DFLA_STRICT_ARCH=FLA_ARCH_AVX2)
        # Normalize to base token so later compiler flag logic still applies
        set(isa_config "avx2")
    elseif("${isa_config}" STREQUAL "avx512-strict")
        message(STATUS "LF_ISA_CONFIG requests a compile-time forced AVX512 ISA")
        add_definitions(-DFLA_STRICT_ARCH=FLA_ARCH_AVX512)
        set(isa_config "avx512")
    endif()
    
    set(${ISA_CONFIG_VAR} ${isa_config} PARENT_SCOPE)
endfunction()

# Function to set ISA-specific compiler flags
function(set_isa_compiler_flags ISA_CONFIG)
    if(WIN32)
        if("${ISA_CONFIG}" STREQUAL "avx512")
            add_compile_options(/arch:AVX512)
        elseif("${ISA_CONFIG}" STREQUAL "avx2" OR "${ISA_CONFIG}" STREQUAL "avx")
            add_compile_options(/arch:AVX2)
        else()
            message(STATUS "No ISA flag set")
        endif()
    elseif(UNIX)
        set(COMPILER_OPTIMIZATION_FLAGS "-mtune=native -O3")
        
        if("${ISA_CONFIG}" STREQUAL "avx512")
            set(COMPILER_OPTIMIZATION_FLAGS "${COMPILER_OPTIMIZATION_FLAGS} -mavx512f -mavx512dq -mfma")
        elseif("${ISA_CONFIG}" STREQUAL "avx2" OR "${ISA_CONFIG}" STREQUAL "avx")
            set(COMPILER_OPTIMIZATION_FLAGS "${COMPILER_OPTIMIZATION_FLAGS} -mavx2 -mfma")
        else()
            message(STATUS "No ISA flag set")
        endif()
        
        # Export to parent scope
        set(COMPILER_OPTIMIZATION_FLAGS ${COMPILER_OPTIMIZATION_FLAGS} PARENT_SCOPE)
    else()
        message(STATUS "OS UNKNOWN CANNOT SET SIMD")
    endif()
endfunction()

# ============================================================================
# Compiler Flag Setup Functions
# ============================================================================

# Function to set security and hardening flags (Unix only)
function(set_security_flags)
    if(UNIX)
        set(COMPILE_FLAGS "-fstack-protector-strong -fpie -Wformat -Wformat-security" PARENT_SCOPE)
        set(LINKER_FLAGS "-Wl,-z,relro -Wl,-z,now" PARENT_SCOPE)
    endif()
endfunction()

# Function to set warning flags
function(set_warning_flags)
    if(UNIX)
        set(GCC_WARNING_FLAGS "-Wall -Wno-comment" PARENT_SCOPE)
    endif()
endfunction()

# Function to set language standard flags
function(set_language_flags)
    if(UNIX)
        set(COMPILER_LANGUAGE_FLAGS "-std=c11 -D_GNU_SOURCE -Wno-unused-function -Wno-parentheses -Wfatal-errors" PARENT_SCOPE)
    endif()
endfunction()

# Function to set debug flags
function(set_debug_flags)
    if(UNIX)
        set(COMPILER_DEBUG_FLAG "-g0" PARENT_SCOPE)
    endif()
endfunction()

# Function to handle Clang-specific flags
function(set_clang_specific_flags)
    # Check C compiler ID for C flags
    if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-parentheses -Wno-deprecated-declarations -Wno-macro-redefined" PARENT_SCOPE)
    endif()
    # Check CXX compiler ID for CXX flags
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-parentheses -Wno-deprecated-declarations -Wno-macro-redefined" PARENT_SCOPE)
    endif()
endfunction()

# Function to set GCOV flags
function(set_gcov_flags)
    if(ENABLE_GCOV)
        set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON PARENT_SCOPE)
        add_compile_options(--coverage)
        add_link_options(--coverage)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage" PARENT_SCOPE)
    endif()
endfunction()

# Function to set Address Sanitizer flags
function(set_asan_flags)
    if(WIN32)
        if(ENABLE_ASAN)
            message(STATUS "Address Sanitizer is not supported on Windows")
        endif()
    elseif(UNIX)
        message(STATUS "Enabled ENABLE_ASAN: ${ENABLE_ASAN}")
        if(ENABLE_ASAN)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fsanitize=address" PARENT_SCOPE)
        endif()
    endif()
endfunction()

# Function to set PIC and other Unix-specific flags
function(set_unix_pic_flags)
    if(UNIX)
        add_compile_options(-fPIC -fmacro-prefix-map=${FLAME_ROOT}=.)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_FORTIFY_SOURCE=2" PARENT_SCOPE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORTIFY_SOURCE=2" PARENT_SCOPE)
    endif()
endfunction()

# Function to set Windows-specific compiler flags
function(set_windows_compiler_flags)
    if(WIN32)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP" PARENT_SCOPE)
    endif()
endfunction()

# ============================================================================
# Main Configuration Function
# ============================================================================

# Main macro to configure all compiler flags (using macro instead of function for proper scope)
macro(configure_compiler_flags)
    # Validate and normalize ISA configuration
    validate_and_normalize_isa_config()
    
    # Handle auto-detection if needed
    if(lf_isa_config_lower STREQUAL "auto")
        auto_detect_isa_config()
    endif()
    
    message(STATUS "LF_ISA_CONFIG selected: ${lf_isa_config_lower}")
    
    # Handle strict ISA variants
    handle_strict_isa_variants(lf_isa_config_lower)
    
    # Set ISA-specific compiler flags
    set_isa_compiler_flags(${lf_isa_config_lower})
    
    # Set security flags
    set_security_flags()
    
    # Set warning flags
    set_warning_flags()
    
    # Set language flags
    set_language_flags()
    
    # Set debug flags
    set_debug_flags()
    
    # Apply all Unix flags
    if(UNIX)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMPILER_OPTIMIZATION_FLAGS} ${COMPILE_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_OPTIMIZATION_FLAGS} ${COMPILE_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LINKER_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LINKER_FLAGS}")
        
        # Apply language and warning flags
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMPILER_DEBUG_FLAG} ${GCC_WARNING_FLAGS} ${COMPILER_LANGUAGE_FLAGS}")
    endif()
    
    # Set Clang-specific flags
    set_clang_specific_flags()
    
    # Set GCOV flags
    set_gcov_flags()
    
    # Set ASAN flags
    set_asan_flags()
    
    # Set Unix PIC flags
    set_unix_pic_flags()
    
    # Set Windows compiler flags
    set_windows_compiler_flags()
endmacro()
