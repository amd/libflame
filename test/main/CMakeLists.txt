###############################################################################
# Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################
if(WIN32)
    project(test_libFLAME_main)
elseif(UNIX)
    project(test_lapack.x)
endif()
message(STATUS "TEST EXE NAME ${PROJECT_NAME}")

option (ENABLE_ILP64 "Enable ILP64" OFF)

option (ENABLE_MEM_UNALIGNED "Allocate unaligned memory" OFF)

option (ENABLE_CPP_TEST "Enable CPP test support in test suite" ON)

if (${ENABLE_CPP_TEST} MATCHES "ON")
add_definitions(-DENABLE_CPP_TEST=1)
else()
add_definitions(-DENABLE_CPP_TEST=0)
endif()
message(STATUS "Enabled ENABLE_CPP_TEST : ${ENABLE_CPP_TEST}")

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin") 
if(WIN32)
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function")
elseif(UNIX)
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

if(WIN32)
    set(OpenMP_libomp_LIBRARY "${OpenMP_libomp_LIB}" CACHE STRING "openmp library path")
    set(TEST_FLAGS "${OpenMP_libomp_LIBRARY}")
elseif(UNIX)
    set(OpenMP_libomp_LIBRARY "-fopenmp" CACHE STRING "openmp library path")
    set(Test_C_Flags " -std=c11 -D_GNU_SOURCE -O2 -Wall -Wno-comment -Wno-unused-function" CACHE STRING "C Flags")
    set(Test_LD_Flags "-lm -lpthread" CACHE STRING "LD FLags")
    if(ENABLE_ASAN)
        set(Test_LD_Flags "${Test_LD_Flags} -fsanitize=address")
    endif()
    if(EXT_OPENMP_PATH)
    	set(TEST_FLAGS "${LIBIOMP} ${OpenMP_libomp_LIBRARY} ${Test_C_Flags} ${Test_LD_Flags}")
    else()
    	set(TEST_FLAGS "${OpenMP_libomp_LIBRARY} ${Test_C_Flags} ${Test_LD_Flags}")
    endif()
    message(STATUS "${TEST_FLAGS}")
else()
message(STATUS "Openmp: OS NOT found")
endif()

set(BLAS_HEADER_PATH "" CACHE STRING "BLAS library header path")

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
remove_definitions(-DBLIS1_FROM_LIBFLAME)

if(ENABLE_ILP64)
    add_definitions(-DFLA_ENABLE_ILP64)
endif()

#
# --- Flag to enable uppercase/lowercase API calling conventions --------------
#     Set API_CALL_CONVENTION = upper  for uppercase without underscore
#                             = upper_ for uppercase with underscore
#                             = lower  for lowercase without underscore
#                             = default for default case (lower_)
#     NOTE: LAPACK default is lowercase with underscore
#
set(API_CALL_CONVENTION "default" CACHE STRING "API calling convention")

#
# -- API_CALLING_CONVENTION set to "upper" enables uppercase w/o underscore ----
#
if (${API_CALL_CONVENTION} MATCHES "upper")
   add_definitions(-DUPPER=1)
endif()
#
# -- API_CALLING_CONVENTION set to "upper_" enables uppercase w/ underscore ----
#
if (${API_CALL_CONVENTION} MATCHES "upper_")
   add_definitions(-DUPPER_=1)
endif()
#
# -- API_CALLING_CONVENTION set to "lower" enables lowercase w/o underscore ----
#
if (${API_CALL_CONVENTION} MATCHES "lower")
   add_definitions(-DLOWER=1)
endif()

if(ENABLE_MEM_UNALIGNED)
    add_definitions(-DFLA_MEM_UNALIGNED)
endif()

add_executable("${PROJECT_NAME}" "")

message("${CMAKE_BUILD_TYPE}")

if(UNIX)
    foreach(BLAS_LIB ${BLAS_LINK_LIST})
      target_link_libraries("${PROJECT_NAME}" PRIVATE "${TEST_FLAGS}" "${CMAKE_PROJECT_NAME}" ${BLAS_LIB})
    endforeach()
    if(NOT ENABLE_EMBED_AOCLUTILS)
      target_link_libraries("${PROJECT_NAME}" PRIVATE "${TEST_FLAGS}" "${CMAKE_PROJECT_NAME}" aoclutils)
      target_link_libraries("${PROJECT_NAME}" PRIVATE "${TEST_FLAGS}" "${CMAKE_PROJECT_NAME}" "libstdc++.so")
    endif()
elseif(WIN32)
    target_link_libraries("${PROJECT_NAME}" debug "${TEST_FLAGS}" "${CMAKE_PROJECT_NAME}" blas)
    target_link_libraries("${PROJECT_NAME}" optimized "${TEST_FLAGS}" "${CMAKE_PROJECT_NAME}" blas)
    if(NOT ENABLE_EMBED_AOCLUTILS)
      target_link_libraries("${PROJECT_NAME}" debug "${TEST_FLAGS}" "${CMAKE_PROJECT_NAME}" aoclutils)
      target_link_libraries("${PROJECT_NAME}" optimized "${TEST_FLAGS}" "${CMAKE_PROJECT_NAME}" aoclutils)
    endif()
endif()

if(BLAS_HEADER_PATH)
   target_include_directories("${PROJECT_NAME}" PRIVATE ${BLAS_HEADER_PATH})
else(BLAS_HEADER_PATH)
   message(STATUS "--> BLAS header file path not provided")
endif()

# Compile options to ignore warnings in Windows.
# TODO: warnings are suppressed temporarily & need to rework.
if(WIN32)
target_compile_options(${PROJECT_NAME} PRIVATE
                        -Wno-comment -Wno-unused-function
                        -Wno-macro-redefined -Wno-parentheses-equality
                        -Wno-format)
endif()

target_include_directories("${PROJECT_NAME}" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories("${PROJECT_NAME}" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/validate_src)
if(ENABLE_CPP_TEST)
target_include_directories("${PROJECT_NAME}" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)
endif()

add_subdirectory(src)
add_subdirectory(validate_src)

include(main_ctest.cmake)
include(conformance_test.cmake)
include(overflow_underflow_test.cmake)
