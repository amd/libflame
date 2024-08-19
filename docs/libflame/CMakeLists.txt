# ##############################################################################
# Copyright (c) 2024 Advanced Micro Devices, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# ##############################################################################

# @file CMakeLists.txt
#
# @brief CMake based sphinx documentation generator.
#
# Checks if necessary packages are installed for generating doxygen & sphinx documentation
# if they are not installed then neccessary steps are shown to install the required packages.
#
#


cmake_minimum_required(VERSION 3.22.0)

project(aocllapack-doc)

find_package(Doxygen 1.9.6)

set(W_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Build Doxygen documentation
if(NOT (DOXYGEN_FOUND))
    message(FATAL_ERROR "Doxygen not found please install doxygen and try again.")
else()
    execute_process(COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
            WORKING_DIRECTORY ${W_DIR}
            COMMAND_ECHO STDOUT)
endif()


function(check_rocm_docs_version)

    execute_process(COMMAND python3 -m pip list
        RESULT_VARIABLE     PIP3_EXIT_CODE
        OUTPUT_VARIABLE     PIP3_LIST
        ERROR_VARIABLE      PIP3_ERROR)

    if(PIP3_EXIT_CODE)
        if(WIN32)
            find_program(PIP3_EXECUTABLE pip3)
            if(NOT PIP3_EXECUTABLE)
                message(FATAL_ERROR "pip3 not found")
            endif()
        endif()
        message(FATAL_ERROR "${PIP3_ERROR}")
    endif()
    # Filter rocm-docs-core among list of python3 libraries
    string(REGEX MATCH "rocm-docs-core[ \t\r]*[0-9]+[\.]?[0-9]*[\.]?[0-9]*" ROCM_DOCS_CORE_FILTERED "${PIP3_LIST}")
    # Extract version of rocm-docs-core and save it in ROCM_DOCS_VERSION
    string(REGEX REPLACE "rocm-docs-core[ \t\r]*" "" ROCM_DOCS_VERSION "${ROCM_DOCS_CORE_FILTERED}")

    if(${ROCM_DOCS_VERSION} VERSION_GREATER "0.33.0")
        message(FATAL_ERROR "Current version of rocm-docs-core is: (${ROCM_DOCS_VERSION}) which is greater than 0.32.0,\
                  please use rocm-docs-core version between 0.19.0 to 0.32.0")
    elseif(${ROCM_DOCS_VERSION} VERSION_LESS "0.19.0")
        message(FATAL_ERROR "Current version of rocm-docs-core is: (${ROCM_DOCS_VERSION}) which is less than 0.19.0,\
                  please use rocm-docs-core version between 0.19.0 to 0.32.0")
    endif()
endfunction()



# Check if Doxygen documentation generated successfully, if yes continue building sphinx documentation
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/html/index.html)
    message(STATUS "Doxygen documentation generated successfully, to view documentation open ${W_DIR}/html/index.html .")


# Check for python and the documentation packages:
# - Sphinx
# - rocm-docs-core
# - breathe

    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import sphinx"
      RESULT_VARIABLE EXIT_CODE
      OUTPUT_QUIET)
    if(NOT ${EXIT_CODE} EQUAL 0)
      message(
        FATAL_ERROR
          "Python3 package sphinx is required to build the documentation\n   Use: ${Python3_EXECUTABLE} -m pip3 install sphinx"
      )
    endif()
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "from rocm_docs import ROCmDocs"
      RESULT_VARIABLE EXIT_CODE
      OUTPUT_QUIET)
    if(NOT ${EXIT_CODE} EQUAL 0)
      message(
        FATAL_ERROR
          "Python3 package rocm_docs is required to build the documentation\n   Use: ${Python3_EXECUTABLE} -m pip3 install rocm-docs-core==0.30.0"
      )
    endif()

    # checks rocm-docs-core version is in between 0.19.0 & 0.32.0
    check_rocm_docs_version()

    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import breathe"
      RESULT_VARIABLE EXIT_CODE
      OUTPUT_QUIET)
    if(NOT ${EXIT_CODE} EQUAL 0)
      message(
        FATAL_ERROR
          "Python3 package breathe is required to build the documentation\n   Use: ${Python3_EXECUTABLE} -m pip3 install breathe"
      )
    endif()

    include (FindPackageHandleStandardArgs)

    ## Look for an executable called sphinx-build
    find_program (SPHINX_EXECUTABLE
                  NAMES sphinx-build sphinx-build.exe
                  DOC "Sphinx Documentation Generator"
    )

    if(NOT SPHINX_EXECUTABLE)
      message(FATAL_ERROR "Failed to find sphinx-build executable.")
    endif()

    message(STATUS "Building Sphinx documentation.")

    execute_process(
      COMMAND ${SPHINX_EXECUTABLE} -M html ./rst/ ./sphinx
      RESULT_VARIABLE EXIT_CODE
      WORKING_DIRECTORY ${W_DIR})

    if(NOT ${EXIT_CODE} EQUAL 0)
      message(
        FATAL_ERROR
          "Sphinx document generation failed: ${EXIT_CODE}"
      )
    else()
      message(STATUS "Sphinx generated successfully at ${W_DIR}/sphinx/html")
      file(REMOVE_RECURSE ${W_DIR}/doxygen/xml)
      file(REMOVE_RECURSE ${W_DIR}/sphinx/doctrees)
    endif()
else() # !(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/docs/html/index.html)
    message(STATUS "Doxygen document generation failed.")
endif() # (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/docs/html/index.html)