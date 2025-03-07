# ########################################################################
#Copyright(c) 2023-2025 Advanced Micro Devices, Inc.
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files(the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following                       conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the                               Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.
#
# ########################################################################

# ============= aocl function ================
function(aocl_libs)

  IF(FLA_ENABLE_ILP64) 
    SET(ILP_DIR "ILP64")
  ELSE(FLA_ENABLE_ILP64)
    SET(ILP_DIR "LP64")
  ENDIF(FLA_ENABLE_ILP64)

  IF(WIN32)
    SET(CMAKE_FIND_LIBRARY_PREFIXES "")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
    
    if(FLA_OPENMP_MULTITHREADING)
      IF(BUILD_SHARED_LIBS)
        SET(BLAS_LIB_NAME "AOCL-LibBlis-Win-MT-dll")
      ELSE(BUILD_SHARED_LIBS)
        SET(BLAS_LIB_NAME "AOCL-LibBlis-Win-MT")
      ENDIF(BUILD_SHARED_LIBS)
    ELSE(FLA_OPENMP_MULTITHREADING)
      IF(BUILD_SHARED_LIBS)
        SET(BLAS_LIB_NAME "AOCL-LibBlis-Win-dll")
      ELSE(BUILD_SHARED_LIBS)
        SET(BLAS_LIB_NAME "AOCL-LibBlis-Win")
      ENDIF(BUILD_SHARED_LIBS)
    ENDIF(FLA_OPENMP_MULTITHREADING)

    # set aocl-utils library name
    IF(BUILD_SHARED_LIBS)
      SET(UTILS_LIB_NAME "libaoclutils")
    ELSE(BUILD_SHARED_LIBS)
      SET(UTILS_LIB_NAME "libaoclutils_static")
    ENDIF(BUILD_SHARED_LIBS)

    find_library(AOCL_BLAS_LIB
    NAMES ${BLAS_LIB_NAME}
    HINTS ${AOCL_ROOT}/blis ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}
    PATH_SUFFIXES "lib/${ILP_DIR}" "lib_${ILP_DIR}" "lib"
    DOC "AOCL-BLAS library"
    )

    # add aocl-utils library
    find_library(AOCL_UTILS_LIB
    NAMES ${UTILS_LIB_NAME}
    HINTS ${AOCL_ROOT}/amd-utils ${AOCL_ROOT}
    PATH_SUFFIXES "lib/${ILP_DIR}" "lib_${ILP_DIR}" "lib"
    DOC "AOCL-UTILS library"
    )

    #====Headers
    find_path(AOCL_BLAS_INCLUDE_DIR
    NAMES blis.h  cblas.h
    HINTS ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}/blis ${AOCL_ROOT}
    PATH_SUFFIXES "include/${ILP_DIR}" "include_${ILP_DIR}" "include" "include/blis"
    DOC "AOCL-BLAS headers"
    )

    # add aocl-utils headers
    find_path(AOCL_UTILS_INCLUDE_DIR
    NAMES alci/alci_c.h  alci/alci.h  alci/arch.h  alci/enum.h  alci/macros.h
    HINTS ${AOCL_ROOT}/amd-utils ${AOCL_ROOT}
    PATH_SUFFIXES "include/${ILP_DIR}" "include_${ILP_DIR}" "include"
    DOC "AOCL-UTILS headers"
    )

  ELSE(WIN32)   
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    IF(BUILD_SHARED_LIBS)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
    ELSE(BUILD_SHARED_LIBS)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    ENDIF(BUILD_SHARED_LIBS) 

    IF(FLA_OPENMP_MULTITHREADING)
      SET(BLAS_LIB_NAME "blis-mt")
    ELSE(FLA_OPENMP_MULTITHREADING)
      SET(BLAS_LIB_NAME "blis")
    ENDIF(FLA_OPENMP_MULTITHREADING)

    set(UTILS_LIB_NAME "aoclutils")

    find_library(AOCL_BLAS_LIB
    NAMES ${BLAS_LIB_NAME}
    HINTS ${AOCL_ROOT}/blis ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}
    PATH_SUFFIXES "lib/${ILP_DIR}" "lib_${ILP_DIR}" "lib"
    DOC "AOCL-BLAS library"
    )
    
    # find aoclutils library
    find_library(AOCL_UTILS_LIB
    NAMES ${UTILS_LIB_NAME}
    HINTS ${AOCL_ROOT}/amd-utils ${AOCL_ROOT}
    PATH_SUFFIXES "lib/${ILP_DIR}" "lib_${ILP_DIR}" "lib"
    DOC "AOCL-UTILS library"
    )

    #====Headers
    find_path(AOCL_BLAS_INCLUDE_DIR
    NAMES blis.h  blis.hh  cblas.h  cblas.hh
    HINTS ${AOCL_ROOT}/blis ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}
    PATH_SUFFIXES "include/${ILP_DIR}" "include_${ILP_DIR}" "include" "include/blis"
    DOC "AOCL-BLAS headers"
    )

    # add aocl-utils headers
    find_path(AOCL_UTILS_INCLUDE_DIR
    NAMES alci/alci_c.h  alci/alci.h  alci/arch.h  alci/enum.h  alci/macros.h
    HINTS ${AOCL_ROOT}/amd-utils ${AOCL_ROOT}
    PATH_SUFFIXES "include/${ILP_DIR}" "include_${ILP_DIR}" "include"
    DOC "AOCL-UTILS headers"
    )

  ENDIF(WIN32)

  #===========
  # Check if the library and headers are found
  if(AOCL_BLAS_LIB AND AOCL_BLAS_INCLUDE_DIR)
    set(BLAS_LIBRARY ${AOCL_BLAS_LIB})

    # Get the path and name of the library
    get_filename_component(BLAS_LIB_PATH ${AOCL_BLAS_LIB} DIRECTORY)
    get_filename_component(BLAS_LIB_NAME ${AOCL_BLAS_LIB} NAME)
  endif()

  # check if aocl-utils library and headers are found
  if(AOCL_UTILS_LIB AND AOCL_UTILS_INCLUDE_DIR)
    set(AOCL_UTILS_LIBRARY ${AOCL_UTILS_LIB})
  endif()

  # Set the variables to the parent scope
set(BLAS_LIBRARY ${BLAS_LIBRARY} PARENT_SCOPE)
set(BLAS_LIB_PATH ${BLAS_LIB_PATH} PARENT_SCOPE)
set(BLAS_LIB_NAME ${BLAS_LIB_NAME} PARENT_SCOPE)
set(AOCL_UTILS_LIBRARY ${AOCL_UTILS_LIBRARY} PARENT_SCOPE)
set(AOCL_UTILS_INCLUDE_DIR ${AOCL_UTILS_INCLUDE_DIR} PARENT_SCOPE)

endfunction(aocl_libs)

#==================main=================
# clear to avoid endless appending on subsequent calls
set(BLAS_LIBRARY)
unset(BLAS_INCLUDE_DIR)
set(AOCL_UTILS_LIBRARY)
unset(AOCL_UTILS_INCLUDE_DIR)

if(DEFINED ENV{AOCL_ROOT})            
    SET(AOCL_ROOT $ENV{AOCL_ROOT})
    message(STATUS "AOCL_ROOT set via environment variable is ${AOCL_ROOT}")
    if(NOT EXISTS ${AOCL_ROOT})
			message(FATAL_ERROR "\n Invalid path to AOCL_ROOT \n")
		endif()
elseif(AOCL_ROOT)
    SET(AOCL_ROOT ${AOCL_ROOT})
    message(STATUS "AOCL_ROOT set from cmake option is ${AOCL_ROOT}")
    if(NOT EXISTS ${AOCL_ROOT})
      message(FATAL_ERROR "\n Invalid path to AOCL_ROOT \n")
    endif()
else()
      # if aocl-root is not set then check using pkg-config by default
      find_package(PkgConfig)
      if (PKG_CONFIG_FOUND)
        pkg_check_modules(PKG_AOCL_BLAS ${REQ_BLAS_PKGNAME})
        if (PKG_AOCL_BLAS_FOUND AND UNIX)
          set (AOCL_BLAS_INCLUDE_DIR ${PKG_AOCL_BLAS_INCLUDE_DIRS})
          get_filename_component(BLAS_LIB_PATH ${PKG_AOCL_BLAS_LINK_LIBRARIES} DIRECTORY)
          get_filename_component(BLAS_LIB_NAME_WE ${PKG_AOCL_BLAS_LINK_LIBRARIES} NAME_WE)
          if (BUILD_SHARED_LIBS)
            set (BLAS_LIB_NAME ${BLAS_LIB_NAME_WE}.so)
          else()
            set (BLAS_LIB_NAME ${BLAS_LIB_NAME_WE}.a)
          endif()
          set (BLAS_LIBRARY ${BLAS_LIB_PATH}/${BLAS_LIB_NAME})
          message(STATUS "Found AOCL-BLAS library using pkg-config: ${BLAS_LIBRARY}")
        endif()

        pkg_check_modules(AOCL_UTILS aocl-utils)
        if (AOCL_UTILS_FOUND AND UNIX)
          set (AOCL_UTILS_INCLUDE_DIR ${AOCL_UTILS_INCLUDE_DIRS})
          get_filename_component(AOCL_UTILS_LIB_PATH ${AOCL_UTILS_LINK_LIBRARIES} DIRECTORY)
          get_filename_component(AOCL_UTILS_LIB_NAME ${AOCL_UTILS_LINK_LIBRARIES} NAME_WE)

          if (BUILD_SHARED_LIBS)
            set (AOCL_UTILS_LIBRARY ${AOCL_UTILS_LIB_PATH}/${AOCL_UTILS_LIB_NAME}.so)
          else()
            set (AOCL_UTILS_LIBRARY ${AOCL_UTILS_LIB_PATH}/${AOCL_UTILS_LIB_NAME}.a)
          endif()
          message(STATUS "Found AOCL-UTILS library using pkg-config: ${AOCL_UTILS_LIBRARY}")
        endif()
      endif()      

      set(AOCL_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

if (NOT AOCL_BLAS_INCLUDE_DIR OR NOT BLAS_LIBRARY
    OR NOT AOCL_UTILS_INCLUDE_DIR OR NOT AOCL_UTILS_LIBRARY)
  aocl_libs()
endif()

set(BLAS_LIBRARY ${BLAS_LIBRARY})
set(BLAS_INCLUDE_DIR ${AOCL_BLAS_INCLUDE_DIR})
set(AOCL_UTILS_LIBRARY ${AOCL_UTILS_LIBRARY})
set(AOCL_UTILS_INCLUDE_DIR ${AOCL_UTILS_INCLUDE_DIR})

if(ENABLE_AOCL_BLAS)
  message(STATUS "AOCL-BLAS LIBRARY= ${BLAS_LIBRARY}")
  message(STATUS "AOCL-BLAS INCLUDE DIRS= ${BLAS_INCLUDE_DIR}")
endif()

# Hide these variables from the CMake GUI
mark_as_advanced(BLAS_LIBRARY BLAS_INCLUDE_DIR BLIS_LIB_NAME BLIS_LIB_PATH AOCL_UTILS_LIBRARY AOCL_UTILS_INCLUDE_DIR)