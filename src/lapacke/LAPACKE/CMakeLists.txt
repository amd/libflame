message(STATUS "LAPACKE enabled")
enable_language(C)
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(utils)

set(LAPACK_INSTALL_EXPORT_NAME ${PROJECT_NAME})

include_directories(include ${LAPACK_BINARY_DIR}/include)

option(LAPACKE_BUILD_SINGLE "Build LAPACKE single precision real" ON)
option(LAPACKE_BUILD_DOUBLE "Build LAPACKE double precision real" ON)
option(LAPACKE_BUILD_COMPLEX "Build LAPACKE single precision complex" ON)
option(LAPACKE_BUILD_COMPLEX16 "Build LAPACKE double precision complex" ON)


macro(append_subdir_files variable dirname)
  get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
  foreach(depfile ${holder})
    list(APPEND ${variable} "${dirname}/${depfile}")
  endforeach()
endmacro()

message(STATUS "Build LAPACKE single precision real: ${LAPACKE_BUILD_SINGLE}")
message(STATUS "Build LAPACKE double precision real: ${LAPACKE_BUILD_DOUBLE}")
message(STATUS "Build LAPACKE single precision complex: ${LAPACKE_BUILD_COMPLEX}")
message(STATUS "Build LAPACKE double precision complex: ${LAPACKE_BUILD_COMPLEX16}")

append_subdir_files(LAPACKE_INCLUDE "include")
append_subdir_files(SOURCES "src")
if (LAPACKE_BUILD_SINGLE)
  append_subdir_files(SOURCES_SINGLE "src")
  list(APPEND SOURCES ${SOURCES_SINGLE})
endif()
if (LAPACKE_BUILD_DOUBLE)
  append_subdir_files(SOURCES_DOUBLE "src")
  list(APPEND SOURCES ${SOURCES_DOUBLE})
endif()
if (LAPACKE_BUILD_COMPLEX)
  append_subdir_files(SOURCES_COMPLEX "src")
  list(APPEND SOURCES ${SOURCES_COMPLEX})
endif()
if (LAPACKE_BUILD_COMPLEX16)
  append_subdir_files(SOURCES_COMPLEX16 "src")
  list(APPEND SOURCES ${SOURCES_COMPLEX16})
endif()
append_subdir_files(DEPRECATED "src")
append_subdir_files(EXTENDED "src")
append_subdir_files(MATGEN "src")
append_subdir_files(UTILS "utils")

list(APPEND SOURCES ${DEPRECATED})

if(USE_XBLAS)
  list(APPEND SOURCES ${EXTENDED})
endif()
if(LAPACKE_WITH_TMG)
  list(APPEND SOURCES ${MATGEN})
endif()
list(APPEND SOURCES ${UTILS})

add_library(${LAPACKELIB}_obj OBJECT ${SOURCES})

if(BUILD_INDEX64_EXT_API)
  # 64bit Integer Extended Interface
  set(SOURCES_64_C)
  list(APPEND SOURCES_64_C ${SOURCES})
  list(REMOVE_ITEM SOURCES_64_C src/lapacke_nancheck.c)
  list(REMOVE_ITEM SOURCES_64_C utils/lapacke_make_complex_float.c)
  list(REMOVE_ITEM SOURCES_64_C utils/lapacke_make_complex_double.c)

  add_library(${LAPACKELIB}_64_obj OBJECT ${SOURCES_64_C})
  set_target_properties(${LAPACKELIB}_64_obj PROPERTIES
          POSITION_INDEPENDENT_CODE ON)
  target_compile_options(${LAPACKELIB}_64_obj PRIVATE
          -DLAPACKE_API64
          -DWeirdNEC)
endif()


target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${LAPACK_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
if(WIN32 AND NOT UNIX)
  target_compile_definitions(${PROJECT_NAME} PUBLIC HAVE_LAPACK_CONFIG_H LAPACK_COMPLEX_STRUCTURE)
  message(STATUS "Windows BUILD")
endif()

if(LAPACKE_WITH_TMG)
  target_link_libraries(${PROJECT_NAME} PRIVATE ${TMGLIB})
endif()
target_link_libraries(${PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES})

install( FILES ${LAPACKE_INCLUDE} ${LAPACK_BINARY_DIR}/include/lapacke_mangling.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  COMPONENT Development
  )

if(BUILD_TESTING)
  add_subdirectory(example)
endif()


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lapacke.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKELIB}.pc @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKELIB}.pc
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig
  COMPONENT Development
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${LAPACKELIB}-config.cmake
  ${LAPACK_BINARY_DIR}/${LAPACKELIB}-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LAPACKELIB}-${LAPACK_VERSION}
  COMPONENT Development
  )
