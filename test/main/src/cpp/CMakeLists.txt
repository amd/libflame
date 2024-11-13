###############################################################################
# Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################

target_include_directories("${PROJECT_NAME}" PRIVATE ${CMAKE_SOURCE_DIR}/src/src_cpp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-comment -Wno-unused-function")

if(WIN32)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

target_sources(${PROJECT_NAME}
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gbtrf.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gbtrs.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_geev.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_geevx.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gehrd.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gelqf.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gels.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gelsd.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gelss.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_geqp3.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_geqrf.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gesdd.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gesv.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gesvd.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gesvdx.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_getrf.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_getri.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_getrs.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_ggev.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_ggevx.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gghrd.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_gtsv.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_hetrf_rook.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_hgeqz.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_hseqr.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_lartg.cc
${CMAKE_CURRENT_SOURCE_DIR}/invoke_common.hh
)
