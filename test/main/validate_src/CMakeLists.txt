###############################################################################
# Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################

target_sources(${PROJECT_NAME}
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/test_common.c
${CMAKE_CURRENT_SOURCE_DIR}/test_overflow_underflow.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gbtrf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_geqrf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gerq2.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gerqf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gelqf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gesdd.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gesvd.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_getrf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_getri.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_getrs.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_orgqr.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_potrf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_potrs.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_syevd.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_geevx.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_geev.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_geqp3.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gesv.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_ggev.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_ggevx.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_stedc.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_hseqr.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_spffrt2.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_spffrtx.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gehrd.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_hgeqz.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gghrd.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_rot.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_lartg.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gels.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_larfg.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gtsv.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_syev.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gesvdx.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gelsd.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_sytrf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_hetrf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_larf.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_sygvd.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_lange.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_gecon.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_getrfnpi.c
${CMAKE_CURRENT_SOURCE_DIR}/validate_hetri_rook.c
)

