##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/bl1_gemm.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_hemm.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_her2k.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_herk.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_symm.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_syr2k.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_syrk.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_trmm.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_trmmsx.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_trsm.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_trsmsx.c
	)