##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR2_UT_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR2_UT_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR2_UT_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR2_UT_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR2_UT_vars.h

	)