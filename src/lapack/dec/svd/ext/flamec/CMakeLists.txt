##Copyright (C) 2020-2026, Advanced Micro Devices, Inc. All rights reserved.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Svd_ext.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Svd_ext_u_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/dgesvd2x2.c
	)
