##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_blk_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_opt_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_opt_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_unb_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Chol_u_unb_var3.c

	)