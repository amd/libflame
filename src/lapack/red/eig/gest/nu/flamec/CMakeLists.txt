##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_blk_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_blk_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_opt_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_opt_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_opt_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_unb_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_unb_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu_unb_var5.c
	)