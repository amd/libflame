##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpy.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpy.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpy_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Axpy.c
	)