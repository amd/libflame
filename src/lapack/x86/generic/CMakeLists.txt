##Copyright(C) 2024, Advanced Micro Devices, Inc.All rights reserved.##
	target_sources("${PROJECT_NAME}" PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/ fla_lapack_lu_small_kernals_common.h
	${CMAKE_CURRENT_SOURCE_DIR}/ fla_lapack_lu_small_kernals_d.h
	${CMAKE_CURRENT_SOURCE_DIR}/ fla_lapack_lu_small_kernals_s.h
	${CMAKE_CURRENT_SOURCE_DIR}/ fla_lapack_lu_getri_small_d.c
	)
