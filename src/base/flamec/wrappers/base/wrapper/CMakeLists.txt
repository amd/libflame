##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/libflame_api_wrapper.c
${CMAKE_CURRENT_SOURCE_DIR}/aocl_blas_wrappers.c
  )
