##Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved. ## 
target_sources("${PROJECT_NAME}"
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_x86_common.h
${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_x86_common.c
)

add_subdirectory(nofma)
