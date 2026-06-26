##Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.## 
target_sources("${PROJECT_NAME}"
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/fla_mem.h
${CMAKE_CURRENT_SOURCE_DIR}/fla_mem.c
)

target_include_directories("${PROJECT_NAME}" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(accum)
add_subdirectory(app)
