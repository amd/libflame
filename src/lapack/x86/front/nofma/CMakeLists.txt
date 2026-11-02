##Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved. ##

add_library(FLA_LAPACK_NOFMA OBJECT 
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_x86_common_nofma.c
)

if(MSVC)
    # /fp:precise drops contractions only from VS2022 on, /fp:strict drops them in every version.
    target_compile_options(FLA_LAPACK_NOFMA PRIVATE /fp:strict)
else()
    target_compile_options(FLA_LAPACK_NOFMA PRIVATE -ffp-contract=off)
endif()
