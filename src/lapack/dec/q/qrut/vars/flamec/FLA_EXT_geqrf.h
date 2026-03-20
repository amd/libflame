/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    May 09, 2021
*/

FLA_Error FLA_EXT_sgeqrf( fla_dim_t m_A, fla_dim_t n_A,
                          float* buff_A, fla_dim_t cs_A,
                          float* buff_t,
                          float* buff_w,
                          fla_dim_t* lwork,
                          fla_dim_t* info );
FLA_Error FLA_EXT_dgeqrf( fla_dim_t m_A, fla_dim_t n_A,
                          double* buff_A, fla_dim_t cs_A,
                          double* buff_t,
                          double* buff_w,
                          fla_dim_t* lwork,
                          fla_dim_t* info );
