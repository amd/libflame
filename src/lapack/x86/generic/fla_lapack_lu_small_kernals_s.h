/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file ffla_lapack_lu_small_kernals_s.h
 *  @brief Common front-end functions
 *         for single precision
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"
#include "fla_lapack_lu_small_kernals_common.h"

#if FLA_ENABLE_AMD_OPT
/*
 * LU 2x2  with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_S_2x2(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_2x2(real, i, n, buff_A, ldim_A, buff_p, info)
/*
 * LU 3x3 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_S_3x3(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_3x3(real, i, n, buff_A, ldim_A, buff_p, info)
/*
 * LU 4x4 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_S_4x4(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_4x4(real, i, n, buff_A, ldim_A, buff_p, info)

#endif