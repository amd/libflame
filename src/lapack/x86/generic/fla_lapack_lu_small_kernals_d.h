/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file ffla_lapack_lu_small_kernals_d.h
 *  @brief Common front-end functions
 *         for double precision
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"
#include "fla_lapack_lu_small_kernals_common.h"

#if FLA_ENABLE_AMD_OPT
/*
 * LU 2x2  with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_D_2x2(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_2x2(doublereal, i, n, buff_A, ldim_A, buff_p, info)
/*
 * LU 3x3 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_D_3x3(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_3x3(doublereal, i, n, buff_A, ldim_A, buff_p, info)
/*
 * LU 4x4 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_D_4x4(i, n, buff_A, ldim_A, buff_p, info) \
    FLA_LU_PIV_SMALL_GEN_4x4(doublereal, i, n, buff_A, ldim_A, buff_p, info)

#define LAPACK_GETRI_SMALL_D_3x3(n, a, a_dim1, ipiv, work)                                  \
    integer jp, a_dim2, a_dim3;                                                             \
    doublereal t_val, *apiv, t_val11, t_val21, t_val31, t_val12, t_val22, t_val32, t_val13, \
        t_val23, t_val33;                                                                   \
    /* trmv loop 1 */                                                                       \
    a_dim2 = 2 * a_dim1;                                                                    \
    a_dim3 = 3 * a_dim1;                                                                    \
    t_val11 = 1. / a[1 + a_dim1];                                                           \
    t_val21 = a[2 + a_dim1];                                                                \
    /* loop 2*/                                                                             \
    t_val22 = 1. / a[2 + a_dim2];                                                           \
    t_val12 = a[1 + a_dim2] * t_val11;                                                      \
    /*  loop 3 */                                                                           \
    t_val12 = -t_val22 * t_val12;                                                           \
    t_val33 = 1. / a[3 + a_dim3];                                                           \
    t_val23 = a[2 + a_dim3];                                                                \
    t_val13 = a[1 + a_dim3];                                                                \
    t_val13 = (t_val13 * t_val11) + (t_val23 * t_val12);                                    \
    t_val23 = t_val23 * t_val22;                                                            \
                                                                                            \
    t_val13 = -t_val33 * t_val13;                                                           \
    t_val23 = -t_val33 * t_val23;                                                           \
    /* gemv inline */                                                                       \
    t_val32 = a[3 + a_dim2];                                                                \
    t_val12 = t_val12 - t_val32 * t_val13;                                                  \
    t_val22 = t_val22 - t_val32 * t_val23;                                                  \
    t_val32 = -t_val32 * t_val33;                                                           \
    t_val31 = a[3 + a_dim1];                                                                \
    work[1] = t_val31;                                                                      \
                                                                                            \
    t_val11 = t_val11 - t_val21 * t_val12;                                                  \
    t_val31 = -t_val21 * t_val32;                                                           \
    t_val21 = -t_val21 * t_val22;                                                           \
    a[1 + a_dim2] = t_val12;                                                                \
    a[2 + a_dim2] = t_val22;                                                                \
    a[3 + a_dim2] = t_val32;                                                                \
    a[1 + a_dim1] = t_val11 - work[1] * t_val13;                                            \
    a[2 + a_dim1] = t_val21 - work[1] * t_val23;                                            \
    a[3 + a_dim1] = t_val31 - work[1] * t_val33;                                            \
    a[1 + a_dim3] = t_val13;                                                                \
    a[2 + a_dim3] = t_val23;                                                                \
    a[3 + a_dim3] = t_val33;                                                                \
    jp = ipiv[2];                                                                           \
    if(jp != 2)                                                                             \
    {                                                                                       \
        /*  swap */                                                                         \
        apiv = &a[jp * a_dim1 + 1];                                                         \
        t_val = apiv[0];                                                                    \
        apiv[0] = a[1 + a_dim2];                                                            \
        a[1 + a_dim2] = t_val;                                                              \
        t_val = apiv[1];                                                                    \
        apiv[1] = a[2 + a_dim2];                                                            \
        a[2 + a_dim2] = t_val;                                                              \
        t_val = apiv[2];                                                                    \
        apiv[2] = a[3 + a_dim2];                                                            \
        a[3 + a_dim2] = t_val;                                                              \
    }                                                                                       \
    jp = ipiv[1];                                                                           \
    if(jp != 1)                                                                             \
    {                                                                                       \
        /*  swap */                                                                         \
        apiv = &a[jp * a_dim1 + 1];                                                         \
        t_val = apiv[0];                                                                    \
        apiv[0] = a[1 + a_dim1];                                                            \
        a[1 + a_dim1] = t_val;                                                              \
        t_val = apiv[1];                                                                    \
        apiv[1] = a[2 + a_dim1];                                                            \
        a[2 + a_dim1] = t_val;                                                              \
        t_val = apiv[2];                                                                    \
        apiv[2] = a[3 + a_dim1];                                                            \
        a[3 + a_dim1] = t_val;                                                              \
    }
#endif