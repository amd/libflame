/******************************************************************************
 * * Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_lu_small_kernels_d.h
 *  @brief Common front-end functions
 *         for double precision
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"
#include "fla_lapack_lu_small_kernels_common.h"

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

#define LAPACK_GETRI_SMALL_D_2x2(n, a, a_dim1, ipiv, work)       \
    aocl_int64_t jp;                                                  \
    doublereal t_val, *apiv, t_val11, t_val21, t_val12, t_val22; \
    /* Phase 1: Triangular inversion U^(-1) */                   \
    t_val11 = 1. / a[1 + a_dim1];                                \
    t_val22 = 1. / a[2 + 2 * a_dim1];                            \
    t_val12 = -a[1 + 2 * a_dim1] * t_val11 * t_val22;            \
    /* Phase 2: Lower solve inv(A)*L = inv(U) - column 1 */      \
    t_val21 = a[2 + a_dim1];                                     \
    work[1] = t_val21;                                           \
    a[1 + a_dim1] = t_val11 - t_val21 * t_val12;                 \
    a[2 + a_dim1] = -t_val21 * t_val22;                          \
    a[1 + 2 * a_dim1] = t_val12;                                 \
    a[2 + 2 * a_dim1] = t_val22;                                 \
    /* Phase 3: Apply column interchanges */                     \
    jp = ipiv[1];                                                \
    if(jp != 1)                                                  \
    {                                                            \
        apiv = &a[jp * a_dim1 + 1];                              \
        t_val = apiv[0];                                         \
        apiv[0] = a[1 + a_dim1];                                 \
        a[1 + a_dim1] = t_val;                                   \
        t_val = apiv[1];                                         \
        apiv[1] = a[2 + a_dim1];                                 \
        a[2 + a_dim1] = t_val;                                   \
    }

#define LAPACK_GETRI_SMALL_D_3x3(n, a, a_dim1, ipiv, work)                                  \
    aocl_int64_t jp, a_dim2, a_dim3;                                                             \
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

/*
 * Specialized triangular inversion for 4x4 upper triangular matrix
 * Computes U^(-1) in-place where U is upper triangular
 * Uses fully unrolled algorithm for optimal performance
 */
#define LAPACK_GETRI_TRIANGULAR_INVERSION_D_4X4(a, a_dim1)                \
{                                                                         \
    /* Invert diagonal elements */                                        \
    doublereal inv11 = 1.0 / a[1 + 1 * a_dim1];                           \
    a[1 + 1 * a_dim1] = inv11;                                            \
    doublereal inv22 = 1.0 / a[2 + 2 * a_dim1];                           \
    doublereal u12 = a[1 + 2 * a_dim1];                                   \
    a[2 + 2 * a_dim1] = inv22;                                            \
    /* Compute column 2: U^(-1)[1,2] = -inv22 * inv11 * u12 */            \
    a[1 + 2 * a_dim1] = -inv22 * (inv11 * u12);                           \
    doublereal inv33 = 1.0 / a[3 + 3 * a_dim1];                           \
    doublereal u13 = a[1 + 3 * a_dim1];                                   \
    doublereal u23 = a[2 + 3 * a_dim1];                                   \
    a[3 + 3 * a_dim1] = inv33;                                            \
    /* Compute column 3 */                                                \
    a[2 + 3 * a_dim1] = -inv33 * (inv22 * u23);                           \
    a[1 + 3 * a_dim1] = -inv33 * (inv11 * u13 + a[1 + 2 * a_dim1] * u23); \
    doublereal inv44 = 1.0 / a[4 + 4 * a_dim1];                           \
    doublereal u14 = a[1 + 4 * a_dim1];                                   \
    doublereal u24 = a[2 + 4 * a_dim1];                                   \
    doublereal u34 = a[3 + 4 * a_dim1];                                   \
    a[4 + 4 * a_dim1] = inv44;                                            \
    /* Compute column 4 */                                                \
    a[3 + 4 * a_dim1] = -inv44 * (inv33 * u34);                           \
    a[2 + 4 * a_dim1] = -inv44 * (inv22 * u24 + a[2 + 3 * a_dim1] * u34); \
    a[1 + 4 * a_dim1] = -inv44 * (inv11 * u14 + a[1 + 2 * a_dim1] * u24 + a[1 + 3 * a_dim1] * u34); \
}

/*
 * Specialized triangular inversion for 5x5 upper triangular matrix
 * Computes U^(-1) in-place where U is upper triangular
 * Uses fully unrolled algorithm for optimal performance on 5x5 matrices
 */
#define LAPACK_GETRI_TRIANGULAR_INVERSION_D_5X5(a, a_dim1)                                 \
{                                                                                          \
    /* Invert diagonal elements and compute column 1 */                                    \
    doublereal inv11 = 1.0 / a[1 + 1 * a_dim1];                                            \
    a[1 + 1 * a_dim1] = inv11;                                                             \
    doublereal inv22 = 1.0 / a[2 + 2 * a_dim1];                                            \
    doublereal u12 = a[1 + 2 * a_dim1];                                                    \
    a[2 + 2 * a_dim1] = inv22;                                                             \
    /* Compute column 2: U^(-1)[1,2] = -inv22 * inv11 * u12 */                             \
    a[1 + 2 * a_dim1] = -inv22 * (inv11 * u12);                                            \
    doublereal inv33 = 1.0 / a[3 + 3 * a_dim1];                                            \
    doublereal u13 = a[1 + 3 * a_dim1];                                                    \
    doublereal u23 = a[2 + 3 * a_dim1];                                                    \
    a[3 + 3 * a_dim1] = inv33;                                                             \
    /* Compute column 3 */                                                                 \
    a[2 + 3 * a_dim1] = -inv33 * (inv22 * u23);                                            \
    a[1 + 3 * a_dim1] = -inv33 * (inv11 * u13 + a[1 + 2 * a_dim1] * u23);                  \
    doublereal inv44 = 1.0 / a[4 + 4 * a_dim1];                                            \
    doublereal u14 = a[1 + 4 * a_dim1];                                                    \
    doublereal u24 = a[2 + 4 * a_dim1];                                                    \
    doublereal u34 = a[3 + 4 * a_dim1];                                                    \
    a[4 + 4 * a_dim1] = inv44;                                                             \
    /* Compute column 4 */                                                                 \
    a[3 + 4 * a_dim1] = -inv44 * (inv33 * u34);                                            \
    a[2 + 4 * a_dim1] = -inv44 * (inv22 * u24 + a[2 + 3 * a_dim1] * u34);                  \
    a[1 + 4 * a_dim1]                                                                      \
        = -inv44 * (inv11 * u14 + a[1 + 2 * a_dim1] * u24 + a[1 + 3 * a_dim1] * u34);      \
    doublereal inv55 = 1.0 / a[5 + 5 * a_dim1];                                            \
    doublereal u15 = a[1 + 5 * a_dim1];                                                    \
    doublereal u25 = a[2 + 5 * a_dim1];                                                    \
    doublereal u35 = a[3 + 5 * a_dim1];                                                    \
    doublereal u45 = a[4 + 5 * a_dim1];                                                    \
    a[5 + 5 * a_dim1] = inv55;                                                             \
    /* Compute column 5 */                                                                 \
    a[4 + 5 * a_dim1] = -inv55 * (inv44 * u45);                                            \
    a[3 + 5 * a_dim1] = -inv55 * (inv33 * u35 + a[3 + 4 * a_dim1] * u45);                  \
    a[2 + 5 * a_dim1]                                                                      \
        = -inv55 * (inv22 * u25 + a[2 + 3 * a_dim1] * u35 + a[2 + 4 * a_dim1] * u45);      \
    a[1 + 5 * a_dim1] = -inv55                                                             \
                        * (inv11 * u15 + a[1 + 2 * a_dim1] * u25 + a[1 + 3 * a_dim1] * u35 \
                           + a[1 + 4 * a_dim1] * u45);                                     \
}

/*
 * Specialized lower solve for 4x4 matrices in GETRI operation
 * Solves inv(U) * L = inv(A) for inv(A) where L is unit lower triangular
 * Uses optimized unrolled loops for 4x4 case
 */
#define LAPACK_GETRI_LOWER_SOLVE_D_4X4(a, a_dim1, work)                                   \
{                                                                                         \
    aocl_int64_t j, i;                                                                         \
    /* Process columns from right to left (j = 4 down to 1) */                            \
    for(j = 4; j >= 1; --j)                                                               \
    {                                                                                     \
        /* Copy lower triangular part of column j to work and zero it */                  \
        for(i = j + 1; i <= 4; ++i)                                                       \
        {                                                                                 \
            work[i] = a[i + j * a_dim1];                                                  \
            a[i + j * a_dim1] = 0.0;                                                      \
        }                                                                                 \
        /* Apply GEMV-like operation: A[:, j] -= L[j+1:4, j] * A[:, j+1:4] */             \
        if(j < 4)                                                                         \
        {                                                                                 \
            if(j == 3)                                                                    \
            {                                                                             \
                doublereal w4 = work[4];                                                  \
                for(i = 1; i <= 4; ++i)                                                   \
                    a[i + 3 * a_dim1] -= w4 * a[i + 4 * a_dim1];                          \
            }                                                                             \
            else if(j == 2)                                                               \
            {                                                                             \
                doublereal w3 = work[3], w4 = work[4];                                    \
                for(i = 1; i <= 4; ++i)                                                   \
                    a[i + 2 * a_dim1] -= w3 * a[i + 3 * a_dim1] + w4 * a[i + 4 * a_dim1]; \
            }                                                                             \
            else /* j == 1 */                                                             \
            {                                                                             \
                doublereal w2 = work[2], w3 = work[3], w4 = work[4];                      \
                for(i = 1; i <= 4; ++i)                                                   \
                    a[i + 1 * a_dim1] -= w2 * a[i + 2 * a_dim1] + w3 * a[i + 3 * a_dim1]  \
                                         + w4 * a[i + 4 * a_dim1];                        \
            }                                                                             \
        }                                                                                 \
    }                                                                                     \
}

/*
 * Specialized lower solve for 5x5 matrices in GETRI operation
 * Solves inv(U) * L = inv(A) for inv(A) where L is unit lower triangular
 * Uses pointer optimization and unrolled loops for 5x5 case
 */
#define LAPACK_GETRI_LOWER_SOLVE_D_5X5(a, a_dim1, work)                                         \
{                                                                                               \
    aocl_int64_t j, i;                                                                               \
    /* Process columns from right to left (j = 5 down to 1) */                                  \
    for(j = 5; j >= 1; --j)                                                                     \
    {                                                                                           \
        /* Copy lower triangular part of column j to work and zero it */                        \
        for(i = j + 1; i <= 5; ++i)                                                             \
        {                                                                                       \
            work[i] = a[i + j * a_dim1];                                                        \
            a[i + j * a_dim1] = 0.0;                                                            \
        }                                                                                       \
        /* Apply GEMV-like operation: A[:, j] -= L[j+1:5, j] * A[:, j+1:5] */                   \
        if(j < 5)                                                                               \
        {                                                                                       \
            doublereal *col_j_ptr = &a[1 + j * a_dim1];                                         \
            if(j == 4)                                                                          \
            {                                                                                   \
                doublereal w5 = work[5];                                                        \
                doublereal *col5 = &a[1 + 5 * a_dim1];                                          \
                for(i = 1; i <= 5; ++i)                                                         \
                    col_j_ptr[i - 1] -= w5 * col5[i - 1];                                       \
            }                                                                                   \
            else if(j == 3)                                                                     \
            {                                                                                   \
                doublereal w4 = work[4], w5 = work[5];                                          \
                doublereal *col4 = &a[1 + 4 * a_dim1], *col5 = &a[1 + 5 * a_dim1];              \
                for(i = 1; i <= 5; ++i)                                                         \
                    col_j_ptr[i - 1] -= w4 * col4[i - 1] + w5 * col5[i - 1];                    \
            }                                                                                   \
            else if(j == 2)                                                                     \
            {                                                                                   \
                doublereal w3 = work[3], w4 = work[4], w5 = work[5];                            \
                doublereal *col3 = &a[1 + 3 * a_dim1], *col4 = &a[1 + 4 * a_dim1],              \
                           *col5 = &a[1 + 5 * a_dim1];                                          \
                for(i = 1; i <= 5; ++i)                                                         \
                    col_j_ptr[i - 1] -= w3 * col3[i - 1] + w4 * col4[i - 1] + w5 * col5[i - 1]; \
            }                                                                                   \
            else /* j == 1 */                                                                   \
            {                                                                                   \
                doublereal w2 = work[2], w3 = work[3], w4 = work[4], w5 = work[5];              \
                doublereal *col2 = &a[1 + 2 * a_dim1], *col3 = &a[1 + 3 * a_dim1],              \
                           *col4 = &a[1 + 4 * a_dim1], *col5 = &a[1 + 5 * a_dim1];              \
                for(i = 1; i <= 5; ++i)                                                         \
                    col_j_ptr[i - 1] -= w2 * col2[i - 1] + w3 * col3[i - 1] + w4 * col4[i - 1]  \
                                        + w5 * col5[i - 1];                                     \
            }                                                                                   \
        }                                                                                       \
    }                                                                                           \
}

/*
 * Specialized lower solve for 6x6 matrices in GETRI operation
 * Solves inv(U) * L = inv(A) for inv(A) where L is unit lower triangular
 * Uses optimized unrolled loops for 6x6 case
 */
#define LAPACK_GETRI_LOWER_SOLVE_D_6X6(a, a_dim1, work)                                          \
{                                                                                                \
    aocl_int64_t j, i;                                                                                \
    /* Process columns from right to left (j = 6 down to 1) */                                   \
    for(j = 6; j >= 1; --j)                                                                      \
    {                                                                                            \
        /* Copy lower triangular part of column j to work and zero it */                         \
        for(i = j + 1; i <= 6; ++i)                                                              \
        {                                                                                        \
            work[i] = a[i + j * a_dim1];                                                         \
            a[i + j * a_dim1] = 0.0;                                                             \
        }                                                                                        \
        /* Apply GEMV-like operation: A[:, j] -= L[j+1:6, j] * A[:, j+1:6] */                    \
        if(j < 6)                                                                                \
        {                                                                                        \
            if(j == 5)                                                                           \
            {                                                                                    \
                doublereal w6 = work[6];                                                         \
                for(i = 1; i <= 6; ++i)                                                          \
                    a[i + 5 * a_dim1] -= w6 * a[i + 6 * a_dim1];                                 \
            }                                                                                    \
            else if(j == 4)                                                                      \
            {                                                                                    \
                doublereal w5 = work[5], w6 = work[6];                                           \
                for(i = 1; i <= 6; ++i)                                                          \
                    a[i + 4 * a_dim1] -= w5 * a[i + 5 * a_dim1] + w6 * a[i + 6 * a_dim1];        \
            }                                                                                    \
            else if(j == 3)                                                                      \
            {                                                                                    \
                doublereal w4 = work[4], w5 = work[5], w6 = work[6];                             \
                for(i = 1; i <= 6; ++i)                                                          \
                    a[i + 3 * a_dim1] -= w4 * a[i + 4 * a_dim1] + w5 * a[i + 5 * a_dim1]         \
                                         + w6 * a[i + 6 * a_dim1];                               \
            }                                                                                    \
            else if(j == 2)                                                                      \
            {                                                                                    \
                doublereal w3 = work[3], w4 = work[4], w5 = work[5], w6 = work[6];               \
                for(i = 1; i <= 6; ++i)                                                          \
                    a[i + 2 * a_dim1] -= w3 * a[i + 3 * a_dim1] + w4 * a[i + 4 * a_dim1]         \
                                         + w5 * a[i + 5 * a_dim1] + w6 * a[i + 6 * a_dim1];      \
            }                                                                                    \
            else /* j == 1 */                                                                    \
            {                                                                                    \
                doublereal w2 = work[2], w3 = work[3], w4 = work[4], w5 = work[5], w6 = work[6]; \
                for(i = 1; i <= 6; ++i)                                                          \
                    a[i + 1 * a_dim1] -= w2 * a[i + 2 * a_dim1] + w3 * a[i + 3 * a_dim1]         \
                                         + w4 * a[i + 4 * a_dim1] + w5 * a[i + 5 * a_dim1]       \
                                         + w6 * a[i + 6 * a_dim1];                               \
            }                                                                                    \
        }                                                                                        \
    }                                                                                            \
}
#endif
