/******************************************************************************
 * * Copyright (C) 2023-24, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/
/*! @file fla_lapack_x86_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_avx512_kernels.h"
#include "fla_lapack_x86_common.h"

#if FLA_ENABLE_AMD_OPT
/* 3x3 Householder Rotation */
int fla_dhrot3(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dhrot3_avx2(n, a, lda, v, tau);
    }
    return 0;
}
/* 2x2 Plane Rotation */
int fla_drot(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy,
             doublereal *c__, doublereal *s)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_drot_avx2(n, dx, incx, dy, incy, c__, s);
    }
    return 0;
}
void fla_zrot(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy,
              doublereal *c__, doublecomplex *s)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_zrot_avx2(n, cx, incx, cy, incy, c__, s);
    }
    return;
}
/* complex vector scaling when increment is 1 and specific threshold */
int fla_zscal(integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx)
{
    /* Initialize global context data */
    aocl_fla_init();
    /* Take AVX path only for increment equal to 1 and particular threshold size*/
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2) && *incx == 1 && *n <= FLA_ZSCAL_INLINE_SMALL)
    {
        fla_zscal_ix1_avx2(n, alpha, x);
    }
    else
    {
        zscal_(n, (dcomplex *)alpha, (dcomplex *)x, incx);
    }
    return 0;
}
/* scales a vector by a constant when threshold <= 128 */
void fla_dscal(integer *n, doublereal *da, doublereal *dx, integer *incx)
{
    if(*incx == 1 && *da != 0 && *n >= 1 && *n <= FLA_DSCAL_INLINE_SMALL)
    {
        /* Initialize global context data */
        aocl_fla_init();
        if(FLA_IS_ARCH_ID(FLA_ARCH_AVX512))
        {
            fla_dscal_ix1_avx512(n, da, dx, incx);
        }
        else if(FLA_IS_ARCH_ID(FLA_ARCH_AVX2))
        {
            fla_dscal_ix1_avx2(n, da, dx, incx);
        }
        else
        {
            dscal_(n, da, dx, incx);
        }
    }
    else
    {
        dscal_(n, da, dx, incx);
    }
    return;
}
/* Double QR (DGEQRF) for small sizes */
int fla_dgeqrf_small(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau,
                     doublereal *work)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgeqrf_small_avx2(m, n, a, lda, tau, work);
    }
    return 0;
}
/* real vector scaling when increment is 1 */
void fla_sscal(integer *n, real *alpha, real *x, integer *incx)
{
    /* Take AVX path only for increment equal to 1 */
    if(*incx == 1 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_sscal_ix1_avx2(n, alpha, x);
    }
    else
    {
        sscal_(n, (real *)alpha, (real *)x, incx);
    }
    return;
}
/* Rank 1 Operation */
void fla_sger(integer *m, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy,
              real *a, integer *lda)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_sger_avx2(m, n, alpha, x, incx, y, incy, a, lda);
    }
    return;
}

/* LU factorization.
 * To be used only when vectorized code via avx2/avx512 is enabled
 * */
int fla_dgetrf_small_simd(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv,
                          integer *info)
{
    if(FLA_IS_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_dgetrf_small_avx512(m, n, a, lda, ipiv, info);
    }
    else if(FLA_IS_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgetrf_small_avx2(m, n, a, lda, ipiv, info);
    }
    return 0;
}

/* Double Complex LU for small sizes,
 * Optimized for AVX2 and AVX512 ISAs
 */
int fla_zgetrf_small_simd(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                          integer *info)
{
    if(FLA_IS_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_zgetrf_small_avx512(m, n, a, lda, ipiv, info);
    }
    else if (FLA_IS_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_zgetrf_small_avx2(m, n, a, lda, ipiv, info);
    }
    else
    {
        lapack_zgetf2(m, n, a, lda, ipiv, info);
    }
    return 0;
}

/* SVD for small tall-matrices in DGESVD
 */
void fla_dgesvd_nn_small10(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                           doublereal *work, integer *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_nn_small10_avx2(m, n, a, lda, s, work, info);
    }
    return;
}

/* SVD for small fat-matrices with LQ factorization
 * already computed
 */
void fla_dgesvd_small6(integer *m, integer *n, doublereal *a, integer *lda, doublereal *qr,
                       integer *ldqr, doublereal *s, doublereal *u, integer *ldu, doublereal *vt,
                       integer *ldvt, doublereal *work, integer *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_small6_avx2(m, n, a, lda, qr, ldqr, s, u, ldu, vt, ldvt, work, info);
    }
    return;
}

/* SVD for small fat-matrices for path 1T in DGESVD
 */
void fla_dgesvd_nn_small1T(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                           doublereal *work, integer *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_nn_small1T_avx2(m, n, a, lda, s, work, info);
    }
    return;
}

/* SVD for small fat-matrices with LQ factorization
 * already computed
 */
void fla_dgesvd_small6T(integer *m, integer *n, doublereal *a, integer *lda, doublereal *ql,
                        integer *ldql, doublereal *s, doublereal *u, integer *ldu, doublereal *vt,
                        integer *ldvt, doublereal *work, integer *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_small6T_avx2(m, n, a, lda, ql, ldql, s, u, ldu, vt, ldvt, work, info);
    }
    return;
}

/* Small DGETRS path (NOTRANS) should only be used for size between 3 to 8 and NRHS <= N */
int fla_dgetrs_small_notrans(char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda,
                             integer *ipiv, doublereal *b, integer *ldb, integer *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgetrs_small_trsm_ll_avx2(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
    }
    return 0;
}
#endif
