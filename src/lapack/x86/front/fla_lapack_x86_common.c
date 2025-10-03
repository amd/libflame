/******************************************************************************
 * * Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/
/*! @file fla_lapack_x86_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "fla_lapack_x86_common.h"
#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Transpose copy from matrix 'a' to 'b'.
 * Memory locations of a and b are assumed to be different
 * without any overlapping memory
 * */
void fla_dtranspose(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *b,
                    aocl_int64_t *ldb)
{
    aocl_int64_t i, j;

    /* Offset adjustments */
    a -= (1 + *lda);
    b -= (1 + *ldb);

    /* Do the transpose copy */
    for(i = 1; i <= *n; i++)
    {
        for(j = 1; j <= *m; j++)
        {
            b[i + j * *ldb] = a[i * *lda + j];
        }
    }
}
/* 3x3 Householder Rotation */
int fla_dhrot3(aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *v, doublereal *tau)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_dhrot3_avx512(n, a, lda, v, tau);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dhrot3_avx2(n, a, lda, v, tau);
    }
    return 0;
}
/* 2x2 Plane Rotation */
int fla_drot(aocl_int64_t *n, doublereal *dx, aocl_int64_t *incx, doublereal *dy, aocl_int64_t *incy,
             doublereal *c__, doublereal *s)
{
#ifndef FLA_ENABLE_AOCL_BLAS
#endif

    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_drot_avx512(n, dx, incx, dy, incy, c__, s);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_drot_avx2(n, dx, incx, dy, incy, c__, s);
    }
    else
    {
        aocl_blas_drot(n, dx, incx, dy, incy, c__, s);
    }
    return 0;
}
void fla_zrot(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy, aocl_int64_t *incy,
              doublereal *c__, dcomplex *s)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_zrot_avx512(n, cx, incx, cy, incy, c__, s);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_zrot_avx2(n, cx, incx, cy, incy, c__, s);
    }
    else
    {
        aocl_lapack_zrot(n, (dcomplex *)cx, incx, (dcomplex *)cy, incy, c__, (dcomplex *)s);
    }
    return;
}
/* scomplex vector scaling when increment is 1 and specific threshold */
int fla_zscal(aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx)
{
    /* Initialize global context data */
    aocl_fla_init();

    /* Take AVX path only for increment equal to 1 and particular threshold size*/
    if(*incx == 1 && *n >= 1 && *n <= FLA_ZSCAL_INLINE_SMALL && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        if(FLA_IS_ARCH_ID(FLA_ARCH_AVX512))
        {
            fla_zscal_ix1_avx512(n, alpha, x);
        }
        else
        {
            fla_zscal_ix1_avx2(n, alpha, x);
        }
    }
    else
    {
        aocl_blas_zscal(n, (dcomplex *)alpha, (dcomplex *)x, incx);
    }
    return 0;
}
/* scales a vector by a constant when threshold <= 128 */
void fla_dscal(aocl_int64_t *n, doublereal *da, doublereal *dx, aocl_int64_t *incx)
{
    /* Initialize global context data */
    aocl_fla_init();

    if(*incx == 1 && *da != 0 && *n >= 1 && *n <= FLA_DSCAL_INLINE_SMALL
       && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
        {
            fla_dscal_ix1_avx512(n, da, dx, incx);
        }
        else
        {
            fla_dscal_ix1_avx2(n, da, dx, incx);
        }
    }
    else
    {
        aocl_blas_dscal(n, da, dx, incx);
    }
    return;
}
/* Double QR (DGEQRF) for small sizes */
int fla_dgeqrf_small(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *tau,
                     doublereal *work)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgeqrf_small_avx2(m, n, a, lda, tau, work);
    }
    return 0;
}
/* Double QR (DGEQRF) for small sizes */
int fla_dgelqf_small(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *tau,
                     doublereal *work)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        doublereal *at;

        /* Allocate transpose matrix */
        at = malloc(*n * *m * sizeof(doublereal));

        /* Do transpose and store it in at */
        fla_dtranspose(m, n, a, lda, at, n);

        /* Call QR for the transposed n x m matrix at */
        fla_dgeqrf_small_avx2(n, m, at, n, tau, work);

        /* Transpose at and store back in a */
        fla_dtranspose(n, m, at, n, a, lda);

        /* Free the transpose matrix */
        free(at);
    }
    return 0;
}
/* real vector scaling when increment is 1 */
void fla_sscal(aocl_int64_t *n, real *alpha, real *x, aocl_int64_t *incx)
{
    /* Initialize global context data */
    aocl_fla_init();

    /* Take AVX path only for increment equal to 1 */
    if(*incx == 1 && *alpha != 0 && *n >= 1 && *n <= FLA_SSCAL_INLINE_SMALL
       && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
        {
            fla_sscal_ix1_avx512(n, alpha, x);
        }
        else
        {
            fla_sscal_ix1_avx2(n, alpha, x);
        }
    }
    else
    {
        aocl_blas_sscal(n, (real *)alpha, (real *)x, incx);
    }
    return;
}
/* Rank 1 Operation */
void fla_sger(aocl_int64_t *m, aocl_int64_t *n, real *alpha, real *x, aocl_int64_t *incx, real *y, aocl_int64_t *incy,
              real *a, aocl_int64_t *lda)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_sger_avx2(m, n, alpha, x, incx, y, incy, a, lda);
    }
    else
    {
        aocl_blas_sger(m, n, alpha, x, incx, y, incy, a, lda);
    }
    return;
}

/* LU factorization.
 * To be used only when vectorized code via avx2/avx512 is enabled
 * */
int fla_dgetrf_small_simd(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                          aocl_int64_t *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_dgetrf_small_avx512(m, n, a, lda, ipiv, info);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgetrf_small_avx2(m, n, a, lda, ipiv, info);
    }
    return 0;
}

/* Double Complex LU for small sizes,
 * Optimized for AVX2 and AVX512 ISAs
 */
int fla_zgetrf_small_simd(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                          aocl_int64_t *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_zgetrf_small_avx512(m, n, a, lda, ipiv, info);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
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
void fla_dgesvd_xx_small10(aocl_int64_t wntus, aocl_int64_t wntvs, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ncu, doublereal *a,
                           aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                           aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_xx_small10_avx2(wntus, wntvs, m, n, ncu, a, lda, s, u, ldu, vt, ldvt, work, info);
    }
    return;
}

/* SVD for small fat-matrices in DGESVD
 */
void fla_dgesvd_xs_small10T(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *s,
                            doublereal *u, aocl_int64_t *ldu, doublereal *vt, aocl_int64_t *ldvt,
                            doublereal *work, aocl_int64_t *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_xs_small10T_avx2(m, n, a, lda, s, u, ldu, vt, ldvt, work, info);
    }
    return;
}

/* SVD for small fat-matrices with LQ factorization
 * already computed
 */
void fla_dgesvd_small6(aocl_int64_t wntus, aocl_int64_t wntvs, aocl_int64_t *m, aocl_int64_t *n, doublereal *a,
                       aocl_int64_t *lda, doublereal *qr, aocl_int64_t *ldqr, doublereal *s, doublereal *u,
                       aocl_int64_t *ldu, doublereal *vt, aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_small6_avx2(wntus, wntvs, m, n, a, lda, qr, ldqr, s, u, ldu, vt, ldvt, work,
                               info);
    }
    return;
}

/* SVD for small fat-matrices for path 1T in DGESVD
 */
void fla_dgesvd_nn_small1T(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *s,
                           doublereal *work, aocl_int64_t *info)
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
void fla_dgesvd_small6T(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, doublereal *ql,
                        aocl_int64_t *ldql, doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                        aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgesvd_small6T_avx2(m, n, a, lda, ql, ldql, s, u, ldu, vt, ldvt, work, info);
    }
    return;
}

/* Small DGETRS path (NOTRANS) should only be used for size between 3 to 8 and NRHS <= N */
int fla_dgetrs_small_notrans(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *a, aocl_int64_t *lda,
                             aocl_int_t *ipiv, doublereal *b, aocl_int64_t *ldb, aocl_int64_t *info)
{
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dgetrs_small_trsm_ll_avx2(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
    }
    return 0;
}

/* Find the maximum element from absolute values of a real vector */
real fla_get_max_sabs_element_vector(aocl_int64_t m, real *a, aocl_int64_t a_diml)
{
    real max_value = 0.0, temp;
    aocl_int64_t i__;
    /* Path when AVX512 ISA is supported */
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        max_value = fla_get_max_sabs_element_vector_avx512(m, a, a_diml);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        /* Path when AVX2 ISA is supported */
        max_value = fla_get_max_sabs_element_vector_avx2(m, a, a_diml);
    }
    else
    {
        for(i__ = 1; i__ <= m; ++i__)
        {
            temp = f2c_abs(a[i__ + a_diml]);
            if(max_value < temp || temp != temp)
            {
                max_value = temp;
            }
        }
    }
    return max_value;
}

/* Find the maximum element from absolute values of a doublereal vector */
doublereal fla_get_max_dabs_element_vector(aocl_int64_t m, doublereal *a, aocl_int64_t a_diml)
{
    doublereal max_value = 0.0, temp;
    aocl_int64_t i__;
    /* Path when AVX512 ISA is supported */
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        max_value = fla_get_max_dabs_element_vector_avx512(m, a, a_diml);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        /* Path when AVX2 ISA is supported */
        max_value = fla_get_max_dabs_element_vector_avx2(m, a, a_diml);
    }
    else
    {
        for(i__ = 1; i__ <= m; ++i__)
        {
            temp = f2c_abs(a[i__ + a_diml]);
            max_value = (max_value < temp || temp != temp) ? temp : max_value;
        }
    }
    return max_value;
}

/* Find the maximum element from absolute values of a scomplex vector */
real fla_get_max_cabs_element_vector(aocl_int64_t m, scomplex *a, aocl_int64_t a_diml)
{
    double c_abs(scomplex *);
    real max_value = 0.0, temp;
    aocl_int64_t i__;
    /* Path when AVX512 ISA is supported */
    if(m > FLA_CLANGEM_SIMD_AVX512_THRESH_M &&  FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        max_value = fla_get_max_cabs_element_vector_avx512(m, a, a_diml);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        /* Path when AVX2 ISA is supported */
        max_value = fla_get_max_cabs_element_vector_avx2(m, a, a_diml);
    }
    else
    {
        for(i__ = 1; i__ <= m; ++i__)
        {
            temp = c_abs(&a[i__ + a_diml]);
            max_value = (max_value < temp || temp != temp) ? temp : max_value;
        }
    }
    return max_value;
}

/* Find the maximum element from absolute values of a dcomplex vector */
doublereal fla_get_max_zabs_element_vector(aocl_int64_t m, dcomplex *a, aocl_int64_t a_diml)
{
    double z_abs(dcomplex *);
    doublereal max_value = 0.0, temp;
    aocl_int64_t i__;
    /* Path when AVX512 ISA is supported */
    if(m > FLA_ZLANGEM_SIMD_AVX512_THRESH_M && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        max_value = fla_get_max_zabs_element_vector_avx512(m, a, a_diml);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        /* Path when AVX2 ISA is supported */
        max_value = fla_get_max_zabs_element_vector_avx2(m, a, a_diml);
    }
    else
    {
        for(i__ = 1; i__ <= m; ++i__)
        {
            temp = z_abs(&a[i__ + a_diml]);
            max_value = (max_value < temp || temp != temp) ? temp : max_value;
        }
    }
    return max_value;
}

/* DLARF for small sizes
 * To be used only when vectorized code via avx2/avx512 is enabled
 * */
void fla_dlarf_small_incv1_simd(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff, aocl_int64_t ldr,
                                doublereal *v, doublereal ntau, doublereal *work)
{
    /* Select AVX512 kernel based on preset threshold and ISA support  */
    if(m > FLA_DGEMV_DGER_SIMD_AVX512_THRESH_M && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_dlarf_left_apply_incv1_avx512(m, n, a_buff, ldr, v, ntau, work);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dlarf_left_apply_incv1_avx2(m, n, a_buff, ldr, v, ntau, work);
    }
    return;
}

/* dnrm2 for small input sizes */
doublereal fla_dnrm2_blas_kernel(aocl_int64_t *sd, doublereal *a, aocl_int64_t *inc)
{
    doublereal value = 0.;
    /* TODO : Call DNRM2 AVX2 and AVX512 kernels using AOCL_BLAS_ENABLE 
       feature and call directly */
    if(*sd > FLA_DNRM2_SMALL_THRESH0 && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        value = fla_dnrm2_blas_avx512(sd, a, inc);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        value = fla_dnrm2_blas_avx2(sd, a, inc);
    }
    else
    {
        value = aocl_blas_dnrm2(sd, a, inc);
    }
    return value;
}

/* ZLARF optimized */
void fla_zlarf_left_invc1_opt(aocl_int64_t m, aocl_int64_t n, dcomplex *a_buff, aocl_int64_t ldr,
                              dcomplex *v, dcomplex *ntau, dcomplex *work)
{
    if(m < FLA_ZGEMV_ZGER_SIMD_AXV2_THRESH_M && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        /* If size is less than AVX2 threshold and AVX2 is available then execute AVX2 kernel */
        fla_zlarf_left_apply_incv1_avx2(m, n, a_buff, ldr, v, ntau, work);
    }
    else if(m < FLA_ZGEMV_ZGER_SIMD_AXV512_THRESH_M && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        /* AVX512 based on size and available ISA */
        fla_zlarf_left_apply_incv1_avx512(m, n, a_buff, ldr, v, ntau, work);
    }
    else
    {
        /* Original code */
        dcomplex c_b1 = {1., 0.};
        dcomplex c_b2 = {0., 0.};
        aocl_int64_t c__1 = 1;
        aocl_int64_t a_offset = 1 + ldr;

        /* w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1) */
        aocl_blas_zgemv("Conjugate transpose", &m, &n, &c_b1, (dcomplex *)&a_buff[a_offset], &ldr,
               (dcomplex *)&v[1], &c__1, (dcomplex *)&c_b2, (dcomplex *)&work[1], &c__1);
        /* C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H */
        aocl_blas_zgerc(&m, &n, (dcomplex *)ntau, (dcomplex *)&v[1], &c__1, (dcomplex *)&work[1], &c__1,
               (dcomplex *)&a_buff[a_offset], &ldr);
    }
}

#endif
