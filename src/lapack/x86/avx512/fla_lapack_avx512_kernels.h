/******************************************************************************
 * Copyright (C) 2023-25, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_LAPACK_AVX512_KERNELS_DEFS_H
#define FLA_LAPACK_AVX512_KERNELS_DEFS_H
/*! @file fla_lapack_avx512_kernels.h
 *  @brief AVX512 Kernel Declarations.
 *  */
#include "FLAME.h"
#include "immintrin.h"

#if FLA_ENABLE_AMD_OPT

int fla_drot_avx512(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy,
                    doublereal *c__, doublereal *s);
int fla_zrot_avx512(integer *n, doublecomplex *dx, integer *incx, doublecomplex *dy, integer *incy,
                    doublereal *c__, doublecomplex *s);
int fla_dhrot3_avx512(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau);
int fla_zgetrf_small_avx512(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                            integer *info);
integer fla_dgetrf_small_avx512(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv,
                                integer *info);
integer fla_sgetrf_small_avx512(integer *m, integer *n, real *a, integer *lda, integer *ipiv,
                                integer *info);
int fla_dscal_ix1_avx512(integer *n, doublereal *da, doublereal *dx, integer *incx);
int fla_sscal_ix1_avx512(integer *n, real *alpha, real *x);
int fla_zscal_ix1_avx512(integer *n, doublecomplex *alpha, doublecomplex *x);
doublereal fla_get_max_dabs_element_vector_avx512(integer m, doublereal *a, integer a_dim);
real fla_get_max_sabs_element_vector_avx512(integer m, real *a, integer a_dim);
doublereal fla_get_max_zabs_element_vector_avx512(integer m, doublecomplex *a, integer a_dim);
real fla_get_max_cabs_element_vector_avx512(integer m, complex *a, integer a_dim);
void fla_dlarf_left_apply_incv1_avx512(integer m, integer n, doublereal *a_buff, integer ldr,
                                        doublereal *v, doublereal ntau, doublereal *work);
doublereal fla_dnrm2_blas_avx512(integer *sd, doublereal *a, integer *inc);
void fla_zlarf_left_apply_incv1_avx512(integer m, integer n, doublecomplex *a_buff, integer ldr,
                                       doublecomplex *v, doublecomplex *ntau, doublecomplex *work);

#endif
#endif
