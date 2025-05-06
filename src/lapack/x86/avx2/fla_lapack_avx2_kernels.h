/******************************************************************************
 * Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_LAPACK_AVX2_KERNELS_DEFS_H
#define FLA_LAPACK_AVX2_KERNELS_DEFS_H

/*! @file fla_lapack_avx2_kernels.h
 *  @brief AVX2 Kernel Declarations.
 *  */

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif
#include "fla_dgesvd_small_avx2.h"
#include "immintrin.h"

#if FLA_ENABLE_AMD_OPT
int fla_dhrot3_avx2(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau);
int fla_drot_avx2(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy,
                  doublereal *c__, doublereal *s);
int fla_zrot_avx2(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy,
                  doublereal *c__, doublecomplex *s);
int fla_sscal_ix1_avx2(integer *n, real *alpha, real *x);
int fla_dscal_ix1_avx2(integer *n, doublereal *da, doublereal *dx, integer *incx);
int fla_zscal_avx2(integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx);
int fla_zscal_ix1_avx2(integer *n, doublecomplex *alpha, doublecomplex *x);
int fla_sger_avx2(integer *m, integer *n, real *alpha, real *x, integer *incx, real *y,
                  integer *incy, real *a, integer *lda);
integer fla_dgetrf_small_avx2(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv,
                              integer *info);
integer fla_sgetrf_small_avx2(integer *m, integer *n, real *a, integer *lda, integer *ipiv,
                              integer *info);
void fla_lu_piv_small_d_update_tr_matrix_avx2(integer i_1, integer mi, integer ni, doublereal *acur,
                                              integer lda_t);
void fla_lu_piv_small_s_update_tr_matrix_avx2(integer i_1, integer mi, integer ni, real *acur,
                                              integer lda_t);
int fla_dgetrs_small_trsm_ll_avx2(char *trans, integer *n, integer *nrhs, doublereal *a,
                                  integer *lda, integer *ipiv, doublereal *b, integer *ldb,
                                  integer *info);
int fla_zgetrf_small_avx2(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                          integer *info);
int fla_dgeqrf_small_avx2(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau,
                          doublereal *work);
void fla_dgesvd_nn_small1T_avx2(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                                doublereal *work, integer *info);
void fla_dgesvd_small6_avx2(integer wntus, integer wntvs, integer *m, integer *n, doublereal *a,
                            integer *lda, doublereal *ql, integer *ldql, doublereal *s,
                            doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
                            doublereal *work, integer *info);
void fla_dgesvd_small6T_avx2(integer *m, integer *n, doublereal *a, integer *lda, doublereal *ql,
                             integer *ldql, doublereal *s, doublereal *u, integer *ldu,
                             doublereal *vt, integer *ldvt, doublereal *work, integer *info);
void fla_dgesvd_xx_small10_avx2(integer wntus, integer wntvs, integer *m, integer *n, doublereal *a,
                                integer *lda, doublereal *s, doublereal *u, integer *ldu,
                                doublereal *vt, integer *ldvt, doublereal *work, integer *info);
void fla_dgesvd_xs_small10T_avx2(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                                 doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
                                 doublereal *work, integer *info);
doublereal fla_get_max_dabs_element_vector_avx2(integer m, doublereal *a, integer a_dim);
real fla_get_max_sabs_element_vector_avx2(integer m, real *a, integer a_dim);
doublereal fla_get_max_zabs_element_vector_avx2(integer m, doublecomplex *a, integer a_dim);
real fla_get_max_cabs_element_vector_avx2(integer m, complex *a, integer a_dim);
void fla_dlarf_left_apply_incv1_avx2(integer m, integer n, doublereal *a_buff, integer ldr,
                                     doublereal *v, doublereal tau, doublereal *work);
doublereal fla_dnrm2_blas_avx2(integer *sd, doublereal *a, integer *inc);
void fla_zlarf_left_apply_incv1_avx2(integer m, integer n, doublecomplex *a_buff, integer ldr,
                                     doublecomplex *v, doublecomplex *ntau, doublecomplex *work);
#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_LAPACK_AVX2_KERNELS_DEFS_H */
