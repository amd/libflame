/******************************************************************************
 * * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_x86_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"

#if FLA_ENABLE_AMD_OPT
void fla_dtranspose(integer *m, integer *n, doublereal *a, integer *lda, doublereal *b,
                    integer *ldb);
int fla_dhrot3(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau);
int fla_drot(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy,
             doublereal *c__, doublereal *s);
int fla_zscal(integer *n, doublecomplex *alpha, doublecomplex *x, integer *incx);
int fla_dgeqrf_small(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau,
                     doublereal *work);
int fla_dgelqf_small(integer *m, integer *n, doublereal *a, integer *lda, doublereal *tau,
                     doublereal *work);
void fla_dscal(integer *n, doublereal *alpha, doublereal *x, integer *incx);
void fla_sscal(integer *n, real *alpha, real *x, integer *incx);
void fla_sger(integer *m, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy,
              real *a, integer *lda);
int fla_zgetrf_small_simd(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                          integer *info);
void fla_dgesvd_small6(integer wntus, integer wntvs, integer *m, integer *n, doublereal *a,
                       integer *lda, doublereal *qr, integer *ldqr, doublereal *s, doublereal *u,
                       integer *ldu, doublereal *vt, integer *ldvt, doublereal *work,
                       integer *info);
void fla_dgesvd_nn_small1T(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                           doublereal *work, integer *info);
void fla_dgesvd_small6T(integer *m, integer *n, doublereal *a, integer *lda, doublereal *ql,
                        integer *ldql, doublereal *s, doublereal *u, integer *ldu, doublereal *vt,
                        integer *ldvt, doublereal *work, integer *info);
void fla_dgesvd_xx_small10(integer wntus, integer wntvs, integer *m, integer *n, doublereal *a,
                           integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt,
                           integer *ldvt, doublereal *work, integer *info);
void fla_dgesvd_xs_small10T(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                            doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
                            doublereal *work, integer *info);
int fla_dgetrs_small_notrans(char *trans, integer *n, integer *nrhs, doublereal *a, integer *lda,
                             integer *ipiv, doublereal *b, integer *ldb, integer *info);
void lapack_getri_small_d(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work,
                          integer *info);
doublereal fla_get_max_abs_element_vector(integer m, doublereal *a, integer a_dim);
#endif
