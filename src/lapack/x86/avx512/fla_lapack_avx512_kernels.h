/******************************************************************************
 * Copyright (C) 2023-24, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_lapack_avx512_kernels.h
 *  @brief AVX512 Kernel Declarations.
 *  */

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
int fla_dscal_ix1_avx512(integer *n, doublereal *da, doublereal *dx, integer *incx);
int fla_sscal_ix1_avx512(integer *n, real *alpha, real *x);
doublereal fla_get_max_abs_element_vector_avx512(integer m, doublereal *a, integer a_dim);

#endif
