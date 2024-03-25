/******************************************************************************
 * Copyright (C) 2023-24, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_lapack_avx512_kernels.h
 *  @brief AVX512 Kernel Declarations.
 *  */

#include "immintrin.h"

#if FLA_ENABLE_AMD_OPT
int fla_zgetrf_small_avx512(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                            integer *info);
integer fla_dgetrf_small_avx512(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv,
                                integer *info);
int fla_dscal_ix1_avx512(integer *n, doublereal *da, doublereal *dx, integer *incx);
#endif

