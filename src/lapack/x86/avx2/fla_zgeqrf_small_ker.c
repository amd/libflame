/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_zgeqrf_small_ker.c
 *  @brief QR for small inputs in AVX2.
 *  */

#include "FLAME.h"
#include "fla_geqrf_small_cmplx_avx2_kernel.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

FLA_DEFINE_GEQRF_CMPLX_SMALL(fla_zgeqrf_small_ker, z, dcomplex, doublereal, dlamch_, fabs, sqrt)

#endif
