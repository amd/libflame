/******************************************************************************
 * * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_fblas_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT

/* IDAMAX for small sizes (<= 128)*/
fla_dim_t fla_idamax(fla_dim_t *n, doublereal *dx, fla_dim_t *incx)
{
    return fla_idamax_small(n, dx, incx);
}
#endif
