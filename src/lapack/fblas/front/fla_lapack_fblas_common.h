/******************************************************************************
 * * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_fblas_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#ifdef FLA_ENABLE_AMD_OPT
fla_dim_t fla_idamax(fla_dim_t *n, doublereal *dx, fla_dim_t *incx);
#endif
