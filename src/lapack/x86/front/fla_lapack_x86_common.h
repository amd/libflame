/******************************************************************************
 * * Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_lapack_x86_common.c
 *  @brief Common front-end functions
 *         to choose optimized paths
 *  *  */

#ifdef FLA_ENABLE_AMD_OPT
int fla_dhrot3(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau);
int fla_drot(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *c__, doublereal *s);
#endif