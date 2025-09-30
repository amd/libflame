/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
*     Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "blis1.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

void bl1_saxpy( fla_dim_t n, float* alpha, float* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_saxpy( n,
	             *alpha,
	             x, incx,
	             y, incy );
#else
	F77_saxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

void bl1_daxpy( fla_dim_t n, double* alpha, double* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_daxpy( n,
	             *alpha,
	             x, incx,
	             y, incy );
#else
	F77_daxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

void bl1_caxpy( fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_caxpy( n,
	             alpha,
	             x, incx,
	             y, incy );
#else
	F77_caxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

void bl1_zaxpy( fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zaxpy( n,
	             alpha,
	             x, incx,
	             y, incy );
#else
	F77_zaxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

