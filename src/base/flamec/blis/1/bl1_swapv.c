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

void bl1_sswapv( fla_dim_t n, float* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_sswap( n,
	           x, incx, 
	           y, incy );
}

void bl1_dswapv( fla_dim_t n, double* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_dswap( n,
	           x, incx, 
	           y, incy );
}

void bl1_cswapv( fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_cswap( n,
	           x, incx, 
	           y, incy );
}

void bl1_zswapv( fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_zswap( n,
	           x, incx, 
	           y, incy );
}

