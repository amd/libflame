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

void bl1_sconjv( fla_dim_t m, float* x, fla_dim_t incx )
{
	return;
}

void bl1_dconjv( fla_dim_t m, double* x, fla_dim_t incx )
{
	return;
}

void bl1_cconjv( fla_dim_t m, scomplex* x, fla_dim_t incx )
{
	float  m1        = bl1_sm1();
	float* x_conj    = ( float* ) x + 1;
	fla_dim_t    incx_conj = 2 * incx;

	bl1_sscal( m,
	           &m1,
	           x_conj, incx_conj );
}

void bl1_zconjv( fla_dim_t m, dcomplex* x, fla_dim_t incx )
{
	double  m1        = bl1_dm1();
	double* x_conj    = ( double* ) x + 1;
	fla_dim_t     incx_conj = 2 * incx;

	bl1_dscal( m,
	           &m1,
	           x_conj, incx_conj );
}

