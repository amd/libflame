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

void bl1_isetv( fla_dim_t n, fla_dim_t* sigma, fla_dim_t* x, fla_dim_t incx )
{
	fla_dim_t*   chi;
	fla_dim_t    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_ssetv( fla_dim_t n, float* sigma, float* x, fla_dim_t incx )
{
	float* chi;
	fla_dim_t    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_dsetv( fla_dim_t n, double* sigma, double* x, fla_dim_t incx )
{
	double* chi;
	fla_dim_t     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_csetv( fla_dim_t n, scomplex* sigma, scomplex* x, fla_dim_t incx )
{
	scomplex* chi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

void bl1_zsetv( fla_dim_t n, dcomplex* sigma, dcomplex* x, fla_dim_t incx )
{
	dcomplex* chi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

