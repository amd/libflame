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

void bl1_srandv( fla_dim_t n, float* x, fla_dim_t incx )
{
	float* chi;
	fla_dim_t    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_srands( chi );
	}
}

void bl1_drandv( fla_dim_t n, double* x, fla_dim_t incx )
{
	double* chi;
	fla_dim_t     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_drands( chi );
	}
}

void bl1_crandv( fla_dim_t n, scomplex* x, fla_dim_t incx )
{
	scomplex* chi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_crands( chi );
	}
}

void bl1_zrandv( fla_dim_t n, dcomplex* x, fla_dim_t incx )
{
	dcomplex* chi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_zrands( chi );
	}
}

