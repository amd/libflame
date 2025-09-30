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

void bl1_isetm( fla_dim_t m, fla_dim_t n, fla_dim_t* sigma, fla_dim_t* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	fla_dim_t*   alpha;
	fla_dim_t    i, j;

	for ( j = 0; j < n; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = *sigma;
		}
	}
}

void bl1_ssetm( fla_dim_t m, fla_dim_t n, float* sigma, float* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	float* alpha;
	fla_dim_t    i, j;

	for ( j = 0; j < n; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = *sigma;
		}
	}
}

void bl1_dsetm( fla_dim_t m, fla_dim_t n, double* sigma, double* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	double* alpha;
	fla_dim_t     i, j;

	for ( j = 0; j < n; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = *sigma;
		}
	}
}

void bl1_csetm( fla_dim_t m, fla_dim_t n, scomplex* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	scomplex* alpha;
	fla_dim_t       i, j;

	for ( j = 0; j < n; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = sigma->real;
			alpha->imag = sigma->imag;
		}
	}
}

void bl1_zsetm( fla_dim_t m, fla_dim_t n, dcomplex* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	dcomplex* alpha;
	fla_dim_t       i, j;

	for ( j = 0; j < n; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = sigma->real;
			alpha->imag = sigma->imag;
		}
	}
}

