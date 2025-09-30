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

void bl1_sident( fla_dim_t m, float* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	float* alpha;
	fla_dim_t    i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = 0.0F;

			if ( i == j )
				*alpha = 1.0F;
		}
	}
}

void bl1_dident( fla_dim_t m, double* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	double* alpha;
	fla_dim_t     i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			*alpha = 0.0;

			if ( i == j )
				*alpha = 1.0;
		}
	}
}

void bl1_cident( fla_dim_t m, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	scomplex* alpha;
	fla_dim_t       i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = 0.0F;
			alpha->imag = 0.0F;

			if ( i == j )
				alpha->real = 1.0F;
		}
	}
}

void bl1_zident( fla_dim_t m, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	dcomplex* alpha;
	fla_dim_t       i, j;

	for ( j = 0; j < m; ++j )
	{
		for ( i = 0; i < m; ++i )
		{
			alpha = a + i*a_rs + j*a_cs;
	
			alpha->real = 0.0;
			alpha->imag = 0.0;

			if ( i == j )
				alpha->real = 1.0;
		}
	}
}

