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

void bl1_sscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, float* sigma, float* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	float* alpha;
	fla_dim_t    i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		*alpha *= *sigma;

		++i;
		++j;
	}
}

void bl1_dscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, double* sigma, double* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	double* alpha;
	fla_dim_t     i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		*alpha *= *sigma;

		++i;
		++j;
	}
}

void bl1_csscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, float* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	scomplex* alpha;
	fla_dim_t       i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real *= *sigma;
		alpha->imag *= *sigma;

		++i;
		++j;
	}
}

void bl1_zdscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, double* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	dcomplex* alpha;
	fla_dim_t       i, j;

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		alpha->real *= *sigma;
		alpha->imag *= *sigma;

		++i;
		++j;
	}
}

void bl1_cscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, scomplex* sigma, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	scomplex* alpha;
	scomplex  sigma_conj;
	fla_dim_t       i, j;

	bl1_ccopys( conj, sigma, &sigma_conj );

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		bl1_cscals( &sigma_conj, alpha );

		++i;
		++j;
	}
}

void bl1_zscalediag( conj1_t conj, fla_dim_t offset, fla_dim_t m, fla_dim_t n, dcomplex* sigma, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	dcomplex* alpha;
	dcomplex  sigma_conj;
	fla_dim_t       i, j;

	bl1_zcopys( conj, sigma, &sigma_conj );

	i = j = 0;

	if      ( offset < 0 ) i = -offset;
	else if ( offset > 0 ) j =  offset;
	
	while ( i < m && j < n )
	{
		alpha = a + i*a_rs + j*a_cs;
	
		bl1_zscals( &sigma_conj, alpha );

		++i;
		++j;
	}
}

