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

void bl1_srandm( fla_dim_t m, fla_dim_t n, float* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	float*    a_begin;
	fla_dim_t       inca, lda;
	fla_dim_t       n_iter;
	fla_dim_t       n_elem;
	fla_dim_t       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_srandv( n_elem,
		            a_begin, inca );
	}
}

void bl1_drandm( fla_dim_t m, fla_dim_t n, double* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	double*   a_begin;
	fla_dim_t       inca, lda;
	fla_dim_t       n_iter;
	fla_dim_t       n_elem;
	fla_dim_t       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_drandv( n_elem,
		            a_begin, inca );
	}
}

void bl1_crandm( fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	scomplex* a_begin;
	fla_dim_t       inca, lda;
	fla_dim_t       n_iter;
	fla_dim_t       n_elem;
	fla_dim_t       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_crandv( n_elem,
		            a_begin, inca );
	}
}

void bl1_zrandm( fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	dcomplex* a_begin;
	fla_dim_t       inca, lda;
	fla_dim_t       n_iter;
	fla_dim_t       n_elem;
	fla_dim_t       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;

		bl1_zrandv( n_elem,
		            a_begin, inca );
	}
}

