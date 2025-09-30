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

void bl1_sconjm( fla_dim_t m, fla_dim_t n, float* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	return;
}

void bl1_dconjm( fla_dim_t m, fla_dim_t n, double* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	return;
}

void bl1_cconjm( fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	float   m1 = bl1_sm1();
	float*  a_conj;
	fla_dim_t     lda, inca;
	fla_dim_t     n_iter;
	fla_dim_t     n_elem;
	fla_dim_t     j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	for ( j = 0; j < n_iter; ++j )
	{
		a_conj = ( float* )( a + j*lda ) + 1;

		bl1_sscal( n_elem,
		           &m1,
		           a_conj, 2*inca );
	}
}

void bl1_zconjm( fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs )
{
	double  m1 = bl1_dm1();
	double* a_conj;
	fla_dim_t     lda, inca;
	fla_dim_t     n_iter;
	fla_dim_t     n_elem;
	fla_dim_t     j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A is a vector to ensure that the underlying axpy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for a vector.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, a_rs, a_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;

		// An optimization: if A is row-major, then let's access the matrix
		// by rows instead of by columns to increase spatial locality.
		if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			bl1_swap_ints( n_iter, n_elem );
			bl1_swap_ints( lda, inca );
		}
	}

	for ( j = 0; j < n_iter; ++j )
	{
		a_conj = ( double* )( a + j*lda ) + 1;

		bl1_dscal( n_elem,
		           &m1,
		           a_conj, 2*inca );
	}
}
