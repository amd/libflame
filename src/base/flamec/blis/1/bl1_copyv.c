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

void bl1_icopyv( conj1_t conj, fla_dim_t m, fla_dim_t* x, fla_dim_t incx, fla_dim_t* y, fla_dim_t incy )
{
	fla_dim_t*      chi;
	fla_dim_t*      psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}

void bl1_scopyv( conj1_t conj, fla_dim_t m, float* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
	bl1_scopy( m,
	           x, incx, 
	           y, incy );
}

void bl1_dcopyv( conj1_t conj, fla_dim_t m, double* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
	bl1_dcopy( m,
	           x, incx, 
	           y, incy );
}

void bl1_ccopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	bl1_ccopy( m,
	           x, incx, 
	           y, incy );

	if ( bl1_is_conj( conj ) )
		bl1_cconjv( m,
	                y, incy );
}

void bl1_zcopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	bl1_zcopy( m,
	           x, incx, 
	           y, incy );

	if ( bl1_is_conj( conj ) )
		bl1_zconjv( m,
	                y, incy );
}

// --- Mixed-datatype and general stride copy routines---------------

// sd ds
void bl1_sdcopyv( conj1_t conj, fla_dim_t m, float* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
	float*    chi;
	double*   psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}
void bl1_dscopyv( conj1_t conj, fla_dim_t m, double* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
	double*   chi;
	float*    psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}

// sc cs
void bl1_sccopyv( conj1_t conj, fla_dim_t m, float* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
	float*    chi;
	scomplex* psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0F;

		chi += incx;
		psi += incy;
	}
}
void bl1_cscopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
	scomplex* chi;
	float*    psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// sz zs
void bl1_szcopyv( conj1_t conj, fla_dim_t m, float* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
	float*    chi;
	dcomplex* psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0;

		chi += incx;
		psi += incy;
	}
}
void bl1_zscopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
	dcomplex* chi;
	float*    psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// dc cd
void bl1_dccopyv( conj1_t conj, fla_dim_t m, double* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
	double*   chi;
	scomplex* psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0F;

		chi += incx;
		psi += incy;
	}
}
void bl1_cdcopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
	scomplex* chi;
	double*   psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// dz zd
void bl1_dzcopyv( conj1_t conj, fla_dim_t m, double* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
	double*   chi;
	dcomplex* psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0;

		chi += incx;
		psi += incy;
	}
}
void bl1_zdcopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
	dcomplex* chi;
	double*   psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// cz zc
void bl1_czcopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
	scomplex* chi;
	dcomplex* psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = chi->real;
		psi->imag = chi->imag;

		chi += incx;
		psi += incy;
	}

	if ( bl1_is_conj( conj ) )
		bl1_zconjv( m,
	                y, incy );
}
void bl1_zccopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
	dcomplex* chi;
	scomplex* psi;
	fla_dim_t       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = chi->real;
		psi->imag = chi->imag;

		chi += incx;
		psi += incy;
	}

	if ( bl1_is_conj( conj ) )
		bl1_cconjv( m,
	                y, incy );
}

