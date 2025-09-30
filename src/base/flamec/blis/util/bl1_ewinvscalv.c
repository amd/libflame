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

void bl1_sewinvscalv( conj1_t conj, fla_dim_t n, float* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
	float*    chi;
	float*    psi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_sinvscals( chi, psi );
	}
}

void bl1_dewinvscalv( conj1_t conj, fla_dim_t n, double* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
	double*   chi;
	double*   psi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_dinvscals( chi, psi );
	}
}

void bl1_csewinvscalv( conj1_t conj, fla_dim_t n, float* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
	float*    chi;
	scomplex* psi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_csinvscals( chi, psi );
	}
}

void bl1_cewinvscalv( conj1_t conj, fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
	scomplex* chi;
	scomplex* psi;
	scomplex  conjchi;
	fla_dim_t       i;

	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bl1_ccopyconj( chi, &conjchi );
			bl1_cinvscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bl1_cinvscals( chi, psi );
		}
	}
}

void bl1_zdewinvscalv( conj1_t conj, fla_dim_t n, double* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
	double*   chi;
	dcomplex* psi;
	fla_dim_t       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_zdinvscals( chi, psi );
	}
}

void bl1_zewinvscalv( conj1_t conj, fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
	dcomplex* chi;
	dcomplex* psi;
	dcomplex  conjchi;
	fla_dim_t       i;

	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bl1_zcopyconj( chi, &conjchi );
			bl1_zinvscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bl1_zinvscals( chi, psi );
		}
	}
}

