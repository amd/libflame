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

void bl1_sinvertv( conj1_t conj, fla_dim_t n, float* x, fla_dim_t incx )
{
	float  one = 1.0F;
	float* chi;
	fla_dim_t    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = one / *chi;
	}
}

void bl1_dinvertv( conj1_t conj, fla_dim_t n, double* x, fla_dim_t incx )
{
	double  one = 1.0;
	double* chi;
	fla_dim_t     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = one / *chi;
	}
}

void bl1_cinvertv( conj1_t conj, fla_dim_t n, scomplex* x, fla_dim_t incx )
{
	float     one = 1.0F;
	float     temp;
	float     s, xr_s, xi_s;
	float     conjsign;
	scomplex* chi;
	fla_dim_t       i;

	if ( bl1_is_conj( conj ) ) conjsign =  one;
	else                       conjsign = -one;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		s         = bl1_fmaxabs( chi->real, chi->imag ); \
		xr_s      = chi->real / s;
		xi_s      = chi->imag / s;
		temp      = xr_s * chi->real + xi_s * chi->imag;

		chi->real =            xr_s / temp;
		chi->imag = conjsign * xi_s / temp;
	}
}

void bl1_zinvertv( conj1_t conj, fla_dim_t n, dcomplex* x, fla_dim_t incx )
{
	double    one = 1.0;
	double    temp;
	double    s, xr_s, xi_s;
	double    conjsign;
	dcomplex* chi;
	fla_dim_t       i;

	if ( bl1_is_conj( conj ) ) conjsign =  one;
	else                       conjsign = -one;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		s         = bl1_fmaxabs( chi->real, chi->imag ); \
		xr_s      = chi->real / s;
		xi_s      = chi->imag / s;
		temp      = xr_s * chi->real + xi_s * chi->imag;

		chi->real =            xr_s / temp;
		chi->imag = conjsign * xi_s / temp;
	}
}

