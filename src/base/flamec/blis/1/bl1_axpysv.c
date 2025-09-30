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

void bl1_saxpysv( fla_dim_t n, float* alpha0, float* alpha1, float* x, fla_dim_t incx, float* beta, float* y, fla_dim_t incy )
{
	float    alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	bl1_sscal( n,
	           beta,
	           y, incy );

	bl1_saxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bl1_daxpysv( fla_dim_t n, double* alpha0, double* alpha1, double* x, fla_dim_t incx, double* beta, double* y, fla_dim_t incy )
{
	double   alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	bl1_dscal( n,
	           beta,
	           y, incy );

	bl1_daxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bl1_caxpysv( fla_dim_t n, scomplex* alpha0, scomplex* alpha1, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy )
{
	scomplex alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	bl1_cscal( n,
	           beta,
	           y, incy );

	bl1_caxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bl1_zaxpysv( fla_dim_t n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy )
{
	dcomplex alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	bl1_zscal( n,
	           beta,
	           y, incy );

	bl1_zaxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

