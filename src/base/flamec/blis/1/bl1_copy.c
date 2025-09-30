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

void bl1_scopy( fla_dim_t m, float* x, fla_dim_t incx, float* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_scopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_scopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_dcopy( fla_dim_t m, double* x, fla_dim_t incx, double* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_dcopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_dcopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_ccopy( fla_dim_t m, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_ccopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_ccopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_zcopy( fla_dim_t m, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zcopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_zcopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

