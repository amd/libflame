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

void bl1_ssymv( uplo1_t uplo, fla_dim_t m, float* alpha, float* a, fla_dim_t a_rs, fla_dim_t a_cs, float* x, fla_dim_t incx, float* beta, float* y, fla_dim_t incy )
{
	float*    a_save    = a;
	fla_dim_t       a_rs_save = a_rs;
	fla_dim_t       a_cs_save = a_cs;
	fla_dim_t       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_screate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	bl1_ssymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_dsymv( uplo1_t uplo, fla_dim_t m, double* alpha, double* a, fla_dim_t a_rs, fla_dim_t a_cs, double* x, fla_dim_t incx, double* beta, double* y, fla_dim_t incy )
{
	double*   a_save    = a;
	fla_dim_t       a_rs_save = a_rs;
	fla_dim_t       a_cs_save = a_cs;
	fla_dim_t       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_dcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	bl1_dsymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_csymv( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy )
{
	scomplex* a_save    = a;
	fla_dim_t       a_rs_save = a_rs;
	fla_dim_t       a_cs_save = a_cs;
	fla_dim_t       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_ccreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	bl1_csymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_zsymv( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy )
{
	dcomplex* a_save    = a;
	fla_dim_t       a_rs_save = a_rs;
	fla_dim_t       a_cs_save = a_cs;
	fla_dim_t       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_zcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	bl1_zsymv_blas( uplo,
	                m,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_ssymv_blas( uplo1_t uplo, fla_dim_t m, float* alpha, float* a, fla_dim_t lda, float* x, fla_dim_t incx, float* beta, float* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_ssymv( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_ssymv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bl1_dsymv_blas( uplo1_t uplo, fla_dim_t m, double* alpha, double* a, fla_dim_t lda, double* x, fla_dim_t incx, double* beta, double* y, fla_dim_t incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_dsymv( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_dsymv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bl1_csymv_blas( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy )
{
	scomplex* x_copy;
	scomplex* y_copy;
	fla_dim_t       n   = 1;
	fla_dim_t       ldx = m;
	fla_dim_t       ldy = m;

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( BLIS1_LEFT, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	x_copy = bl1_callocv( m );
	y_copy = bl1_callocv( m );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	cblas_csymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             alpha,
	             a,      lda,
	             x_copy, ldx,
	             beta,
	             y_copy, ldy );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bl1_cfree( x_copy );
	bl1_cfree( y_copy );

#else
	char blas_side;
	char blas_uplo;

	bl1_param_map_to_netlib_side( BLIS1_LEFT, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	x_copy = bl1_callocv( m );
	y_copy = bl1_callocv( m );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	F77_csymm ( &blas_side,
	            &blas_uplo,
	            &m,
	            &n,
	            alpha,
	            a,      &lda,
	            x_copy, &ldx,
	            beta,
	            y_copy, &ldy );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bl1_cfree( x_copy );
	bl1_cfree( y_copy );
#endif
}

void bl1_zsymv_blas( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy )
{
	dcomplex*        x_copy;
	dcomplex*        y_copy;
	fla_dim_t              n   = 1;
	fla_dim_t              ldx = m;
	fla_dim_t              ldy = m;

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( BLIS1_LEFT, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	x_copy = bl1_zallocv( m );
	y_copy = bl1_zallocv( m );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	cblas_zsymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             alpha,
	             a,      lda,
	             x_copy, ldx,
	             beta,
	             y_copy, ldy );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bl1_zfree( x_copy );
	bl1_zfree( y_copy );

#else
	char blas_side;
	char blas_uplo;

	bl1_param_map_to_netlib_side( BLIS1_LEFT, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	x_copy = bl1_zallocv( m );
	y_copy = bl1_zallocv( m );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	F77_zsymm ( &blas_side,
	            &blas_uplo,
	            &m,
	            &n,
	            alpha,
	            a,      &lda,
	            x_copy, &ldx,
	            beta,
	            y_copy, &ldy );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y_copy, 1,
	            y,      incy );

	bl1_zfree( x_copy );
	bl1_zfree( y_copy );
#endif
}

