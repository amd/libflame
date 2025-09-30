/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-1 BLAS-like prototypes --------------------------------------------

// --- amax ---

void bl1_samax( fla_dim_t n, float*    x, fla_dim_t incx, fla_dim_t* index );
void bl1_damax( fla_dim_t n, double*   x, fla_dim_t incx, fla_dim_t* index );
void bl1_camax( fla_dim_t n, scomplex* x, fla_dim_t incx, fla_dim_t* index );
void bl1_zamax( fla_dim_t n, dcomplex* x, fla_dim_t incx, fla_dim_t* index );

// --- asum ---

void bl1_sasum( fla_dim_t n, float*    x, fla_dim_t incx, float*  norm );
void bl1_dasum( fla_dim_t n, double*   x, fla_dim_t incx, double* norm );
void bl1_casum( fla_dim_t n, scomplex* x, fla_dim_t incx, float*  norm );
void bl1_zasum( fla_dim_t n, dcomplex* x, fla_dim_t incx, double* norm );

// --- axpy ---

void bl1_saxpy( fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_daxpy( fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_caxpy( fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zaxpy( fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

// --- axpyv ---

void bl1_saxpyv( conj1_t conj, fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_daxpyv( conj1_t conj, fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_caxpyv( conj1_t conj, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zaxpyv( conj1_t conj, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

// --- axpymt ---

void bl1_saxpymt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_daxpymt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_caxpymt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zaxpymt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- axpymrt ---

void bl1_saxpymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_daxpymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_caxpymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zaxpymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- axpysv ---

void bl1_saxpysv( fla_dim_t n, float*    alpha0, float*    alpha1, float*    x, fla_dim_t incx, float*    beta, float*    y, fla_dim_t incy );
void bl1_daxpysv( fla_dim_t n, double*   alpha0, double*   alpha1, double*   x, fla_dim_t incx, double*   beta, double*   y, fla_dim_t incy );
void bl1_caxpysv( fla_dim_t n, scomplex* alpha0, scomplex* alpha1, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_zaxpysv( fla_dim_t n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

// --- axpysmt ---

void bl1_saxpysmt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    alpha0, float*    alpha1, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    beta, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_daxpysmt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   alpha0, double*   alpha1, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   beta, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_caxpysmt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* alpha0, scomplex* alpha1, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* beta, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zaxpysmt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* beta, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- conjv ---

void bl1_sconjv( fla_dim_t m, float* x, fla_dim_t incx );
void bl1_dconjv( fla_dim_t m, double* x, fla_dim_t incx );
void bl1_cconjv( fla_dim_t m, scomplex* x, fla_dim_t incx );
void bl1_zconjv( fla_dim_t m, dcomplex* x, fla_dim_t incx );

// --- conjm ---

void bl1_sconjm( fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dconjm( fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cconjm( fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zconjm( fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- conjmr ---

void bl1_sconjmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dconjmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cconjmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zconjmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- copy ---

void bl1_scopy( fla_dim_t m, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_dcopy( fla_dim_t m, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_ccopy( fla_dim_t m, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zcopy( fla_dim_t m, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

// --- copyv ---

void bl1_icopyv( conj1_t conj, fla_dim_t m, fla_dim_t*      x, fla_dim_t incx, fla_dim_t*      y, fla_dim_t incy );
void bl1_scopyv( conj1_t conj, fla_dim_t m, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_dcopyv( conj1_t conj, fla_dim_t m, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_ccopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zcopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

void bl1_sdcopyv( conj1_t conj, fla_dim_t m, float*    x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_dscopyv( conj1_t conj, fla_dim_t m, double*   x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_sccopyv( conj1_t conj, fla_dim_t m, float*    x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_cscopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_szcopyv( conj1_t conj, fla_dim_t m, float*    x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );
void bl1_zscopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_dccopyv( conj1_t conj, fla_dim_t m, double*   x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_cdcopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_dzcopyv( conj1_t conj, fla_dim_t m, double*   x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );
void bl1_zdcopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_czcopyv( conj1_t conj, fla_dim_t m, scomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );
void bl1_zccopyv( conj1_t conj, fla_dim_t m, dcomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );

// --- copymr ---

void bl1_scopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ccopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

void bl1_sscopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_sdcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dscopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_sccopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cscopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_szcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zscopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ddcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dccopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cdcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dzcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zdcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cccopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_czcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zccopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zzcopymr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- copymrt ---

void bl1_scopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ccopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

void bl1_sscopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_sdcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_sccopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_szcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dscopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ddcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dccopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dzcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cscopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cdcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cccopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_czcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zscopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zdcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zccopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zzcopymrt( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- copymt ---

void bl1_icopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, fla_dim_t*      a, fla_dim_t a_rs, fla_dim_t a_cs, fla_dim_t*      b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_scopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ccopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

void bl1_sscopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_sdcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dscopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_sccopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cscopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_szcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zscopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ddcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dccopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cdcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dzcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zdcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cccopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_czcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zccopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zzcopymt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

// --- dot ---

void bl1_cdot_in( conj1_t conj, fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* rho );
void bl1_zdot_in( conj1_t conj, fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* rho );

void bl1_sdot( conj1_t conj, fla_dim_t n, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    rho );
void bl1_ddot( conj1_t conj, fla_dim_t n, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   rho );
void bl1_cdot( conj1_t conj, fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* rho );
void bl1_zdot( conj1_t conj, fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* rho );

// --- dots ---

void bl1_sdots( conj1_t conj, fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    beta, float*    rho );
void bl1_ddots( conj1_t conj, fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   beta, double*   rho );
void bl1_cdots( conj1_t conj, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* beta, scomplex* rho );
void bl1_zdots( conj1_t conj, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* beta, dcomplex* rho );

// --- dot2s ---

void bl1_sdot2s( conj1_t conj, fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    beta, float*    rho );
void bl1_ddot2s( conj1_t conj, fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   beta, double*   rho );
void bl1_cdot2s( conj1_t conj, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* beta, scomplex* rho );
void bl1_zdot2s( conj1_t conj, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* beta, dcomplex* rho );

// --- fnorm ---

void bl1_sfnorm( fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*  norm );
void bl1_dfnorm( fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double* norm );
void bl1_cfnorm( fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*  norm );
void bl1_zfnorm( fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double* norm );

// --- invscalv ---

void bl1_sinvscalv(  conj1_t conj, fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx );
void bl1_dinvscalv(  conj1_t conj, fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx );
void bl1_csinvscalv( conj1_t conj, fla_dim_t n, float*    alpha, scomplex* x, fla_dim_t incx );
void bl1_cinvscalv(  conj1_t conj, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx );
void bl1_zdinvscalv( conj1_t conj, fla_dim_t n, double*   alpha, dcomplex* x, fla_dim_t incx );
void bl1_zinvscalv(  conj1_t conj, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx );

// --- invscalm ---

void bl1_sinvscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dinvscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csinvscalm( conj1_t conj, fla_dim_t m, fla_dim_t n, float*    alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cinvscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zdinvscalm( conj1_t conj, fla_dim_t m, fla_dim_t n, double*   alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zinvscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- nrm2 ---

void bl1_snrm2( fla_dim_t n, float*    x, fla_dim_t incx, float*  norm );
void bl1_dnrm2( fla_dim_t n, double*   x, fla_dim_t incx, double* norm );
void bl1_cnrm2( fla_dim_t n, scomplex* x, fla_dim_t incx, float*  norm );
void bl1_znrm2( fla_dim_t n, dcomplex* x, fla_dim_t incx, double* norm );

// --- scal ---

void bl1_sscal(  fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx );
void bl1_dscal(  fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx );
void bl1_csscal( fla_dim_t n, float*    alpha, scomplex* x, fla_dim_t incx );
void bl1_cscal(  fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx );
void bl1_zdscal( fla_dim_t n, double*   alpha, dcomplex* x, fla_dim_t incx );
void bl1_zscal(  fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx );

// --- scalv ---

void bl1_sscalv(  conj1_t conj, fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx );
void bl1_dscalv(  conj1_t conj, fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx );
void bl1_csscalv( conj1_t conj, fla_dim_t n, float*    alpha, scomplex* x, fla_dim_t incx );
void bl1_cscalv(  conj1_t conj, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx );
void bl1_zdscalv( conj1_t conj, fla_dim_t n, double*   alpha, dcomplex* x, fla_dim_t incx );
void bl1_zscalv(  conj1_t conj, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx );

// --- scalm ---

void bl1_sscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csscalm( conj1_t conj, fla_dim_t m, fla_dim_t n, float*    alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zdscalm( conj1_t conj, fla_dim_t m, fla_dim_t n, double*   alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zscalm(  conj1_t conj, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- scalmr ---

void bl1_sscalmr(  uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dscalmr(  uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csscalmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cscalmr(  uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zdscalmr( uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zscalmr(  uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

// --- swap ---

void bl1_sswap( fla_dim_t n, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_dswap( fla_dim_t n, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_cswap( fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zswap( fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

// --- swapv ---

void bl1_sswapv( fla_dim_t n, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy );
void bl1_dswapv( fla_dim_t n, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy );
void bl1_cswapv( fla_dim_t n, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy );
void bl1_zswapv( fla_dim_t n, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy );

// --- swapmt ---

void bl1_sswapmt( trans1_t trans, fla_dim_t m, fla_dim_t n, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dswapmt( trans1_t trans, fla_dim_t m, fla_dim_t n, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_cswapmt( trans1_t trans, fla_dim_t m, fla_dim_t n, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_zswapmt( trans1_t trans, fla_dim_t m, fla_dim_t n, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

