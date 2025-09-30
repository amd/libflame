/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-2 BLAS-like prototypes --------------------------------------------

// --- gemv ---

void bl1_sgemv( trans1_t transa, conj1_t conjx, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    x, fla_dim_t incx, float*    beta, float*    y, fla_dim_t incy );
void bl1_dgemv( trans1_t transa, conj1_t conjx, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   x, fla_dim_t incx, double*   beta, double*   y, fla_dim_t incy );
void bl1_cgemv( trans1_t transa, conj1_t conjx, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_zgemv( trans1_t transa, conj1_t conjx, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

void bl1_sgemv_blas( trans1_t transa, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t lda, float*    x, fla_dim_t incx, float*    beta, float*    y, fla_dim_t incy );
void bl1_dgemv_blas( trans1_t transa, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t lda, double*   x, fla_dim_t incx, double*   beta, double*   y, fla_dim_t incy );
void bl1_cgemv_blas( trans1_t transa, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_zgemv_blas( trans1_t transa, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

// --- ger ---

void bl1_sger( conj1_t conjx, conj1_t conjy, fla_dim_t m, fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dger( conj1_t conjx, conj1_t conjy, fla_dim_t m, fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cger( conj1_t conjx, conj1_t conjy, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zger( conj1_t conjx, conj1_t conjy, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

void bl1_sger_blas(  fla_dim_t m, fla_dim_t n, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    a, fla_dim_t lda );
void bl1_dger_blas(  fla_dim_t m, fla_dim_t n, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   a, fla_dim_t lda );
void bl1_cgerc_blas( fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* a, fla_dim_t lda );
void bl1_cgeru_blas( fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* a, fla_dim_t lda );
void bl1_zgerc_blas( fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* a, fla_dim_t lda );
void bl1_zgeru_blas( fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* a, fla_dim_t lda );

// --- hemv ---

void bl1_shemv( uplo1_t uplo, conj1_t conj, fla_dim_t m, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    x, fla_dim_t incx, float*    beta, float*    y, fla_dim_t incy );
void bl1_dhemv( uplo1_t uplo, conj1_t conj, fla_dim_t m, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   x, fla_dim_t incx, double*   beta, double*   y, fla_dim_t incy );
void bl1_chemv( uplo1_t uplo, conj1_t conj, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_zhemv( uplo1_t uplo, conj1_t conj, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

void bl1_chemv_blas( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_zhemv_blas( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

// --- her ---

void bl1_sher( uplo1_t uplo, conj1_t conj, fla_dim_t m, float*  alpha, float*    x, fla_dim_t incx, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dher( uplo1_t uplo, conj1_t conj, fla_dim_t m, double* alpha, double*   x, fla_dim_t incx, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cher( uplo1_t uplo, conj1_t conj, fla_dim_t m, float*  alpha, scomplex* x, fla_dim_t incx, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zher( uplo1_t uplo, conj1_t conj, fla_dim_t m, double* alpha, dcomplex* x, fla_dim_t incx, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

void bl1_cher_blas( uplo1_t uplo, fla_dim_t m, float*  alpha, scomplex* x, fla_dim_t incx, scomplex* a, fla_dim_t lda );
void bl1_zher_blas( uplo1_t uplo, fla_dim_t m, double* alpha, dcomplex* x, fla_dim_t incx, dcomplex* a, fla_dim_t lda );

// --- her2 ---

void bl1_sher2( uplo1_t uplo, conj1_t conj, fla_dim_t m, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dher2( uplo1_t uplo, conj1_t conj, fla_dim_t m, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_cher2( uplo1_t uplo, conj1_t conj, fla_dim_t m, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zher2( uplo1_t uplo, conj1_t conj, fla_dim_t m, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

void bl1_cher2_blas( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* a, fla_dim_t lda );
void bl1_zher2_blas( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* a, fla_dim_t lda );

// --- symv ---

void bl1_ssymv( uplo1_t uplo, fla_dim_t m, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    x, fla_dim_t incx, float*    beta, float*    y, fla_dim_t incy );
void bl1_dsymv( uplo1_t uplo, fla_dim_t m, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   x, fla_dim_t incx, double*   beta, double*   y, fla_dim_t incy );
void bl1_csymv( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_zsymv( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

void bl1_ssymv_blas( uplo1_t uplo, fla_dim_t m, float*    alpha, float*    a, fla_dim_t lda, float*    x, fla_dim_t incx, float*    beta, float*    y, fla_dim_t incy );
void bl1_dsymv_blas( uplo1_t uplo, fla_dim_t m, double*   alpha, double*   a, fla_dim_t lda, double*   x, fla_dim_t incx, double*   beta, double*   y, fla_dim_t incy );
void bl1_csymv_blas( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_zsymv_blas( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

// --- syr ---

void bl1_ssyr( uplo1_t uplo, fla_dim_t m, float*    alpha, float*    x, fla_dim_t incx, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dsyr( uplo1_t uplo, fla_dim_t m, double*   alpha, double*   x, fla_dim_t incx, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csyr( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zsyr( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

void bl1_ssyr_blas( uplo1_t uplo, fla_dim_t m, float*    alpha, float*    x, fla_dim_t incx, float*    a, fla_dim_t lda );
void bl1_dsyr_blas( uplo1_t uplo, fla_dim_t m, double*   alpha, double*   x, fla_dim_t incx, double*   a, fla_dim_t lda );
void bl1_csyr_blas( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* a, fla_dim_t lda );
void bl1_zsyr_blas( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* a, fla_dim_t lda );

// --- syr2 ---

void bl1_ssyr2( uplo1_t uplo, fla_dim_t m, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_dsyr2( uplo1_t uplo, fla_dim_t m, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_csyr2( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );
void bl1_zsyr2( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs );

void bl1_ssyr2_blas( uplo1_t uplo, fla_dim_t m, float*    alpha, float*    x, fla_dim_t incx, float*    y, fla_dim_t incy, float*    a, fla_dim_t lda );
void bl1_dsyr2_blas( uplo1_t uplo, fla_dim_t m, double*   alpha, double*   x, fla_dim_t incx, double*   y, fla_dim_t incy, double*   a, fla_dim_t lda );
void bl1_csyr2_blas( uplo1_t uplo, fla_dim_t m, scomplex* alpha, scomplex* x, fla_dim_t incx, scomplex* y, fla_dim_t incy, scomplex* a, fla_dim_t lda );
void bl1_zsyr2_blas( uplo1_t uplo, fla_dim_t m, dcomplex* alpha, dcomplex* x, fla_dim_t incx, dcomplex* y, fla_dim_t incy, dcomplex* a, fla_dim_t lda );

// --- trmv ---

void bl1_strmv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    x, fla_dim_t incx );
void bl1_dtrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   x, fla_dim_t incx );
void bl1_ctrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx );
void bl1_ztrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx );

void bl1_strmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, float*    a, fla_dim_t lda, float*    x, fla_dim_t incx );
void bl1_dtrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, double*   a, fla_dim_t lda, double*   x, fla_dim_t incx );
void bl1_ctrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, scomplex* a, fla_dim_t lda, scomplex* x, fla_dim_t incx );
void bl1_ztrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, dcomplex* a, fla_dim_t lda, dcomplex* x, fla_dim_t incx );

// --- trsv ---

void bl1_strsv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    x, fla_dim_t incx );
void bl1_dtrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   x, fla_dim_t incx );
void bl1_ctrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx );
void bl1_ztrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx );

void bl1_strsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, float*    a, fla_dim_t lda, float*    x, fla_dim_t incx );
void bl1_dtrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, double*   a, fla_dim_t lda, double*   x, fla_dim_t incx );
void bl1_ctrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, scomplex* a, fla_dim_t lda, scomplex* x, fla_dim_t incx );
void bl1_ztrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, dcomplex* a, fla_dim_t lda, dcomplex* x, fla_dim_t incx );

// --- trmvsx ---

void bl1_strmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, float* alpha, float* a, fla_dim_t a_rs, fla_dim_t a_cs, float* x, fla_dim_t incx, float* beta, float* y, fla_dim_t incy );
void bl1_dtrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, double* alpha, double* a, fla_dim_t a_rs, fla_dim_t a_cs, double* x, fla_dim_t incx, double* beta, double* y, fla_dim_t incy );
void bl1_ctrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_ztrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

// --- trsvsx ---

void bl1_strsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, float* alpha, float* a, fla_dim_t a_rs, fla_dim_t a_cs, float* x, fla_dim_t incx, float* beta, float* y, fla_dim_t incy );
void bl1_dtrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, double* alpha, double* a, fla_dim_t a_rs, fla_dim_t a_cs, double* x, fla_dim_t incx, double* beta, double* y, fla_dim_t incy );
void bl1_ctrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* x, fla_dim_t incx, scomplex* beta, scomplex* y, fla_dim_t incy );
void bl1_ztrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* x, fla_dim_t incx, dcomplex* beta, dcomplex* y, fla_dim_t incy );

