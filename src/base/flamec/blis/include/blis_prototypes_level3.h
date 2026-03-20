/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-3 BLAS-like prototypes --------------------------------------------

// --- gemm ---

void bl1_sgemm( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t k, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs, float*    beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dgemm( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t k, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs, double*   beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_cgemm( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t k, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, scomplex* beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_zgemm( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t k, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, dcomplex* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

void bl1_sgemm_blas( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t n, fla_dim_t k, float*    alpha, float*    a, fla_dim_t lda, float*    b, fla_dim_t ldb, float*    beta, float*    c, fla_dim_t ldc );
void bl1_dgemm_blas( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t n, fla_dim_t k, double*   alpha, double*   a, fla_dim_t lda, double*   b, fla_dim_t ldb, double*   beta, double*   c, fla_dim_t ldc );
void bl1_cgemm_blas( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t n, fla_dim_t k, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* b, fla_dim_t ldb, scomplex* beta, scomplex* c, fla_dim_t ldc );
void bl1_zgemm_blas( trans1_t transa, trans1_t transb, fla_dim_t m, fla_dim_t n, fla_dim_t k, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* b, fla_dim_t ldb, dcomplex* beta, dcomplex* c, fla_dim_t ldc );

// --- hemm ---

void bl1_shemm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs, float*    beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dhemm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs, double*   beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_chemm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, scomplex* beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_zhemm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, dcomplex* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

void bl1_chemm_blas( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* b, fla_dim_t ldb, scomplex* beta, scomplex* c, fla_dim_t ldc );
void bl1_zhemm_blas( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* b, fla_dim_t ldb, dcomplex* beta, dcomplex* c, fla_dim_t ldc );

// --- herk ---

void bl1_sherk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*  alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*  beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dherk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double* alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double* beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_cherk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*  alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, float*  beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_zherk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, double* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

void bl1_cherk_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*  alpha, scomplex* a, fla_dim_t lda, float*  beta, scomplex* c, fla_dim_t ldc );
void bl1_zherk_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double* alpha, dcomplex* a, fla_dim_t lda, double* beta, dcomplex* c, fla_dim_t ldc );

// --- her2k ---

void bl1_sher2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs, float*  beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dher2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs, double* beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_cher2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, float*  beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_zher2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, double* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

void bl1_cher2k_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* b, fla_dim_t ldb, float*  beta, scomplex* c, fla_dim_t ldc );
void bl1_zher2k_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* b, fla_dim_t ldb, double* beta, dcomplex* c, fla_dim_t ldc );

// --- symm ---

void bl1_ssymm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs, float*    beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dsymm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs, double*   beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_csymm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, scomplex* beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_zsymm( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, dcomplex* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

void bl1_ssymm_blas( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t lda, float*    b, fla_dim_t ldb, float*    beta, float*    c, fla_dim_t ldc );
void bl1_dsymm_blas( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t lda, double*   b, fla_dim_t ldb, double*   beta, double*   c, fla_dim_t ldc );
void bl1_csymm_blas( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* b, fla_dim_t ldb, scomplex* beta, scomplex* c, fla_dim_t ldc );
void bl1_zsymm_blas( side1_t side, uplo1_t uplo, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* b, fla_dim_t ldb, dcomplex* beta, dcomplex* c, fla_dim_t ldc );

// --- syrk ---

void bl1_ssyrk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dsyrk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_csyrk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_zsyrk( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

void bl1_ssyrk_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*    alpha, float*    a, fla_dim_t lda, float*    beta, float*    c, fla_dim_t ldc );
void bl1_dsyrk_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double*   alpha, double*   a, fla_dim_t lda, double*   beta, double*   c, fla_dim_t ldc );
void bl1_csyrk_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* beta, scomplex* c, fla_dim_t ldc );
void bl1_zsyrk_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* beta, dcomplex* c, fla_dim_t ldc );

// --- syr2k ---

void bl1_ssyr2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs, float*    beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dsyr2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs, double*   beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_csyr2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, scomplex* beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_zsyr2k( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, dcomplex* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

void bl1_ssyr2k_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, float*    alpha, float*    a, fla_dim_t lda, float*    b, fla_dim_t ldb, float*    beta, float*    c, fla_dim_t ldc );
void bl1_dsyr2k_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, double*   alpha, double*   a, fla_dim_t lda, double*   b, fla_dim_t ldb, double*   beta, double*   c, fla_dim_t ldc );
void bl1_csyr2k_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* b, fla_dim_t ldb, scomplex* beta, scomplex* c, fla_dim_t ldc );
void bl1_zsyr2k_blas( uplo1_t uplo, trans1_t trans, fla_dim_t m, fla_dim_t k, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* b, fla_dim_t ldb, dcomplex* beta, dcomplex* c, fla_dim_t ldc );

// --- trmm ---

void bl1_strmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dtrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ctrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ztrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

void bl1_strmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t lda, float*    b, fla_dim_t ldb );
void bl1_dtrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t lda, double*   b, fla_dim_t ldb );
void bl1_ctrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* b, fla_dim_t ldb );
void bl1_ztrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* b, fla_dim_t ldb );

// --- trsm ---

void bl1_strsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_dtrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ctrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );
void bl1_ztrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs );

void bl1_strsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t lda, float*    b, fla_dim_t ldb );
void bl1_dtrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t lda, double*   b, fla_dim_t ldb );
void bl1_ctrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t lda, scomplex* b, fla_dim_t ldb );
void bl1_ztrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t lda, dcomplex* b, fla_dim_t ldb );

// --- trmmsx ---

void bl1_strmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs, float*    beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dtrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs, double*   beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_ctrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, scomplex* beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_ztrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, dcomplex* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

// --- trsmsx ---

void bl1_strsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, float*    alpha, float*    a, fla_dim_t a_rs, fla_dim_t a_cs, float*    b, fla_dim_t b_rs, fla_dim_t b_cs, float*    beta, float*    c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_dtrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, double*   alpha, double*   a, fla_dim_t a_rs, fla_dim_t a_cs, double*   b, fla_dim_t b_rs, fla_dim_t b_cs, double*   beta, double*   c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_ctrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, scomplex* alpha, scomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, scomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, scomplex* beta, scomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );
void bl1_ztrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, fla_dim_t m, fla_dim_t n, dcomplex* alpha, dcomplex* a, fla_dim_t a_rs, fla_dim_t a_cs, dcomplex* b, fla_dim_t b_rs, fla_dim_t b_cs, dcomplex* beta, dcomplex* c, fla_dim_t c_rs, fla_dim_t c_cs );

