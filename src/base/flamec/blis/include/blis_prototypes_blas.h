/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/**
 * Modifications Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 */

// --- Name-mangling macro definitions -----------------------------------------

// --- Name-mangle level-1 BLAS routines ---------------------------

// Allow C++ users to include this header file in their source code. However,
// we make the extern "C" conditional on whether we're using a C++ compiler,
// since regular C compilers don't understand the extern "C" construct.
#ifdef __cplusplus
extern "C" {
#endif

// #define F77_isamax F77_FUNC( isamax , ISAMAX )
#define F77_isamax aocl_blas_isamax
#define F77_idamax aocl_blas_idamax
#define F77_icamax aocl_blas_icamax
#define F77_izamax aocl_blas_izamax
#define F77_sasum  aocl_blas_sasum
#define F77_dasum  aocl_blas_dasum
#define F77_scasum aocl_blas_scasum
#define F77_dzasum aocl_blas_dzasum
#define F77_saxpy  aocl_blas_saxpy
#define F77_daxpy  aocl_blas_daxpy
#define F77_caxpy  aocl_blas_caxpy
#define F77_zaxpy  aocl_blas_zaxpy
#define F77_scopy  aocl_blas_scopy
#define F77_dcopy  aocl_blas_dcopy
#define F77_ccopy  aocl_blas_ccopy
#define F77_zcopy  aocl_blas_zcopy
#define F77_sdot   aocl_blas_sdot
#define F77_ddot   aocl_blas_ddot
#define F77_cdotu  aocl_blas_cdotu
#define F77_cdotc  aocl_blas_cdotc
#define F77_zdotu  aocl_blas_zdotu
#define F77_zdotc  aocl_blas_zdotc
#define F77_snrm2  aocl_blas_snrm2
#define F77_dnrm2  aocl_blas_dnrm2
#define F77_scnrm2 aocl_blas_scnrm2
#define F77_dznrm2 aocl_blas_dznrm2
#define F77_sscal  aocl_blas_sscal
#define F77_dscal  aocl_blas_dscal
#define F77_cscal  aocl_blas_cscal
#define F77_csscal aocl_blas_csscal
#define F77_zscal  aocl_blas_zscal
#define F77_zdscal aocl_blas_zdscal
#define F77_sswap  aocl_blas_sswap
#define F77_dswap  aocl_blas_dswap
#define F77_cswap  aocl_blas_cswap
#define F77_zswap  aocl_blas_zswap  

// --- Name-mangle level-2 BLAS routines ---------------------------

#define F77_sgemv  aocl_blas_sgemv
#define F77_dgemv  aocl_blas_dgemv
#define F77_cgemv  aocl_blas_cgemv
#define F77_zgemv  aocl_blas_zgemv
#define F77_sger   aocl_blas_sger
#define F77_dger   aocl_blas_dger
#define F77_cgerc  aocl_blas_cgerc
#define F77_cgeru  aocl_blas_cgeru
#define F77_zgerc  aocl_blas_zgerc
#define F77_zgeru  aocl_blas_zgeru
#define F77_chemv  aocl_blas_chemv
#define F77_zhemv  aocl_blas_zhemv
#define F77_cher   aocl_blas_cher
#define F77_zher   aocl_blas_zher
#define F77_cher2  aocl_blas_cher2
#define F77_zher2  aocl_blas_zher2
#define F77_ssymv  aocl_blas_ssymv
#define F77_dsymv  aocl_blas_dsymv
#define F77_ssyr   aocl_blas_ssyr
#define F77_dsyr   aocl_blas_dsyr
#define F77_ssyr2  aocl_blas_ssyr2
#define F77_dsyr2  aocl_blas_dsyr2
#define F77_strmv  aocl_blas_strmv
#define F77_dtrmv  aocl_blas_dtrmv
#define F77_ctrmv  aocl_blas_ctrmv
#define F77_ztrmv  aocl_blas_ztrmv
#define F77_strsv  aocl_blas_strsv
#define F77_dtrsv  aocl_blas_dtrsv
#define F77_ctrsv  aocl_blas_ctrsv
#define F77_ztrsv  aocl_blas_ztrsv

// --- Name-mangle level-3 BLAS routines ---------------------------

#define F77_sgemm  aocl_blas_sgemm
#define F77_dgemm  aocl_blas_dgemm
#define F77_cgemm  aocl_blas_cgemm
#define F77_zgemm  aocl_blas_zgemm
#define F77_chemm  aocl_blas_chemm
#define F77_zhemm  aocl_blas_zhemm
#define F77_cherk  aocl_blas_cherk
#define F77_zherk  aocl_blas_zherk
#define F77_cher2k aocl_blas_cher2k
#define F77_zher2k aocl_blas_zher2k
#define F77_ssymm  aocl_blas_ssymm
#define F77_dsymm  aocl_blas_dsymm
#define F77_csymm  aocl_blas_csymm
#define F77_zsymm  aocl_blas_zsymm
#define F77_ssyrk  aocl_blas_ssyrk
#define F77_dsyrk  aocl_blas_dsyrk
#define F77_csyrk  aocl_blas_csyrk
#define F77_zsyrk  aocl_blas_zsyrk
#define F77_ssyr2k aocl_blas_ssyr2k
#define F77_dsyr2k aocl_blas_dsyr2k
#define F77_csyr2k aocl_blas_csyr2k
#define F77_zsyr2k aocl_blas_zsyr2k
#define F77_strmm  aocl_blas_strmm
#define F77_dtrmm  aocl_blas_dtrmm
#define F77_ctrmm  aocl_blas_ctrmm
#define F77_ztrmm  aocl_blas_ztrmm
#define F77_strsm  aocl_blas_strsm
#define F77_dtrsm  aocl_blas_dtrsm
#define F77_ctrsm  aocl_blas_ctrsm
#define F77_ztrsm  aocl_blas_ztrsm  



#ifdef BLIS1_FROM_LIBFLAME

#define BL1_EXPORT_isamax F77_FUNC( isamax,  ISAMAX )
#define BL1_EXPORT_idamax F77_FUNC( idamax,  IDAMAX )
#define BL1_EXPORT_icamax F77_FUNC( icamax,  ICAMAX )
#define BL1_EXPORT_izamax F77_FUNC( izamax,  IZAMAX )
#define BL1_EXPORT_sasum  F77_FUNC( sasum,  SASUM )
#define BL1_EXPORT_dasum  F77_FUNC( dasum,  DASUM )
#define BL1_EXPORT_scasum F77_FUNC( scasum,  SCASUM )
#define BL1_EXPORT_dzasum F77_FUNC( dzasum,  DZASUM )
#define BL1_EXPORT_saxpy  F77_FUNC( saxpy,  SAXPY )
#define BL1_EXPORT_daxpy  F77_FUNC( daxpy,  DAXPY )
#define BL1_EXPORT_caxpy  F77_FUNC( caxpy,  CAXPY )
#define BL1_EXPORT_zaxpy  F77_FUNC( zaxpy,  ZAXPY )
#define BL1_EXPORT_scopy  F77_FUNC( scopy,  SCOPY )
#define BL1_EXPORT_dcopy  F77_FUNC( dcopy,  DCOPY )
#define BL1_EXPORT_ccopy  F77_FUNC( ccopy,  CCOPY )
#define BL1_EXPORT_zcopy  F77_FUNC( zcopy,  ZCOPY )
#define BL1_EXPORT_sdot   F77_FUNC( sdot,  SDOT )
#define BL1_EXPORT_ddot   F77_FUNC( ddot,  DDOT )
#define BL1_EXPORT_cdotu  F77_FUNC( cdotu,  CDOTU )
#define BL1_EXPORT_cdotc  F77_FUNC( cdotc,  CDOTC )
#define BL1_EXPORT_zdotu  F77_FUNC( zdotu,  ZDOTU )
#define BL1_EXPORT_zdotc  F77_FUNC( zdotc,  ZDOTC )
#define BL1_EXPORT_snrm2  F77_FUNC( snrm2,  SNRM2 )
#define BL1_EXPORT_dnrm2  F77_FUNC( dnrm2,  DNRM2 )
#define BL1_EXPORT_scnrm2 F77_FUNC( scnrm2,  SCNRM2 )
#define BL1_EXPORT_dznrm2 F77_FUNC( dznrm2,  DZNRM2 )
#define BL1_EXPORT_sscal  F77_FUNC( sscal,  SSCAL )
#define BL1_EXPORT_dscal  F77_FUNC( dscal,  DSCAL )
#define BL1_EXPORT_cscal  F77_FUNC( cscal,  CSCAL )
#define BL1_EXPORT_csscal F77_FUNC( csscal,  CSSCAL )
#define BL1_EXPORT_zscal  F77_FUNC( zscal,  ZSCAL )
#define BL1_EXPORT_zdscal F77_FUNC( zdscal,  ZDSCAL )
#define BL1_EXPORT_sswap  F77_FUNC( sswap,  SSWAP )
#define BL1_EXPORT_dswap  F77_FUNC( dswap,  DSWAP )
#define BL1_EXPORT_cswap  F77_FUNC( cswap,  CSWAP )
#define BL1_EXPORT_zswap  F77_FUNC( zswap,  ZSWAP )  

// --- Name-mangle level-2 BLAS routines ---------------------------

#define BL1_EXPORT_sgemv  F77_FUNC( sgemv,  SGEMV )
#define BL1_EXPORT_dgemv  F77_FUNC( dgemv,  DGEMV )
#define BL1_EXPORT_cgemv  F77_FUNC( cgemv,  CGEMV )
#define BL1_EXPORT_zgemv  F77_FUNC( zgemv,  ZGEMV )
#define BL1_EXPORT_sger   F77_FUNC( sger,  SGER )
#define BL1_EXPORT_dger   F77_FUNC( dger,  DGER )
#define BL1_EXPORT_cgerc  F77_FUNC( cgerc,  CGERC )
#define BL1_EXPORT_cgeru  F77_FUNC( cgeru,  CGERU )
#define BL1_EXPORT_zgerc  F77_FUNC( zgerc,  ZGERC )
#define BL1_EXPORT_zgeru  F77_FUNC( zgeru,  ZGERU )
#define BL1_EXPORT_chemv  F77_FUNC( chemv,  CHEMV )
#define BL1_EXPORT_zhemv  F77_FUNC( zhemv,  ZHEMV )
#define BL1_EXPORT_cher   F77_FUNC( cher,  CHER )
#define BL1_EXPORT_zher   F77_FUNC( zher,  ZHER )
#define BL1_EXPORT_cher2  F77_FUNC( cher2,  CHER2 )
#define BL1_EXPORT_zher2  F77_FUNC( zher2,  ZHER2 )
#define BL1_EXPORT_ssymv  F77_FUNC( ssymv,  SSYMV )
#define BL1_EXPORT_dsymv  F77_FUNC( dsymv,  DSYMV )
#define BL1_EXPORT_ssyr   F77_FUNC( ssyr,  SSYR )
#define BL1_EXPORT_dsyr   F77_FUNC( dsyr,  DSYR )
#define BL1_EXPORT_ssyr2  F77_FUNC( ssyr2,  SSYR2 )
#define BL1_EXPORT_dsyr2  F77_FUNC( dsyr2,  DSYR2 )
#define BL1_EXPORT_strmv  F77_FUNC( strmv,  STRMV )
#define BL1_EXPORT_dtrmv  F77_FUNC( dtrmv,  DTRMV )
#define BL1_EXPORT_ctrmv  F77_FUNC( ctrmv,  CTRMV )
#define BL1_EXPORT_ztrmv  F77_FUNC( ztrmv,  ZTRMV )
#define BL1_EXPORT_strsv  F77_FUNC( strsv,  STRSV )
#define BL1_EXPORT_dtrsv  F77_FUNC( dtrsv,  DTRSV )
#define BL1_EXPORT_ctrsv  F77_FUNC( ctrsv,  CTRSV )
#define BL1_EXPORT_ztrsv  F77_FUNC( ztrsv,  ZTRSV )

// --- Name-mangle level-3 BLAS routines ---------------------------

#define BL1_EXPORT_sgemm  F77_FUNC( sgemm,  SGEMM )
#define BL1_EXPORT_dgemm  F77_FUNC( dgemm,  DGEMM )
#define BL1_EXPORT_cgemm  F77_FUNC( cgemm,  CGEMM )
#define BL1_EXPORT_zgemm  F77_FUNC( zgemm,  ZGEMM )
#define BL1_EXPORT_chemm  F77_FUNC( chemm,  CHEMM )
#define BL1_EXPORT_zhemm  F77_FUNC( zhemm,  ZHEMM )
#define BL1_EXPORT_cherk  F77_FUNC( cherk,  CHERK )
#define BL1_EXPORT_zherk  F77_FUNC( zherk,  ZHERK )
#define BL1_EXPORT_cher2k F77_FUNC( cher2k,  CHER2K )
#define BL1_EXPORT_zher2k F77_FUNC( zher2k,  ZHER2K )
#define BL1_EXPORT_ssymm  F77_FUNC( ssymm,  SSYMM )
#define BL1_EXPORT_dsymm  F77_FUNC( dsymm,  DSYMM )
#define BL1_EXPORT_csymm  F77_FUNC( csymm,  CSYMM )
#define BL1_EXPORT_zsymm  F77_FUNC( zsymm,  ZSYMM )
#define BL1_EXPORT_ssyrk  F77_FUNC( ssyrk,  SSYRK )
#define BL1_EXPORT_dsyrk  F77_FUNC( dsyrk,  DSYRK )
#define BL1_EXPORT_csyrk  F77_FUNC( csyrk,  CSYRK )
#define BL1_EXPORT_zsyrk  F77_FUNC( zsyrk,  ZSYRK )
#define BL1_EXPORT_ssyr2k F77_FUNC( ssyr2k,  SSYR2K )
#define BL1_EXPORT_dsyr2k F77_FUNC( dsyr2k,  DSYR2K )
#define BL1_EXPORT_csyr2k F77_FUNC( csyr2k,  CSYR2K )
#define BL1_EXPORT_zsyr2k F77_FUNC( zsyr2k,  ZSYR2K )
#define BL1_EXPORT_strmm  F77_FUNC( strmm,  STRMM )
#define BL1_EXPORT_dtrmm  F77_FUNC( dtrmm,  DTRMM )
#define BL1_EXPORT_ctrmm  F77_FUNC( ctrmm,  CTRMM )
#define BL1_EXPORT_ztrmm  F77_FUNC( ztrmm,  ZTRMM )
#define BL1_EXPORT_strsm  F77_FUNC( strsm,  STRSM )
#define BL1_EXPORT_dtrsm  F77_FUNC( dtrsm,  DTRSM )
#define BL1_EXPORT_ctrsm  F77_FUNC( ctrsm,  CTRSM )
#define BL1_EXPORT_ztrsm  F77_FUNC( ztrsm,  ZTRSM )  

// --- Prototypes --------------------------------------------------------------

// --- Level-1 BLAS prototypes -------------------

// --- amax ---
integer  F77_isamax ( integer* n, float*    x, integer* incx );
integer  F77_idamax ( integer* n, double*   x, integer* incx );
integer  F77_icamax ( integer* n, scomplex* x, integer* incx );
integer  F77_izamax ( integer* n, dcomplex* x, integer* incx );
// --- asum ---
float    F77_sasum  ( integer* n, float*    x, integer* incx );
double   F77_dasum  ( integer* n, double*   x, integer* incx );
float    F77_scasum ( integer* n, scomplex* x, integer* incx );
double   F77_dzasum ( integer* n, dcomplex* x, integer* incx );
// --- axpy ---
void     F77_saxpy  ( integer* n, float*    alpha, float*    x, integer* incx,  float*    y, integer* incy );
void     F77_daxpy  ( integer* n, double*   alpha, double*   x, integer* incx,  double*   y, integer* incy );
void     F77_caxpy  ( integer* n, scomplex* alpha, scomplex* x, integer* incx,  scomplex* y, integer* incy );
void     F77_zaxpy  ( integer* n, dcomplex* alpha, dcomplex* x, integer* incx,  dcomplex* y, integer* incy );
// --- copy ---
void     F77_scopy  ( integer* n, float*    x, integer* incx, float*    y, integer* incy );
void     F77_dcopy  ( integer* n, double*   x, integer* incx, double*   y, integer* incy );
void     F77_ccopy  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
void     F77_zcopy  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );
// --- dot ---
float    F77_sdot   ( integer* n, float*    x, integer* incx, float*    y, integer* incy );
double   F77_ddot   ( integer* n, double*   x, integer* incx, double*   y, integer* incy );
scomplex F77_cdotu  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
scomplex F77_cdotc  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
dcomplex F77_zdotu  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );
dcomplex F77_zdotc  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );
// --- nrm2 ---
float    F77_snrm2  ( integer* n, float*    x, integer* incx );
double   F77_dnrm2  ( integer* n, double*   x, integer* incx );
float    F77_scnrm2 ( integer* n, scomplex* x, integer* incx );
double   F77_dznrm2 ( integer* n, dcomplex* x, integer* incx );
// --- scal ---
void     F77_sscal  ( integer* n, float*    alpha, float*    y, integer* incy );
void     F77_dscal  ( integer* n, double*   alpha, double*   y, integer* incy );
void     F77_cscal  ( integer* n, scomplex* alpha, scomplex* y, integer* incy );
void     F77_csscal ( integer* n, float*    alpha, scomplex* y, integer* incy );
void     F77_zscal  ( integer* n, dcomplex* alpha, dcomplex* y, integer* incy );
void     F77_zdscal ( integer* n, double*   alpha, dcomplex* y, integer* incy );
// --- swap ---
void     F77_sswap  ( integer* n, float*    x, integer* incx, float*    y, integer* incy );
void     F77_dswap  ( integer* n, double*   x, integer* incx, double*   y, integer* incy );
void     F77_cswap  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
void     F77_zswap  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );

// --- Level-2 BLAS prototypes -------------------

// --- gemv ---
void     F77_sgemv  ( char* transa, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    x, integer* incx, float*    beta, float*    y, integer* incy );
void     F77_dgemv  ( char* transa, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   x, integer* incx, double*   beta, double*   y, integer* incy );
void     F77_cgemv  ( char* transa, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* x, integer* incx, scomplex* beta, scomplex* y, integer* incy );
void     F77_zgemv  ( char* transa, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* x, integer* incx, dcomplex* beta, dcomplex* y, integer* incy );
// --- ger ---
void     F77_sger   ( integer* m, integer* n, float*    alpha, float*    x, integer* incx, float*    y, integer* incy, float*    a, integer* lda );
void     F77_dger   ( integer* m, integer* n, double*   alpha, double*   x, integer* incx, double*   y, integer* incy, double*   a, integer* lda );
void     F77_cgerc  ( integer* m, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* y, integer* incy, scomplex* a, integer* lda );
void     F77_cgeru  ( integer* m, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* y, integer* incy, scomplex* a, integer* lda );
void     F77_zgerc  ( integer* m, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* y, integer* incy, dcomplex* a, integer* lda );
void     F77_zgeru  ( integer* m, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* y, integer* incy, dcomplex* a, integer* lda );
// --- hemv ---
void     F77_chemv  ( char* uplo, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* x, integer* incx, scomplex* beta, scomplex* y, integer* incy );
void     F77_zhemv  ( char* uplo, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* x, integer* incx, dcomplex* beta, dcomplex* y, integer* incy );
// --- her ---
void     F77_cher   ( char* uplo, integer* n, float*    alpha, scomplex* x, integer* incx, scomplex* a, integer* lda );
void     F77_zher   ( char* uplo, integer* n, double*   alpha, dcomplex* x, integer* incx, dcomplex* a, integer* lda );
// --- her2 ---
void     F77_cher2  ( char* uplo, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* y, integer* incy, scomplex* a, integer* lda );
void     F77_zher2  ( char* uplo, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* y, integer* incy, dcomplex* a, integer* lda );
// --- symv ---
void     F77_ssymv  ( char* uplo, integer* n, float*    alpha, float*    a, integer* lda, float*    x, integer* incx, float*    beta, float*    y, integer* incy );
void     F77_dsymv  ( char* uplo, integer* n, double*   alpha, double*   a, integer* lda, double*   x, integer* incx, double*   beta, double*   y, integer* incy );
// --- syr ---
void     F77_ssyr   ( char* uplo, integer* n, float*    alpha, float*    x, integer* incx, float*    a, integer* lda );
void     F77_dsyr   ( char* uplo, integer* n, double*   alpha, double*   x, integer* incx, double*   a, integer* lda );
// --- syr2 ---
void     F77_ssyr2  ( char* uplo, integer* n, float*    alpha, float*    x, integer* incx, float*    y, integer* incy, float*    a, integer* lda );
void     F77_dsyr2  ( char* uplo, integer* n, double*   alpha, double*   x, integer* incx, double*   y, integer* incy, double*   a, integer* lda );
// --- trmv ---
void     F77_strmv  ( char* uplo, char* transa, char* diag, integer* n,  float*    a, integer* lda, float*    y, integer* incy );
void     F77_dtrmv  ( char* uplo, char* transa, char* diag, integer* n,  double*   a, integer* lda, double*   y, integer* incy );
void     F77_ctrmv  ( char* uplo, char* transa, char* diag, integer* n,  scomplex* a, integer* lda, scomplex* y, integer* incy );
void     F77_ztrmv  ( char* uplo, char* transa, char* diag, integer* n,  dcomplex* a, integer* lda, dcomplex* y, integer* incy );
// --- trsv ---
void     F77_strsv  ( char* uplo, char* transa, char* diag, integer* n,  float*    a, integer* lda, float*    y, integer* incy );
void     F77_dtrsv  ( char* uplo, char* transa, char* diag, integer* n,  double*   a, integer* lda, double*   y, integer* incy );
void     F77_ctrsv  ( char* uplo, char* transa, char* diag, integer* n,  scomplex* a, integer* lda, scomplex* y, integer* incy );
void     F77_ztrsv  ( char* uplo, char* transa, char* diag, integer* n,  dcomplex* a, integer* lda, dcomplex* y, integer* incy );

// --- Level-3 BLAS prototypes -------------------

// --- gemm ---
void     F77_sgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb, float*    beta, float*    c, integer* ldc );
void     F77_dgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb, double*   beta, double*   c, integer* ldc );
void     F77_cgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- hemm ---
void     F77_chemm  ( char* side, char* uplo, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zhemm  ( char* side, char* uplo, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- herk ---
void     F77_cherk  ( char* uplo, char* transa, integer* n, integer* k, float*  alpha, scomplex* a, integer* lda, float*  beta, scomplex* c, integer* ldc );
void     F77_zherk  ( char* uplo, char* transa, integer* n, integer* k, double* alpha, dcomplex* a, integer* lda, double* beta, dcomplex* c, integer* ldc );
// --- her2k ---
void     F77_cher2k ( char* uplo, char* transa, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, float*  beta, scomplex* c, integer* ldc );
void     F77_zher2k ( char* uplo, char* transa, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* beta, dcomplex* c, integer* ldc );
// --- symm ---
void     F77_ssymm  ( char* side, char* uplo, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb, float*    beta, float*    c, integer* ldc );
void     F77_dsymm  ( char* side, char* uplo, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb, double*   beta, double*   c, integer* ldc );
void     F77_csymm  ( char* side, char* uplo, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zsymm  ( char* side, char* uplo, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- syrk ---
void     F77_ssyrk  ( char* uplo, char* transa, integer* n, integer* k, float*    alpha, float*    a, integer* lda, float*    beta, float*    c, integer* ldc );
void     F77_dsyrk  ( char* uplo, char* transa, integer* n, integer* k, double*   alpha, double*   a, integer* lda, double*   beta, double*   c, integer* ldc );
void     F77_csyrk  ( char* uplo, char* transa, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zsyrk  ( char* uplo, char* transa, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* beta, dcomplex* c, integer* ldc );
// --- syr2k ---
void     F77_ssyr2k ( char* uplo, char* transa, integer* n, integer* k, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb, float*    beta, float*    c, integer* ldc );
void     F77_dsyr2k ( char* uplo, char* transa, integer* n, integer* k, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb, double*   beta, double*   c, integer* ldc );
void     F77_csyr2k ( char* uplo, char* transa, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zsyr2k ( char* uplo, char* transa, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- trmm ---
void     F77_strmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb );
void     F77_dtrmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb );
void     F77_ctrmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb );
void     F77_ztrmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb );
// --- trsm ---
void     F77_strsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb );
void     F77_dtrsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb );
void     F77_ctrsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb );
void     F77_ztrsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb );

