#ifndef AOCL_BLAS_H
#define AOCL_BLAS_H

#include "FLA_type_defs.h"

void aocl_blas_caxpy(const aocl_int64_t *n, const scomplex *ca, const scomplex *cx,
                     const aocl_int64_t *incx, scomplex *cy, const aocl_int64_t *incy);
void aocl_blas_ccopy(const aocl_int64_t *n, const scomplex *cx, const aocl_int64_t *incx,
                     scomplex *cy, const aocl_int64_t *incy);
void aocl_blas_cdotc(scomplex *r, const aocl_int64_t *n, const scomplex *cx, const aocl_int64_t *incx,
                        const scomplex *cy, const aocl_int64_t *incy);
void aocl_blas_cdotu(scomplex *r, const aocl_int64_t *n, const scomplex *cx, const aocl_int64_t *incx,
                        const scomplex *cy, const aocl_int64_t *incy);
void aocl_blas_cgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *x,
                     const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy);
void aocl_blas_cgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *b,
                     const aocl_int64_t *ldb, const scomplex *beta, scomplex *c,
                     const aocl_int64_t *ldc);
void aocl_blas_cgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const scomplex *alpha,
                       const scomplex *a, const aocl_int64_t *lda, const scomplex *b,
                       const aocl_int64_t *ldb, const scomplex *beta, scomplex *c,
                       const aocl_int64_t *ldc);
void aocl_blas_cgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const scomplex *alpha, const scomplex *a, const aocl_int64_t *lda,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy);
void aocl_blas_cgerc(const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *a, const aocl_int64_t *lda);
void aocl_blas_cgeru(const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *a, const aocl_int64_t *lda);
void aocl_blas_chbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const scomplex *alpha, const scomplex *a, const aocl_int64_t *lda,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy);
void aocl_blas_chemm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const scomplex *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                     const scomplex *beta, scomplex *c, const aocl_int64_t *ldc);
void aocl_blas_chemv(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *x,
                     const aocl_int64_t *incx, const scomplex *beta, scomplex *y,
                     const aocl_int64_t *incy);
void aocl_blas_cher(const char *uplo, const aocl_int64_t *n, const real *alpha, const scomplex *x,
                    const aocl_int64_t *incx, scomplex *a, const aocl_int64_t *lda);
void aocl_blas_cher2(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *a, const aocl_int64_t *lda);
void aocl_blas_cher2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const scomplex *alpha, const scomplex *a,
                      const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                      const real *beta, scomplex *c, const aocl_int64_t *ldc);
void aocl_blas_cherk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const real *beta, scomplex *c,
                     const aocl_int64_t *ldc);
void aocl_blas_chpmv(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *ap, const scomplex *x, const aocl_int64_t *incx,
                     const scomplex *beta, scomplex *y, const aocl_int64_t *incy);
void aocl_blas_chpr(const char *uplo, const aocl_int64_t *n, const real *alpha, const scomplex *x,
                    const aocl_int64_t *incx, scomplex *ap);
void aocl_blas_chpr2(const char *uplo, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *x, const aocl_int64_t *incx, const scomplex *y,
                     const aocl_int64_t *incy, scomplex *ap);
void aocl_blas_cscal(const aocl_int64_t *n, const scomplex *ca, scomplex *cx,
                     const aocl_int64_t *incx);
void aocl_blas_csrot(const aocl_int64_t *n, scomplex *cx, const aocl_int64_t *incx, scomplex *cy,
                     const aocl_int64_t *incy, const real *c, const real *s);
void aocl_blas_csscal(const aocl_int64_t *n, const real *sa, scomplex *cx, const aocl_int64_t *incx);
void aocl_blas_cswap(const aocl_int64_t *n, scomplex *cx, const aocl_int64_t *incx, scomplex *cy,
                     const aocl_int64_t *incy);
void aocl_blas_csymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const scomplex *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                     const scomplex *beta, scomplex *c, const aocl_int64_t *ldc);
void aocl_blas_csyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const scomplex *alpha, const scomplex *a,
                      const aocl_int64_t *lda, const scomplex *b, const aocl_int64_t *ldb,
                      const scomplex *beta, scomplex *c, const aocl_int64_t *ldc);
void aocl_blas_csyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const scomplex *alpha, const scomplex *a,
                     const aocl_int64_t *lda, const scomplex *beta, scomplex *c,
                     const aocl_int64_t *ldc);
void aocl_blas_ctbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx);
void aocl_blas_ctbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx);
void aocl_blas_ctpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *ap, scomplex *x, const aocl_int64_t *incx);
void aocl_blas_ctpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *ap, scomplex *x, const aocl_int64_t *incx);
void aocl_blas_ctrmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *b,
                     const aocl_int64_t *ldb);
void aocl_blas_ctrmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx);
void aocl_blas_ctrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *b,
                     const aocl_int64_t *ldb);
void aocl_blas_ctrsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const scomplex *a, const aocl_int64_t *lda, scomplex *x,
                     const aocl_int64_t *incx);
doublereal aocl_blas_dasum(const aocl_int64_t *n, const doublereal *dx, const aocl_int64_t *incx);
void aocl_blas_daxpy(const aocl_int64_t *n, const doublereal *da, const doublereal *dx,
                     const aocl_int64_t *incx, doublereal *dy, const aocl_int64_t *incy);
void aocl_blas_dcopy(const aocl_int64_t *n, const doublereal *dx, const aocl_int64_t *incx,
                     doublereal *dy, const aocl_int64_t *incy);
doublereal aocl_blas_ddot(const aocl_int64_t *n, const doublereal *dx, const aocl_int64_t *incx,
                          const doublereal *dy, const aocl_int64_t *incy);
void aocl_blas_dgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *x,
                     const aocl_int64_t *incx, const doublereal *beta, doublereal *y,
                     const aocl_int64_t *incy);
void aocl_blas_dgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *b,
                     const aocl_int64_t *ldb, const doublereal *beta, doublereal *c,
                     const aocl_int64_t *ldc);
void aocl_blas_dgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const doublereal *alpha,
                       const doublereal *a, const aocl_int64_t *lda, const doublereal *b,
                       const aocl_int64_t *ldb, const doublereal *beta, doublereal *c,
                       const aocl_int64_t *ldc);
void aocl_blas_dgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const doublereal *alpha, const doublereal *a, const aocl_int64_t *lda,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *beta,
                     doublereal *y, const aocl_int64_t *incy);
void aocl_blas_dger(const aocl_int64_t *m, const aocl_int64_t *n, const doublereal *alpha,
                    const doublereal *x, const aocl_int64_t *incx, const doublereal *y,
                    const aocl_int64_t *incy, doublereal *a, const aocl_int64_t *lda);
doublereal aocl_blas_dnrm2(const aocl_int64_t *n, const doublereal *x, const aocl_int64_t *incx);
void aocl_blas_drot(const aocl_int64_t *n, doublereal *dx, const aocl_int64_t *incx, doublereal *dy,
                    const aocl_int64_t *incy, const doublereal *c, const doublereal *s);
void aocl_blas_drotm(const aocl_int64_t *n, doublereal *dx, const aocl_int64_t *incx,
                     doublereal *dy, const aocl_int64_t *incy, const doublereal *dparam);
void aocl_blas_dsbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const doublereal *alpha, const doublereal *a, const aocl_int64_t *lda,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *beta,
                     doublereal *y, const aocl_int64_t *incy);
void aocl_blas_dscal(const aocl_int64_t *n, const doublereal *da, doublereal *dx,
                     const aocl_int64_t *incx);
doublereal aocl_blas_dsdot(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx,
                           const real *sy, const aocl_int64_t *incy);
void aocl_blas_dspmv(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *ap, const doublereal *x, const aocl_int64_t *incx,
                     const doublereal *beta, doublereal *y, const aocl_int64_t *incy);
void aocl_blas_dspr(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const doublereal *x, const aocl_int64_t *incx, doublereal *ap);
void aocl_blas_dspr2(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *y,
                     const aocl_int64_t *incy, doublereal *ap);
void aocl_blas_dswap(const aocl_int64_t *n, doublereal *dx, const aocl_int64_t *incx,
                     doublereal *dy, const aocl_int64_t *incy);
void aocl_blas_dsymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const doublereal *alpha, const doublereal *a,
                     const aocl_int64_t *lda, const doublereal *b, const aocl_int64_t *ldb,
                     const doublereal *beta, doublereal *c, const aocl_int64_t *ldc);
void aocl_blas_dsymv(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *x,
                     const aocl_int64_t *incx, const doublereal *beta, doublereal *y,
                     const aocl_int64_t *incy);
void aocl_blas_dsyr(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const doublereal *x, const aocl_int64_t *incx, doublereal *a,
                    const aocl_int64_t *lda);
void aocl_blas_dsyr2(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *x, const aocl_int64_t *incx, const doublereal *y,
                     const aocl_int64_t *incy, doublereal *a, const aocl_int64_t *lda);
void aocl_blas_dsyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const doublereal *alpha, const doublereal *a,
                      const aocl_int64_t *lda, const doublereal *b, const aocl_int64_t *ldb,
                      const doublereal *beta, doublereal *c, const aocl_int64_t *ldc);
void aocl_blas_dsyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *alpha, const doublereal *a,
                     const aocl_int64_t *lda, const doublereal *beta, doublereal *c,
                     const aocl_int64_t *ldc);
void aocl_blas_dtbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *a, const aocl_int64_t *lda,
                     doublereal *x, const aocl_int64_t *incx);
void aocl_blas_dtbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *a, const aocl_int64_t *lda,
                     doublereal *x, const aocl_int64_t *incx);
void aocl_blas_dtpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *ap, doublereal *x, const aocl_int64_t *incx);
void aocl_blas_dtpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *ap, doublereal *x, const aocl_int64_t *incx);
void aocl_blas_dtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *b,
                     const aocl_int64_t *ldb);
void aocl_blas_dtrmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *x,
                     const aocl_int64_t *incx);
void aocl_blas_dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *b,
                     const aocl_int64_t *ldb);
void aocl_blas_dtrsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const doublereal *a, const aocl_int64_t *lda, doublereal *x,
                     const aocl_int64_t *incx);
doublereal aocl_blas_dzasum(const aocl_int64_t *n, dcomplex *zx, const aocl_int64_t *incx);
doublereal aocl_blas_dznrm2(const aocl_int64_t *n, const dcomplex *x,
                            const aocl_int64_t *incx);
aocl_int64_t aocl_blas_icamax(const aocl_int64_t *n, const scomplex *x, const aocl_int64_t *incx);
aocl_int64_t aocl_blas_idamax(const aocl_int64_t *n, const doublereal *dx,
                              const aocl_int64_t *incx);
aocl_int64_t aocl_blas_isamax(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx);
aocl_int64_t aocl_blas_izamax(const aocl_int64_t *n, const dcomplex *x,
                              const aocl_int64_t *incx);
real aocl_blas_sasum(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx);
void aocl_blas_saxpy(const aocl_int64_t *n, const real *sa, const real *sx,
                     const aocl_int64_t *incx, real *sy, const aocl_int64_t *incy);
real aocl_blas_scasum(const aocl_int64_t *n, scomplex *cx, const aocl_int64_t *incx);
real aocl_blas_scnrm2(const aocl_int64_t *n, const scomplex *x, const aocl_int64_t *incx);
void aocl_blas_scopy(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx, real *sy,
                     const aocl_int64_t *incy);
real aocl_blas_sdot(const aocl_int64_t *n, const real *sx, const aocl_int64_t *incx, const real *sy,
                    const aocl_int64_t *incy);
real aocl_blas_sdsdot(const aocl_int64_t *n, const real *sb, const real *sx,
                      const aocl_int64_t *incx, const real *sy, const aocl_int64_t *incy);
void aocl_blas_sgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const real *alpha,
                     const real *a, const aocl_int64_t *lda, const real *x,
                     const aocl_int64_t *incx, const real *beta, real *y, const aocl_int64_t *incy);
void aocl_blas_sgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *b, const aocl_int64_t *ldb,
                     const real *beta, real *c, const aocl_int64_t *ldc);
void aocl_blas_sgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const real *alpha,
                       const real *a, const aocl_int64_t *lda, const real *b,
                       const aocl_int64_t *ldb, const real *beta, real *c, const aocl_int64_t *ldc);
void aocl_blas_sgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const real *alpha, const real *a, const aocl_int64_t *lda, const real *x,
                     const aocl_int64_t *incx, const real *beta, real *y, const aocl_int64_t *incy);
void aocl_blas_sger(const aocl_int64_t *m, const aocl_int64_t *n, const real *alpha, const real *x,
                    const aocl_int64_t *incx, const real *y, const aocl_int64_t *incy, real *a,
                    const aocl_int64_t *lda);
real aocl_blas_snrm2(const aocl_int64_t *n, const real *x, const aocl_int64_t *incx);
void aocl_blas_srot(const aocl_int64_t *n, real *sx, const aocl_int64_t *incx, real *sy,
                    const aocl_int64_t *incy, const real *c, const real *s);
void aocl_blas_srotm(const aocl_int64_t *n, real *sx, const aocl_int64_t *incx, real *sy,
                     const aocl_int64_t *incy, const real *sparam);
void aocl_blas_ssbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const real *alpha, const real *a, const aocl_int64_t *lda, const real *x,
                     const aocl_int64_t *incx, const real *beta, real *y, const aocl_int64_t *incy);
void aocl_blas_sscal(const aocl_int64_t *n, const real *sa, real *sx, const aocl_int64_t *incx);
void aocl_blas_sspmv(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *ap,
                     const real *x, const aocl_int64_t *incx, const real *beta, real *y,
                     const aocl_int64_t *incy);
void aocl_blas_sspr(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                    const aocl_int64_t *incx, real *ap);
void aocl_blas_sspr2(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                     const aocl_int64_t *incx, const real *y, const aocl_int64_t *incy, real *ap);
void aocl_blas_sswap(const aocl_int64_t *n, real *sx, const aocl_int64_t *incx, real *sy,
                     const aocl_int64_t *incy);
void aocl_blas_ssymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *b, const aocl_int64_t *ldb,
                     const real *beta, real *c, const aocl_int64_t *ldc);
void aocl_blas_ssymv(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *x, const aocl_int64_t *incx,
                     const real *beta, real *y, const aocl_int64_t *incy);
void aocl_blas_ssyr(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                    const aocl_int64_t *incx, real *a, const aocl_int64_t *lda);
void aocl_blas_ssyr2(const char *uplo, const aocl_int64_t *n, const real *alpha, const real *x,
                     const aocl_int64_t *incx, const real *y, const aocl_int64_t *incy, real *a,
                     const aocl_int64_t *lda);
void aocl_blas_ssyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const real *alpha, const real *a,
                      const aocl_int64_t *lda, const real *b, const aocl_int64_t *ldb,
                      const real *beta, real *c, const aocl_int64_t *ldc);
void aocl_blas_ssyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *alpha, const real *a,
                     const aocl_int64_t *lda, const real *beta, real *c, const aocl_int64_t *ldc);
void aocl_blas_stbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *a, const aocl_int64_t *lda, real *x,
                     const aocl_int64_t *incx);
void aocl_blas_stbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const real *a, const aocl_int64_t *lda, real *x,
                     const aocl_int64_t *incx);
void aocl_blas_stpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *ap, real *x, const aocl_int64_t *incx);
void aocl_blas_stpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *ap, real *x, const aocl_int64_t *incx);
void aocl_blas_strmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, real *b, const aocl_int64_t *ldb);
void aocl_blas_strmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *a, const aocl_int64_t *lda, real *x, const aocl_int64_t *incx);
void aocl_blas_strsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const real *alpha, const real *a,
                     const aocl_int64_t *lda, real *b, const aocl_int64_t *ldb);
void aocl_blas_strsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const real *a, const aocl_int64_t *lda, real *x, const aocl_int64_t *incx);
void aocl_blas_xerbla(const char *srname, const aocl_int64_t *info, aocl_int64_t srname_len);
void aocl_blas_xerbla_array(const char *srname_array, const aocl_int64_t *srname_len,
                            const aocl_int64_t *info);
void aocl_blas_zaxpy(const aocl_int64_t *n, const dcomplex *za, const dcomplex *zx,
                     const aocl_int64_t *incx, dcomplex *zy, const aocl_int64_t *incy);
void aocl_blas_zcopy(const aocl_int64_t *n, const dcomplex *zx, const aocl_int64_t *incx,
                     dcomplex *zy, const aocl_int64_t *incy);
void aocl_blas_zdotc(dcomplex *r, const aocl_int64_t *n, const dcomplex *zx,
                              const aocl_int64_t *incx, const dcomplex *zy,
                              const aocl_int64_t *incy);
void aocl_blas_zdotu(dcomplex *r, const aocl_int64_t *n, const dcomplex *zx,
                              const aocl_int64_t *incx, const dcomplex *zy,
                              const aocl_int64_t *incy);
void aocl_blas_zdrot(const aocl_int64_t *n, dcomplex *zx, const aocl_int64_t *incx,
                     dcomplex *zy, const aocl_int64_t *incy, const doublereal *c,
                     const doublereal *s);
void aocl_blas_zdscal(const aocl_int64_t *n, const doublereal *da, dcomplex *zx,
                      const aocl_int64_t *incx);
void aocl_blas_zgbmv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const aocl_int64_t *kl, const aocl_int64_t *ku, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *x,
                     const aocl_int64_t *incx, const dcomplex *beta, dcomplex *y,
                     const aocl_int64_t *incy);
void aocl_blas_zgemm(const char *transa, const char *transb, const aocl_int64_t *m,
                     const aocl_int64_t *n, const aocl_int64_t *k, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *b,
                     const aocl_int64_t *ldb, const dcomplex *beta, dcomplex *c,
                     const aocl_int64_t *ldc);
void aocl_blas_zgemmtr(const char *uplo, const char *transa, const char *transb,
                       const aocl_int64_t *n, const aocl_int64_t *k, const dcomplex *alpha,
                       const dcomplex *a, const aocl_int64_t *lda, const dcomplex *b,
                       const aocl_int64_t *ldb, const dcomplex *beta, dcomplex *c,
                       const aocl_int64_t *ldc);
void aocl_blas_zgemv(const char *trans, const aocl_int64_t *m, const aocl_int64_t *n,
                     const dcomplex *alpha, const dcomplex *a, const aocl_int64_t *lda,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *beta,
                     dcomplex *y, const aocl_int64_t *incy);
void aocl_blas_zgerc(const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *a, const aocl_int64_t *lda);
void aocl_blas_zgeru(const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *a, const aocl_int64_t *lda);
void aocl_blas_zhbmv(const char *uplo, const aocl_int64_t *n, const aocl_int64_t *k,
                     const dcomplex *alpha, const dcomplex *a, const aocl_int64_t *lda,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *beta,
                     dcomplex *y, const aocl_int64_t *incy);
void aocl_blas_zhemm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const dcomplex *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                     const dcomplex *beta, dcomplex *c, const aocl_int64_t *ldc);
void aocl_blas_zhemv(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *x,
                     const aocl_int64_t *incx, const dcomplex *beta, dcomplex *y,
                     const aocl_int64_t *incy);
void aocl_blas_zher(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const dcomplex *x, const aocl_int64_t *incx, dcomplex *a,
                    const aocl_int64_t *lda);
void aocl_blas_zher2(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *a, const aocl_int64_t *lda);
void aocl_blas_zher2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const dcomplex *alpha, const dcomplex *a,
                      const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                      const doublereal *beta, dcomplex *c, const aocl_int64_t *ldc);
void aocl_blas_zherk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const doublereal *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const doublereal *beta, dcomplex *c,
                     const aocl_int64_t *ldc);
void aocl_blas_zhpmv(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *ap, const dcomplex *x, const aocl_int64_t *incx,
                     const dcomplex *beta, dcomplex *y, const aocl_int64_t *incy);
void aocl_blas_zhpr(const char *uplo, const aocl_int64_t *n, const doublereal *alpha,
                    const dcomplex *x, const aocl_int64_t *incx, dcomplex *ap);
void aocl_blas_zhpr2(const char *uplo, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *x, const aocl_int64_t *incx, const dcomplex *y,
                     const aocl_int64_t *incy, dcomplex *ap);
void aocl_blas_zscal(const aocl_int64_t *n, const dcomplex *za, dcomplex *zx,
                     const aocl_int64_t *incx);
void aocl_blas_zswap(const aocl_int64_t *n, dcomplex *zx, const aocl_int64_t *incx,
                     dcomplex *zy, const aocl_int64_t *incy);
void aocl_blas_zsymm(const char *side, const char *uplo, const aocl_int64_t *m,
                     const aocl_int64_t *n, const dcomplex *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                     const dcomplex *beta, dcomplex *c, const aocl_int64_t *ldc);
void aocl_blas_zsyr2k(const char *uplo, const char *trans, const aocl_int64_t *n,
                      const aocl_int64_t *k, const dcomplex *alpha, const dcomplex *a,
                      const aocl_int64_t *lda, const dcomplex *b, const aocl_int64_t *ldb,
                      const dcomplex *beta, dcomplex *c, const aocl_int64_t *ldc);
void aocl_blas_zsyrk(const char *uplo, const char *trans, const aocl_int64_t *n,
                     const aocl_int64_t *k, const dcomplex *alpha, const dcomplex *a,
                     const aocl_int64_t *lda, const dcomplex *beta, dcomplex *c,
                     const aocl_int64_t *ldc);
void aocl_blas_ztbmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const dcomplex *a, const aocl_int64_t *lda,
                     dcomplex *x, const aocl_int64_t *incx);
void aocl_blas_ztbsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const aocl_int64_t *k, const dcomplex *a, const aocl_int64_t *lda,
                     dcomplex *x, const aocl_int64_t *incx);
void aocl_blas_ztpmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *ap, dcomplex *x, const aocl_int64_t *incx);
void aocl_blas_ztpsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *ap, dcomplex *x, const aocl_int64_t *incx);
void aocl_blas_ztrmm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *b,
                     const aocl_int64_t *ldb);
void aocl_blas_ztrmv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *x,
                     const aocl_int64_t *incx);
void aocl_blas_ztrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                     const aocl_int64_t *m, const aocl_int64_t *n, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *b,
                     const aocl_int64_t *ldb);
void aocl_blas_ztrsv(const char *uplo, const char *trans, const char *diag, const aocl_int64_t *n,
                     const dcomplex *a, const aocl_int64_t *lda, dcomplex *x,
                     const aocl_int64_t *incx);
#ifdef FLA_ENABLE_BLAS_EXT_GEMMT
void aocl_blas_cgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const scomplex *alpha,
                     const scomplex *a, const aocl_int64_t *lda, const scomplex *b,
                     const aocl_int64_t *ldb, const scomplex *beta, scomplex *c,
                     const aocl_int64_t *ldc);
void aocl_blas_dgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const doublereal *alpha,
                     const doublereal *a, const aocl_int64_t *lda, const doublereal *b,
                     const aocl_int64_t *ldb, const doublereal *beta, doublereal *c,
                     const aocl_int64_t *ldc);
void aocl_blas_sgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const real *alpha,
                     const real *a, const aocl_int64_t *lda, const real *b,
                     const aocl_int64_t *ldb, const real *beta, real *c,
                     const aocl_int64_t *ldc);
void aocl_blas_zgemmt(const char* uplo, const char *transa, const char *transb,
                     const aocl_int64_t *n, const aocl_int64_t *k, const dcomplex *alpha,
                     const dcomplex *a, const aocl_int64_t *lda, const dcomplex *b,
                     const aocl_int64_t *ldb, const dcomplex *beta, dcomplex *c,
                     const aocl_int64_t *ldc);                                                              
#endif /* FLA_ENABLE_BLAS_EXT_GEMMT */

#endif /* AOCL_BLAS_H */