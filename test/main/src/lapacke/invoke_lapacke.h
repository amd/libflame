/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#ifndef INVOKE_LAPACKE_H
#define INVOKE_LAPACKE_H

integer invoke_lapacke_orgqr(integer datatype, int layout, integer m, integer n, integer k, void *a,
                             integer lda, const void *tau);
integer invoke_lapacke_potrf(integer datatype, int layout, char uplo, integer n, void *a,
                             integer lda);
integer invoke_lapacke_potrs(integer datatype, int layout, char uplo, integer n, integer nrhs,
                             const void *A, integer lda, void *B, integer ldb);
integer invoke_lapacke_stedc(integer datatype, int layout, char compz, integer n, void *D, void *E,
                             void *Z, integer ldz);
integer invoke_lapacke_steqr(integer datatype, int layout, char compz, integer n, void *d, void *e,
                             void *z, integer ldz);
integer invoke_lapacke_stevd(integer datatype, int layout, char jobz, integer n, void *d, void *e,
                             void *z, integer ldz);
integer invoke_lapacke_syev(integer datatype, int layout, char jobz, char uplo, integer n, void *a,
                            integer lda, void *w);
integer invoke_lapacke_syevd(integer datatype, int layout, char jobz, char uplo, integer n, void *a,
                             integer lda, void *w);
integer invoke_lapacke_syevx(integer datatype, int layout, char jobz, char range, char uplo,
                             integer n, void *a, integer lda, void *vl, void *vu, integer il,
                             integer iu, void *abstol, integer *m, void *w, void *z, integer ldz,
                             integer *ifail);
integer invoke_lapacke_sygvd(integer datatype, int layout, int itype, char jobz, char uplo,
                             integer n, void *a, integer lda, void *b, integer ldb, void *w);
integer invoke_lapacke_sytrf(integer datatype, integer layout, char uplo, integer n, void *a,
                             integer lda, integer *ipiv);
integer invoke_lapacke_gbtrf(integer datatype, int layout, integer m, integer n, integer kl,
                             integer ku, void *ab, integer ldab, integer *ipiv);
integer invoke_lapacke_gbtrs(integer datatype, int layout, char trans, integer n, integer kl,
                             integer ku, integer nrhs, void *ab, integer ldab, integer *ipiv,
                             void *b, integer ldb);
integer invoke_lapacke_gecon(integer datatype, integer layout, char norm, integer n, void *A,
                             integer lda, void *anorm, void *rcond);
integer invoke_lapacke_geev(integer datatype, int layout, char jobvl, char jobvr, integer n,
                            void *a, integer lda, void *wr, void *wi, void *w, void *vl,
                            integer ldvl, void *vr, integer ldvr);
integer invoke_lapacke_geevx(integer datatype, int layout, char balanc, char jobvl, char jobvr,
                             char sense, integer n, void *a, integer lda, void *wr, void *wi,
                             void *w, void *vl, integer ldvl, void *vr, integer ldvr, integer *ilo,
                             integer *ihi, void *scale, void *abnrm, void *rconde, void *rcondv);
integer invoke_lapacke_gehrd(integer datatype, int layout, integer n, integer ilo, integer ihi,
                             void *a, integer lda, void *tau);
integer invoke_lapacke_gelqf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, void *tau);
integer invoke_lapacke_gels(integer datatype, integer layout, char trans, integer m, integer n,
                            integer nrhs, void *A, integer lda, void *B, integer ldb);
integer invoke_lapacke_gelsd(integer datatype, integer layout, integer m, integer n, integer nrhs,
                             void *A, integer lda, void *B, integer ldb, void *s, void *rcond,
                             integer *rank);
integer invoke_lapacke_gelss(integer datatype, integer layout, integer m, integer n, integer nrhs,
                             void *A, integer lda, void *B, integer ldb, void *s, void *rcond,
                             integer *rank);
integer invoke_lapacke_geqp3(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, integer *jpvt, void *tau);
integer invoke_lapacke_geqrf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, void *tau);
integer invoke_lapacke_gerqf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, void *tau);
integer invoke_lapacke_gesdd(integer datatype, int layout, char jobz, integer m, integer n, void *a,
                             integer lda, void *s, void *u, integer ldu, void *vt, integer ldvt);
integer invoke_lapacke_gesv(integer datatype, int layout, integer n, integer nrhs, void *a,
                            integer lda, integer *ipiv, void *b, integer ldb);
integer invoke_lapacke_gesvd(integer datatype, int layout, char jobu, char jobvt, integer m,
                             integer n, void *a, integer lda, void *s, void *u, integer ldu,
                             void *vt, integer ldvt, void *work, void *rwork);
integer invoke_lapacke_gesvdx(integer datatype, int layout, char jobu, char jobvt, char range,
                              integer m, integer n, void *a, integer lda, void *vl, void *vu,
                              integer il, integer iu, integer *ns, void *s, void *u, integer ldu,
                              void *vt, integer ldvt, void *superb);
integer invoke_lapacke_getrf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, integer *ipiv);
integer invoke_lapacke_getri(integer datatype, int layout, integer n, void *a, integer lda,
                             const integer *ipiv);
integer invoke_lapacke_getrs(integer datatype, int layout, char trans, integer n, integer nrhs,
                             const void *a, integer lda, const integer *ipiv, void *b, integer ldb);
integer invoke_lapacke_ggev(integer datatype, int layout, char jobvl, char jobvr, integer n,
                            void *a, integer lda, void *b, integer ldb, void *alpha, void *alphar,
                            void *alphai, void *beta, void *vl, integer ldvl, void *vr,
                            integer ldvr);
integer invoke_lapacke_ggevx(integer datatype, int layout, char balanc, char jobvl, char jobvr,
                             char sense, integer n, void *a, integer lda, void *b, integer ldb,
                             void *alpha, void *alphar, void *alphai, void *beta, void *vl,
                             integer ldvl, void *vr, integer ldvr, integer *ilo, integer *ihi,
                             void *lscale, void *rscale, void *abnrm, void *bbnrm, void *rconde,
                             void *rcondv);
integer invoke_lapacke_gghrd(integer datatype, int layout, char compq, char compz, integer n,
                             integer ilo, integer ihi, void *a, integer lda, void *b, integer ldb,
                             void *q, integer ldq, void *z, integer ldz);
integer invoke_lapacke_gtsv(integer datatype, integer layout, integer n, integer nrhs, void *dl,
                            void *d, void *du, void *B, integer ldb);
integer invoke_lapacke_hetrf(integer datatype, integer layout, char uplo, integer n, void *a,
                             integer lda, integer *ipiv);
integer invoke_lapacke_hgeqz(integer datatype, int layout, char job, char compq, char compz,
                             integer n, integer ilo, integer ihi, void *h, integer ldh, void *t,
                             integer ldt, void *alpha, void *alphar, void *alphai, void *beta,
                             void *q, integer ldq, void *z, integer ldz);
integer invoke_lapacke_hseqr(integer datatype, int layout, char job, char compz, integer n,
                             integer ilo, integer ihi, void *h, integer ldh, void *w, void *wr,
                             void *wi, void *z, integer ldz);
integer invoke_lapacke_larfg(integer datatype, integer *n, void *x, integer *incx,
                             integer *abs_incx, void *tau);
integer invoke_lapacke_sytrf_rook(integer datatype, integer layout, char uplo, integer n, void *a,
                                  integer lda, integer *ipiv);
integer invoke_lapacke_hetrf_rook(integer datatype, integer layout, char uplo, integer n, void *a,
                                  integer lda, integer *ipiv);
integer invoke_lapacke_ormqr(integer datatype, int matrix_layout, char side, char trans, integer m,
                             integer n, integer k, void *a, integer lda, const void *tau, void *c,
                             integer ldc);
#endif