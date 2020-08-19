/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

    Copyright (c) 2025 Advanced Micro Devices, Inc. All rights reserved.
*/
int cbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, int *m, int *p, int *q, float *theta, float *phi, scomplex *u1, int *ldu1, scomplex *u2, int *ldu2, scomplex *v1t, int *ldv1t, scomplex *v2t, int *ldv2t, float * b11d, float *b11e, float *b12d, float *b12e, float *b21d, float *b21e, float *b22d, float *b22e, float *rwork, int *lrwork, int *info);
int cbdsqr_check(char *uplo, int *n, int *ncvt, int * nru, int *ncc, float *d__, float *e, scomplex *vt, int *ldvt, scomplex *u, int *ldu, scomplex *c__, int *ldc, float *rwork, int *info);
int cgbbrd_check(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, scomplex *ab, int *ldab, float *d__, float *e, scomplex *q, int *ldq, scomplex *pt, int *ldpt, scomplex *c__, int *ldc, scomplex *work, float *rwork, int *info);
int cgbcon_check(char *norm, int *n, int *kl, int *ku, scomplex *ab, int *ldab, int *ipiv, float *anorm, float *rcond, scomplex *work, float *rwork, int *info);
int cgbequ_check(int *m, int *n, int *kl, int *ku, scomplex *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);
int cgbequb_check(int *m, int *n, int *kl, int * ku, scomplex *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);
int cgbrfs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, scomplex *ab, int *ldab, scomplex *afb, int * ldafb, int *ipiv, scomplex *b, int *ldb, scomplex *x, int * ldx, float *ferr, float *berr, scomplex *work, float *rwork, int * info);
int cgbrfsx_check(char *trans, char *equed, int *n, int * kl, int *ku, int *nrhs, scomplex *ab, int *ldab, scomplex * afb, int *ldafb, int *ipiv, float *r__, float *c__, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int *nparams, float *params, scomplex *work, float *rwork, int * info);
int cgbsv_check(int *n, int *kl, int *ku, int * nrhs, scomplex *ab, int *ldab, int *ipiv, scomplex *b, int * ldb, int *info);
int cgbsvx_check(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, scomplex *ab, int *ldab, scomplex *afb, int *ldafb, int *ipiv, char *equed, float *r__, float *c__, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cgbsvxx_check(char *fact, char *trans, int *n, int * kl, int *ku, int *nrhs, scomplex *ab, int *ldab, scomplex * afb, int *ldafb, int *ipiv, char *equed, float *r__, float *c__, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int *n_err_bnds__, float * err_bnds_norm__, float *err_bnds_comp__, int *nparams, float * params, scomplex *work, float *rwork, int *info);
int cgbtf2_check(int *m, int *n, int *kl, int *ku, scomplex *ab, int *ldab, int *ipiv, int *info);
int cgbtrf_check(int *m, int *n, int *kl, int *ku, scomplex *ab, int *ldab, int *ipiv, int *info);
int cgbtrs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, scomplex *ab, int *ldab, int *ipiv, scomplex *b, int *ldb, int *info);
int cgebak_check(char *job, char *side, int *n, int *ilo, int *ihi, float *scale, int *m, scomplex *v, int *ldv, int *info);
int cgebal_check(char *job, int *n, scomplex *a, int *lda, int *ilo, int *ihi, float *scale, int *info);
int cgebd2_check(int *m, int *n, scomplex *a, int *lda, float *d__, float *e, scomplex *tauq, scomplex *taup, scomplex *work, int *info);
int cgebrd_check(int *m, int *n, scomplex *a, int *lda, float *d__, float *e, scomplex *tauq, scomplex *taup, scomplex *work, int *lwork, int *info);
int cgecon_check(char *norm, int *n, scomplex *a, int *lda, float *anorm, float *rcond, scomplex *work, float *rwork, int *info);
int cgeequ_check(int *m, int *n, scomplex *a, int *lda, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);
int cgeequb_check(int *m, int *n, scomplex *a, int * lda, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);
int cgees_check(char *jobvs, char *sort, L_fp select, int *n, scomplex *a, int *lda, int *sdim, scomplex *w, scomplex *vs, int *ldvs, scomplex *work, int *lwork, float *rwork, logical * bwork, int *info);
int cgeesx_check(char *jobvs, char *sort, L_fp select, char * sense, int *n, scomplex *a, int *lda, int *sdim, scomplex * w, scomplex *vs, int *ldvs, float *rconde, float *rcondv, scomplex * work, int *lwork, float *rwork, logical *bwork, int *info);
int cgeev_check(char *jobvl, char *jobvr, int *n, scomplex *a, int *lda, scomplex *w, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, scomplex *work, int *lwork, float *rwork, int * info);
int cgeevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, scomplex *a, int *lda, scomplex *w, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, int *ilo, int *ihi, float *scale, float *abnrm, float *rconde, float *rcondv, scomplex *work, int *lwork, float *rwork, int *info);
int cgegs_check(char *jobvsl, char *jobvsr, int *n, scomplex * a, int *lda, scomplex *b, int *ldb, scomplex *alpha, scomplex * beta, scomplex *vsl, int *ldvsl, scomplex *vsr, int *ldvsr, scomplex *work, int *lwork, float *rwork, int *info);
int cgegv_check(char *jobvl, char *jobvr, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, scomplex * work, int *lwork, float *rwork, int *info);
int cgehd2_check(int *n, int *ilo, int *ihi, scomplex * a, int *lda, scomplex *tau, scomplex *work, int *info);
int cgehrd_check(int *n, int *ilo, int *ihi, scomplex * a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cgelq2_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cgelqf_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cgels_check(char *trans, int *m, int *n, int * nrhs, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex * work, int *lwork, int *info);
int cgelsd_check(int *m, int *n, int *nrhs, scomplex * a, int *lda, scomplex *b, int *ldb, float *s, float *rcond, int *rank, scomplex *work, int *lwork, float *rwork, int * iwork, int *info);
int cgelss_check(int *m, int *n, int *nrhs, scomplex * a, int *lda, scomplex *b, int *ldb, float *s, float *rcond, int *rank, scomplex *work, int *lwork, float *rwork, int * info);
int cgelsx_check(int *m, int *n, int *nrhs, scomplex * a, int *lda, scomplex *b, int *ldb, int *jpvt, float *rcond, int *rank, scomplex *work, float *rwork, int *info);
int cgelsy_check(int *m, int *n, int *nrhs, scomplex * a, int *lda, scomplex *b, int *ldb, int *jpvt, float *rcond, int *rank, scomplex *work, int *lwork, float *rwork, int * info);
int cgemqrt_check(char *side, char *trans, int *m, int *n, int *k, int *nb, scomplex *v, int *ldv, scomplex *t, int *ldt, scomplex *c__, int *ldc, scomplex *work, int * info);
int cgeql2_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cgeqlf_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cgeqp3_check(int *m, int *n, scomplex *a, int *lda, int *jpvt, scomplex *tau, scomplex *work, int *lwork, float * rwork, int *info);
int cgeqpf_check(int *m, int *n, scomplex *a, int *lda, int *jpvt, scomplex *tau, scomplex *work, float *rwork, int * info);
int cgeqr2_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cgeqr2p_check(int *m, int *n, scomplex *a, int * lda, scomplex *tau, scomplex *work, int *info);
int cgeqrf_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cgeqrfp_check(int *m, int *n, scomplex *a, int * lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cgeqrt_check(int *m, int *n, int *nb, scomplex *a, int *lda, scomplex *t, int *ldt, scomplex *work, int *info);
int cgeqrt2_check(int *m, int *n, scomplex *a, int * lda, scomplex *t, int *ldt, int *info);
int cgeqrt3_check(int *m, int *n, scomplex *a, int * lda, scomplex *t, int *ldt, int *info);
int cgerfs_check(char *trans, int *n, int *nrhs, scomplex * a, int *lda, scomplex *af, int *ldaf, int *ipiv, scomplex * b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cgerfsx_check(char *trans, char *equed, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, float *r__, float *c__, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float * err_bnds_norm__, float *err_bnds_comp__, int *nparams, float * params, scomplex *work, float *rwork, int *info);
int cgerq2_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cgerqf_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cgesc2_check(int *n, scomplex *a, int *lda, scomplex * rhs, int *ipiv, int *jpiv, float *scale);
int cgesdd_check(char *jobz, int *m, int *n, scomplex *a, int *lda, float *s, scomplex *u, int *ldu, scomplex *vt, int *ldvt, scomplex *work, int *lwork, float *rwork, int *iwork, int *info);
int cgesv_check(int *n, int *nrhs, scomplex *a, int * lda, int *ipiv, scomplex *b, int *ldb, int *info);
int cgesvd_check(char *jobu, char *jobvt, int *m, int *n, scomplex *a, int *lda, float *s, scomplex *u, int *ldu, scomplex * vt, int *ldvt, scomplex *work, int *lwork, float *rwork, int *info);
int cgesvx_check(char *fact, char *trans, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, char *equed, float *r__, float *c__, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cgesvxx_check(char *fact, char *trans, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, char *equed, float *r__, float *c__, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int *nparams, float *params, scomplex *work, float *rwork, int * info);
int cgetc2_check(int *n, scomplex *a, int *lda, int * ipiv, int *jpiv, int *info);
int cgetf2_check(int *m, int *n, scomplex *a, int *lda, int *ipiv, int *info);
int cgetrf_check(int *m, int *n, scomplex *a, int *lda, int *ipiv, int *info);
int cgetrfnp_check(int *m, int *n, scomplex *a, int *lda, int *info);
int cgetri_check(int *n, scomplex *a, int *lda, int * ipiv, scomplex *work, int *lwork, int *info);
int cgetrs_check(char *trans, int *n, int *nrhs, scomplex * a, int *lda, int *ipiv, scomplex *b, int *ldb, int * info);
int cggbak_check(char *job, char *side, int *n, int *ilo, int *ihi, float *lscale, float *rscale, int *m, scomplex *v, int *ldv, int *info);
int cggbal_check(char *job, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, int *ilo, int *ihi, float *lscale, float *rscale, float *work, int *info);
int cgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, int *n, scomplex *a, int *lda, scomplex *b, int * ldb, int *sdim, scomplex *alpha, scomplex *beta, scomplex *vsl, int *ldvsl, scomplex *vsr, int *ldvsr, scomplex *work, int * lwork, float *rwork, logical *bwork, int *info);
int cggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, int *sdim, scomplex *alpha, scomplex *beta, scomplex * vsl, int *ldvsl, scomplex *vsr, int *ldvsr, float *rconde, float *rcondv, scomplex *work, int *lwork, float *rwork, int *iwork, int *liwork, logical *bwork, int *info);
int cggev_check(char *jobvl, char *jobvr, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, scomplex * work, int *lwork, float *rwork, int *info);
int cggevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, int *ldvl, scomplex * vr, int *ldvr, int *ilo, int *ihi, float *lscale, float * rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, scomplex *work, int *lwork, float *rwork, int *iwork, logical *bwork, int *info);
int cggglm_check(int *n, int *m, int *p, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *d__, scomplex *x, scomplex *y, scomplex *work, int *lwork, int *info);
int cgghrd_check(char *compq, char *compz, int *n, int * ilo, int *ihi, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *q, int *ldq, scomplex *z__, int *ldz, int *info);
int cgglse_check(int *m, int *n, int *p, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *c__, scomplex *d__, scomplex *x, scomplex *work, int *lwork, int *info);
int cggqrf_check(int *n, int *m, int *p, scomplex *a, int *lda, scomplex *taua, scomplex *b, int *ldb, scomplex *taub, scomplex *work, int *lwork, int *info);
int cggrqf_check(int *m, int *p, int *n, scomplex *a, int *lda, scomplex *taua, scomplex *b, int *ldb, scomplex *taub, scomplex *work, int *lwork, int *info);
int cggsvd_check(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, scomplex *a, int * lda, scomplex *b, int *ldb, float *alpha, float *beta, scomplex *u, int *ldu, scomplex *v, int *ldv, scomplex *q, int *ldq, scomplex *work, float *rwork, int *iwork, int *info);
int cggsvp_check(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, float *tola, float *tolb, int *k, int *l, scomplex *u, int *ldu, scomplex *v, int *ldv, scomplex *q, int *ldq, int *iwork, float *rwork, scomplex *tau, scomplex *work, int * info);
int cgtcon_check(char *norm, int *n, scomplex *dl, scomplex * d__, scomplex *du, scomplex *du2, int *ipiv, float *anorm, float * rcond, scomplex *work, int *info);
int cgtrfs_check(char *trans, int *n, int *nrhs, scomplex * dl, scomplex *d__, scomplex *du, scomplex *dlf, scomplex *df, scomplex * duf, scomplex *du2, int *ipiv, scomplex *b, int *ldb, scomplex * x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cgtsv_check(int *n, int *nrhs, scomplex *dl, scomplex * d__, scomplex *du, scomplex *b, int *ldb, int *info);
int cgtsvx_check(char *fact, char *trans, int *n, int * nrhs, scomplex *dl, scomplex *d__, scomplex *du, scomplex *dlf, scomplex * df, scomplex *duf, scomplex *du2, int *ipiv, scomplex *b, int * ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cgttrf_check(int *n, scomplex *dl, scomplex *d__, scomplex * du, scomplex *du2, int *ipiv, int *info);
int cgttrs_check(char *trans, int *n, int *nrhs, scomplex * dl, scomplex *d__, scomplex *du, scomplex *du2, int *ipiv, scomplex * b, int *ldb, int *info);
int cgtts2_check(int *itrans, int *n, int *nrhs, scomplex *dl, scomplex *d__, scomplex *du, scomplex *du2, int *ipiv, scomplex *b, int *ldb);
int chbev_check(char *jobz, char *uplo, int *n, int *kd, scomplex *ab, int *ldab, float *w, scomplex *z__, int *ldz, scomplex *work, float *rwork, int *info);
int chbevd_check(char *jobz, char *uplo, int *n, int *kd, scomplex *ab, int *ldab, float *w, scomplex *z__, int *ldz, scomplex *work, int *lwork, float *rwork, int *lrwork, int * iwork, int *liwork, int *info);
int chbevx_check(char *jobz, char *range, char *uplo, int *n, int *kd, scomplex *ab, int *ldab, scomplex *q, int *ldq, float *vl, float *vu, int *il, int *iu, float *abstol, int * m, float *w, scomplex *z__, int *ldz, scomplex *work, float *rwork, int *iwork, int *ifail, int *info);
int chbgst_check(char *vect, char *uplo, int *n, int *ka, int *kb, scomplex *ab, int *ldab, scomplex *bb, int *ldbb, scomplex *x, int *ldx, scomplex *work, float *rwork, int *info);
int chbgv_check(char *jobz, char *uplo, int *n, int *ka, int *kb, scomplex *ab, int *ldab, scomplex *bb, int *ldbb, float *w, scomplex *z__, int *ldz, scomplex *work, float *rwork, int *info);
int chbgvd_check(char *jobz, char *uplo, int *n, int *ka, int *kb, scomplex *ab, int *ldab, scomplex *bb, int *ldbb, float *w, scomplex *z__, int *ldz, scomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);
int chbgvx_check(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, scomplex *ab, int *ldab, scomplex *bb, int *ldbb, scomplex *q, int *ldq, float *vl, float *vu, int * il, int *iu, float *abstol, int *m, float *w, scomplex *z__, int *ldz, scomplex *work, float *rwork, int *iwork, int * ifail, int *info);
int chbtrd_check(char *vect, char *uplo, int *n, int *kd, scomplex *ab, int *ldab, float *d__, float *e, scomplex *q, int * ldq, scomplex *work, int *info);
int checon_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, float *anorm, float *rcond, scomplex *work, int * info);
int checon_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, float *anorm, float *rcond, scomplex *work, int *info);
int cheequb_check(char *uplo, int *n, scomplex *a, int * lda, float *s, float *scond, float *amax, scomplex *work, int *info);
int cheev_check(char *jobz, char *uplo, int *n, scomplex *a, int *lda, float *w, scomplex *work, int *lwork, float *rwork, int *info);
int cheevd_check(char *jobz, char *uplo, int *n, scomplex *a, int *lda, float *w, scomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);
int cheevr_check(char *jobz, char *range, char *uplo, int *n, scomplex *a, int *lda, float *vl, float *vu, int *il, int * iu, float *abstol, int *m, float *w, scomplex *z__, int *ldz, int *isuppz, scomplex *work, int *lwork, float *rwork, int * lrwork, int *iwork, int *liwork, int *info);
int cheevx_check(char *jobz, char *range, char *uplo, int *n, scomplex *a, int *lda, float *vl, float *vu, int *il, int * iu, float *abstol, int *m, float *w, scomplex *z__, int *ldz, scomplex *work, int *lwork, float *rwork, int *iwork, int * ifail, int *info);
int chegs2_check(int *itype, char *uplo, int *n, scomplex * a, int *lda, scomplex *b, int *ldb, int *info);
int chegst_check(int *itype, char *uplo, int *n, scomplex * a, int *lda, scomplex *b, int *ldb, int *info);
int chegv_check(int *itype, char *jobz, char *uplo, int * n, scomplex *a, int *lda, scomplex *b, int *ldb, float *w, scomplex *work, int *lwork, float *rwork, int *info);
int chegvd_check(int *itype, char *jobz, char *uplo, int * n, scomplex *a, int *lda, scomplex *b, int *ldb, float *w, scomplex *work, int *lwork, float *rwork, int *lrwork, int * iwork, int *liwork, int *info);
int chegvx_check(int *itype, char *jobz, char *range, char * uplo, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, float *vl, float *vu, int *il, int *iu, float *abstol, int * m, float *w, scomplex *z__, int *ldz, scomplex *work, int *lwork, float *rwork, int *iwork, int *ifail, int *info);
int cherfs_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, scomplex *af, int *ldaf, int *ipiv, scomplex * b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cherfsx_check(char *uplo, char *equed, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int *nparams, float *params, scomplex *work, float *rwork, int *info);
int chesv_check(char *uplo, int *n, int *nrhs, scomplex *a, int *lda, int *ipiv, scomplex *b, int *ldb, scomplex *work, int *lwork, int *info);
int chesv_rook_check(char *uplo, int *n, int *nrhs, scomplex *a, int *lda, int *ipiv, scomplex *b, int *ldb, scomplex *work, int *lwork, int *info);
int chesvx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, int *lwork, float *rwork, int *info);
int chesvxx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, char *equed, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int * n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int * nparams, float *params, scomplex *work, float *rwork, int *info);
int cheswapr_check(char *uplo, int *n, scomplex *a, int * lda, int *i1, int *i2);
int chetd2_check(char *uplo, int *n, scomplex *a, int *lda, float *d__, float *e, scomplex *tau, int *info);
int chetf2_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, int *info);
int chetf2_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, int *info);
int chetrd_check(char *uplo, int *n, scomplex *a, int *lda, float *d__, float *e, scomplex *tau, scomplex *work, int *lwork, int *info);
int chetrf_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *lwork, int *info);
int chetrf_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *lwork, int * info);
int chetri_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *info);
int chetri2_check(char *uplo, int *n, scomplex *a, int * lda, int *ipiv, scomplex *work, int *lwork, int *info);
int chetri2x_check(char *uplo, int *n, scomplex *a, int * lda, int *ipiv, scomplex *work, int *nb, int *info);
int chetri_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *info);
int chetrs_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, int *ipiv, scomplex *b, int *ldb, int * info);
int chetrs2_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, int *ipiv, scomplex *b, int *ldb, scomplex * work, int *info);
int chetrs_rook_check(char *uplo, int *n, int *nrhs, scomplex *a, int *lda, int *ipiv, scomplex *b, int *ldb, int *info);
int chfrk_check(char *transr, char *uplo, char *trans, int *n, int *k, float *alpha, scomplex *a, int *lda, float *beta, scomplex *c__);
int chgeqz_check(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, scomplex *h__, int *ldh, scomplex *t, int *ldt, scomplex *alpha, scomplex *beta, scomplex *q, int *ldq, scomplex *z__, int *ldz, scomplex *work, int *lwork, float * rwork, int *info);
VOID chla_transtype_check(char *ret_val, int *trans);
int chpcon_check(char *uplo, int *n, scomplex *ap, int * ipiv, float *anorm, float *rcond, scomplex *work, int *info);
int chpev_check(char *jobz, char *uplo, int *n, scomplex *ap, float *w, scomplex *z__, int *ldz, scomplex *work, float *rwork, int *info);
int chpevd_check(char *jobz, char *uplo, int *n, scomplex *ap, float *w, scomplex *z__, int *ldz, scomplex *work, int *lwork, float *rwork, int *lrwork, int *iwork, int *liwork, int *info);
int chpevx_check(char *jobz, char *range, char *uplo, int *n, scomplex *ap, float *vl, float *vu, int *il, int *iu, float * abstol, int *m, float *w, scomplex *z__, int *ldz, scomplex * work, float *rwork, int *iwork, int *ifail, int *info);
int chpgst_check(int *itype, char *uplo, int *n, scomplex * ap, scomplex *bp, int *info);
int chpgv_check(int *itype, char *jobz, char *uplo, int * n, scomplex *ap, scomplex *bp, float *w, scomplex *z__, int *ldz, scomplex *work, float *rwork, int *info);
int chpgvd_check(int *itype, char *jobz, char *uplo, int * n, scomplex *ap, scomplex *bp, float *w, scomplex *z__, int *ldz, scomplex *work, int *lwork, float *rwork, int *lrwork, int * iwork, int *liwork, int *info);
int chpgvx_check(int *itype, char *jobz, char *range, char * uplo, int *n, scomplex *ap, scomplex *bp, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, scomplex * z__, int *ldz, scomplex *work, float *rwork, int *iwork, int *ifail, int *info);
int chprfs_check(char *uplo, int *n, int *nrhs, scomplex * ap, scomplex *afp, int *ipiv, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int chpsv_check(char *uplo, int *n, int *nrhs, scomplex * ap, int *ipiv, scomplex *b, int *ldb, int *info);
int chpsvx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *ap, scomplex *afp, int *ipiv, scomplex *b, int * ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int chptrd_check(char *uplo, int *n, scomplex *ap, float *d__, float *e, scomplex *tau, int *info);
int chptrf_check(char *uplo, int *n, scomplex *ap, int * ipiv, int *info);
int chptri_check(char *uplo, int *n, scomplex *ap, int * ipiv, scomplex *work, int *info);
int chptrs_check(char *uplo, int *n, int *nrhs, scomplex * ap, int *ipiv, scomplex *b, int *ldb, int *info);
int chsein_check(char *side, char *eigsrc, char *initv, logical * select, int *n, scomplex *h__, int *ldh, scomplex *w, scomplex * vl, int *ldvl, scomplex *vr, int *ldvr, int *mm, int * m, scomplex *work, float *rwork, int *ifaill, int *ifailr, int *info);
int chseqr_check(char *job, char *compz, int *n, int *ilo, int *ihi, scomplex *h__, int *ldh, scomplex *w, scomplex *z__, int *ldz, scomplex *work, int *lwork, int *info);
int cla_gbamv_check(int *trans, int *m, int *n, int *kl, int *ku, float *alpha, scomplex *ab, int *ldab, scomplex *x, int *incx, float *beta, float *y, int *incy);
float cla_gbrcond_c_check(char *trans, int *n, int *kl, int *ku, scomplex *ab, int *ldab, scomplex *afb, int *ldafb, int * ipiv, float *c__, logical *capply, int *info, scomplex *work, float * rwork);
float cla_gbrcond_x_check(char *trans, int *n, int *kl, int *ku, scomplex *ab, int *ldab, scomplex *afb, int *ldafb, int * ipiv, scomplex *x, int *info, scomplex *work, float *rwork);
int cla_gbrfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *kl, int *ku, int *nrhs, scomplex *ab, int *ldab, scomplex *afb, int *ldafb, int * ipiv, logical *colequ, float *c__, scomplex *b, int *ldb, scomplex * y, int *ldy, float *berr_out__, int *n_norms__, float * err_bnds_norm__, float *err_bnds_comp__, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond, int *ithresh, float * rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float cla_gbrpvgrw_check(int *n, int *kl, int *ku, int *ncols, scomplex *ab, int *ldab, scomplex *afb, int *ldafb);
int cla_geamv_check(int *trans, int *m, int *n, float *alpha, scomplex *a, int *lda, scomplex *x, int *incx, float * beta, float *y, int *incy);
float cla_gercond_c_check(char *trans, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, float *c__, logical *capply, int *info, scomplex *work, float *rwork);
float cla_gercond_x_check(char *trans, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, scomplex *x, int *info, scomplex *work, float *rwork);
int cla_gerfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, logical *colequ, float *c__, scomplex *b, int *ldb, scomplex *y, int *ldy, float *berr_out__, int *n_norms__, float *errs_n__, float *errs_c__, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond, int * ithresh, float *rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float cla_gerpvgrw_check(int *n, int *ncols, scomplex *a, int *lda, scomplex *af, int *ldaf);
int cla_heamv_check(int *uplo, int *n, float *alpha, scomplex *a, int *lda, scomplex *x, int *incx, float *beta, float *y, int *incy);
float cla_hercond_c_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, float *c__, logical *capply, int *info, scomplex *work, float *rwork);
float cla_hercond_x_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, scomplex *x, int *info, scomplex *work, float *rwork);
int cla_herfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, logical *colequ, float *c__, scomplex *b, int *ldb, scomplex *y, int *ldy, float *berr_out__, int * n_norms__, float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond, int * ithresh, float *rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float cla_herpvgrw_check(char *uplo, int *n, int *info, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, float *work);
int cla_lin_berr_check(int *n, int *nz, int *nrhs, scomplex *res, float *ayb, float *berr);
float cla_porcond_c_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, float *c__, logical *capply, int *info, scomplex *work, float *rwork);
float cla_porcond_x_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, scomplex *x, int *info, scomplex *work, float *rwork);
int cla_porfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, logical *colequ, float *c__, scomplex *b, int *ldb, scomplex *y, int *ldy, float *berr_out__, int *n_norms__, float * err_bnds_norm__, float *err_bnds_comp__, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond, int *ithresh, float * rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float cla_porpvgrw_check(char *uplo, int *ncols, scomplex *a, int *lda, scomplex *af, int *ldaf, float *work);
int cla_syamv_check(int *uplo, int *n, float *alpha, scomplex *a, int *lda, scomplex *x, int *incx, float *beta, float *y, int *incy);
float cla_syrcond_c_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, float *c__, logical *capply, int *info, scomplex *work, float *rwork);
float cla_syrcond_x_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, scomplex *x, int *info, scomplex *work, float *rwork);
int cla_syrfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, logical *colequ, float *c__, scomplex *b, int *ldb, scomplex *y, int *ldy, float *berr_out__, int * n_norms__, float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res, float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond, int * ithresh, float *rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float cla_syrpvgrw_check(char *uplo, int *n, int *info, scomplex *a, int *lda, scomplex *af, int *ldaf, int *ipiv, float *work);
int cla_wwaddw_check(int *n, scomplex *x, scomplex *y, scomplex *w);
int clabrd_check(int *m, int *n, int *nb, scomplex *a, int *lda, float *d__, float *e, scomplex *tauq, scomplex *taup, scomplex *x, int *ldx, scomplex *y, int *ldy);
int clacgv_check(int *n, scomplex *x, int *incx);
int clacn2_check(int *n, scomplex *v, scomplex *x, float *est, int *kase, int *isave);
int clacon_check(int *n, scomplex *v, scomplex *x, float *est, int *kase);
int clacp2_check(char *uplo, int *m, int *n, float *a, int *lda, scomplex *b, int *ldb);
int clacpy_check(char *uplo, int *m, int *n, scomplex *a, int *lda, scomplex *b, int *ldb);
int clacrm_check(int *m, int *n, scomplex *a, int *lda, float *b, int *ldb, scomplex *c__, int *ldc, float *rwork);
int clacrt_check(int *n, scomplex *cx, int *incx, scomplex * cy, int *incy, scomplex *c__, scomplex *s);
VOID cladiv_check(scomplex * ret_val, scomplex *x, scomplex *y);
int claed0_check(int *qsiz, int *n, float *d__, float *e, scomplex *q, int *ldq, scomplex *qstore, int *ldqs, float *rwork, int *iwork, int *info);
int claed7_check(int *n, int *cutpnt, int *qsiz, int *tlvls, int *curlvl, int *curpbm, float *d__, scomplex * q, int *ldq, float *rho, int *indxq, float *qstore, int * qptr, int *prmptr, int *perm, int *givptr, int * givcol, float *givnum, scomplex *work, float *rwork, int *iwork, int *info);
int claed8_check(int *k, int *n, int *qsiz, scomplex * q, int *ldq, float *d__, float *rho, int *cutpnt, float *z__, float *dlamda, scomplex *q2, int *ldq2, float *w, int *indxp, int *indx, int *indxq, int *perm, int *givptr, int *givcol, float *givnum, int *info);
int claein_check(logical *rightv, logical *noinit, int *n, scomplex *h__, int *ldh, scomplex *w, scomplex *v, scomplex *b, int *ldb, float *rwork, float *eps3, float *smlnum, int *info);
int claesy_check(scomplex *a, scomplex *b, scomplex *c__, scomplex * rt1, scomplex *rt2, scomplex *evscal, scomplex *cs1, scomplex *sn1);
int claev2_check(scomplex *a, scomplex *b, scomplex *c__, float *rt1, float *rt2, float *cs1, scomplex *sn1);
int clag2z_check(int *m, int *n, scomplex *sa, int * ldsa, dcomplex *a, int *lda, int *info);
int clags2_check(logical *upper, float *a1, scomplex *a2, float *a3, float *b1, scomplex *b2, float *b3, float *csu, scomplex *snu, float *csv, scomplex *snv, float *csq, scomplex *snq);
int clagtm_check(char *trans, int *n, int *nrhs, float * alpha, scomplex *dl, scomplex *d__, scomplex *du, scomplex *x, int * ldx, float *beta, scomplex *b, int *ldb);
int clahef_check(char *uplo, int *n, int *nb, int *kb, scomplex *a, int *lda, int *ipiv, scomplex *w, int *ldw, int *info);
int clahef_rook_check(char *uplo, int *n, int *nb, int *kb, scomplex *a, int *lda, int *ipiv, scomplex *w, int *ldw, int *info);
int clahqr_check(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, scomplex *h__, int *ldh, scomplex *w, int *iloz, int *ihiz, scomplex *z__, int *ldz, int * info);
int clahr2_check(int *n, int *k, int *nb, scomplex *a, int *lda, scomplex *tau, scomplex *t, int *ldt, scomplex *y, int *ldy);
int clahrd_check(int *n, int *k, int *nb, scomplex *a, int *lda, scomplex *tau, scomplex *t, int *ldt, scomplex *y, int *ldy);
int claic1_check(int *job, int *j, scomplex *x, float *sest, scomplex *w, scomplex *gamma, float *sestpr, scomplex *s, scomplex *c__);
int clals0_check(int *icompq, int *nl, int *nr, int *sqre, int *nrhs, scomplex *b, int *ldb, scomplex *bx, int *ldbx, int *perm, int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum, float *poles, float * difl, float *difr, float *z__, int *k, float *c__, float *s, float * rwork, int *info);
int clalsa_check(int *icompq, int *smlsiz, int *n, int *nrhs, scomplex *b, int *ldb, scomplex *bx, int *ldbx, float *u, int *ldu, float *vt, int *k, float *difl, float *difr, float *z__, float *poles, int *givptr, int *givcol, int * ldgcol, int *perm, float *givnum, float *c__, float *s, float *rwork, int *iwork, int *info);
int clalsd_check(char *uplo, int *smlsiz, int *n, int *nrhs, float *d__, float *e, scomplex *b, int *ldb, float *rcond, int *rank, scomplex *work, float *rwork, int *iwork, int * info);
float clangb_check(char *norm, int *n, int *kl, int *ku, scomplex *ab, int *ldab, float *work);
float clange_check(char *norm, int *m, int *n, scomplex *a, int *lda, float *work);
float clangt_check(char *norm, int *n, scomplex *dl, scomplex *d__, scomplex *du);
float clanhb_check(char *norm, char *uplo, int *n, int *k, scomplex *ab, int *ldab, float *work);
float clanhe_check(char *norm, char *uplo, int *n, scomplex *a, int *lda, float *work);
float clanhf_check(char *norm, char *transr, char *uplo, int *n, scomplex *a, float *work);
float clanhp_check(char *norm, char *uplo, int *n, scomplex *ap, float *work);
float clanhs_check(char *norm, int *n, scomplex *a, int *lda, float *work);
float clanht_check(char *norm, int *n, float *d__, scomplex *e);
float clansb_check(char *norm, char *uplo, int *n, int *k, scomplex *ab, int *ldab, float *work);
float clansp_check(char *norm, char *uplo, int *n, scomplex *ap, float *work);
float clansy_check(char *norm, char *uplo, int *n, scomplex *a, int *lda, float *work);
float clantb_check(char *norm, char *uplo, char *diag, int *n, int *k, scomplex *ab, int *ldab, float *work);
float clantp_check(char *norm, char *uplo, char *diag, int *n, scomplex *ap, float *work);
float clantr_check(char *norm, char *uplo, char *diag, int *m, int *n, scomplex *a, int *lda, float *work);
int clapll_check(int *n, scomplex *x, int *incx, scomplex * y, int *incy, float *ssmin);
int clapmr_check(logical *forwrd, int *m, int *n, scomplex *x, int *ldx, int *k);
int clapmt_check(logical *forwrd, int *m, int *n, scomplex *x, int *ldx, int *k);
int claqgb_check(int *m, int *n, int *kl, int *ku, scomplex *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *equed);
int claqge_check(int *m, int *n, scomplex *a, int *lda, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char * equed);
int claqhb_check(char *uplo, int *n, int *kd, scomplex *ab, int *ldab, float *s, float *scond, float *amax, char *equed);
int claqhe_check(char *uplo, int *n, scomplex *a, int *lda, float *s, float *scond, float *amax, char *equed);
int claqhp_check(char *uplo, int *n, scomplex *ap, float *s, float *scond, float *amax, char *equed);
int claqp2_check(int *m, int *n, int *offset, scomplex *a, int *lda, int *jpvt, scomplex *tau, float *vn1, float *vn2, scomplex *work);
int claqps_check(int *m, int *n, int *offset, int *nb, int *kb, scomplex *a, int *lda, int *jpvt, scomplex * tau, float *vn1, float *vn2, scomplex *auxv, scomplex *f, int *ldf);
int claqr0_check(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, scomplex *h__, int *ldh, scomplex *w, int *iloz, int *ihiz, scomplex *z__, int *ldz, scomplex * work, int *lwork, int *info);
int claqr1_check(int *n, scomplex *h__, int *ldh, scomplex * s1, scomplex *s2, scomplex *v);
int claqr2_check(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, scomplex *h__, int *ldh, int *iloz, int *ihiz, scomplex *z__, int *ldz, int * ns, int *nd, scomplex *sh, scomplex *v, int *ldv, int *nh, scomplex *t, int *ldt, int *nv, scomplex *wv, int *ldwv, scomplex *work, int *lwork);
int claqr3_check(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, scomplex *h__, int *ldh, int *iloz, int *ihiz, scomplex *z__, int *ldz, int * ns, int *nd, scomplex *sh, scomplex *v, int *ldv, int *nh, scomplex *t, int *ldt, int *nv, scomplex *wv, int *ldwv, scomplex *work, int *lwork);
int claqr4_check(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, scomplex *h__, int *ldh, scomplex *w, int *iloz, int *ihiz, scomplex *z__, int *ldz, scomplex * work, int *lwork, int *info);
int claqr5_check(logical *wantt, logical *wantz, int *kacc22, int *n, int *ktop, int *kbot, int *nshfts, scomplex *s, scomplex *h__, int *ldh, int *iloz, int *ihiz, scomplex * z__, int *ldz, scomplex *v, int *ldv, scomplex *u, int *ldu, int *nv, scomplex *wv, int *ldwv, int *nh, scomplex *wh, int *ldwh);
int claqsb_check(char *uplo, int *n, int *kd, scomplex *ab, int *ldab, float *s, float *scond, float *amax, char *equed);
int claqsp_check(char *uplo, int *n, scomplex *ap, float *s, float *scond, float *amax, char *equed);
int claqsy_check(char *uplo, int *n, scomplex *a, int *lda, float *s, float *scond, float *amax, char *equed);
int clar1v_check(int *n, int *b1, int *bn, float * lambda, float *d__, float *l, float *ld, float *lld, float *pivmin, float * gaptol, scomplex *z__, logical *wantnc, int *negcnt, float *ztz, float *mingma, int *r__, int *isuppz, float *nrminv, float * resid, float *rqcorr, float *work);
int clar2v_check(int *n, scomplex *x, scomplex *y, scomplex *z__, int *incx, float *c__, scomplex *s, int *incc);
int clarcm_check(int *m, int *n, float *a, int *lda, scomplex *b, int *ldb, scomplex *c__, int *ldc, float *rwork);
int clarf_check(char *side, int *m, int *n, scomplex *v, int *incv, scomplex *tau, scomplex *c__, int *ldc, scomplex * work);
int clarfb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, scomplex *v, int *ldv, scomplex *t, int *ldt, scomplex *c__, int *ldc, scomplex *work, int *ldwork);
int clarfg_check(int *n, scomplex *alpha, scomplex *x, int * incx, scomplex *tau);
int clarfgp_check(int *n, scomplex *alpha, scomplex *x, int *incx, scomplex *tau);
int clarft_check(char *direct, char *storev, int *n, int * k, scomplex *v, int *ldv, scomplex *tau, scomplex *t, int *ldt);
int clarfx_check(char *side, int *m, int *n, scomplex *v, scomplex *tau, scomplex *c__, int *ldc, scomplex *work);
int clargv_check(int *n, scomplex *x, int *incx, scomplex * y, int *incy, float *c__, int *incc);
int clarnv_check(int *idist, int *iseed, int *n, scomplex *x);
int clarrv_check(int *n, float *vl, float *vu, float *d__, float * l, float *pivmin, int *isplit, int *m, int *dol, int * dou, float *minrgp, float *rtol1, float *rtol2, float *w, float *werr, float *wgap, int *iblock, int *indexw, float *gers, scomplex * z__, int *ldz, int *isuppz, float *work, int *iwork, int *info);
int clarscl2_check(int *m, int *n, float *d__, scomplex *x, int *ldx);
int clartg_check(scomplex *f, scomplex *g, float *cs, scomplex *sn, scomplex *r__);
int clartv_check(int *n, scomplex *x, int *incx, scomplex * y, int *incy, float *c__, scomplex *s, int *incc);
int clarz_check(char *side, int *m, int *n, int *l, scomplex *v, int *incv, scomplex *tau, scomplex *c__, int *ldc, scomplex *work);
int clarzb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, int *l, scomplex *v, int *ldv, scomplex *t, int *ldt, scomplex *c__, int *ldc, scomplex *work, int *ldwork);
int clarzt_check(char *direct, char *storev, int *n, int * k, scomplex *v, int *ldv, scomplex *tau, scomplex *t, int *ldt);
int clascl_check(char *type__, int *kl, int *ku, float * cfrom, float *cto, int *m, int *n, scomplex *a, int *lda, int *info);
int clascl2_check(int *m, int *n, float *d__, scomplex *x, int *ldx);
int claset_check(char *uplo, int *m, int *n, scomplex * alpha, scomplex *beta, scomplex *a, int *lda);
int clasr_check(char *side, char *pivot, char *direct, int *m, int *n, float *c__, float *s, scomplex *a, int *lda);
int classq_check(int *n, scomplex *x, int *incx, float * scale, float *sumsq);
int claswp_check(int *n, scomplex *a, int *lda, int * k1, int *k2, int *ipiv, int *incx);
int clasyf_check(char *uplo, int *n, int *nb, int *kb, scomplex *a, int *lda, int *ipiv, scomplex *w, int *ldw, int *info);
int clasyf_rook_check(char *uplo, int *n, int *nb, int *kb, scomplex *a, int *lda, int *ipiv, scomplex *w, int *ldw, int *info);
int clatbs_check(char *uplo, char *trans, char *diag, char * normin, int *n, int *kd, scomplex *ab, int *ldab, scomplex * x, float *scale, float *cnorm, int *info);
int clatdf_check(int *ijob, int *n, scomplex *z__, int *ldz, scomplex *rhs, float *rdsum, float *rdscal, int *ipiv, int *jpiv);
int clatps_check(char *uplo, char *trans, char *diag, char * normin, int *n, scomplex *ap, scomplex *x, float *scale, float *cnorm, int *info);
int clatrd_check(char *uplo, int *n, int *nb, scomplex *a, int *lda, float *e, scomplex *tau, scomplex *w, int *ldw);
int clatrs_check(char *uplo, char *trans, char *diag, char * normin, int *n, scomplex *a, int *lda, scomplex *x, float *scale, float *cnorm, int *info);
int clatrz_check(int *m, int *n, int *l, scomplex *a, int *lda, scomplex *tau, scomplex *work);
int clatzm_check(char *side, int *m, int *n, scomplex *v, int *incv, scomplex *tau, scomplex *c1, scomplex *c2, int *ldc, scomplex *work);
int clauu2_check(char *uplo, int *n, scomplex *a, int *lda, int *info);
int clauum_check(char *uplo, int *n, scomplex *a, int *lda, int *info);
int cpbcon_check(char *uplo, int *n, int *kd, scomplex *ab, int *ldab, float *anorm, float *rcond, scomplex *work, float *rwork, int *info);
int cpbequ_check(char *uplo, int *n, int *kd, scomplex *ab, int *ldab, float *s, float *scond, float *amax, int *info);
int cpbrfs_check(char *uplo, int *n, int *kd, int * nrhs, scomplex *ab, int *ldab, scomplex *afb, int *ldafb, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float * berr, scomplex *work, float *rwork, int *info);
int cpbstf_check(char *uplo, int *n, int *kd, scomplex *ab, int *ldab, int *info);
int cpbsv_check(char *uplo, int *n, int *kd, int * nrhs, scomplex *ab, int *ldab, scomplex *b, int *ldb, int * info);
int cpbsvx_check(char *fact, char *uplo, int *n, int *kd, int *nrhs, scomplex *ab, int *ldab, scomplex *afb, int * ldafb, char *equed, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cpbtf2_check(char *uplo, int *n, int *kd, scomplex *ab, int *ldab, int *info);
int cpbtrf_check(char *uplo, int *n, int *kd, scomplex *ab, int *ldab, int *info);
int cpbtrs_check(char *uplo, int *n, int *kd, int * nrhs, scomplex *ab, int *ldab, scomplex *b, int *ldb, int * info);
int cpftrf_check(char *transr, char *uplo, int *n, scomplex *a, int *info);
int cpftri_check(char *transr, char *uplo, int *n, scomplex *a, int *info);
int cpftrs_check(char *transr, char *uplo, int *n, int * nrhs, scomplex *a, scomplex *b, int *ldb, int *info);
int cpocon_check(char *uplo, int *n, scomplex *a, int *lda, float *anorm, float *rcond, scomplex *work, float *rwork, int *info);
int cpoequ_check(int *n, scomplex *a, int *lda, float *s, float *scond, float *amax, int *info);
int cpoequb_check(int *n, scomplex *a, int *lda, float *s, float *scond, float *amax, int *info);
int cporfs_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, scomplex *af, int *ldaf, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cporfsx_check(char *uplo, char *equed, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float * err_bnds_comp__, int *nparams, float *params, scomplex *work, float * rwork, int *info);
int cposv_check(char *uplo, int *n, int *nrhs, scomplex *a, int *lda, scomplex *b, int *ldb, int *info);
int cposvx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, char * equed, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cposvxx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, char * equed, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int *n_err_bnds__, float * err_bnds_norm__, float *err_bnds_comp__, int *nparams, float * params, scomplex *work, float *rwork, int *info);
int cpotf2_check(char *uplo, int *n, scomplex *a, int *lda, int *info);
int cpotrf_check(char *uplo, int *n, scomplex *a, int *lda, int *info);
int cpotri_check(char *uplo, int *n, scomplex *a, int *lda, int *info);
int cpotrs_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, scomplex *b, int *ldb, int *info);
int cppcon_check(char *uplo, int *n, scomplex *ap, float *anorm, float *rcond, scomplex *work, float *rwork, int *info);
int cppequ_check(char *uplo, int *n, scomplex *ap, float *s, float *scond, float *amax, int *info);
int cpprfs_check(char *uplo, int *n, int *nrhs, scomplex * ap, scomplex *afp, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cppsv_check(char *uplo, int *n, int *nrhs, scomplex * ap, scomplex *b, int *ldb, int *info);
int cppsvx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *ap, scomplex *afp, char *equed, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cpptrf_check(char *uplo, int *n, scomplex *ap, int * info);
int cpptri_check(char *uplo, int *n, scomplex *ap, int * info);
int cpptrs_check(char *uplo, int *n, int *nrhs, scomplex * ap, scomplex *b, int *ldb, int *info);
int cpstf2_check(char *uplo, int *n, scomplex *a, int *lda, int *piv, int *rank, float *tol, float *work, int *info);
int cpstrf_check(char *uplo, int *n, scomplex *a, int *lda, int *piv, int *rank, float *tol, float *work, int *info);
int cptcon_check(int *n, float *d__, scomplex *e, float *anorm, float *rcond, float *rwork, int *info);
int cpteqr_check(char *compz, int *n, float *d__, float *e, scomplex *z__, int *ldz, float *work, int *info);
int cptrfs_check(char *uplo, int *n, int *nrhs, float *d__, scomplex *e, float *df, scomplex *ef, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cptsv_check(int *n, int *nrhs, float *d__, scomplex *e, scomplex *b, int *ldb, int *info);
int cptsvx_check(char *fact, int *n, int *nrhs, float *d__, scomplex *e, float *df, scomplex *ef, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cpttrf_check(int *n, float *d__, scomplex *e, int *info);
int cpttrs_check(char *uplo, int *n, int *nrhs, float *d__, scomplex *e, scomplex *b, int *ldb, int *info);
int cptts2_check(int *iuplo, int *n, int *nrhs, float * d__, scomplex *e, scomplex *b, int *ldb);
int crot_check(int *n, scomplex *cx, int *incx, scomplex * cy, int *incy, float *c__, scomplex *s);
int cspcon_check(char *uplo, int *n, scomplex *ap, int * ipiv, float *anorm, float *rcond, scomplex *work, int *info);
int cspmv_check(char *uplo, int *n, scomplex *alpha, scomplex * ap, scomplex *x, int *incx, scomplex *beta, scomplex *y, int * incy);
int cspr_check(char *uplo, int *n, scomplex *alpha, scomplex *x, int *incx, scomplex *ap);
int csprfs_check(char *uplo, int *n, int *nrhs, scomplex * ap, scomplex *afp, int *ipiv, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int cspsv_check(char *uplo, int *n, int *nrhs, scomplex * ap, int *ipiv, scomplex *b, int *ldb, int *info);
int cspsvx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *ap, scomplex *afp, int *ipiv, scomplex *b, int * ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int csptrf_check(char *uplo, int *n, scomplex *ap, int * ipiv, int *info);
int csptri_check(char *uplo, int *n, scomplex *ap, int * ipiv, scomplex *work, int *info);
int csptrs_check(char *uplo, int *n, int *nrhs, scomplex * ap, int *ipiv, scomplex *b, int *ldb, int *info);
int csrscl_check(int *n, float *sa, scomplex *sx, int *incx);
int cstedc_check(char *compz, int *n, float *d__, float *e, scomplex *z__, int *ldz, scomplex *work, int *lwork, float * rwork, int *lrwork, int *iwork, int *liwork, int * info);
int cstegr_check(char *jobz, char *range, int *n, float *d__, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, scomplex *z__, int *ldz, int *isuppz, float *work, int *lwork, int *iwork, int *liwork, int * info);
int cstein_check(int *n, float *d__, float *e, int *m, float *w, int *iblock, int *isplit, scomplex *z__, int *ldz, float *work, int *iwork, int *ifail, int *info);
int cstemr_check(char *jobz, char *range, int *n, float *d__, float *e, float *vl, float *vu, int *il, int *iu, int *m, float *w, scomplex *z__, int *ldz, int *nzc, int *isuppz, logical *tryrac, float *work, int *lwork, int *iwork, int * liwork, int *info);
int csteqr_check(char *compz, int *n, float *d__, float *e, scomplex *z__, int *ldz, float *work, int *info);
int csycon_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, float *anorm, float *rcond, scomplex *work, int * info);
int csycon_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, float *anorm, float *rcond, scomplex *work, int *info);
int csyconv_check(char *uplo, char *way, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *info);
int csyequb_check(char *uplo, int *n, scomplex *a, int * lda, float *s, float *scond, float *amax, scomplex *work, int *info);
int csymv_check(char *uplo, int *n, scomplex *alpha, scomplex * a, int *lda, scomplex *x, int *incx, scomplex *beta, scomplex *y, int *incy);
int csyr_check(char *uplo, int *n, scomplex *alpha, scomplex *x, int *incx, scomplex *a, int *lda);
int csyrfs_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, scomplex *af, int *ldaf, int *ipiv, scomplex * b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int csyrfsx_check(char *uplo, char *equed, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int *nparams, float *params, scomplex *work, float *rwork, int *info);
int csysv_check(char *uplo, int *n, int *nrhs, scomplex *a, int *lda, int *ipiv, scomplex *b, int *ldb, scomplex *work, int *lwork, int *info);
int csysv_rook_check(char *uplo, int *n, int *nrhs, scomplex *a, int *lda, int *ipiv, scomplex *b, int *ldb, scomplex *work, int *lwork, int *info);
int csysvx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *ferr, float *berr, scomplex *work, int *lwork, float *rwork, int *info);
int csysvxx_check(char *fact, char *uplo, int *n, int * nrhs, scomplex *a, int *lda, scomplex *af, int *ldaf, int * ipiv, char *equed, float *s, scomplex *b, int *ldb, scomplex *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int * n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int * nparams, float *params, scomplex *work, float *rwork, int *info);
int csyswapr_check(char *uplo, int *n, scomplex *a, int * lda, int *i1, int *i2);
int csytf2_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, int *info);
int csytf2_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, int *info);
int csytrf_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *lwork, int *info);
int csytrf_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *lwork, int * info);
int csytri_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *info);
int csytri2_check(char *uplo, int *n, scomplex *a, int * lda, int *ipiv, scomplex *work, int *lwork, int *info);
int csytri2x_check(char *uplo, int *n, scomplex *a, int * lda, int *ipiv, scomplex *work, int *nb, int *info);
int csytri_rook_check(char *uplo, int *n, scomplex *a, int *lda, int *ipiv, scomplex *work, int *info);
int csytrs_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, int *ipiv, scomplex *b, int *ldb, int * info);
int csytrs2_check(char *uplo, int *n, int *nrhs, scomplex * a, int *lda, int *ipiv, scomplex *b, int *ldb, scomplex * work, int *info);
int csytrs_rook_check(char *uplo, int *n, int *nrhs, scomplex *a, int *lda, int *ipiv, scomplex *b, int *ldb, int *info);
int ctbcon_check(char *norm, char *uplo, char *diag, int *n, int *kd, scomplex *ab, int *ldab, float *rcond, scomplex *work, float *rwork, int *info);
int ctbrfs_check(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, scomplex *ab, int *ldab, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int ctbtrs_check(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, scomplex *ab, int *ldab, scomplex *b, int *ldb, int *info);
int ctfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, int *m, int *n, scomplex *alpha, scomplex *a, scomplex *b, int *ldb);
int ctftri_check(char *transr, char *uplo, char *diag, int *n, scomplex *a, int *info);
int ctfttp_check(char *transr, char *uplo, int *n, scomplex * arf, scomplex *ap, int *info);
int ctfttr_check(char *transr, char *uplo, int *n, scomplex * arf, scomplex *a, int *lda, int *info);
int ctgevc_check(char *side, char *howmny, logical *select, int *n, scomplex *s, int *lds, scomplex *p, int *ldp, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, int *mm, int *m, scomplex *work, float *rwork, int *info);
int ctgex2_check(logical *wantq, logical *wantz, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *q, int *ldq, scomplex *z__, int *ldz, int *j1, int *info);
int ctgexc_check(logical *wantq, logical *wantz, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *q, int *ldq, scomplex *z__, int *ldz, int *ifst, int * ilst, int *info);
int ctgsen_check(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *alpha, scomplex *beta, scomplex *q, int *ldq, scomplex *z__, int *ldz, int *m, float *pl, float *pr, float * dif, scomplex *work, int *lwork, int *iwork, int *liwork, int *info);
int ctgsja_check(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, scomplex *a, int * lda, scomplex *b, int *ldb, float *tola, float *tolb, float *alpha, float *beta, scomplex *u, int *ldu, scomplex *v, int *ldv, scomplex *q, int *ldq, scomplex *work, int *ncycle, int * info);
int ctgsna_check(char *job, char *howmny, logical *select, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, float *s, float *dif, int *mm, int *m, scomplex *work, int *lwork, int *iwork, int *info);
int ctgsy2_check(char *trans, int *ijob, int *m, int * n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *c__, int *ldc, scomplex *d__, int *ldd, scomplex *e, int *lde, scomplex *f, int *ldf, float *scale, float *rdsum, float *rdscal, int *info);
int ctgsyl_check(char *trans, int *ijob, int *m, int * n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *c__, int *ldc, scomplex *d__, int *ldd, scomplex *e, int *lde, scomplex *f, int *ldf, float *scale, float *dif, scomplex *work, int *lwork, int *iwork, int *info);
int ctpcon_check(char *norm, char *uplo, char *diag, int *n, scomplex *ap, float *rcond, scomplex *work, float *rwork, int *info);
int ctpmqrt_check(char *side, char *trans, int *m, int *n, int *k, int *l, int *nb, scomplex *v, int *ldv, scomplex *t, int *ldt, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *work, int *info);
int ctpqrt_check(int *m, int *n, int *l, int *nb, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *t, int *ldt, scomplex *work, int *info);
int ctpqrt2_check(int *m, int *n, int *l, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *t, int *ldt, int *info);
int ctprfb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, int *l, scomplex *v, int *ldv, scomplex *t, int *ldt, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *work, int *ldwork);
int ctprfs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, scomplex *ap, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int ctptri_check(char *uplo, char *diag, int *n, scomplex *ap, int *info);
int ctptrs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, scomplex *ap, scomplex *b, int *ldb, int *info);
int ctpttf_check(char *transr, char *uplo, int *n, scomplex * ap, scomplex *arf, int *info);
int ctpttr_check(char *uplo, int *n, scomplex *ap, scomplex *a, int *lda, int *info);
int ctrcon_check(char *norm, char *uplo, char *diag, int *n, scomplex *a, int *lda, float *rcond, scomplex *work, float *rwork, int *info);
int ctrevc_check(char *side, char *howmny, logical *select, int *n, scomplex *t, int *ldt, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, int *mm, int *m, scomplex *work, float *rwork, int *info);
int ctrexc_check(char *compq, int *n, scomplex *t, int * ldt, scomplex *q, int *ldq, int *ifst, int *ilst, int * info);
int ctrrfs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *x, int *ldx, float *ferr, float *berr, scomplex *work, float *rwork, int *info);
int ctrsen_check(char *job, char *compq, logical *select, int *n, scomplex *t, int *ldt, scomplex *q, int *ldq, scomplex *w, int *m, float *s, float *sep, scomplex *work, int *lwork, int *info);
int ctrsna_check(char *job, char *howmny, logical *select, int *n, scomplex *t, int *ldt, scomplex *vl, int *ldvl, scomplex *vr, int *ldvr, float *s, float *sep, int *mm, int * m, scomplex *work, int *ldwork, float *rwork, int *info);
int ctrsyl_check(char *trana, char *tranb, int *isgn, int *m, int *n, scomplex *a, int *lda, scomplex *b, int *ldb, scomplex *c__, int *ldc, float *scale, int *info);
int ctrti2_check(char *uplo, char *diag, int *n, scomplex *a, int *lda, int *info);
int ctrtri_check(char *uplo, char *diag, int *n, scomplex *a, int *lda, int *info);
int ctrtrs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, scomplex *a, int *lda, scomplex *b, int *ldb, int *info);
int ctrttf_check(char *transr, char *uplo, int *n, scomplex *a, int *lda, scomplex *arf, int *info);
int ctrttp_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *ap, int *info);
int ctzrqf_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, int *info);
int ctzrzf_check(int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cunbdb_check(char *trans, char *signs, int *m, int *p, int *q, scomplex *x11, int *ldx11, scomplex *x12, int * ldx12, scomplex *x21, int *ldx21, scomplex *x22, int *ldx22, float *theta, float *phi, scomplex *taup1, scomplex *taup2, scomplex * tauq1, scomplex *tauq2, scomplex *work, int *lwork, int *info);
int cunbdb1_check(int *m, int *p, int *q, scomplex * x11, int *ldx11, scomplex *x21, int *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, int *lwork, int *info);
int cunbdb2_check(int *m, int *p, int *q, scomplex * x11, int *ldx11, scomplex *x21, int *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, int *lwork, int *info);
int cunbdb3_check(int *m, int *p, int *q, scomplex * x11, int *ldx11, scomplex *x21, int *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *work, int *lwork, int *info);
int cunbdb4_check(int *m, int *p, int *q, scomplex * x11, int *ldx11, scomplex *x21, int *ldx21, float *theta, float * phi, scomplex *taup1, scomplex *taup2, scomplex *tauq1, scomplex *phantom, scomplex *work, int *lwork, int *info);
int cunbdb5_check(int *m1, int *m2, int *n, scomplex * x1, int *incx1, scomplex *x2, int *incx2, scomplex *q1, int *ldq1, scomplex *q2, int *ldq2, scomplex *work, int *lwork, int *info);
int cunbdb6_check(int *m1, int *m2, int *n, scomplex * x1, int *incx1, scomplex *x2, int *incx2, scomplex *q1, int *ldq1, scomplex *q2, int *ldq2, scomplex *work, int *lwork, int *info);
int cuncsd_check(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, char *signs, int *m, int *p, int *q, scomplex *x11, int *ldx11, scomplex *x12, int *ldx12, scomplex * x21, int *ldx21, scomplex *x22, int *ldx22, float *theta, scomplex *u1, int *ldu1, scomplex *u2, int *ldu2, scomplex *v1t, int *ldv1t, scomplex *v2t, int *ldv2t, scomplex *work, int * lwork, float *rwork, int *lrwork, int *iwork, int *info);
int cuncsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, int *m, int *p, int *q, scomplex *x11, int *ldx11, scomplex *x21, int *ldx21, float *theta, scomplex *u1, int *ldu1, scomplex *u2, int *ldu2, scomplex *v1t, int *ldv1t, scomplex * work, int *lwork, float *rwork, int *lrwork, int *iwork, int *info);
int cung2l_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cung2r_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cungbr_check(char *vect, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cunghr_check(int *n, int *ilo, int *ihi, scomplex * a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cungl2_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cunglq_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int * info);
int cungql_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int * info);
int cungqr_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int * info);
int cungr2_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *info);
int cungrq_check(int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int * info);
int cungtr_check(char *uplo, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *work, int *lwork, int *info);
int cunm2l_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *info);
int cunm2r_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *info);
int cunmbr_check(char *vect, char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int * info);
int cunmhr_check(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int * info);
int cunml2_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *info);
int cunmlq_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int *info);
int cunmql_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int *info);
int cunmqr_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int *info);
int cunmr2_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *info);
int cunmr3_check(char *side, char *trans, int *m, int *n, int *k, int *l, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *info);
int cunmrq_check(char *side, char *trans, int *m, int *n, int *k, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int *info);
int cunmrz_check(char *side, char *trans, int *m, int *n, int *k, int *l, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int * info);
int cunmtr_check(char *side, char *uplo, char *trans, int *m, int *n, scomplex *a, int *lda, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *lwork, int *info);
int cupgtr_check(char *uplo, int *n, scomplex *ap, scomplex * tau, scomplex *q, int *ldq, scomplex *work, int *info);
int cupmtr_check(char *side, char *uplo, char *trans, int *m, int *n, scomplex *ap, scomplex *tau, scomplex *c__, int *ldc, scomplex *work, int *info);
int dbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, int *m, int *p, int *q, double * theta, double *phi, double *u1, int *ldu1, double *u2, int *ldu2, double *v1t, int *ldv1t, double *v2t, int *ldv2t, double *b11d, double *b11e, double *b12d, double *b12e, double *b21d, double *b21e, double * b22d, double *b22e, double *work, int *lwork, int * info);
int dbdsdc_check(char *uplo, char *compq, int *n, double * d__, double *e, double *u, int *ldu, double *vt, int *ldvt, double *q, int *iq, double *work, int * iwork, int *info);
int dbdsqr_check(char *uplo, int *n, int *ncvt, int * nru, int *ncc, double *d__, double *e, double *vt, int *ldvt, double *u, int *ldu, double *c__, int * ldc, double *work, int *info);
int ddisna_check(char *job, int *m, int *n, double * d__, double *sep, int *info);
int dgbbrd_check(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, double *ab, int *ldab, double * d__, double *e, double *q, int *ldq, double *pt, int *ldpt, double *c__, int *ldc, double *work, int *info);
int dgbcon_check(char *norm, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);
int dgbequ_check(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int * info);
int dgbequb_check(int *m, int *n, int *kl, int * ku, double *ab, int *ldab, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int * info);
int dgbrfs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
int dgbrfsx_check(char *trans, char *equed, int *n, int * kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, double *r__, double *c__, double *b, int *ldb, double *x, int * ldx, double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int * nparams, double *params, double *work, int *iwork, int *info);
int dgbsv_check(int *n, int *kl, int *ku, int * nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);
int dgbsvx_check(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, char *equed, double *r__, double *c__, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
int dgbsvxx_check(char *fact, char *trans, int *n, int * kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, char *equed, double *r__, double *c__, double *b, int *ldb, double *x, int *ldx, double *rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *params, double *work, int *iwork, int *info);
int dgbtf2_check(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);
int dgbtrf_check(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);
int dgbtrs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);
int dgebak_check(char *job, char *side, int *n, int *ilo, int *ihi, double *scale, int *m, double *v, int * ldv, int *info);
int dgebal_check(char *job, int *n, double *a, int * lda, int *ilo, int *ihi, double *scale, int *info);
int dgebd2_check(int *m, int *n, double *a, int * lda, double *d__, double *e, double *tauq, double * taup, double *work, int *info);
int dgebrd_check(int *m, int *n, double *a, int * lda, double *d__, double *e, double *tauq, double * taup, double *work, int *lwork, int *info);
int dgecon_check(char *norm, int *n, double *a, int * lda, double *anorm, double *rcond, double *work, int * iwork, int *info);
int dgeequ_check(int *m, int *n, double *a, int * lda, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int *info);
int dgeequb_check(int *m, int *n, double *a, int * lda, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int *info);
int dgees_check(char *jobvs, char *sort, L_fp select, int *n, double *a, int *lda, int *sdim, double *wr, double *wi, double *vs, int *ldvs, double *work, int *lwork, logical *bwork, int *info);
int dgeesx_check(char *jobvs, char *sort, L_fp select, char * sense, int *n, double *a, int *lda, int *sdim, double *wr, double *wi, double *vs, int *ldvs, double *rconde, double *rcondv, double *work, int * lwork, int *iwork, int *liwork, logical *bwork, int *info);
int dgeev_check(char *jobvl, char *jobvr, int *n, double * a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
int dgeevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, int *ilo, int *ihi, double *scale, double *abnrm, double *rconde, double *rcondv, double *work, int *lwork, int *iwork, int *info);
int dgegs_check(char *jobvsl, char *jobvsr, int *n, double *a, int *lda, double *b, int *ldb, double * alphar, double *alphai, double *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *work, int *lwork, int *info);
int dgegv_check(char *jobvl, char *jobvr, int *n, double * a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
int dgehd2_check(int *n, int *ilo, int *ihi, double *a, int *lda, double *tau, double *work, int *info);
int dgehrd_check(int *n, int *ilo, int *ihi, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
int dgejsv_check(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, int *m, int *n, double *a, int *lda, double *sva, double *u, int *ldu, double *v, int *ldv, double *work, int *lwork, int *iwork, int *info);
int dgelq2_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int dgelqf_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int dgels_check(char *trans, int *m, int *n, int * nrhs, double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);
int dgelsd_check(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double * s, double *rcond, int *rank, double *work, int *lwork, int *iwork, int *info);
int dgelss_check(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double * s, double *rcond, int *rank, double *work, int *lwork, int *info);
int dgelsx_check(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int * jpvt, double *rcond, int *rank, double *work, int * info);
int dgelsy_check(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int * jpvt, double *rcond, int *rank, double *work, int * lwork, int *info);
int dgemqrt_check(char *side, char *trans, int *m, int *n, int *k, int *nb, double *v, int *ldv, double *t, int *ldt, double *c__, int *ldc, double *work, int *info);
int dgeql2_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int dgeqlf_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int dgeqp3_check(int *m, int *n, double *a, int * lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
int dgeqpf_check(int *m, int *n, double *a, int * lda, int *jpvt, double *tau, double *work, int *info);
int dgeqr2_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int dgeqr2p_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int dgeqrf_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int dgeqrfp_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int dgeqrt_check(int *m, int *n, int *nb, double * a, int *lda, double *t, int *ldt, double *work, int *info);
int dgeqrt2_check(int *m, int *n, double *a, int * lda, double *t, int *ldt, int *info);
int dgeqrt3_check(int *m, int *n, double *a, int * lda, double *t, int *ldt, int *info);
int dgerfs_check(char *trans, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int * ipiv, double *b, int *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *iwork, int *info);
int dgerfsx_check(char *trans, char *equed, int *n, int * nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, double *r__, double *c__, double *b, int *ldb, double *x, int *ldx, double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *params, double *work, int *iwork, int *info);
int dgerq2_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int dgerqf_check(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int dgesc2_check(int *n, double *a, int *lda, double *rhs, int *ipiv, int *jpiv, double *scale);
int dgesdd_check(char *jobz, int *m, int *n, double * a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *iwork, int *info);
int dgesv_check(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
int dgesvd_check(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int * ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);
int dgesvj_check(char *joba, char *jobu, char *jobv, int *m, int *n, double *a, int *lda, double *sva, int *mv, double *v, int *ldv, double *work, int *lwork, int *info);
int dgesvx_check(char *fact, char *trans, int *n, int * nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, char *equed, double *r__, double *c__, double *b, int *ldb, double *x, int *ldx, double * rcond, double *ferr, double *berr, double *work, int * iwork, int *info);
int dgesvxx_check(char *fact, char *trans, int *n, int * nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, char *equed, double *r__, double *c__, double *b, int *ldb, double *x, int *ldx, double * rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int * nparams, double *params, double *work, int *iwork, int *info);
int dgetc2_check(int *n, double *a, int *lda, int *ipiv, int *jpiv, int *info);
int dgetf2_check(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
int dgetrf_check(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
int dgetrfnp_check(int *m, int *n, double *a, int * lda, int *info);
int dgetri_check(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
int dgetrs_check(char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int * ldb, int *info);
int dggbak_check(char *job, char *side, int *n, int *ilo, int *ihi, double *lscale, double *rscale, int *m, double *v, int *ldv, int *info);
int dggbal_check(char *job, int *n, double *a, int * lda, double *b, int *ldb, int *ilo, int *ihi, double *lscale, double *rscale, double *work, int * info);
int dgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, int *n, double *a, int *lda, double *b, int *ldb, int *sdim, double *alphar, double *alphai, double *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *work, int *lwork, logical *bwork, int *info);
int dggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, int *n, double *a, int *lda, double *b, int *ldb, int *sdim, double *alphar, double *alphai, double *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *rconde, double * rcondv, double *work, int *lwork, int *iwork, int * liwork, logical *bwork, int *info);
int dggev_check(char *jobvl, char *jobvr, int *n, double * a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double *beta, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
int dggevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, double *a, int *lda, double *b, int *ldb, double *alphar, double *alphai, double * beta, double *vl, int *ldvl, double *vr, int *ldvr, int *ilo, int *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm, double *rconde, double * rcondv, double *work, int *lwork, int *iwork, logical * bwork, int *info);
int dggglm_check(int *n, int *m, int *p, double * a, int *lda, double *b, int *ldb, double *d__, double *x, double *y, double *work, int *lwork, int *info);
int dgghrd_check(char *compq, char *compz, int *n, int * ilo, int *ihi, double *a, int *lda, double *b, int *ldb, double *q, int *ldq, double *z__, int * ldz, int *info);
int dgglse_check(int *m, int *n, int *p, double * a, int *lda, double *b, int *ldb, double *c__, double *d__, double *x, double *work, int *lwork, int *info);
int dggqrf_check(int *n, int *m, int *p, double * a, int *lda, double *taua, double *b, int *ldb, double *taub, double *work, int *lwork, int *info);
int dggrqf_check(int *m, int *p, int *n, double * a, int *lda, double *taua, double *b, int *ldb, double *taub, double *work, int *lwork, int *info);
int dggsvd_check(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, double *a, int *lda, double *b, int *ldb, double *alpha, double *beta, double *u, int *ldu, double *v, int *ldv, double *q, int *ldq, double *work, int *iwork, int *info);
int dggsvp_check(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, double *a, int *lda, double *b, int *ldb, double *tola, double *tolb, int *k, int *l, double *u, int *ldu, double *v, int *ldv, double *q, int *ldq, int *iwork, double *tau, double *work, int *info);
int dgsvj0_check(char *jobv, int *m, int *n, double * a, int *lda, double *d__, double *sva, int *mv, double *v, int *ldv, double *eps, double *sfmin, double *tol, int *nsweep, double *work, int *lwork, int *info);
int dgsvj1_check(char *jobv, int *m, int *n, int *n1, double *a, int *lda, double *d__, double *sva, int *mv, double *v, int *ldv, double *eps, double *sfmin, double *tol, int *nsweep, double *work, int * lwork, int *info);
int dgtcon_check(char *norm, int *n, double *dl, double *d__, double *du, double *du2, int *ipiv, double *anorm, double *rcond, double *work, int * iwork, int *info);
int dgtrfs_check(char *trans, int *n, int *nrhs, double *dl, double *d__, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, int *ldb, double *x, int *ldx, double * ferr, double *berr, double *work, int *iwork, int * info);
int dgtsv_check(int *n, int *nrhs, double *dl, double *d__, double *du, double *b, int *ldb, int *info);
int dgtsvx_check(char *fact, char *trans, int *n, int * nrhs, double *dl, double *d__, double *du, double * dlf, double *df, double *duf, double *du2, int *ipiv, double *b, int *ldb, double *x, int *ldx, double * rcond, double *ferr, double *berr, double *work, int * iwork, int *info);
int dgttrf_check(int *n, double *dl, double *d__, double *du, double *du2, int *ipiv, int *info);
int dgttrs_check(char *trans, int *n, int *nrhs, double *dl, double *d__, double *du, double *du2, int *ipiv, double *b, int *ldb, int *info);
int dgtts2_check(int *itrans, int *n, int *nrhs, double *dl, double *d__, double *du, double *du2, int *ipiv, double *b, int *ldb);
int dhgeqz_check(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, double *h__, int *ldh, double *t, int *ldt, double *alphar, double *alphai, double * beta, double *q, int *ldq, double *z__, int *ldz, double *work, int *lwork, int *info);
int dhsein_check(char *side, char *eigsrc, char *initv, logical * select, int *n, double *h__, int *ldh, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, int *mm, int *m, double *work, int * ifaill, int *ifailr, int *info);
int dhseqr_check(char *job, char *compz, int *n, int *ilo, int *ihi, double *h__, int *ldh, double *wr, double *wi, double *z__, int *ldz, double *work, int *lwork, int *info);
logical disnan_check(double *din);
int dla_gbamv_check(int *trans, int *m, int *n, int *kl, int *ku, double *alpha, double *ab, int * ldab, double *x, int *incx, double *beta, double *y, int *incy);
double dla_gbrcond_check(char *trans, int *n, int *kl, int *ku, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, int *cmode, double *c__, int *info, double *work, int *iwork);
int dla_gbrfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, double *afb, int *ldafb, int *ipiv, logical *colequ, double *c__, double *b, int *ldb, double *y, int *ldy, double *berr_out__, int *n_norms__, double *err_bnds_norm__, double * err_bnds_comp__, double *res, double *ayb, double *dy, double *y_tail__, double *rcond, int *ithresh, double *rthresh, double *dz_ub__, logical *ignore_cwise__, int *info);
double dla_gbrpvgrw_check(int *n, int *kl, int *ku, int * ncols, double *ab, int *ldab, double *afb, int *ldafb);
int dla_geamv_check(int *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
double dla_gercond_check(char *trans, int *n, double *a, int *lda, double *af, int *ldaf, int *ipiv, int *cmode, double *c__, int *info, double *work, int *iwork);
int dla_gerfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf, int *ipiv, logical *colequ, double *c__, double *b, int *ldb, double *y, int * ldy, double *berr_out__, int *n_norms__, double *errs_n__, double *errs_c__, double *res, double *ayb, double * dy, double *y_tail__, double *rcond, int *ithresh, double *rthresh, double *dz_ub__, logical *ignore_cwise__, int *info);
double dla_gerpvgrw_check(int *n, int *ncols, double *a, int * lda, double *af, int *ldaf);
int dla_lin_berr_check(int *n, int *nz, int *nrhs, double *res, double *ayb, double *berr);
double dla_porcond_check(char *uplo, int *n, double *a, int *lda, double *af, int *ldaf, int *cmode, double *c__, int *info, double *work, int *iwork);
int dla_porfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, double *a, int *lda, double * af, int *ldaf, logical *colequ, double *c__, double *b, int *ldb, double *y, int *ldy, double *berr_out__, int *n_norms__, double *err_bnds_norm__, double * err_bnds_comp__, double *res, double *ayb, double *dy, double *y_tail__, double *rcond, int *ithresh, double *rthresh, double *dz_ub__, logical *ignore_cwise__, int *info);
double dla_porpvgrw_check(char *uplo, int *ncols, double *a, int * lda, double *af, int *ldaf, double *work);
int dla_syamv_check(int *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
double dla_syrcond_check(char *uplo, int *n, double *a, int *lda, double *af, int *ldaf, int *ipiv, int *cmode, double *c__, int *info, double *work, int *iwork);
int dla_syrfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, double *a, int *lda, double * af, int *ldaf, int *ipiv, logical *colequ, double *c__, double *b, int *ldb, double *y, int *ldy, double * berr_out__, int *n_norms__, double *err_bnds_norm__, double *err_bnds_comp__, double *res, double *ayb, double *dy, double *y_tail__, double *rcond, int * ithresh, double *rthresh, double *dz_ub__, logical * ignore_cwise__, int *info);
double dla_syrpvgrw_check(char *uplo, int *n, int *info, double * a, int *lda, double *af, int *ldaf, int *ipiv, double *work);
int dla_wwaddw_check(int *n, double *x, double *y, double *w);
int dlabad_check(double *small_, double *large);
int dlabrd_check(int *m, int *n, int *nb, double * a, int *lda, double *d__, double *e, double *tauq, double *taup, double *x, int *ldx, double *y, int *ldy);
int dlacn2_check(int *n, double *v, double *x, int *isgn, double *est, int *kase, int *isave);
int dlacon_check(int *n, double *v, double *x, int *isgn, double *est, int *kase);
int dlacpy_check(char *uplo, int *m, int *n, double * a, int *lda, double *b, int *ldb);
int dladiv_check(double *a, double *b, double *c__, double *d__, double *p, double *q);
int dlae2_check(double *a, double *b, double *c__, double *rt1, double *rt2);
int dlaebz_check(aocl_int64_t *ijob, aocl_int64_t *nitmax, aocl_int64_t *n, aocl_int64_t *mmax, aocl_int64_t *minp,
                 aocl_int64_t *nbmin, double *abstol, double *reltol, double *pivmin, double *d__,
                 double *e, double *e2, aocl_int64_t *nval, double *ab, double *c__, aocl_int64_t *mout,
                 aocl_int64_t *nab, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlaed0_check(aocl_int64_t *icompq, aocl_int64_t *qsiz, aocl_int64_t *n, double *d__, double *e, double *q,
                 aocl_int64_t *ldq, double *qstore, aocl_int64_t *ldqs, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dlaed1_check(aocl_int64_t *n, double *d__, double *q, aocl_int64_t *ldq, aocl_int64_t *indxq, double *rho,
                 aocl_int64_t *cutpnt, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlaed2_check(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *n1, double *d__, double *q, aocl_int64_t *ldq,
                 aocl_int64_t *indxq, double *rho, double *z__, double *dlamda, double *w, double *q2,
                 aocl_int64_t *indx, aocl_int64_t *indxc, aocl_int64_t *indxp, aocl_int64_t *coltyp, aocl_int64_t *info);
int dlaed3_check(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *n1, double *d__, double *q, aocl_int64_t *ldq,
                 double *rho, double *dlamda, double *q2, aocl_int64_t *indx, aocl_int64_t *ctot, double *w,
                 double *s, aocl_int64_t *info);
int dlaed4_check(aocl_int64_t *n, aocl_int64_t *i__, double *d__, double *z__, double *delta, double *rho,
                 double *dlam, aocl_int64_t *info);
int dlaed5_check(aocl_int64_t *i__, double *d__, double *z__, double *delta, double *rho, double *dlam);
int dlaed6_check(aocl_int64_t *kniter, logical *orgati, double *rho, double *d__, double *z__,
                 double *finit, double *tau, aocl_int64_t *info);
int dlaed7_check(aocl_int64_t *icompq, aocl_int64_t *n, aocl_int64_t *qsiz, aocl_int64_t *tlvls, aocl_int64_t *curlvl,
                 aocl_int64_t *curpbm, double *d__, double *q, aocl_int64_t *ldq, aocl_int64_t *indxq, double *rho,
                 aocl_int64_t *cutpnt, double *qstore, aocl_int64_t *qptr, aocl_int64_t *prmptr, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, double *givnum, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dlaed8_check(aocl_int64_t *icompq, aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *qsiz, double *d__, double *q,
                 aocl_int64_t *ldq, aocl_int64_t *indxq, double *rho, aocl_int64_t *cutpnt, double *z__,
                 double *dlamda, double *q2, aocl_int64_t *ldq2, double *w, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, double *givnum, aocl_int64_t *indxp, aocl_int64_t *indx,
                 aocl_int64_t *info);
int dlaed9_check(aocl_int64_t *k, aocl_int64_t *kstart, aocl_int64_t *kstop, aocl_int64_t *n, double *d__, double *q,
                 aocl_int64_t *ldq, double *rho, double *dlamda, double *w, double *s, aocl_int64_t *lds,
                 aocl_int64_t *info);
int dlaeda_check(aocl_int64_t *n, aocl_int64_t *tlvls, aocl_int64_t *curlvl, aocl_int64_t *curpbm, aocl_int64_t *prmptr,
                 aocl_int64_t *perm, aocl_int64_t *givptr, aocl_int64_t *givcol, double *givnum, double *q,
                 aocl_int64_t *qptr, double *z__, double *ztemp, aocl_int64_t *info);
int dlaein_check(logical *rightv, logical *noinit, aocl_int64_t *n, double *h__, aocl_int64_t *ldh,
                 double *wr, double *wi, double *vr, double *vi, double *b, aocl_int64_t *ldb,
                 double *work, double *eps3, double *smlnum, double *bignum, aocl_int64_t *info);
int dlaev2_check(double *a, double *b, double *c__, double *rt1, double *rt2, double *cs1,
                 double *sn1);
int dlaexc_check(logical *wantq, aocl_int64_t *n, double *t, aocl_int64_t *ldt, double *q, aocl_int64_t *ldq,
                 aocl_int64_t *j1, aocl_int64_t *n1, aocl_int64_t *n2, double *work, aocl_int64_t *info);
int dlag2_check(double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *safmin, double *scale1,
                double *scale2, double *wr1, double *wr2, double *wi);
int dlag2s_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, float *sa, aocl_int64_t *ldsa,
                 aocl_int64_t *info);
int dlags2_check(logical *upper, double *a1, double *a2, double *a3, double *b1, double *b2,
                 double *b3, double *csu, double *snu, double *csv, double *snv, double *csq,
                 double *snq);
int dlagtf_check(aocl_int64_t *n, double *a, double *lambda, double *b, double *c__, double *tol,
                 double *d__, aocl_int64_t *in, aocl_int64_t *info);
int dlagtm_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *alpha, double *dl, double *d__,
                 double *du, double *x, aocl_int64_t *ldx, double *beta, double *b, aocl_int64_t *ldb);
int dlagts_check(aocl_int64_t *job, aocl_int64_t *n, double *a, double *b, double *c__, double *d__,
                 aocl_int64_t *in, double *y, double *tol, aocl_int64_t *info);
int dlagv2_check(double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *alphar, double *alphai,
                 double *beta, double *csl, double *snl, double *csr, double *snr);
int dlahqr_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 double *h__, aocl_int64_t *ldh, double *wr, double *wi, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 double *z__, aocl_int64_t *ldz, aocl_int64_t *info);
int dlahr2_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, double *a, aocl_int64_t *lda, double *tau,
                 double *t, aocl_int64_t *ldt, double *y, aocl_int64_t *ldy);
int dlahrd_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, double *a, aocl_int64_t *lda, double *tau,
                 double *t, aocl_int64_t *ldt, double *y, aocl_int64_t *ldy);
int dlaic1_check(aocl_int64_t *job, aocl_int64_t *j, double *x, double *sest, double *w, double *gamma,
                 double *sestpr, double *s, double *c__);
logical dlaisnan_check(double *din1, double *din2);
int dlaln2_check(logical *ltrans, aocl_int64_t *na, aocl_int64_t *nw, double *smin, double *ca, double *a,
                 aocl_int64_t *lda, double *d1, double *d2, double *b, aocl_int64_t *ldb, double *wr,
                 double *wi, double *x, aocl_int64_t *ldx, double *scale, double *xnorm, aocl_int64_t *info);
int dlals0_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *nrhs, double *b,
                 aocl_int64_t *ldb, double *bx, aocl_int64_t *ldbx, aocl_int64_t *perm, aocl_int64_t *givptr,
                 aocl_int64_t *givcol, aocl_int64_t *ldgcol, double *givnum, aocl_int64_t *ldgnum, double *poles,
                 double *difl, double *difr, double *z__, aocl_int64_t *k, double *c__, double *s,
                 double *work, aocl_int64_t *info);
int dlalsa_check(aocl_int64_t *icompq, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, double *b,
                 aocl_int64_t *ldb, double *bx, aocl_int64_t *ldbx, double *u, aocl_int64_t *ldu, double *vt,
                 aocl_int64_t *k, double *difl, double *difr, double *z__, double *poles,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol, aocl_int64_t *perm, double *givnum,
                 double *c__, double *s, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlalsd_check(char *uplo, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, double *e,
                 double *b, aocl_int64_t *ldb, double *rcond, aocl_int64_t *rank, double *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dlamrg_check(aocl_int64_t *n1, aocl_int64_t *n2, double *a, aocl_int64_t *dtrd1, aocl_int64_t *dtrd2,
                 aocl_int64_t *index);
int dlaneg_check(aocl_int64_t *n, double *d__, double *lld, double *sigma, double *pivmin, aocl_int64_t *r__);
double dlangb_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab, aocl_int64_t *ldab,
                    double *work);
double dlange_check(char *norm, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *work);
double dlangt_check(char *norm, aocl_int64_t *n, double *dl, double *d__, double *du);
double dlanhs_check(char *norm, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *work);
double dlansb_check(char *norm, char *uplo, aocl_int64_t *n, aocl_int64_t *k, double *ab, aocl_int64_t *ldab,
                    double *work);
double dlansf_check(char *norm, char *transr, char *uplo, aocl_int64_t *n, double *a, double *work);
double dlansp_check(char *norm, char *uplo, aocl_int64_t *n, double *ap, double *work);
double dlanst_check(char *norm, aocl_int64_t *n, double *d__, double *e);
double dlansy_check(char *norm, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *work);
double dlantb_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *k, double *ab,
                    aocl_int64_t *ldab, double *work);
double dlantp_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, double *ap, double *work);
double dlantr_check(char *norm, char *uplo, char *diag, aocl_int64_t *m, aocl_int64_t *n, double *a,
                    aocl_int64_t *lda, double *work);
int dlanv2_check(double *a, double *b, double *c__, double *d__, double *rt1r, double *rt1i,
                 double *rt2r, double *rt2i, double *cs, double *sn);
int dlapll_check(aocl_int64_t *n, double *x, aocl_int64_t *incx, double *y, aocl_int64_t *incy, double *ssmin);
int dlapmr_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, double *x, aocl_int64_t *ldx, aocl_int64_t *k);
int dlapmt_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, double *x, aocl_int64_t *ldx, aocl_int64_t *k);
double dlapy2_check(double *x, double *y);
double dlapy3_check(double *x, double *y, double *z__);
int dlaqgb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab, aocl_int64_t *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 char *equed);
int dlaqge_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, char *equed);
int dlaqp2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, double *a, aocl_int64_t *lda, aocl_int64_t *jpvt,
                 double *tau, double *vn1, double *vn2, double *work);
int dlaqps_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, aocl_int64_t *nb, aocl_int64_t *kb, double *a,
                 aocl_int64_t *lda, aocl_int64_t *jpvt, double *tau, double *vn1, double *vn2, double *auxv,
                 double *f, aocl_int64_t *ldf);
int dlaqr0_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 double *h__, aocl_int64_t *ldh, double *wr, double *wi, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dlaqr1_check(aocl_int64_t *n, double *h__, aocl_int64_t *ldh, double *sr1, double *si1, double *sr2,
                 double *si2, double *v);
int dlaqr2_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, double *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz, double *z__,
                 aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, double *sr, double *si, double *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, double *t, aocl_int64_t *ldt, aocl_int64_t *nv, double *wv,
                 aocl_int64_t *ldwv, double *work, aocl_int64_t *lwork);
int dlaqr3_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, double *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz, double *z__,
                 aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, double *sr, double *si, double *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, double *t, aocl_int64_t *ldt, aocl_int64_t *nv, double *wv,
                 aocl_int64_t *ldwv, double *work, aocl_int64_t *lwork);
int dlaqr4_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 double *h__, aocl_int64_t *ldh, double *wr, double *wi, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dlaqr5_check(logical *wantt, logical *wantz, aocl_int64_t *kacc22, aocl_int64_t *n, aocl_int64_t *ktop,
                 aocl_int64_t *kbot, aocl_int64_t *nshfts, double *sr, double *si, double *h__, aocl_int64_t *ldh,
                 aocl_int64_t *iloz, aocl_int64_t *ihiz, double *z__, aocl_int64_t *ldz, double *v, aocl_int64_t *ldv,
                 double *u, aocl_int64_t *ldu, aocl_int64_t *nv, double *wv, aocl_int64_t *ldwv, aocl_int64_t *nh,
                 double *wh, aocl_int64_t *ldwh);
int dlaqsb_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab, double *s,
                 double *scond, double *amax, char *equed);
int dlaqsp_check(char *uplo, aocl_int64_t *n, double *ap, double *s, double *scond, double *amax,
                 char *equed);
int dlaqsy_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *s, double *scond,
                 double *amax, char *equed);
int dlaqtr_check(logical *ltran, logical *lfloat, aocl_int64_t *n, double *t, aocl_int64_t *ldt, double *b,
                 double *w, double *scale, double *x, double *work, aocl_int64_t *info);
int dlar1v_check(aocl_int64_t *n, aocl_int64_t *b1, aocl_int64_t *bn, double *lambda, double *d__, double *l,
                 double *ld, double *lld, double *pivmin, double *gaptol, double *z__,
                 logical *wantnc, aocl_int64_t *negcnt, double *ztz, double *mingma, aocl_int64_t *r__,
                 aocl_int64_t *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);
int dlar2v_check(aocl_int64_t *n, double *x, double *y, double *z__, aocl_int64_t *incx, double *c__,
                 double *s, aocl_int64_t *incc);
int dlarf_check(char *side, aocl_int64_t *m, aocl_int64_t *n, double *v, aocl_int64_t *incv, double *tau,
                double *c__, aocl_int64_t *ldc, double *work);
int dlarfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, double *v, aocl_int64_t *ldv, double *t, aocl_int64_t *ldt, double *c__,
                 aocl_int64_t *ldc, double *work, aocl_int64_t *ldwork);
int dlarfg_check(aocl_int64_t *n, double *alpha, double *x, aocl_int64_t *incx, double *tau);
int dlarfgp_check(aocl_int64_t *n, double *alpha, double *x, aocl_int64_t *incx, double *tau);
int dlarft_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, double *v, aocl_int64_t *ldv,
                 double *tau, double *t, aocl_int64_t *ldt);
int dlarfx_check(char *side, aocl_int64_t *m, aocl_int64_t *n, double *v, double *tau, double *c__,
                 aocl_int64_t *ldc, double *work);
int dlargv_check(aocl_int64_t *n, double *x, aocl_int64_t *incx, double *y, aocl_int64_t *incy, double *c__,
                 aocl_int64_t *incc);
int dlarnv_check(aocl_int64_t *idist, aocl_int64_t *iseed, aocl_int64_t *n, double *x);
int dlarra_check(aocl_int64_t *n, double *d__, double *e, double *e2, double *spltol, double *tnrm,
                 aocl_int64_t *nsplit, aocl_int64_t *isplit, aocl_int64_t *info);
int dlarrb_check(aocl_int64_t *n, double *d__, double *lld, aocl_int64_t *ifirst, aocl_int64_t *ilast,
                 double *rtol1, double *rtol2, aocl_int64_t *offset, double *w, double *wgap,
                 double *werr, double *work, aocl_int64_t *iwork, double *pivmin, double *spdiam,
                 aocl_int64_t *twist, aocl_int64_t *info);
int dlarrc_check(char *jobt, aocl_int64_t *n, double *vl, double *vu, double *d__, double *e,
                 double *pivmin, aocl_int64_t *eigcnt, aocl_int64_t *lcnt, aocl_int64_t *rcnt, aocl_int64_t *info);
int dlarrd_check(char *range, char *order, aocl_int64_t *n, double *vl, double *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, double *gers, double *reltol, double *d__, double *e, double *e2,
                 double *pivmin, aocl_int64_t *nsplit, aocl_int64_t *isplit, aocl_int64_t *m, double *w,
                 double *werr, double *wl, double *wu, aocl_int64_t *iblock, aocl_int64_t *indexw,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlarre_check(char *range, aocl_int64_t *n, double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu,
                 double *d__, double *e, double *e2, double *rtol1, double *rtol2, double *spltol,
                 aocl_int64_t *nsplit, aocl_int64_t *isplit, aocl_int64_t *m, double *w, double *werr,
                 double *wgap, aocl_int64_t *iblock, aocl_int64_t *indexw, double *gers, double *pivmin,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlarrf_check(aocl_int64_t *n, double *d__, double *l, double *ld, aocl_int64_t *clstrt, aocl_int64_t *clend,
                 double *w, double *wgap, double *werr, double *spdiam, double *clgapl,
                 double *clgapr, double *pivmin, double *sigma, double *dplus, double *lplus,
                 double *work, aocl_int64_t *info);
int dlarrj_check(aocl_int64_t *n, double *d__, double *e2, aocl_int64_t *ifirst, aocl_int64_t *ilast, double *rtol,
                 aocl_int64_t *offset, double *w, double *werr, double *work, aocl_int64_t *iwork,
                 double *pivmin, double *spdiam, aocl_int64_t *info);
int dlarrk_check(aocl_int64_t *n, aocl_int64_t *iw, double *gl, double *gu, double *d__, double *e2,
                 double *pivmin, double *reltol, double *w, double *werr, aocl_int64_t *info);
int dlarrr_check(aocl_int64_t *n, double *d__, double *e, aocl_int64_t *info);
int dlarrv_check(aocl_int64_t *n, double *vl, double *vu, double *d__, double *l, double *pivmin,
                 aocl_int64_t *isplit, aocl_int64_t *m, aocl_int64_t *dol, aocl_int64_t *dou, double *minrgp,
                 double *rtol1, double *rtol2, double *w, double *werr, double *wgap,
                 aocl_int64_t *iblock, aocl_int64_t *indexw, double *gers, double *z__, aocl_int64_t *ldz,
                 aocl_int64_t *isuppz, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlarscl2_check(aocl_int64_t *m, aocl_int64_t *n, double *d__, double *x, aocl_int64_t *ldx);
int dlartg_check(double *f, double *g, double *cs, double *sn, double *r__);
int dlartgp_check(double *f, double *g, double *cs, double *sn, double *r__);
int dlartgs_check(double *x, double *y, double *sigma, double *cs, double *sn);
int dlartv_check(aocl_int64_t *n, double *x, aocl_int64_t *incx, double *y, aocl_int64_t *incy, double *c__,
                 double *s, aocl_int64_t *incc);
int dlaruv_check(aocl_int64_t *iseed, aocl_int64_t *n, double *x);
int dlarz_check(char *side, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, double *v, aocl_int64_t *incv,
                double *tau, double *c__, aocl_int64_t *ldc, double *work);
int dlarzb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, double *v, aocl_int64_t *ldv, double *t, aocl_int64_t *ldt,
                 double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *ldwork);
int dlarzt_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, double *v, aocl_int64_t *ldv,
                 double *tau, double *t, aocl_int64_t *ldt);
int dlas2_check(double *f, double *g, double *h__, double *ssmin, double *ssmax);
int dlascl_check(char *type__, aocl_int64_t *kl, aocl_int64_t *ku, double *cfrom, double *cto, aocl_int64_t *m,
                 aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dlascl2_check(aocl_int64_t *m, aocl_int64_t *n, double *d__, double *x, aocl_int64_t *ldx);
int dlasd0_check(aocl_int64_t *n, aocl_int64_t *sqre, double *d__, double *e, double *u, aocl_int64_t *ldu,
                 double *vt, aocl_int64_t *ldvt, aocl_int64_t *smlsiz, aocl_int64_t *iwork, double *work,
                 aocl_int64_t *info);
int dlasd1_check(aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, double *d__, double *alpha, double *beta,
                 double *u, aocl_int64_t *ldu, double *vt, aocl_int64_t *ldvt, aocl_int64_t *idxq, aocl_int64_t *iwork,
                 double *work, aocl_int64_t *info);
int dlasd2_check(aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *k, double *d__, double *z__,
                 double *alpha, double *beta, double *u, aocl_int64_t *ldu, double *vt, aocl_int64_t *ldvt,
                 double *dsigma, double *u2, aocl_int64_t *ldu2, double *vt2, aocl_int64_t *ldvt2,
                 aocl_int64_t *idxp, aocl_int64_t *idx, aocl_int64_t *idxc, aocl_int64_t *idxq, aocl_int64_t *coltyp,
                 aocl_int64_t *info);
int dlasd3_check(aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *k, double *d__, double *q,
                 aocl_int64_t *ldq, double *dsigma, double *u, aocl_int64_t *ldu, double *u2, aocl_int64_t *ldu2,
                 double *vt, aocl_int64_t *ldvt, double *vt2, aocl_int64_t *ldvt2, aocl_int64_t *idxc,
                 aocl_int64_t *ctot, double *z__, aocl_int64_t *info);
int dlasd4_check(aocl_int64_t *n, aocl_int64_t *i__, double *d__, double *z__, double *delta, double *rho,
                 double *sigma, double *work, aocl_int64_t *info);
int dlasd5_check(aocl_int64_t *i__, double *d__, double *z__, double *delta, double *rho, double *dsigma,
                 double *work);
int dlasd6_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, double *d__, double *vf,
                 double *vl, double *alpha, double *beta, aocl_int64_t *idxq, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol, double *givnum, aocl_int64_t *ldgnum,
                 double *poles, double *difl, double *difr, double *z__, aocl_int64_t *k, double *c__,
                 double *s, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlasd7_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *k, double *d__,
                 double *z__, double *zw, double *vf, double *vfw, double *vl, double *vlw,
                 double *alpha, double *beta, double *dsigma, aocl_int64_t *idx, aocl_int64_t *idxp,
                 aocl_int64_t *idxq, aocl_int64_t *perm, aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol,
                 double *givnum, aocl_int64_t *ldgnum, double *c__, double *s, aocl_int64_t *info);
int dlasd8_check(aocl_int64_t *icompq, aocl_int64_t *k, double *d__, double *z__, double *vf, double *vl,
                 double *difl, double *difr, aocl_int64_t *lddifr, double *dsigma, double *work,
                 aocl_int64_t *info);
int dlasda_check(aocl_int64_t *icompq, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *sqre, double *d__,
                 double *e, double *u, aocl_int64_t *ldu, double *vt, aocl_int64_t *k, double *difl,
                 double *difr, double *z__, double *poles, aocl_int64_t *givptr, aocl_int64_t *givcol,
                 aocl_int64_t *ldgcol, aocl_int64_t *perm, double *givnum, double *c__, double *s,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dlasdq_check(char *uplo, aocl_int64_t *sqre, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru, aocl_int64_t *ncc,
                 double *d__, double *e, double *vt, aocl_int64_t *ldvt, double *u, aocl_int64_t *ldu,
                 double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dlasdt_check(aocl_int64_t *n, aocl_int64_t *lvl, aocl_int64_t *nd, aocl_int64_t *inode, aocl_int64_t *ndiml,
                 aocl_int64_t *ndimr, aocl_int64_t *msub);
int dlaset_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, double *alpha, double *beta, double *a,
                 aocl_int64_t *lda);
int dlasq1_check(aocl_int64_t *n, double *d__, double *e, double *work, aocl_int64_t *info);
int dlasq2_check(aocl_int64_t *n, double *z__, aocl_int64_t *info);
int dlasq3_check(aocl_int64_t *i0, aocl_int64_t *n0, double *z__, aocl_int64_t *pp, double *dmin__, double *sigma,
                 double *desig, double *qmax, aocl_int64_t *nfail, aocl_int64_t *iter, aocl_int64_t *ndiv,
                 logical *ieee, aocl_int64_t *ttype, double *dmin1, double *dmin2, double *dn,
                 double *dn1, double *dn2, double *g, double *tau);
int dlasq4_check(aocl_int64_t *i0, aocl_int64_t *n0, double *z__, aocl_int64_t *pp, aocl_int64_t *n0in, double *dmin__,
                 double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau,
                 aocl_int64_t *ttype, double *g);
int dlasq5_check(aocl_int64_t *i0, aocl_int64_t *n0, double *z__, aocl_int64_t *pp, double *tau, double *sigma,
                 double *dmin__, double *dmin1, double *dmin2, double *dn, double *dnm1,
                 double *dnm2, logical *ieee, double *eps);
int dlasq6_check(aocl_int64_t *i0, aocl_int64_t *n0, double *z__, aocl_int64_t *pp, double *dmin__, double *dmin1,
                 double *dmin2, double *dn, double *dnm1, double *dnm2);
int dlasr_check(char *side, char *pivot, char *direct, aocl_int64_t *m, aocl_int64_t *n, double *c__,
                double *s, double *a, aocl_int64_t *lda);
int dlasrt_check(char *id, aocl_int64_t *n, double *d__, aocl_int64_t *info);
int dlassq_check(aocl_int64_t *n, double *x, aocl_int64_t *incx, double *scale, double *sumsq);
int dlasv2_check(double *f, double *g, double *h__, double *ssmin, double *ssmax, double *snr,
                 double *csr, double *snl, double *csl);
int dlaswp_check(aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *k1, aocl_int64_t *k2, aocl_int_t *ipiv,
                 aocl_int64_t *incx);
int dlasy2_check(logical *ltranl, logical *ltranr, aocl_int64_t *isgn, aocl_int64_t *n1, aocl_int64_t *n2,
                 double *tl, aocl_int64_t *ldtl, double *tr, aocl_int64_t *ldtr, double *b, aocl_int64_t *ldb,
                 double *scale, double *x, aocl_int64_t *ldx, double *xnorm, aocl_int64_t *info);
int dlasyf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, double *a, aocl_int64_t *lda,
                 aocl_int_t *ipiv, double *w, aocl_int64_t *ldw, aocl_int64_t *info);
int dlasyf_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, double *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, double *w, aocl_int64_t *ldw, aocl_int64_t *info);
int dlat2s_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, float *sa, aocl_int64_t *ldsa,
                 aocl_int64_t *info);
int dlatbs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, aocl_int64_t *kd,
                 double *ab, aocl_int64_t *ldab, double *x, double *scale, double *cnorm, aocl_int64_t *info);
int dlatdf_check(aocl_int64_t *ijob, aocl_int64_t *n, double *z__, aocl_int64_t *ldz, double *rhs, double *rdsum,
                 double *rdscal, aocl_int_t *ipiv, aocl_int64_t *jpiv);
int dlatps_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, double *ap,
                 double *x, double *scale, double *cnorm, aocl_int64_t *info);
int dlatrd_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, double *a, aocl_int64_t *lda, double *e,
                 double *tau, double *w, aocl_int64_t *ldw);
int dlatrs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *x, double *scale, double *cnorm, aocl_int64_t *info);
int dlatrz_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, double *a, aocl_int64_t *lda, double *tau,
                 double *work);
int dlatzm_check(char *side, aocl_int64_t *m, aocl_int64_t *n, double *v, aocl_int64_t *incv, double *tau,
                 double *c1, double *c2, aocl_int64_t *ldc, double *work);
int dlauu2_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dlauum_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dopgtr_check(char *uplo, aocl_int64_t *n, double *ap, double *tau, double *q, aocl_int64_t *ldq,
                 double *work, aocl_int64_t *info);
int dopmtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, double *ap,
                 double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dorbdb_check(char *trans, char *signs, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, double *x11,
                 aocl_int64_t *ldx11, double *x12, aocl_int64_t *ldx12, double *x21, aocl_int64_t *ldx21,
                 double *x22, aocl_int64_t *ldx22, double *theta, double *phi, double *taup1,
                 double *taup2, double *tauq1, double *tauq2, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dorbdb1_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, double *x11, aocl_int64_t *ldx11, double *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorbdb2_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, double *x11, aocl_int64_t *ldx11, double *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorbdb3_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, double *x11, aocl_int64_t *ldx11, double *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorbdb4_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, double *x11, aocl_int64_t *ldx11, double *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *phantom, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorbdb5_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, double *x1, aocl_int64_t *incx1, double *x2,
                  aocl_int64_t *incx2, double *q1, aocl_int64_t *ldq1, double *q2, aocl_int64_t *ldq2,
                  double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorbdb6_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, double *x1, aocl_int64_t *incx1, double *x2,
                  aocl_int64_t *incx2, double *q1, aocl_int64_t *ldq1, double *q2, aocl_int64_t *ldq2,
                  double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, double *x11, aocl_int64_t *ldx11, double *x12,
                 aocl_int64_t *ldx12, double *x21, aocl_int64_t *ldx21, double *x22, aocl_int64_t *ldx22,
                 double *theta, double *u1, aocl_int64_t *ldu1, double *u2, aocl_int64_t *ldu2, double *v1t,
                 aocl_int64_t *ldv1t, double *v2t, aocl_int64_t *ldv2t, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dorcsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q,
                     double *x11, aocl_int64_t *ldx11, double *x21, aocl_int64_t *ldx21, double *theta,
                     double *u1, aocl_int64_t *ldu1, double *u2, aocl_int64_t *ldu2, double *v1t,
                     aocl_int64_t *ldv1t, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int dorg2l_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *info);
int dorg2r_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *info);
int dorgbr_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda,
                 double *tau, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorghr_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorgl2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *info);
int dorglq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorgql_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorgqr_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorgr2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *info);
int dorgrq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dorgtr_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dorm2l_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dorm2r_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dormbr_check(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dormhr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 double *a, aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dorml2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dormlq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dormql_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dormqr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dormr2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dormr3_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dormrq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dormrz_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dormtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dpbcon_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab, double *anorm,
                 double *rcond, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dpbequ_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab, double *s,
                 double *scond, double *amax, aocl_int64_t *info);
int dpbrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, double *ab, aocl_int64_t *ldab,
                 double *afb, aocl_int64_t *ldafb, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dpbstf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int dpbsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, double *ab, aocl_int64_t *ldab,
                double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dpbsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, double *ab,
                 aocl_int64_t *ldab, double *afb, aocl_int64_t *ldafb, char *equed, double *s, double *b,
                 aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond, double *ferr, double *berr,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dpbtf2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int dpbtrf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int dpbtrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, double *ab, aocl_int64_t *ldab,
                 double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dpftrf_check(char *transr, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *info);
int dpftri_check(char *transr, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *info);
int dpftrs_check(char *transr, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, double *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int dpocon_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *anorm, double *rcond,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dpoequ_check(aocl_int64_t *n, double *a, aocl_int64_t *lda, double *s, double *scond, double *amax,
                 aocl_int64_t *info);
int dpoequb_check(aocl_int64_t *n, double *a, aocl_int64_t *lda, double *s, double *scond, double *amax,
                  aocl_int64_t *info);
int dporfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *af,
                 aocl_int64_t *ldaf, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *ferr,
                 double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dporfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                  double *af, aocl_int64_t *ldaf, double *s, double *b, aocl_int64_t *ldb, double *x,
                  aocl_int64_t *ldx, double *rcond, double *berr, aocl_int64_t *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, aocl_int64_t *nparams,
                  double *params, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dposv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int dposvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                 double *af, aocl_int64_t *ldaf, char *equed, double *s, double *b, aocl_int64_t *ldb,
                 double *x, aocl_int64_t *ldx, double *rcond, double *ferr, double *berr, double *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dposvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                  double *af, aocl_int64_t *ldaf, char *equed, double *s, double *b, aocl_int64_t *ldb,
                  double *x, aocl_int64_t *ldx, double *rcond, double *rpvgrw, double *berr,
                  aocl_int64_t *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  aocl_int64_t *nparams, double *params, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dpotf2_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dpotrf_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dpotri_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dpotrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int dppcon_check(char *uplo, aocl_int64_t *n, double *ap, double *anorm, double *rcond, double *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dppequ_check(char *uplo, aocl_int64_t *n, double *ap, double *s, double *scond, double *amax,
                 aocl_int64_t *info);
int dpprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, double *afp, double *b,
                 aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *ferr, double *berr, double *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dppsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, double *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int dppsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, double *afp,
                 char *equed, double *s, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx,
                 double *rcond, double *ferr, double *berr, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dpptrf_check(char *uplo, aocl_int64_t *n, double *ap, aocl_int64_t *info);
int dpptri_check(char *uplo, aocl_int64_t *n, double *ap, aocl_int64_t *info);
int dpptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, double *b, aocl_int64_t *ldb,
                 aocl_int64_t *info);
int dpstf2_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 double *tol, double *work, aocl_int64_t *info);
int dpstrf_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 double *tol, double *work, aocl_int64_t *info);
int dptcon_check(aocl_int64_t *n, double *d__, double *e, double *anorm, double *rcond, double *work,
                 aocl_int64_t *info);
int dpteqr_check(char *compz, aocl_int64_t *n, double *d__, double *e, double *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *info);
int dptrfs_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, double *e, double *df, double *ef,
                 double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 double *work, aocl_int64_t *info);
int dptsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, double *e, double *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int dptsvx_check(char *fact, aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, double *e, double *df,
                 double *ef, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond,
                 double *ferr, double *berr, double *work, aocl_int64_t *info);
int dpttrf_check(aocl_int64_t *n, double *d__, double *e, aocl_int64_t *info);
int dpttrs_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, double *e, double *b, aocl_int64_t *ldb,
                 aocl_int64_t *info);
int dptts2_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, double *e, double *b, aocl_int64_t *ldb);
int drscl_check(aocl_int64_t *n, double *sa, double *sx, aocl_int64_t *incx);
int dsbev_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab,
                double *w, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *info);
int dsbevd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab,
                 double *w, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int dsbevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab,
                 aocl_int64_t *ldab, double *q, aocl_int64_t *ldq, double *vl, double *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w, double *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int dsbgst_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, double *ab,
                 aocl_int64_t *ldab, double *bb, aocl_int64_t *ldbb, double *x, aocl_int64_t *ldx, double *work,
                 aocl_int64_t *info);
int dsbgv_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, double *ab,
                aocl_int64_t *ldab, double *bb, aocl_int64_t *ldbb, double *w, double *z__, aocl_int64_t *ldz,
                double *work, aocl_int64_t *info);
int dsbgvd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, double *ab,
                 aocl_int64_t *ldab, double *bb, aocl_int64_t *ldbb, double *w, double *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dsbgvx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb,
                 double *ab, aocl_int64_t *ldab, double *bb, aocl_int64_t *ldbb, double *q, aocl_int64_t *ldq,
                 double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m,
                 double *w, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *iwork, aocl_int64_t *ifail,
                 aocl_int64_t *info);
int dsbtrd_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, double *ab, aocl_int64_t *ldab,
                 double *d__, double *e, double *q, aocl_int64_t *ldq, double *work, aocl_int64_t *info);
int dsfrk_check(char *transr, char *uplo, char *trans, aocl_int64_t *n, aocl_int64_t *k, double *alpha,
                double *a, aocl_int64_t *lda, double *beta, double *c__);
int dsgesv_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *b,
                 aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *work, float *swork, aocl_int64_t *iter,
                 aocl_int64_t *info);
int dspcon_check(char *uplo, aocl_int64_t *n, double *ap, aocl_int_t *ipiv, double *anorm, double *rcond,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dspev_check(char *jobz, char *uplo, aocl_int64_t *n, double *ap, double *w, double *z__,
                aocl_int64_t *ldz, double *work, aocl_int64_t *info);
int dspevd_check(char *jobz, char *uplo, aocl_int64_t *n, double *ap, double *w, double *z__,
                 aocl_int64_t *ldz, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork,
                 aocl_int64_t *info);
int dspevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, double *ap, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w,
                 double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *iwork, aocl_int64_t *ifail,
                 aocl_int64_t *info);
int dspgst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, double *ap, double *bp, aocl_int64_t *info);
int dspgv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, double *ap, double *bp,
                double *w, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *info);
int dspgvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, double *ap, double *bp,
                 double *w, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int dspgvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, double *ap,
                 double *bp, double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol,
                 aocl_int64_t *m, double *w, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int dsposv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *work, float *swork, aocl_int64_t *iter,
                 aocl_int64_t *info);
int dsprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, double *afp, aocl_int_t *ipiv,
                 double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dspsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, aocl_int_t *ipiv, double *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int dspsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, double *afp,
                 aocl_int_t *ipiv, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond,
                 double *ferr, double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dsptrd_check(char *uplo, aocl_int64_t *n, double *ap, double *d__, double *e, double *tau,
                 aocl_int64_t *info);
int dsptrf_check(char *uplo, aocl_int64_t *n, double *ap, aocl_int_t *ipiv, aocl_int64_t *info);
int dsptri_check(char *uplo, aocl_int64_t *n, double *ap, aocl_int_t *ipiv, double *work, aocl_int64_t *info);
int dsptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap, aocl_int_t *ipiv, double *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int dstebz_check(char *range, char *order, aocl_int64_t *n, double *vl, double *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, double *abstol, double *d__, double *e, aocl_int64_t *m, aocl_int64_t *nsplit,
                 double *w, aocl_int64_t *iblock, aocl_int64_t *isplit, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dstedc_check(char *compz, aocl_int64_t *n, double *d__, double *e, double *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dstegr_check(char *jobz, char *range, aocl_int64_t *n, double *d__, double *e, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w,
                 double *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dstein_check(aocl_int64_t *n, double *d__, double *e, aocl_int64_t *m, double *w, aocl_int64_t *iblock,
                 aocl_int64_t *isplit, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int dstemr_check(char *jobz, char *range, aocl_int64_t *n, double *d__, double *e, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, aocl_int64_t *m, double *w, double *z__,
                 aocl_int64_t *ldz, aocl_int64_t *nzc, aocl_int64_t *isuppz, logical *tryrac, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dsteqr_check(char *compz, aocl_int64_t *n, double *d__, double *e, double *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *info);
int dsterf_check(aocl_int64_t *n, double *d__, double *e, aocl_int64_t *info);
int dstev_check(char *jobz, aocl_int64_t *n, double *d__, double *e, double *z__, aocl_int64_t *ldz,
                double *work, aocl_int64_t *info);
int dstevd_check(char *jobz, aocl_int64_t *n, double *d__, double *e, double *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dstevr_check(char *jobz, char *range, aocl_int64_t *n, double *d__, double *e, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w,
                 double *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dstevx_check(char *jobz, char *range, aocl_int64_t *n, double *d__, double *e, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w,
                 double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *iwork, aocl_int64_t *ifail,
                 aocl_int64_t *info);
int dsycon_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *anorm,
                 double *rcond, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dsycon_rook_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *anorm,
                      double *rcond, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dsyconv_check(char *uplo, char *way, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  double *work, aocl_int64_t *info);
int dsyequb_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *s, double *scond,
                  double *amax, double *work, aocl_int64_t *info);
int dsyev_check(char *jobz, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *w,
                double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dsyevd_check(char *jobz, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *w,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dsyevr_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m,
                 double *w, double *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int dsyevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m,
                 double *w, double *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int dsygs2_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int dsygst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int dsygv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                double *b, aocl_int64_t *ldb, double *w, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dsygvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *b, aocl_int64_t *ldb, double *w, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int dsygvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *vl, double *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w, double *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int dsyrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dsyrfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                  double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *s, double *b, aocl_int64_t *ldb,
                  double *x, aocl_int64_t *ldx, double *rcond, double *berr, aocl_int64_t *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, aocl_int64_t *nparams,
                  double *params, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dsysv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                double *b, aocl_int64_t *ldb, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dsysv_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                     double *b, aocl_int64_t *ldb, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dsysvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                 double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *b, aocl_int64_t *ldb, double *x,
                 aocl_int64_t *ldx, double *rcond, double *ferr, double *berr, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int dsysvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                  double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, double *s, double *b,
                  aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond, double *rpvgrw,
                  double *berr, aocl_int64_t *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, aocl_int64_t *nparams, double *params, double *work,
                  aocl_int64_t *iwork, aocl_int64_t *info);
int dsyswapr_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *i1, aocl_int64_t *i2);
int dsytd2_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *d__, double *e,
                 double *tau, aocl_int64_t *info);
int dsytf2_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int dsytf2_rook_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      aocl_int64_t *info);
int dsytrd_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *d__, double *e,
                 double *tau, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dsytrf_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dsytrf_rook_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *work,
                      aocl_int64_t *lwork, aocl_int64_t *info);
int dsytri_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *work,
                 aocl_int64_t *info);
int dsytri2_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int dsytri2x_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *work,
                   aocl_int64_t *nb, aocl_int64_t *info);
int dsytri_rook_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *work,
                      aocl_int64_t *info);
int dsytrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dsytrs2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  double *b, aocl_int64_t *ldb, double *work, aocl_int64_t *info);
int dsytrs_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dtbcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *kd, double *ab,
                 aocl_int64_t *ldab, double *rcond, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dtbrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 double *ab, aocl_int64_t *ldab, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dtbtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 double *ab, aocl_int64_t *ldab, double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dtfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, aocl_int64_t *m,
                aocl_int64_t *n, double *alpha, double *a, double *b, aocl_int64_t *ldb);
int dtftri_check(char *transr, char *uplo, char *diag, aocl_int64_t *n, double *a, aocl_int64_t *info);
int dtfttp_check(char *transr, char *uplo, aocl_int64_t *n, double *arf, double *ap, aocl_int64_t *info);
int dtfttr_check(char *transr, char *uplo, aocl_int64_t *n, double *arf, double *a, aocl_int64_t *lda,
                 aocl_int64_t *info);
int dtgevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, double *s, aocl_int64_t *lds,
                 double *p, aocl_int64_t *ldp, double *vl, aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr,
                 aocl_int64_t *mm, aocl_int64_t *m, double *work, aocl_int64_t *info);
int dtgex2_check(logical *wantq, logical *wantz, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, double *q, aocl_int64_t *ldq, double *z__, aocl_int64_t *ldz, aocl_int64_t *j1,
                 aocl_int64_t *n1, aocl_int64_t *n2, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dtgexc_check(logical *wantq, logical *wantz, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, double *q, aocl_int64_t *ldq, double *z__, aocl_int64_t *ldz, aocl_int64_t *ifst,
                 aocl_int64_t *ilst, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dtgsen_check(aocl_int64_t *ijob, logical *wantq, logical *wantz, logical *select, aocl_int64_t *n,
                 double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *alphar, double *alphai,
                 double *beta, double *q, aocl_int64_t *ldq, double *z__, aocl_int64_t *ldz, aocl_int64_t *m,
                 double *pl, double *pr, double *dif, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int dtgsja_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, aocl_int64_t *k,
                 aocl_int64_t *l, double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *tola,
                 double *tolb, double *alpha, double *beta, double *u, aocl_int64_t *ldu, double *v,
                 aocl_int64_t *ldv, double *q, aocl_int64_t *ldq, double *work, aocl_int64_t *ncycle,
                 aocl_int64_t *info);
int dtgsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *b, aocl_int64_t *ldb, double *vl, aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr,
                 double *s, double *dif, aocl_int64_t *mm, aocl_int64_t *m, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dtgsy2_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *b, aocl_int64_t *ldb, double *c__, aocl_int64_t *ldc, double *d__, aocl_int64_t *ldd,
                 double *e, aocl_int64_t *lde, double *f, aocl_int64_t *ldf, double *scale, double *rdsum,
                 double *rdscal, aocl_int64_t *iwork, aocl_int64_t *pq, aocl_int64_t *info);
int dtgsyl_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *b, aocl_int64_t *ldb, double *c__, aocl_int64_t *ldc, double *d__, aocl_int64_t *ldd,
                 double *e, aocl_int64_t *lde, double *f, aocl_int64_t *ldf, double *scale, double *dif,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int dtpcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, double *ap, double *rcond,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dtpmqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                  aocl_int64_t *nb, double *v, aocl_int64_t *ldv, double *t, aocl_int64_t *ldt, double *a,
                  aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *work, aocl_int64_t *info);
int dtpqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, aocl_int64_t *nb, double *a, aocl_int64_t *lda,
                 double *b, aocl_int64_t *ldb, double *t, aocl_int64_t *ldt, double *work, aocl_int64_t *info);
int dtpqrt2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, double *a, aocl_int64_t *lda, double *b,
                  aocl_int64_t *ldb, double *t, aocl_int64_t *ldt, aocl_int64_t *info);
int dtprfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, double *v, aocl_int64_t *ldv, double *t, aocl_int64_t *ldt,
                 double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *work, aocl_int64_t *ldwork);
int dtprfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap,
                 double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dtptri_check(char *uplo, char *diag, aocl_int64_t *n, double *ap, aocl_int64_t *info);
int dtptrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, double *ap,
                 double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dtpttf_check(char *transr, char *uplo, aocl_int64_t *n, double *ap, double *arf, aocl_int64_t *info);
int dtpttr_check(char *uplo, aocl_int64_t *n, double *ap, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dtrcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *rcond, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dtrevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, double *t, aocl_int64_t *ldt,
                 double *vl, aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m,
                 double *work, aocl_int64_t *info);
int dtrexc_check(char *compq, aocl_int64_t *n, double *t, aocl_int64_t *ldt, double *q, aocl_int64_t *ldq,
                 aocl_int64_t *ifst, aocl_int64_t *ilst, double *work, aocl_int64_t *info);
int dtrrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, double *a,
                 aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *ferr,
                 double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dtrsen_check(char *job, char *compq, logical *select, aocl_int64_t *n, double *t, aocl_int64_t *ldt,
                 double *q, aocl_int64_t *ldq, double *wr, double *wi, aocl_int64_t *m, double *s,
                 double *sep, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork,
                 aocl_int64_t *info);
int dtrsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, double *t, aocl_int64_t *ldt,
                 double *vl, aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr, double *s, double *sep,
                 aocl_int64_t *mm, aocl_int64_t *m, double *work, aocl_int64_t *ldwork, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dtrsyl_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *c__, aocl_int64_t *ldc, double *scale,
                 aocl_int64_t *info);
int dtrsyl3_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, doublereal *a,
                  aocl_int64_t *lda, doublereal *b, aocl_int64_t *ldb, doublereal *c__, aocl_int64_t *ldc,
                  doublereal *scale, aocl_int64_t *iwork, aocl_int64_t *liwork, doublereal *swork,
                  aocl_int64_t *ldswork, aocl_int64_t *info);
int dtrti2_check(char *uplo, char *diag, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dtrtri_check(char *uplo, char *diag, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dtrtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, double *a,
                 aocl_int64_t *lda, double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dtrttf_check(char *transr, char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *arf,
                 aocl_int64_t *info);
int dtrttp_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *ap, aocl_int64_t *info);
int dtzrqf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, aocl_int64_t *info);
int dtzrzf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
double dzsum1_check(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx);
int icmax1_check(aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx);
int ilaclc_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda);
int ilaclr_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda);
int iladiag_check(char *diag);
int iladlc_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda);
int iladlr_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda);
int ilaprec_check(char *prec);
int ilaslc_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda);
int ilaslr_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda);
int ilatrans_check(char *trans);
int ilauplo_check(char *uplo);
int ilaver_check(int *vers_major__, int *vers_minor__, int *vers_patch__);
int ilazlc_check(int *m, int *n, dcomplex *a, int *lda);
int ilazlr_check(int *m, int *n, dcomplex *a, int *lda);
int izmax1_check(int *n, dcomplex *cx, int *incx);
int sbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, int *m, int *p, int *q, float *theta, float *phi, float *u1, int *ldu1, float *u2, int *ldu2, float * v1t, int *ldv1t, float *v2t, int *ldv2t, float *b11d, float * b11e, float *b12d, float *b12e, float *b21d, float *b21e, float *b22d, float *b22e, float *work, int *lwork, int *info);
int sbdsdc_check(char *uplo, char *compq, int *n, float *d__, float *e, float *u, int *ldu, float *vt, int *ldvt, float *q, int *iq, float *work, int *iwork, int *info);
int sbdsqr_check(char *uplo, int *n, int *ncvt, int * nru, int *ncc, float *d__, float *e, float *vt, int *ldvt, float * u, int *ldu, float *c__, int *ldc, float *work, int *info);
float scsum1_check(int *n, scomplex *cx, int *incx);
int sdisna_check(char *job, int *m, int *n, float *d__, float *sep, int *info);
int sgbbrd_check(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, float *ab, int *ldab, float *d__, float * e, float *q, int *ldq, float *pt, int *ldpt, float *c__, int *ldc, float *work, int *info);
int sgbcon_check(char *norm, int *n, int *kl, int *ku, float *ab, int *ldab, int *ipiv, float *anorm, float *rcond, float *work, int *iwork, int *info);
int sgbequ_check(int *m, int *n, int *kl, int *ku, float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float * colcnd, float *amax, int *info);
int sgbequb_check(int *m, int *n, int *kl, int * ku, float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);
int sgbrfs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, int *ipiv, float *b, int *ldb, float *x, int *ldx, float * ferr, float *berr, float *work, int *iwork, int *info);
int sgbrfsx_check(char *trans, char *equed, int *n, int * kl, int *ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, int *ipiv, float *r__, float *c__, float *b, int *ldb, float *x, int *ldx, float *rcond, float *berr, int * n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int * nparams, float *params, float *work, int *iwork, int *info);
int sgbsv_check(int *n, int *kl, int *ku, int * nrhs, float *ab, int *ldab, int *ipiv, float *b, int *ldb, int *info);
int sgbsvx_check(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, int *ipiv, char *equed, float *r__, float *c__, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
int sgbsvxx_check(char *fact, char *trans, int *n, int * kl, int *ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, int *ipiv, char *equed, float *r__, float *c__, float *b, int *ldb, float *x, int *ldx, float *rcond, float * rpvgrw, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int *nparams, float *params, float *work, int *iwork, int *info);
int sgbtf2_check(int *m, int *n, int *kl, int *ku, float *ab, int *ldab, int *ipiv, int *info);
int sgbtrf_check(int *m, int *n, int *kl, int *ku, float *ab, int *ldab, int *ipiv, int *info);
int sgbtrs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, float *ab, int *ldab, int *ipiv, float *b, int *ldb, int *info);
int sgebak_check(char *job, char *side, int *n, int *ilo, int *ihi, float *scale, int *m, float *v, int *ldv, int *info);
int sgebal_check(char *job, int *n, float *a, int *lda, int *ilo, int *ihi, float *scale, int *info);
int sgebd2_check(int *m, int *n, float *a, int *lda, float *d__, float *e, float *tauq, float *taup, float *work, int *info);
int sgebrd_check(int *m, int *n, float *a, int *lda, float *d__, float *e, float *tauq, float *taup, float *work, int * lwork, int *info);
int sgecon_check(char *norm, int *n, float *a, int *lda, float *anorm, float *rcond, float *work, int *iwork, int *info);
int sgeequ_check(int *m, int *n, float *a, int *lda, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);
int sgeequb_check(int *m, int *n, float *a, int *lda, float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int *info);
int sgees_check(char *jobvs, char *sort, L_fp select, int *n, float *a, int *lda, int *sdim, float *wr, float *wi, float *vs, int *ldvs, float *work, int *lwork, logical *bwork, int * info);
int sgeesx_check(char *jobvs, char *sort, L_fp select, char * sense, int *n, float *a, int *lda, int *sdim, float *wr, float *wi, float *vs, int *ldvs, float *rconde, float *rcondv, float * work, int *lwork, int *iwork, int *liwork, logical *bwork, int *info);
int sgeev_check(char *jobvl, char *jobvr, int *n, float *a, int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);
int sgeevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, float *a, int *lda, float *wr, float *wi, float * vl, int *ldvl, float *vr, int *ldvr, int *ilo, int * ihi, float *scale, float *abnrm, float *rconde, float *rcondv, float *work, int *lwork, int *iwork, int *info);
int sgegs_check(char *jobvsl, char *jobvsr, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float * work, int *lwork, int *info);
int sgegv_check(char *jobvl, char *jobvr, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);
int sgehd2_check(int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *work, int *info);
int sgehrd_check(int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sgejsv_check(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, int *m, int *n, float *a, int *lda, float *sva, float *u, int *ldu, float *v, int *ldv, float *work, int *lwork, int *iwork, int *info);
int sgelq2_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);
int sgelqf_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sgels_check(char *trans, int *m, int *n, int * nrhs, float *a, int *lda, float *b, int *ldb, float *work, int *lwork, int *info);
int sgelsd_check(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *s, float *rcond, int * rank, float *work, int *lwork, int *iwork, int *info);
int sgelss_check(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *s, float *rcond, int * rank, float *work, int *lwork, int *info);
int sgelsx_check(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *jpvt, float *rcond, int *rank, float *work, int *info);
int sgelsy_check(int *m, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *jpvt, float *rcond, int *rank, float *work, int *lwork, int *info);
int sgemqrt_check(char *side, char *trans, int *m, int *n, int *k, int *nb, float *v, int *ldv, float *t, int * ldt, float *c__, int *ldc, float *work, int *info);
int sgeql2_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);
int sgeqlf_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sgeqp3_check(int *m, int *n, float *a, int *lda, int *jpvt, float *tau, float *work, int *lwork, int *info);
int sgeqpf_check(int *m, int *n, float *a, int *lda, int *jpvt, float *tau, float *work, int *info);
int sgeqr2_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);
int sgeqr2p_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);
int sgeqrf_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sgeqrfp_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sgeqrt_check(int *m, int *n, int *nb, float *a, int *lda, float *t, int *ldt, float *work, int *info);
int sgeqrt2_check(int *m, int *n, float *a, int *lda, float *t, int *ldt, int *info);
int sgeqrt3_check(int *m, int *n, float *a, int *lda, float *t, int *ldt, int *info);
int sgerfs_check(char *trans, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float * work, int *iwork, int *info);
int sgerfsx_check(char *trans, char *equed, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *r__, float *c__, float *b, int *ldb, float *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int *nparams, float *params, float *work, int *iwork, int *info);
int sgerq2_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);
int sgerqf_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sgesc2_check(int *n, float *a, int *lda, float *rhs, int *ipiv, int *jpiv, float *scale);
int sgesdd_check(char *jobz, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *iwork, int *info);
int sgesv_check(int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);
int sgesvd_check(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info);
int sgesvj_check(char *joba, char *jobu, char *jobv, int *m, int *n, float *a, int *lda, float *sva, int *mv, float *v, int *ldv, float *work, int *lwork, int *info);
int sgesvx_check(char *fact, char *trans, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, char *equed, float *r__, float *c__, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
int sgesvxx_check(char *fact, char *trans, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, char *equed, float *r__, float *c__, float *b, int *ldb, float *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int * n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int * nparams, float *params, float *work, int *iwork, int *info);
int sgetc2_check(int *n, float *a, int *lda, int *ipiv, int *jpiv, int *info);
int sgetf2_check(int *m, int *n, float *a, int *lda, int *ipiv, int *info);
int sgetrf_check(int *m, int *n, float *a, int *lda, int *ipiv, int *info);
int sgetrfnp_check(int *m, int *n, float *a, int *lda, int *info);
int sgetri_check(int *n, float *a, int *lda, int *ipiv, float *work, int *lwork, int *info);
int sgetrs_check(char *trans, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);
int sggbak_check(char *job, char *side, int *n, int *ilo, int *ihi, float *lscale, float *rscale, int *m, float *v, int *ldv, int *info);
int sggbal_check(char *job, int *n, float *a, int *lda, float *b, int *ldb, int *ilo, int *ihi, float *lscale, float *rscale, float *work, int *info);
int sgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, int *n, float *a, int *lda, float *b, int *ldb, int *sdim, float *alphar, float *alphai, float *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *work, int *lwork, logical *bwork, int *info);
int sggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, int *n, float *a, int *lda, float *b, int *ldb, int *sdim, float *alphar, float *alphai, float *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *rconde, float *rcondv, float *work, int *lwork, int *iwork, int * liwork, logical *bwork, int *info);
int sggev_check(char *jobvl, char *jobvr, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, int *lwork, int *info);
int sggevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, int *ldvr, int *ilo, int *ihi, float *lscale, float *rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, float *work, int *lwork, int *iwork, logical *bwork, int *info);
int sggglm_check(int *n, int *m, int *p, float *a, int *lda, float *b, int *ldb, float *d__, float *x, float *y, float *work, int *lwork, int *info);
int sgghrd_check(char *compq, char *compz, int *n, int * ilo, int *ihi, float *a, int *lda, float *b, int *ldb, float *q, int *ldq, float *z__, int *ldz, int *info);
int sgglse_check(int *m, int *n, int *p, float *a, int *lda, float *b, int *ldb, float *c__, float *d__, float *x, float *work, int *lwork, int *info);
int sggqrf_check(int *n, int *m, int *p, float *a, int *lda, float *taua, float *b, int *ldb, float *taub, float * work, int *lwork, int *info);
int sggrqf_check(int *m, int *p, int *n, float *a, int *lda, float *taua, float *b, int *ldb, float *taub, float * work, int *lwork, int *info);
int sggsvd_check(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, float *a, int *lda, float *b, int *ldb, float *alpha, float *beta, float *u, int * ldu, float *v, int *ldv, float *q, int *ldq, float *work, int *iwork, int *info);
int sggsvp_check(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, float *a, int *lda, float *b, int *ldb, float *tola, float *tolb, int *k, int *l, float *u, int *ldu, float *v, int *ldv, float *q, int *ldq, int *iwork, float * tau, float *work, int *info);
int sgsvj0_check(char *jobv, int *m, int *n, float *a, int *lda, float *d__, float *sva, int *mv, float *v, int * ldv, float *eps, float *sfmin, float *tol, int *nsweep, float *work, int *lwork, int *info);
int sgsvj1_check(char *jobv, int *m, int *n, int *n1, float *a, int *lda, float *d__, float *sva, int *mv, float *v, int *ldv, float *eps, float *sfmin, float *tol, int *nsweep, float *work, int *lwork, int *info);
int sgtcon_check(char *norm, int *n, float *dl, float *d__, float *du, float *du2, int *ipiv, float *anorm, float *rcond, float * work, int *iwork, int *info);
int sgtrfs_check(char *trans, int *n, int *nrhs, float *dl, float *d__, float *du, float *dlf, float *df, float *duf, float *du2, int *ipiv, float *b, int *ldb, float *x, int *ldx, float * ferr, float *berr, float *work, int *iwork, int *info);
int sgtsv_check(int *n, int *nrhs, float *dl, float *d__, float *du, float *b, int *ldb, int *info);
int sgtsvx_check(char *fact, char *trans, int *n, int * nrhs, float *dl, float *d__, float *du, float *dlf, float *df, float *duf, float *du2, int *ipiv, float *b, int *ldb, float *x, int * ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
int sgttrf_check(int *n, float *dl, float *d__, float *du, float * du2, int *ipiv, int *info);
int sgttrs_check(char *trans, int *n, int *nrhs, float *dl, float *d__, float *du, float *du2, int *ipiv, float *b, int *ldb, int *info);
int sgtts2_check(int *itrans, int *n, int *nrhs, float *dl, float *d__, float *du, float *du2, int *ipiv, float *b, int * ldb);
int shgeqz_check(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, float *h__, int *ldh, float *t, int *ldt, float *alphar, float *alphai, float *beta, float *q, int *ldq, float *z__, int *ldz, float *work, int *lwork, int *info);
int shsein_check(char *side, char *eigsrc, char *initv, logical * select, int *n, float *h__, int *ldh, float *wr, float *wi, float *vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, float *work, int *ifaill, int *ifailr, int *info);
int shseqr_check(char *job, char *compz, int *n, int *ilo, int *ihi, float *h__, int *ldh, float *wr, float *wi, float *z__, int *ldz, float *work, int *lwork, int *info);
logical sisnan_check(float *sin__);
int sla_gbamv_check(int *trans, int *m, int *n, int *kl, int *ku, float *alpha, float *ab, int *ldab, float * x, int *incx, float *beta, float *y, int *incy);
float sla_gbrcond_check(char *trans, int *n, int *kl, int *ku, float * ab, int *ldab, float *afb, int *ldafb, int *ipiv, int * cmode, float *c__, int *info, float *work, int *iwork);
int sla_gbrfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, int *ipiv, logical *colequ, float *c__, float *b, int *ldb, float *y, int * ldy, float *berr_out__, int *n_norms__, float *err_bnds_norm__, float *err_bnds_comp__, float *res, float *ayb, float *dy, float *y_tail__, float *rcond, int *ithresh, float *rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float sla_gbrpvgrw_check(int *n, int *kl, int *ku, int *ncols, float *ab, int *ldab, float *afb, int *ldafb);
int sla_geamv_check(int *trans, int *m, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);
float sla_gercond_check(char *trans, int *n, float *a, int *lda, float *af, int *ldaf, int *ipiv, int *cmode, float *c__, int * info, float *work, int *iwork);
int sla_gerfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *nrhs, float *a, int *lda, float * af, int *ldaf, int *ipiv, logical *colequ, float *c__, float *b, int *ldb, float *y, int *ldy, float *berr_out__, int * n_norms__, float *errs_n__, float *errs_c__, float *res, float *ayb, float *dy, float *y_tail__, float *rcond, int *ithresh, float *rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float sla_gerpvgrw_check(int *n, int *ncols, float *a, int *lda, float * af, int *ldaf);
int sla_lin_berr_check(int *n, int *nz, int *nrhs, float *res, float *ayb, float *berr);
float sla_porcond_check(char *uplo, int *n, float *a, int *lda, float *af, int *ldaf, int *cmode, float *c__, int *info, float *work, int *iwork);
int sla_porfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int * ldaf, logical *colequ, float *c__, float *b, int *ldb, float *y, int *ldy, float *berr_out__, int *n_norms__, float * err_bnds_norm__, float *err_bnds_comp__, float *res, float *ayb, float * dy, float *y_tail__, float *rcond, int *ithresh, float *rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float sla_porpvgrw_check(char *uplo, int *ncols, float *a, int *lda, float * af, int *ldaf, float *work);
int sla_syamv_check(int *uplo, int *n, float *alpha, float *a, int *lda, float *x, int *incx, float *beta, float *y, int *incy);
float sla_syrcond_check(char *uplo, int *n, float *a, int *lda, float *af, int *ldaf, int *ipiv, int *cmode, float *c__, int * info, float *work, int *iwork);
int sla_syrfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int * ldaf, int *ipiv, logical *colequ, float *c__, float *b, int * ldb, float *y, int *ldy, float *berr_out__, int *n_norms__, float *err_bnds_norm__, float *err_bnds_comp__, float *res, float *ayb, float *dy, float *y_tail__, float *rcond, int *ithresh, float * rthresh, float *dz_ub__, logical *ignore_cwise__, int *info);
float sla_syrpvgrw_check(char *uplo, int *n, int *info, float *a, int * lda, float *af, int *ldaf, int *ipiv, float *work);
int sla_wwaddw_check(int *n, float *x, float *y, float *w);
int slabad_check(float *small_, float *large);
int slabrd_check(int *m, int *n, int *nb, float *a, int *lda, float *d__, float *e, float *tauq, float *taup, float *x, int *ldx, float *y, int *ldy);
int slacn2_check(int *n, float *v, float *x, int *isgn, float *est, int *kase, int *isave);
int slacon_check(int *n, float *v, float *x, int *isgn, float *est, int *kase);
int slacpy_check(char *uplo, int *m, int *n, float *a, int *lda, float *b, int *ldb);
int sladiv_check(float *a, float *b, float *c__, float *d__, float *p, float *q);
int slae2_check(float *a, float *b, float *c__, float *rt1, float *rt2);
int slaebz_check(aocl_int64_t *ijob, aocl_int64_t *nitmax, aocl_int64_t *n, aocl_int64_t *mmax, aocl_int64_t *minp,
                 aocl_int64_t *nbmin, float *abstol, float *reltol, float *pivmin, float *d__, float *e,
                 float *e2, aocl_int64_t *nval, float *ab, float *c__, aocl_int64_t *mout, aocl_int64_t *nab,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int slaed0_check(aocl_int64_t *icompq, aocl_int64_t *qsiz, aocl_int64_t *n, float *d__, float *e, float *q,
                 aocl_int64_t *ldq, float *qstore, aocl_int64_t *ldqs, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int slaed1_check(aocl_int64_t *n, float *d__, float *q, aocl_int64_t *ldq, aocl_int64_t *indxq, float *rho,
                 aocl_int64_t *cutpnt, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int slaed2_check(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *n1, float *d__, float *q, aocl_int64_t *ldq,
                 aocl_int64_t *indxq, float *rho, float *z__, float *dlamda, float *w, float *q2,
                 aocl_int64_t *indx, aocl_int64_t *indxc, aocl_int64_t *indxp, aocl_int64_t *coltyp, aocl_int64_t *info);
int slaed3_check(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *n1, float *d__, float *q, aocl_int64_t *ldq,
                 float *rho, float *dlamda, float *q2, aocl_int64_t *indx, aocl_int64_t *ctot, float *w,
                 float *s, aocl_int64_t *info);
int slaed4_check(aocl_int64_t *n, aocl_int64_t *i__, float *d__, float *z__, float *delta, float *rho,
                 float *dlam, aocl_int64_t *info);
int slaed5_check(aocl_int64_t *i__, float *d__, float *z__, float *delta, float *rho, float *dlam);
int slaed6_check(aocl_int64_t *kniter, logical *orgati, float *rho, float *d__, float *z__, float *finit,
                 float *tau, aocl_int64_t *info);
int slaed7_check(aocl_int64_t *icompq, aocl_int64_t *n, aocl_int64_t *qsiz, aocl_int64_t *tlvls, aocl_int64_t *curlvl,
                 aocl_int64_t *curpbm, float *d__, float *q, aocl_int64_t *ldq, aocl_int64_t *indxq, float *rho,
                 aocl_int64_t *cutpnt, float *qstore, aocl_int64_t *qptr, aocl_int64_t *prmptr, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, float *givnum, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int slaed8_check(aocl_int64_t *icompq, aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *qsiz, float *d__, float *q,
                 aocl_int64_t *ldq, aocl_int64_t *indxq, float *rho, aocl_int64_t *cutpnt, float *z__,
                 float *dlamda, float *q2, aocl_int64_t *ldq2, float *w, aocl_int64_t *perm, aocl_int64_t *givptr,
                 aocl_int64_t *givcol, float *givnum, aocl_int64_t *indxp, aocl_int64_t *indx, aocl_int64_t *info);
int slaed9_check(aocl_int64_t *k, aocl_int64_t *kstart, aocl_int64_t *kstop, aocl_int64_t *n, float *d__, float *q,
                 aocl_int64_t *ldq, float *rho, float *dlamda, float *w, float *s, aocl_int64_t *lds,
                 aocl_int64_t *info);
int slaeda_check(aocl_int64_t *n, aocl_int64_t *tlvls, aocl_int64_t *curlvl, aocl_int64_t *curpbm, aocl_int64_t *prmptr,
                 aocl_int64_t *perm, aocl_int64_t *givptr, aocl_int64_t *givcol, float *givnum, float *q,
                 aocl_int64_t *qptr, float *z__, float *ztemp, aocl_int64_t *info);
int slaein_check(logical *rightv, logical *noinit, aocl_int64_t *n, float *h__, aocl_int64_t *ldh, float *wr,
                 float *wi, float *vr, float *vi, float *b, aocl_int64_t *ldb, float *work, float *eps3,
                 float *smlnum, float *bignum, aocl_int64_t *info);
int slaev2_check(float *a, float *b, float *c__, float *rt1, float *rt2, float *cs1, float *sn1);
int slaexc_check(logical *wantq, aocl_int64_t *n, float *t, aocl_int64_t *ldt, float *q, aocl_int64_t *ldq,
                 aocl_int64_t *j1, aocl_int64_t *n1, aocl_int64_t *n2, float *work, aocl_int64_t *info);
int slag2_check(float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *safmin, float *scale1,
                float *scale2, float *wr1, float *wr2, float *wi);
int slag2d_check(aocl_int64_t *m, aocl_int64_t *n, float *sa, aocl_int64_t *ldsa, double *a, aocl_int64_t *lda,
                 aocl_int64_t *info);
int slags2_check(logical *upper, float *a1, float *a2, float *a3, float *b1, float *b2, float *b3,
                 float *csu, float *snu, float *csv, float *snv, float *csq, float *snq);
int slagtf_check(aocl_int64_t *n, float *a, float *lambda, float *b, float *c__, float *tol, float *d__,
                 aocl_int64_t *in, aocl_int64_t *info);
int slagtm_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *alpha, float *dl, float *d__,
                 float *du, float *x, aocl_int64_t *ldx, float *beta, float *b, aocl_int64_t *ldb);
int slagts_check(aocl_int64_t *job, aocl_int64_t *n, float *a, float *b, float *c__, float *d__, aocl_int64_t *in,
                 float *y, float *tol, aocl_int64_t *info);
int slagv2_check(float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *alphar, float *alphai,
                 float *beta, float *csl, float *snl, float *csr, float *snr);
int slahqr_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *h__,
                 aocl_int64_t *ldh, float *wr, float *wi, aocl_int64_t *iloz, aocl_int64_t *ihiz, float *z__,
                 aocl_int64_t *ldz, aocl_int64_t *info);
int slahr2_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, float *a, aocl_int64_t *lda, float *tau, float *t,
                 aocl_int64_t *ldt, float *y, aocl_int64_t *ldy);
int slahrd_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, float *a, aocl_int64_t *lda, float *tau, float *t,
                 aocl_int64_t *ldt, float *y, aocl_int64_t *ldy);
int slaic1_check(aocl_int64_t *job, aocl_int64_t *j, float *x, float *sest, float *w, float *gamma,
                 float *sestpr, float *s, float *c__);
logical slaisnan_check(float *sin1, float *sin2);
int slaln2_check(logical *ltrans, aocl_int64_t *na, aocl_int64_t *nw, float *smin, float *ca, float *a,
                 aocl_int64_t *lda, float *d1, float *d2, float *b, aocl_int64_t *ldb, float *wr, float *wi,
                 float *x, aocl_int64_t *ldx, float *scale, float *xnorm, aocl_int64_t *info);
int slals0_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *nrhs, float *b,
                 aocl_int64_t *ldb, float *bx, aocl_int64_t *ldbx, aocl_int64_t *perm, aocl_int64_t *givptr,
                 aocl_int64_t *givcol, aocl_int64_t *ldgcol, float *givnum, aocl_int64_t *ldgnum, float *poles,
                 float *difl, float *difr, float *z__, aocl_int64_t *k, float *c__, float *s,
                 float *work, aocl_int64_t *info);
int slalsa_check(aocl_int64_t *icompq, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, float *b,
                 aocl_int64_t *ldb, float *bx, aocl_int64_t *ldbx, float *u, aocl_int64_t *ldu, float *vt,
                 aocl_int64_t *k, float *difl, float *difr, float *z__, float *poles, aocl_int64_t *givptr,
                 aocl_int64_t *givcol, aocl_int64_t *ldgcol, aocl_int64_t *perm, float *givnum, float *c__,
                 float *s, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int slalsd_check(char *uplo, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, float *e,
                 float *b, aocl_int64_t *ldb, float *rcond, aocl_int64_t *rank, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int slamrg_check(aocl_int64_t *n1, aocl_int64_t *n2, float *a, aocl_int64_t *strd1, aocl_int64_t *strd2,
                 aocl_int64_t *index);
int slaneg_check(aocl_int64_t *n, float *d__, float *lld, float *sigma, float *pivmin, aocl_int64_t *r__);
float slangb_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                   float *work);
float slange_check(char *norm, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *work);
float slangt_check(char *norm, aocl_int64_t *n, float *dl, float *d__, float *du);
float slanhs_check(char *norm, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *work);
float slansb_check(char *norm, char *uplo, aocl_int64_t *n, aocl_int64_t *k, float *ab, aocl_int64_t *ldab,
                   float *work);
float slansf_check(char *norm, char *transr, char *uplo, aocl_int64_t *n, float *a, float *work);
float slansp_check(char *norm, char *uplo, aocl_int64_t *n, float *ap, float *work);
float slanst_check(char *norm, aocl_int64_t *n, float *d__, float *e);
float slansy_check(char *norm, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *work);
float slantb_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *k, float *ab,
                   aocl_int64_t *ldab, float *work);
float slantp_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, float *ap, float *work);
float slantr_check(char *norm, char *uplo, char *diag, aocl_int64_t *m, aocl_int64_t *n, float *a,
                   aocl_int64_t *lda, float *work);
int slanv2_check(float *a, float *b, float *c__, float *d__, float *rt1r, float *rt1i, float *rt2r,
                 float *rt2i, float *cs, float *sn);
int slapll_check(aocl_int64_t *n, float *x, aocl_int64_t *incx, float *y, aocl_int64_t *incy, float *ssmin);
int slapmr_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, float *x, aocl_int64_t *ldx, aocl_int64_t *k);
int slapmt_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, float *x, aocl_int64_t *ldx, aocl_int64_t *k);
float slapy2_check(float *x, float *y);
float slapy3_check(float *x, float *y, float *z__);
int slaqgb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *equed);
int slaqge_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, char *equed);
int slaqp2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, float *a, aocl_int64_t *lda, aocl_int64_t *jpvt,
                 float *tau, float *vn1, float *vn2, float *work);
int slaqps_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, aocl_int64_t *nb, aocl_int64_t *kb, float *a,
                 aocl_int64_t *lda, aocl_int64_t *jpvt, float *tau, float *vn1, float *vn2, float *auxv,
                 float *f, aocl_int64_t *ldf);
int slaqr0_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *h__,
                 aocl_int64_t *ldh, float *wr, float *wi, aocl_int64_t *iloz, aocl_int64_t *ihiz, float *z__,
                 aocl_int64_t *ldz, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int slaqr1_check(aocl_int64_t *n, float *h__, aocl_int64_t *ldh, float *sr1, float *si1, float *sr2,
                 float *si2, float *v);
int slaqr2_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, float *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz, float *z__,
                 aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, float *sr, float *si, float *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, float *t, aocl_int64_t *ldt, aocl_int64_t *nv, float *wv,
                 aocl_int64_t *ldwv, float *work, aocl_int64_t *lwork);
int slaqr3_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, float *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz, float *z__,
                 aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, float *sr, float *si, float *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, float *t, aocl_int64_t *ldt, aocl_int64_t *nv, float *wv,
                 aocl_int64_t *ldwv, float *work, aocl_int64_t *lwork);
int slaqr4_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *h__,
                 aocl_int64_t *ldh, float *wr, float *wi, aocl_int64_t *iloz, aocl_int64_t *ihiz, float *z__,
                 aocl_int64_t *ldz, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int slaqr5_check(logical *wantt, logical *wantz, aocl_int64_t *kacc22, aocl_int64_t *n, aocl_int64_t *ktop,
                 aocl_int64_t *kbot, aocl_int64_t *nshfts, float *sr, float *si, float *h__, aocl_int64_t *ldh,
                 aocl_int64_t *iloz, aocl_int64_t *ihiz, float *z__, aocl_int64_t *ldz, float *v, aocl_int64_t *ldv,
                 float *u, aocl_int64_t *ldu, aocl_int64_t *nv, float *wv, aocl_int64_t *ldwv, aocl_int64_t *nh,
                 float *wh, aocl_int64_t *ldwh);
int slaqsb_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab, float *s,
                 float *scond, float *amax, char *equed);
int slaqsp_check(char *uplo, aocl_int64_t *n, float *ap, float *s, float *scond, float *amax,
                 char *equed);
int slaqsy_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *s, float *scond,
                 float *amax, char *equed);
int slaqtr_check(logical *ltran, logical *lfloat, aocl_int64_t *n, float *t, aocl_int64_t *ldt, float *b,
                 float *w, float *scale, float *x, float *work, aocl_int64_t *info);
int slar1v_check(aocl_int64_t *n, aocl_int64_t *b1, aocl_int64_t *bn, float *lambda, float *d__, float *l,
                 float *ld, float *lld, float *pivmin, float *gaptol, float *z__, logical *wantnc,
                 aocl_int64_t *negcnt, float *ztz, float *mingma, aocl_int64_t *r__, aocl_int64_t *isuppz,
                 float *nrminv, float *resid, float *rqcorr, float *work);
int slar2v_check(aocl_int64_t *n, float *x, float *y, float *z__, aocl_int64_t *incx, float *c__, float *s,
                 aocl_int64_t *incc);
int slarf_check(char *side, aocl_int64_t *m, aocl_int64_t *n, float *v, aocl_int64_t *incv, float *tau, float *c__,
                aocl_int64_t *ldc, float *work);
int slarfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, float *v, aocl_int64_t *ldv, float *t, aocl_int64_t *ldt, float *c__,
                 aocl_int64_t *ldc, float *work, aocl_int64_t *ldwork);
int slarfg_check(aocl_int64_t *n, float *alpha, float *x, aocl_int64_t *incx, float *tau);
int slarfgp_check(aocl_int64_t *n, float *alpha, float *x, aocl_int64_t *incx, float *tau);
int slarft_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, float *v, aocl_int64_t *ldv,
                 float *tau, float *t, aocl_int64_t *ldt);
int slarfx_check(char *side, aocl_int64_t *m, aocl_int64_t *n, float *v, float *tau, float *c__, aocl_int64_t *ldc,
                 float *work);
int slargv_check(aocl_int64_t *n, float *x, aocl_int64_t *incx, float *y, aocl_int64_t *incy, float *c__,
                 aocl_int64_t *incc);
int slarnv_check(aocl_int64_t *idist, aocl_int64_t *iseed, aocl_int64_t *n, float *x);
int slarra_check(aocl_int64_t *n, float *d__, float *e, float *e2, float *spltol, float *tnrm,
                 aocl_int64_t *nsplit, aocl_int64_t *isplit, aocl_int64_t *info);
int slarrb_check(aocl_int64_t *n, float *d__, float *lld, aocl_int64_t *ifirst, aocl_int64_t *ilast, float *rtol1,
                 float *rtol2, aocl_int64_t *offset, float *w, float *wgap, float *werr, float *work,
                 aocl_int64_t *iwork, float *pivmin, float *spdiam, aocl_int64_t *twist, aocl_int64_t *info);
int slarrc_check(char *jobt, aocl_int64_t *n, float *vl, float *vu, float *d__, float *e, float *pivmin,
                 aocl_int64_t *eigcnt, aocl_int64_t *lcnt, aocl_int64_t *rcnt, aocl_int64_t *info);
int slarrd_check(char *range, char *order, aocl_int64_t *n, float *vl, float *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, float *gers, float *reltol, float *d__, float *e, float *e2,
                 float *pivmin, aocl_int64_t *nsplit, aocl_int64_t *isplit, aocl_int64_t *m, float *w, float *werr,
                 float *wl, float *wu, aocl_int64_t *iblock, aocl_int64_t *indexw, float *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int slarre_check(char *range, aocl_int64_t *n, float *vl, float *vu, aocl_int64_t *il, aocl_int64_t *iu,
                 float *d__, float *e, float *e2, float *rtol1, float *rtol2, float *spltol,
                 aocl_int64_t *nsplit, aocl_int64_t *isplit, aocl_int64_t *m, float *w, float *werr, float *wgap,
                 aocl_int64_t *iblock, aocl_int64_t *indexw, float *gers, float *pivmin, float *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int slarrf_check(aocl_int64_t *n, float *d__, float *l, float *ld, aocl_int64_t *clstrt, aocl_int64_t *clend,
                 float *w, float *wgap, float *werr, float *spdiam, float *clgapl, float *clgapr,
                 float *pivmin, float *sigma, float *dplus, float *lplus, float *work,
                 aocl_int64_t *info);
int slarrj_check(aocl_int64_t *n, float *d__, float *e2, aocl_int64_t *ifirst, aocl_int64_t *ilast, float *rtol,
                 aocl_int64_t *offset, float *w, float *werr, float *work, aocl_int64_t *iwork, float *pivmin,
                 float *spdiam, aocl_int64_t *info);
int slarrk_check(aocl_int64_t *n, aocl_int64_t *iw, float *gl, float *gu, float *d__, float *e2,
                 float *pivmin, float *reltol, float *w, float *werr, aocl_int64_t *info);
int slarrr_check(aocl_int64_t *n, float *d__, float *e, aocl_int64_t *info);
int slarrv_check(aocl_int64_t *n, float *vl, float *vu, float *d__, float *l, float *pivmin,
                 aocl_int64_t *isplit, aocl_int64_t *m, aocl_int64_t *dol, aocl_int64_t *dou, float *minrgp,
                 float *rtol1, float *rtol2, float *w, float *werr, float *wgap, aocl_int64_t *iblock,
                 aocl_int64_t *indexw, float *gers, float *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int slarscl2_check(aocl_int64_t *m, aocl_int64_t *n, float *d__, float *x, aocl_int64_t *ldx);
int slartg_check(float *f, float *g, float *cs, float *sn, float *r__);
int slartgp_check(float *f, float *g, float *cs, float *sn, float *r__);
int slartgs_check(float *x, float *y, float *sigma, float *cs, float * sn);
int slartv_check(int *n, float *x, int *incx, float *y, int *incy, float *c__, float *s, int *incc);
int slaruv_check(int *iseed, int *n, float *x);
int slarz_check(char *side, int *m, int *n, int *l, float *v, int *incv, float *tau, float *c__, int *ldc, float * work);
int slarzb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, int *l, float *v, int *ldv, float *t, int *ldt, float *c__, int *ldc, float * work, int *ldwork);
int slarzt_check(char *direct, char *storev, int *n, int * k, float *v, int *ldv, float *tau, float *t, int *ldt);
int slas2_check(float *f, float *g, float *h__, float *ssmin, float * ssmax);
int slascl_check(char *type__, int *kl, int *ku, float * cfrom, float *cto, int *m, int *n, float *a, int *lda, int *info);
int slascl2_check(int *m, int *n, float *d__, float *x, int *ldx);
int slasd0_check(int *n, int *sqre, float *d__, float *e, float *u, int *ldu, float *vt, int *ldvt, int *smlsiz, int *iwork, float *work, int *info);
int slasd1_check(int *nl, int *nr, int *sqre, float * d__, float *alpha, float *beta, float *u, int *ldu, float *vt, int *ldvt, int *idxq, int *iwork, float *work, int * info);
int slasd2_check(int *nl, int *nr, int *sqre, int *k, float *d__, float *z__, float *alpha, float *beta, float *u, int * ldu, float *vt, int *ldvt, float *dsigma, float *u2, int *ldu2, float *vt2, int *ldvt2, int *idxp, int *idx, int *idxc, int *idxq, int *coltyp, int *info);
int slasd3_check(int *nl, int *nr, int *sqre, int *k, float *d__, float *q, int *ldq, float *dsigma, float *u, int * ldu, float *u2, int *ldu2, float *vt, int *ldvt, float *vt2, int *ldvt2, int *idxc, int *ctot, float *z__, int * info);
int slasd4_check(int *n, int *i__, float *d__, float *z__, float *delta, float *rho, float *sigma, float *work, int *info);
int slasd5_check(int *i__, float *d__, float *z__, float *delta, float *rho, float *dsigma, float *work);
int slasd6_check(int *icompq, int *nl, int *nr, int *sqre, float *d__, float *vf, float *vl, float *alpha, float *beta, int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum, float *poles, float * difl, float *difr, float *z__, int *k, float *c__, float *s, float * work, int *iwork, int *info);
int slasd7_check(int *icompq, int *nl, int *nr, int *sqre, int *k, float *d__, float *z__, float *zw, float *vf, float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma, int *idx, int *idxp, int *idxq, int *perm, int * givptr, int *givcol, int *ldgcol, float *givnum, int * ldgnum, float *c__, float *s, int *info);
int slasd8_check(int *icompq, int *k, float *d__, float * z__, float *vf, float *vl, float *difl, float *difr, int *lddifr, float *dsigma, float *work, int *info);
int slasda_check(int *icompq, int *smlsiz, int *n, int *sqre, float *d__, float *e, float *u, int *ldu, float *vt, int *k, float *difl, float *difr, float *z__, float *poles, int * givptr, int *givcol, int *ldgcol, int *perm, float *givnum, float *c__, float *s, float *work, int *iwork, int *info);
int slasdq_check(char *uplo, int *sqre, int *n, int * ncvt, int *nru, int *ncc, float *d__, float *e, float *vt, int *ldvt, float *u, int *ldu, float *c__, int *ldc, float * work, int *info);
int slasdt_check(int *n, int *lvl, int *nd, int * inode, int *ndiml, int *ndimr, int *msub);
int slaset_check(char *uplo, int *m, int *n, float *alpha, float *beta, float *a, int *lda);
int slasq1_check(int *n, float *d__, float *e, float *work, int *info);
int slasq2_check(int *n, float *z__, int *info);
int slasq3_check(int *i0, int *n0, float *z__, int *pp, float *dmin__, float *sigma, float *desig, float *qmax, int *nfail, int *iter, int *ndiv, logical *ieee, int *ttype, float * dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *g, float * tau);
int slasq4_check(int *i0, int *n0, float *z__, int *pp, int *n0in, float *dmin__, float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau, int *ttype, float *g);
int slasq5_check(int *i0, int *n0, float *z__, int *pp, float *tau, float *sigma, float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float *dnm2, logical *ieee, float *eps);
int slasq6_check(int *i0, int *n0, float *z__, int *pp, float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float * dnm2);
int slasr_check(char *side, char *pivot, char *direct, int *m, int *n, float *c__, float *s, float *a, int *lda);
int slasrt_check(char *id, int *n, float *d__, int *info);
int slassq_check(int *n, float *x, int *incx, float *scale, float *sumsq);
int slasv2_check(float *f, float *g, float *h__, float *ssmin, float * ssmax, float *snr, float *csr, float *snl, float *csl);
int slaswp_check(int *n, float *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);
int slasy2_check(logical *ltranl, logical *ltranr, int *isgn, int *n1, int *n2, float *tl, int *ldtl, float *tr, int * ldtr, float *b, int *ldb, float *scale, float *x, int *ldx, float *xnorm, int *info);
int slasyf_check(char *uplo, int *n, int *nb, int *kb, float *a, int *lda, int *ipiv, float *w, int *ldw, int *info);
int slasyf_rook_check(char *uplo, int *n, int *nb, int *kb, float *a, int *lda, int *ipiv, float *w, int * ldw, int *info);
int slatbs_check(char *uplo, char *trans, char *diag, char * normin, int *n, int *kd, float *ab, int *ldab, float *x, float *scale, float *cnorm, int *info);
int slatdf_check(int *ijob, int *n, float *z__, int * ldz, float *rhs, float *rdsum, float *rdscal, int *ipiv, int * jpiv);
int slatps_check(char *uplo, char *trans, char *diag, char * normin, int *n, float *ap, float *x, float *scale, float *cnorm, int *info);
int slatrd_check(char *uplo, int *n, int *nb, float *a, int *lda, float *e, float *tau, float *w, int *ldw);
int slatrs_check(char *uplo, char *trans, char *diag, char * normin, int *n, float *a, int *lda, float *x, float *scale, float *cnorm, int *info);
int slatrz_check(int *m, int *n, int *l, float *a, int *lda, float *tau, float *work);
int slatzm_check(char *side, int *m, int *n, float *v, int *incv, float *tau, float *c1, float *c2, int *ldc, float * work);
int slauu2_check(char *uplo, int *n, float *a, int *lda, int *info);
int slauum_check(char *uplo, int *n, float *a, int *lda, int *info);
int sopgtr_check(char *uplo, int *n, float *ap, float *tau, float *q, int *ldq, float *work, int *info);
int sopmtr_check(char *side, char *uplo, char *trans, int *m, int *n, float *ap, float *tau, float *c__, int *ldc, float *work, int *info);
int sorbdb_check(char *trans, char *signs, int *m, int *p, int *q, float *x11, int *ldx11, float *x12, int *ldx12, float *x21, int *ldx21, float *x22, int *ldx22, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *tauq2, float * work, int *lwork, int *info);
int sorbdb1_check(int *m, int *p, int *q, float *x11, int *ldx11, float *x21, int *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, int *lwork, int *info);
int sorbdb2_check(int *m, int *p, int *q, float *x11, int *ldx11, float *x21, int *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, int *lwork, int *info);
int sorbdb3_check(int *m, int *p, int *q, float *x11, int *ldx11, float *x21, int *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *work, int *lwork, int *info);
int sorbdb4_check(int *m, int *p, int *q, float *x11, int *ldx11, float *x21, int *ldx21, float *theta, float *phi, float *taup1, float *taup2, float *tauq1, float *phantom, float *work, int *lwork, int *info);
int sorbdb5_check(int *m1, int *m2, int *n, float *x1, int *incx1, float *x2, int *incx2, float *q1, int *ldq1, float *q2, int *ldq2, float *work, int *lwork, int *info);
int sorbdb6_check(int *m1, int *m2, int *n, float *x1, int *incx1, float *x2, int *incx2, float *q1, int *ldq1, float *q2, int *ldq2, float *work, int *lwork, int *info);
int sorcsd_check(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, char *signs, int *m, int *p, int *q, float *x11, int *ldx11, float *x12, int *ldx12, float *x21, int *ldx21, float *x22, int *ldx22, float *theta, float *u1, int *ldu1, float *u2, int *ldu2, float *v1t, int *ldv1t, float *v2t, int *ldv2t, float *work, int *lwork, int *iwork, int *info);
int sorcsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, int *m, int *p, int *q, float *x11, int *ldx11, float * x21, int *ldx21, float *theta, float *u1, int *ldu1, float *u2, int *ldu2, float *v1t, int *ldv1t, float *work, int *lwork, int *iwork, int *info);
int sorg2l_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);
int sorg2r_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);
int sorgbr_check(char *vect, int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sorghr_check(int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sorgl2_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);
int sorglq_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sorgql_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sorgqr_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sorgr2_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *info);
int sorgrq_check(int *m, int *n, int *k, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sorgtr_check(char *uplo, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int sorm2l_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *info);
int sorm2r_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *info);
int sormbr_check(char *vect, char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *lwork, int *info);
int sormhr_check(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, float *a, int *lda, float *tau, float * c__, int *ldc, float *work, int *lwork, int *info);
int sorml2_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *info);
int sormlq_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *lwork, int *info);
int sormql_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *lwork, int *info);
int sormqr_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *lwork, int *info);
int sormr2_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *info);
int sormr3_check(char *side, char *trans, int *m, int *n, int *k, int *l, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *info);
int sormrq_check(char *side, char *trans, int *m, int *n, int *k, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *lwork, int *info);
int sormrz_check(char *side, char *trans, int *m, int *n, int *k, int *l, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *lwork, int *info);
int sormtr_check(char *side, char *uplo, char *trans, int *m, int *n, float *a, int *lda, float *tau, float *c__, int *ldc, float *work, int *lwork, int *info);
int spbcon_check(char *uplo, int *n, int *kd, float *ab, int *ldab, float *anorm, float *rcond, float *work, int *iwork, int *info);
int spbequ_check(char *uplo, int *n, int *kd, float *ab, int *ldab, float *s, float *scond, float *amax, int *info);
int spbrfs_check(char *uplo, int *n, int *kd, int * nrhs, float *ab, int *ldab, float *afb, int *ldafb, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float * work, int *iwork, int *info);
int spbstf_check(char *uplo, int *n, int *kd, float *ab, int *ldab, int *info);
int spbsv_check(char *uplo, int *n, int *kd, int * nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);
int spbsvx_check(char *fact, char *uplo, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, char *equed, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
int spbtf2_check(char *uplo, int *n, int *kd, float *ab, int *ldab, int *info);
int spbtrf_check(char *uplo, int *n, int *kd, float *ab, int *ldab, int *info);
int spbtrs_check(char *uplo, int *n, int *kd, int * nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);
int spftrf_check(char *transr, char *uplo, int *n, float *a, int *info);
int spftri_check(char *transr, char *uplo, int *n, float *a, int *info);
int spftrs_check(char *transr, char *uplo, int *n, int * nrhs, float *a, float *b, int *ldb, int *info);
int spocon_check(char *uplo, int *n, float *a, int *lda, float *anorm, float *rcond, float *work, int *iwork, int *info);
int spoequ_check(int *n, float *a, int *lda, float *s, float *scond, float *amax, int *info);
int spoequb_check(int *n, float *a, int *lda, float *s, float *scond, float *amax, int *info);
int sporfs_check(char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);
int sporfsx_check(char *uplo, char *equed, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, float *s, float * b, int *ldb, float *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__, int *nparams, float *params, float *work, int *iwork, int * info);
int sposv_check(char *uplo, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *info);
int sposvx_check(char *fact, char *uplo, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, char *equed, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
int sposvxx_check(char *fact, char *uplo, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, char *equed, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int *n_err_bnds__, float * err_bnds_norm__, float *err_bnds_comp__, int *nparams, float * params, float *work, int *iwork, int *info);
int spotf2_check(char *uplo, int *n, float *a, int *lda, int *info);
int spotrf_check(char *uplo, int *n, float *a, int *lda, int *info);
int spotri_check(char *uplo, int *n, float *a, int *lda, int *info);
int spotrs_check(char *uplo, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int *info);
int sppcon_check(char *uplo, int *n, float *ap, float *anorm, float *rcond, float *work, int *iwork, int *info);
int sppequ_check(char *uplo, int *n, float *ap, float *s, float * scond, float *amax, int *info);
int spprfs_check(char *uplo, int *n, int *nrhs, float *ap, float *afp, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);
int sppsv_check(char *uplo, int *n, int *nrhs, float *ap, float *b, int *ldb, int *info);
int sppsvx_check(char *fact, char *uplo, int *n, int * nrhs, float *ap, float *afp, char *equed, float *s, float *b, int * ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
int spptrf_check(char *uplo, int *n, float *ap, int *info);
int spptri_check(char *uplo, int *n, float *ap, int *info);
int spptrs_check(char *uplo, int *n, int *nrhs, float *ap, float *b, int *ldb, int *info);
int spstf2_check(char *uplo, int *n, float *a, int *lda, int *piv, int *rank, float *tol, float *work, int *info);
int spstrf_check(char *uplo, int *n, float *a, int *lda, int *piv, int *rank, float *tol, float *work, int *info);
int sptcon_check(int *n, float *d__, float *e, float *anorm, float *rcond, float *work, int *info);
int spteqr_check(char *compz, int *n, float *d__, float *e, float *z__, int *ldz, float *work, int *info);
int sptrfs_check(int *n, int *nrhs, float *d__, float *e, float *df, float *ef, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *info);
int sptsv_check(int *n, int *nrhs, float *d__, float *e, float *b, int *ldb, int *info);
int sptsvx_check(char *fact, int *n, int *nrhs, float *d__, float *e, float *df, float *ef, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *info);
int spttrf_check(int *n, float *d__, float *e, int *info);
int spttrs_check(int *n, int *nrhs, float *d__, float *e, float *b, int *ldb, int *info);
int sptts2_check(int *n, int *nrhs, float *d__, float *e, float *b, int *ldb);
int srscl_check(int *n, float *sa, float *sx, int *incx);
int ssbev_check(char *jobz, char *uplo, int *n, int *kd, float *ab, int *ldab, float *w, float *z__, int *ldz, float *work, int *info);
int ssbevd_check(char *jobz, char *uplo, int *n, int *kd, float *ab, int *ldab, float *w, float *z__, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);
int ssbevx_check(char *jobz, char *range, char *uplo, int *n, int *kd, float *ab, int *ldab, float *q, int *ldq, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float * w, float *z__, int *ldz, float *work, int *iwork, int * ifail, int *info);
int ssbgst_check(char *vect, char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int *ldbb, float * x, int *ldx, float *work, int *info);
int ssbgv_check(char *jobz, char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int *ldbb, float * w, float *z__, int *ldz, float *work, int *info);
int ssbgvd_check(char *jobz, char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int *ldbb, float * w, float *z__, int *ldz, float *work, int *lwork, int * iwork, int *liwork, int *info);
int ssbgvx_check(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, float *ab, int *ldab, float *bb, int * ldbb, float *q, int *ldq, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, float *work, int *iwork, int *ifail, int *info);
int ssbtrd_check(char *vect, char *uplo, int *n, int *kd, float *ab, int *ldab, float *d__, float *e, float *q, int *ldq, float *work, int *info);
int ssfrk_check(char *transr, char *uplo, char *trans, int *n, int *k, float *alpha, float *a, int *lda, float *beta, float * c__);
int sspcon_check(char *uplo, int *n, float *ap, int *ipiv, float *anorm, float *rcond, float *work, int *iwork, int *info);
int sspev_check(char *jobz, char *uplo, int *n, float *ap, float *w, float *z__, int *ldz, float *work, int *info);
int sspevd_check(char *jobz, char *uplo, int *n, float *ap, float *w, float *z__, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);
int sspevx_check(char *jobz, char *range, char *uplo, int *n, float *ap, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, float *work, int * iwork, int *ifail, int *info);
int sspgst_check(int *itype, char *uplo, int *n, float *ap, float *bp, int *info);
int sspgv_check(int *itype, char *jobz, char *uplo, int * n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, int *info);
int sspgvd_check(int *itype, char *jobz, char *uplo, int * n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);
int sspgvx_check(int *itype, char *jobz, char *range, char * uplo, int *n, float *ap, float *bp, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int * ldz, float *work, int *iwork, int *ifail, int *info);
int ssprfs_check(char *uplo, int *n, int *nrhs, float *ap, float *afp, int *ipiv, float *b, int *ldb, float *x, int * ldx, float *ferr, float *berr, float *work, int *iwork, int * info);
int sspsv_check(char *uplo, int *n, int *nrhs, float *ap, int *ipiv, float *b, int *ldb, int *info);
int sspsvx_check(char *fact, char *uplo, int *n, int * nrhs, float *ap, float *afp, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, int *info);
int ssptrd_check(char *uplo, int *n, float *ap, float *d__, float *e, float *tau, int *info);
int ssptrf_check(char *uplo, int *n, float *ap, int *ipiv, int *info);
int ssptri_check(char *uplo, int *n, float *ap, int *ipiv, float *work, int *info);
int ssptrs_check(char *uplo, int *n, int *nrhs, float *ap, int *ipiv, float *b, int *ldb, int *info);
int sstebz_check(char *range, char *order, int *n, float *vl, float *vu, int *il, int *iu, float *abstol, float *d__, float *e, int *m, int *nsplit, float *w, int *iblock, int * isplit, float *work, int *iwork, int *info);
int sstedc_check(char *compz, int *n, float *d__, float *e, float *z__, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);
int sstegr_check(char *jobz, char *range, int *n, float *d__, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, int *isuppz, float * work, int *lwork, int *iwork, int *liwork, int *info);
int sstein_check(int *n, float *d__, float *e, int *m, float *w, int *iblock, int *isplit, float *z__, int *ldz, float * work, int *iwork, int *ifail, int *info);
int sstemr_check(char *jobz, char *range, int *n, float *d__, float *e, float *vl, float *vu, int *il, int *iu, int *m, float *w, float *z__, int *ldz, int *nzc, int *isuppz, logical *tryrac, float *work, int *lwork, int *iwork, int * liwork, int *info);
int ssteqr_check(char *compz, int *n, float *d__, float *e, float *z__, int *ldz, float *work, int *info);
int ssterf_check(int *n, float *d__, float *e, int *info);
int sstev_check(char *jobz, int *n, float *d__, float *e, float * z__, int *ldz, float *work, int *info);
int sstevd_check(char *jobz, int *n, float *d__, float *e, float *z__, int *ldz, float *work, int *lwork, int *iwork, int *liwork, int *info);
int sstevr_check(char *jobz, char *range, int *n, float *d__, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, int *isuppz, float * work, int *lwork, int *iwork, int *liwork, int *info);
int sstevx_check(char *jobz, char *range, int *n, float *d__, float *e, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, float *work, int * iwork, int *ifail, int *info);
int ssycon_check(char *uplo, int *n, float *a, int *lda, int *ipiv, float *anorm, float *rcond, float *work, int *iwork, int *info);
int ssycon_rook_check(char *uplo, int *n, float *a, int * lda, int *ipiv, float *anorm, float *rcond, float *work, int * iwork, int *info);
int ssyconv_check(char *uplo, char *way, int *n, float *a, int *lda, int *ipiv, float *work, int *info);
int ssyequb_check(char *uplo, int *n, float *a, int *lda, float *s, float *scond, float *amax, float *work, int *info);
int ssyev_check(char *jobz, char *uplo, int *n, float *a, int *lda, float *w, float *work, int *lwork, int *info);
int ssyevd_check(char *jobz, char *uplo, int *n, float *a, int *lda, float *w, float *work, int *lwork, int *iwork, int *liwork, int *info);
int ssyevr_check(char *jobz, char *range, char *uplo, int *n, float *a, int *lda, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, int * isuppz, float *work, int *lwork, int *iwork, int *liwork, int *info);
int ssyevx_check(char *jobz, char *range, char *uplo, int *n, float *a, int *lda, float *vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, float * work, int *lwork, int *iwork, int *ifail, int *info);
int ssygs2_check(int *itype, char *uplo, int *n, float *a, int *lda, float *b, int *ldb, int *info);
int ssygst_check(int *itype, char *uplo, int *n, float *a, int *lda, float *b, int *ldb, int *info);
int ssygv_check(int *itype, char *jobz, char *uplo, int * n, float *a, int *lda, float *b, int *ldb, float *w, float *work, int *lwork, int *info);
int ssygvd_check(int *itype, char *jobz, char *uplo, int * n, float *a, int *lda, float *b, int *ldb, float *w, float *work, int *lwork, int *iwork, int *liwork, int *info);
int ssygvx_check(int *itype, char *jobz, char *range, char * uplo, int *n, float *a, int *lda, float *b, int *ldb, float * vl, float *vu, int *il, int *iu, float *abstol, int *m, float *w, float *z__, int *ldz, float *work, int *lwork, int *iwork, int *ifail, int *info);
int ssyrfs_check(char *uplo, int *n, int *nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float * work, int *iwork, int *info);
int ssyrfsx_check(char *uplo, char *equed, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *berr, int *n_err_bnds__, float *err_bnds_norm__, float * err_bnds_comp__, int *nparams, float *params, float *work, int * iwork, int *info);
int ssysv_check(char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, float *work, int *lwork, int *info);
int ssysv_rook_check(char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, float *work, int *lwork, int *info);
int ssysvx_check(char *fact, char *uplo, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, int *lwork, int *iwork, int * info);
int ssysvxx_check(char *fact, char *uplo, int *n, int * nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, char *equed, float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, float *rpvgrw, float *berr, int *n_err_bnds__, float * err_bnds_norm__, float *err_bnds_comp__, int *nparams, float * params, float *work, int *iwork, int *info);
int ssyswapr_check(char *uplo, int *n, float *a, int *lda, int *i1, int *i2);
int ssytd2_check(char *uplo, int *n, float *a, int *lda, float *d__, float *e, float *tau, int *info);
int ssytf2_check(char *uplo, int *n, float *a, int *lda, int *ipiv, int *info);
int ssytf2_rook_check(char *uplo, int *n, float *a, int * lda, int *ipiv, int *info);
int ssytrd_check(char *uplo, int *n, float *a, int *lda, float *d__, float *e, float *tau, float *work, int *lwork, int * info);
int ssytrf_check(char *uplo, int *n, float *a, int *lda, int *ipiv, float *work, int *lwork, int *info);
int ssytrf_rook_check(char *uplo, int *n, float *a, int * lda, int *ipiv, float *work, int *lwork, int *info);
int ssytri_check(char *uplo, int *n, float *a, int *lda, int *ipiv, float *work, int *info);
int ssytri2_check(char *uplo, int *n, float *a, int *lda, int *ipiv, float *work, int *lwork, int *info);
int ssytri2x_check(char *uplo, int *n, float *a, int *lda, int *ipiv, float *work, int *nb, int *info);
int ssytri_rook_check(char *uplo, int *n, float *a, int * lda, int *ipiv, float *work, int *info);
int ssytrs_check(char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);
int ssytrs2_check(char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, float *work, int *info);
int ssytrs_rook_check(char *uplo, int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int * info);
int stbcon_check(char *norm, char *uplo, char *diag, int *n, int *kd, float *ab, int *ldab, float *rcond, float *work, int *iwork, int *info);
int stbrfs_check(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);
int stbtrs_check(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);
int stfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, int *m, int *n, float *alpha, float *a, float *b, int *ldb);
int stftri_check(char *transr, char *uplo, char *diag, int *n, float *a, int *info);
int stfttp_check(char *transr, char *uplo, int *n, float *arf, float *ap, int *info);
int stfttr_check(char *transr, char *uplo, int *n, float *arf, float *a, int *lda, int *info);
int stgevc_check(char *side, char *howmny, logical *select, int *n, float *s, int *lds, float *p, int *ldp, float *vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, float *work, int *info);
int stgex2_check(logical *wantq, logical *wantz, int *n, float *a, int *lda, float *b, int *ldb, float *q, int *ldq, float * z__, int *ldz, int *j1, int *n1, int *n2, float *work, int *lwork, int *info);
int stgexc_check(logical *wantq, logical *wantz, int *n, float *a, int *lda, float *b, int *ldb, float *q, int *ldq, float * z__, int *ldz, int *ifst, int *ilst, float *work, int * lwork, int *info);
int stgsen_check(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, float *a, int *lda, float *b, int * ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, float *z__, int *ldz, int *m, float *pl, float *pr, float *dif, float *work, int *lwork, int *iwork, int *liwork, int * info);
int stgsja_check(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, float *a, int *lda, float *b, int *ldb, float *tola, float *tolb, float *alpha, float * beta, float *u, int *ldu, float *v, int *ldv, float *q, int * ldq, float *work, int *ncycle, int *info);
int stgsna_check(char *job, char *howmny, logical *select, int *n, float *a, int *lda, float *b, int *ldb, float *vl, int *ldvl, float *vr, int *ldvr, float *s, float *dif, int * mm, int *m, float *work, int *lwork, int *iwork, int * info);
int stgsy2_check(char *trans, int *ijob, int *m, int * n, float *a, int *lda, float *b, int *ldb, float *c__, int * ldc, float *d__, int *ldd, float *e, int *lde, float *f, int *ldf, float *scale, float *rdsum, float *rdscal, int *iwork, int *pq, int *info);
int stgsyl_check(char *trans, int *ijob, int *m, int * n, float *a, int *lda, float *b, int *ldb, float *c__, int * ldc, float *d__, int *ldd, float *e, int *lde, float *f, int *ldf, float *scale, float *dif, float *work, int *lwork, int * iwork, int *info);
int stpcon_check(char *norm, char *uplo, char *diag, int *n, float *ap, float *rcond, float *work, int *iwork, int *info);
int stpmqrt_check(char *side, char *trans, int *m, int *n, int *k, int *l, int *nb, float *v, int *ldv, float *t, int *ldt, float *a, int *lda, float *b, int *ldb, float * work, int *info);
int stpqrt_check(int *m, int *n, int *l, int *nb, float *a, int *lda, float *b, int *ldb, float *t, int *ldt, float *work, int *info);
int stpqrt2_check(int *m, int *n, int *l, float *a, int *lda, float *b, int *ldb, float *t, int *ldt, int * info);
int stprfb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, int *l, float *v, int *ldv, float *t, int *ldt, float *a, int *lda, float *b, int *ldb, float *work, int *ldwork);
int stprfs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, float *ap, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);
int stptri_check(char *uplo, char *diag, int *n, float *ap, int *info);
int stptrs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, float *ap, float *b, int *ldb, int *info);
int stpttf_check(char *transr, char *uplo, int *n, float *ap, float *arf, int *info);
int stpttr_check(char *uplo, int *n, float *ap, float *a, int *lda, int *info);
int strcon_check(char *norm, char *uplo, char *diag, int *n, float *a, int *lda, float *rcond, float *work, int *iwork, int *info);
int strevc_check(char *side, char *howmny, logical *select, int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, float *work, int *info);
int strexc_check(char *compq, int *n, float *t, int *ldt, float *q, int *ldq, int *ifst, int *ilst, float *work, int *info);
int strrfs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, int *iwork, int *info);
int strsen_check(char *job, char *compq, logical *select, int *n, float *t, int *ldt, float *q, int *ldq, float *wr, float *wi, int *m, float *s, float *sep, float *work, int *lwork, int * iwork, int *liwork, int *info);
int strsna_check(char *job, char *howmny, logical *select, int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, int *ldvr, float *s, float *sep, int *mm, int *m, float * work, int *ldwork, int *iwork, int *info);
int strsyl_check(char *trana, char *tranb, int *isgn, int *m, int *n, float *a, int *lda, float *b, int *ldb, float * c__, int *ldc, float *scale, int *info);
int strti2_check(char *uplo, char *diag, int *n, float *a, int *lda, int *info);
int strtri_check(char *uplo, char *diag, int *n, float *a, int *lda, int *info);
int strtrs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, float *a, int *lda, float *b, int *ldb, int * info);
int strttf_check(char *transr, char *uplo, int *n, float *a, int *lda, float *arf, int *info);
int strttp_check(char *uplo, int *n, float *a, int *lda, float *ap, int *info);
int stzrqf_check(int *m, int *n, float *a, int *lda, float *tau, int *info);
int stzrzf_check(int *m, int *n, float *a, int *lda, float *tau, float *work, int *lwork, int *info);
int zbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, int *m, int *p, int *q, double * theta, double *phi, dcomplex *u1, int *ldu1, dcomplex *u2, int *ldu2, dcomplex *v1t, int *ldv1t, dcomplex *v2t, int *ldv2t, double *b11d, double * b11e, double *b12d, double *b12e, double *b21d, double *b21e, double *b22d, double *b22e, double * rwork, int *lrwork, int *info);
int zbdsqr_check(char *uplo, int *n, int *ncvt, int * nru, int *ncc, double *d__, double *e, dcomplex *vt, int *ldvt, dcomplex *u, int *ldu, dcomplex *c__, int *ldc, double *rwork, int *info);
int zcgesv_check(int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, dcomplex *work, scomplex *swork, double *rwork, int *iter, int *info);
int zcposv_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *x, int *ldx, dcomplex *work, scomplex *swork, double *rwork, int *iter, int *info);
int zdrscl_check(int *n, double *sa, dcomplex *sx, int *incx);
int zgbbrd_check(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, dcomplex *ab, int *ldab, double *d__, double *e, dcomplex *q, int *ldq, dcomplex *pt, int *ldpt, dcomplex *c__, int *ldc, dcomplex *work, double *rwork, int *info);
int zgbcon_check(char *norm, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, int *ipiv, double *anorm, double *rcond, dcomplex *work, double *rwork, int * info);
int zgbequ_check(int *m, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int * info);
int zgbequb_check(int *m, int *n, int *kl, int * ku, dcomplex *ab, int *ldab, double *r__, double * c__, double *rowcnd, double *colcnd, double *amax, int *info);
int zgbrfs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex * afb, int *ldafb, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *info);
int zgbrfsx_check(char *trans, char *equed, int *n, int * kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *ipiv, double *r__, double *c__, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *berr, int * n_err_bnds__, double *err_bnds_norm__, double * err_bnds_comp__, int *nparams, double *params, dcomplex * work, double *rwork, int *info);
int zgbsv_check(int *n, int *kl, int *ku, int * nrhs, dcomplex *ab, int *ldab, int *ipiv, dcomplex * b, int *ldb, int *info);
int zgbsvx_check(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *ipiv, char *equed, double *r__, double *c__, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex *work, double *rwork, int * info);
int zgbsvxx_check(char *fact, char *trans, int *n, int * kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *ipiv, char *equed, double *r__, double *c__, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *params, dcomplex *work, double *rwork, int *info);
int zgbtf2_check(int *m, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, int *ipiv, int *info);
int zgbtrf_check(int *m, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, int *ipiv, int *info);
int zgbtrs_check(char *trans, int *n, int *kl, int * ku, int *nrhs, dcomplex *ab, int *ldab, int *ipiv, dcomplex *b, int *ldb, int *info);
int zgebak_check(char *job, char *side, int *n, int *ilo, int *ihi, double *scale, int *m, dcomplex *v, int *ldv, int *info);
int zgebal_check(char *job, int *n, dcomplex *a, int *lda, int *ilo, int *ihi, double *scale, int *info);
int zgebd2_check(int *m, int *n, dcomplex *a, int *lda, double *d__, double *e, dcomplex *tauq, dcomplex *taup, dcomplex *work, int *info);
int zgebrd_check(int *m, int *n, dcomplex *a, int *lda, double *d__, double *e, dcomplex *tauq, dcomplex *taup, dcomplex *work, int *lwork, int * info);
int zgecon_check(char *norm, int *n, dcomplex *a, int *lda, double *anorm, double *rcond, dcomplex * work, double *rwork, int *info);
int zgeequ_check(int *m, int *n, dcomplex *a, int *lda, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int *info);
int zgeequb_check(int *m, int *n, dcomplex *a, int *lda, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, int *info);
int zgees_check(char *jobvs, char *sort, L_fp select, int *n, dcomplex *a, int *lda, int *sdim, dcomplex *w, dcomplex *vs, int *ldvs, dcomplex *work, int *lwork, double *rwork, logical *bwork, int *info);
int zgeesx_check(char *jobvs, char *sort, L_fp select, char * sense, int *n, dcomplex *a, int *lda, int *sdim, dcomplex *w, dcomplex *vs, int *ldvs, double * rconde, double *rcondv, dcomplex *work, int *lwork, double *rwork, logical *bwork, int *info);
int zgeev_check(char *jobvl, char *jobvr, int *n, dcomplex *a, int *lda, dcomplex *w, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, dcomplex *work, int *lwork, double *rwork, int *info);
int zgeevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, dcomplex *a, int *lda, dcomplex *w, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, int *ilo, int *ihi, double *scale, double *abnrm, double *rconde, double *rcondv, dcomplex *work, int * lwork, double *rwork, int *info);
int zgegs_check(char *jobvsl, char *jobvsr, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vsl, int *ldvsl, dcomplex *vsr, int *ldvsr, dcomplex * work, int *lwork, double *rwork, int *info);
int zgegv_check(char *jobvl, char *jobvr, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, dcomplex *work, int *lwork, double *rwork, int *info);
int zgehd2_check(int *n, int *ilo, int *ihi, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *info);
int zgehrd_check(int *n, int *ilo, int *ihi, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info);
int zgelq2_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *info);
int zgelqf_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zgels_check(char *trans, int *m, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *work, int *lwork, int *info);
int zgelsd_check(int *m, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *s, double *rcond, int *rank, dcomplex *work, int *lwork, double *rwork, int *iwork, int *info);
int zgelss_check(int *m, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *s, double *rcond, int *rank, dcomplex *work, int *lwork, double *rwork, int *info);
int zgelsx_check(int *m, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *jpvt, double *rcond, int *rank, dcomplex *work, double *rwork, int *info);
int zgelsy_check(int *m, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *jpvt, double *rcond, int *rank, dcomplex *work, int *lwork, double *rwork, int *info);
int zgemqrt_check(char *side, char *trans, int *m, int *n, int *k, int *nb, dcomplex *v, int *ldv, dcomplex *t, int *ldt, dcomplex *c__, int *ldc, dcomplex *work, int *info);
int zgeql2_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *info);
int zgeqlf_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zgeqp3_check(int *m, int *n, dcomplex *a, int *lda, int *jpvt, dcomplex *tau, dcomplex *work, int *lwork, double *rwork, int *info);
int zgeqpf_check(int *m, int *n, dcomplex *a, int *lda, int *jpvt, dcomplex *tau, dcomplex *work, double *rwork, int *info);
int zgeqr2_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *info);
int zgeqr2p_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *info);
int zgeqrf_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zgeqrfp_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zgeqrt_check(int *m, int *n, int *nb, dcomplex *a, int *lda, dcomplex *t, int *ldt, dcomplex *work, int *info);
int zgeqrt2_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *t, int *ldt, int *info);
int zgeqrt3_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *t, int *ldt, int *info);
int zgerfs_check(char *trans, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *info);
int zgerfsx_check(char *trans, char *equed, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, double *r__, double *c__, dcomplex * b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *params, dcomplex *work, double *rwork, int *info);
int zgerq2_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *info);
int zgerqf_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zgesc2_check(int *n, dcomplex *a, int *lda, dcomplex *rhs, int *ipiv, int *jpiv, double *scale);
int zgesdd_check(char *jobz, int *m, int *n, dcomplex *a, int *lda, double *s, dcomplex *u, int *ldu, dcomplex *vt, int *ldvt, dcomplex *work, int *lwork, double *rwork, int *iwork, int *info);
int zgesv_check(int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, int * info);
int zgesvd_check(char *jobu, char *jobvt, int *m, int *n, dcomplex *a, int *lda, double *s, dcomplex *u, int *ldu, dcomplex *vt, int *ldvt, dcomplex *work, int *lwork, double *rwork, int *info);
int zgesvx_check(char *fact, char *trans, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, char *equed, double *r__, double *c__, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex * work, double *rwork, int *info);
int zgesvxx_check(char *fact, char *trans, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, char *equed, double *r__, double *c__, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *rpvgrw, double *berr, int * n_err_bnds__, double *err_bnds_norm__, double * err_bnds_comp__, int *nparams, double *params, dcomplex * work, double *rwork, int *info);
int zgetc2_check(int *n, dcomplex *a, int *lda, int *ipiv, int *jpiv, int *info);
int zgetf2_check(int *m, int *n, dcomplex *a, int *lda, int *ipiv, int *info);
int zgetrf_check(int *m, int *n, dcomplex *a, int *lda, int *ipiv, int *info);
int zgetrfnp_check(int *m, int *n, dcomplex *a, int *lda, int *info);
int zgetri_check(int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *lwork, int *info);
int zgetrs_check(char *trans, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *info);
int zggbak_check(char *job, char *side, int *n, int *ilo, int *ihi, double *lscale, double *rscale, int *m, dcomplex *v, int *ldv, int *info);
int zggbal_check(char *job, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *ilo, int *ihi, double *lscale, double *rscale, double *work, int * info);
int zgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *sdim, dcomplex *alpha, dcomplex * beta, dcomplex *vsl, int *ldvsl, dcomplex *vsr, int *ldvsr, dcomplex *work, int *lwork, double *rwork, logical *bwork, int *info);
int zggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *sdim, dcomplex *alpha, dcomplex *beta, dcomplex *vsl, int *ldvsl, dcomplex *vsr, int *ldvsr, double *rconde, double * rcondv, dcomplex *work, int *lwork, double *rwork, int *iwork, int *liwork, logical *bwork, int *info);
int zggev_check(char *jobvl, char *jobvr, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, dcomplex *work, int *lwork, double *rwork, int *info);
int zggevx_check(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, int *ilo, int *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm, double *rconde, double * rcondv, dcomplex *work, int *lwork, double *rwork, int *iwork, logical *bwork, int *info);
int zggglm_check(int *n, int *m, int *p, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *d__, dcomplex *x, dcomplex *y, dcomplex *work, int *lwork, int *info);
int zgghrd_check(char *compq, char *compz, int *n, int * ilo, int *ihi, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *q, int *ldq, dcomplex *z__, int *ldz, int *info);
int zgglse_check(int *m, int *n, int *p, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *c__, dcomplex *d__, dcomplex *x, dcomplex *work, int *lwork, int *info);
int zggqrf_check(int *n, int *m, int *p, dcomplex *a, int *lda, dcomplex *taua, dcomplex *b, int *ldb, dcomplex *taub, dcomplex *work, int * lwork, int *info);
int zggrqf_check(int *m, int *p, int *n, dcomplex *a, int *lda, dcomplex *taua, dcomplex *b, int *ldb, dcomplex *taub, dcomplex *work, int * lwork, int *info);
int zggsvd_check(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *alpha, double *beta, dcomplex *u, int *ldu, dcomplex *v, int *ldv, dcomplex *q, int *ldq, dcomplex *work, double *rwork, int *iwork, int *info);
int zggsvp_check(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *tola, double *tolb, int *k, int *l, dcomplex *u, int *ldu, dcomplex *v, int *ldv, dcomplex *q, int *ldq, int *iwork, double * rwork, dcomplex *tau, dcomplex *work, int *info);
int zgtcon_check(char *norm, int *n, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2, int * ipiv, double *anorm, double *rcond, dcomplex *work, int *info);
int zgtrfs_check(char *trans, int *n, int *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *dlf, dcomplex *df, dcomplex *duf, dcomplex *du2, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *info);
int zgtsv_check(int *n, int *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *b, int *ldb, int *info);
int zgtsvx_check(char *fact, char *trans, int *n, int * nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *dlf, dcomplex *df, dcomplex *duf, dcomplex *du2, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex *work, double *rwork, int * info);
int zgttrf_check(int *n, dcomplex *dl, dcomplex * d__, dcomplex *du, dcomplex *du2, int *ipiv, int * info);
int zgttrs_check(char *trans, int *n, int *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2, int *ipiv, dcomplex *b, int *ldb, int *info);
int zgtts2_check(int *itrans, int *n, int *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2, int *ipiv, dcomplex *b, int *ldb);
int zhbev_check(char *jobz, char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, int *info);
int zhbevd_check(char *jobz, char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
int zhbevx_check(char *jobz, char *range, char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, dcomplex *q, int *ldq, double *vl, double *vu, int *il, int * iu, double *abstol, int *m, double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, int *iwork, int *ifail, int *info);
int zhbgst_check(char *vect, char *uplo, int *n, int *ka, int *kb, dcomplex *ab, int *ldab, dcomplex *bb, int *ldbb, dcomplex *x, int *ldx, dcomplex *work, double *rwork, int *info);
int zhbgv_check(char *jobz, char *uplo, int *n, int *ka, int *kb, dcomplex *ab, int *ldab, dcomplex *bb, int *ldbb, double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, int *info);
int zhbgvd_check(char *jobz, char *uplo, int *n, int *ka, int *kb, dcomplex *ab, int *ldab, dcomplex *bb, int *ldbb, double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double *rwork, int * lrwork, int *iwork, int *liwork, int *info);
int zhbgvx_check(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, dcomplex *ab, int *ldab, dcomplex *bb, int *ldbb, dcomplex *q, int *ldq, double *vl, double *vu, int *il, int *iu, double * abstol, int *m, double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, int *iwork, int * ifail, int *info);
int zhbtrd_check(char *vect, char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, double *d__, double *e, dcomplex *q, int *ldq, dcomplex *work, int *info);
int zhecon_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, double *anorm, double *rcond, dcomplex *work, int *info);
int zhecon_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, double *anorm, double *rcond, dcomplex *work, int *info);
int zheequb_check(char *uplo, int *n, dcomplex *a, int *lda, double *s, double *scond, double *amax, dcomplex *work, int *info);
int zheev_check(char *jobz, char *uplo, int *n, dcomplex *a, int *lda, double *w, dcomplex *work, int *lwork, double *rwork, int *info);
int zheevd_check(char *jobz, char *uplo, int *n, dcomplex *a, int *lda, double *w, dcomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
int zheevr_check(char *jobz, char *range, char *uplo, int *n, dcomplex *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double * w, dcomplex *z__, int *ldz, int *isuppz, dcomplex * work, int *lwork, double *rwork, int *lrwork, int * iwork, int *liwork, int *info);
int zheevx_check(char *jobz, char *range, char *uplo, int *n, dcomplex *a, int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double * w, dcomplex *z__, int *ldz, dcomplex *work, int * lwork, double *rwork, int *iwork, int *ifail, int * info);
int zhegs2_check(int *itype, char *uplo, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *info);
int zhegst_check(int *itype, char *uplo, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *info);
int zhegv_check(int *itype, char *jobz, char *uplo, int * n, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *w, dcomplex *work, int *lwork, double *rwork, int *info);
int zhegvd_check(int *itype, char *jobz, char *uplo, int * n, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *w, dcomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
int zhegvx_check(int *itype, char *jobz, char *range, char * uplo, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *vl, double *vu, int *il, int * iu, double *abstol, int *m, double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double *rwork, int *iwork, int *ifail, int *info);
int zherfs_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *info);
int zherfsx_check(char *uplo, char *equed, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double * err_bnds_comp__, int *nparams, double *params, dcomplex * work, double *rwork, int *info);
int zhesv_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, int *lwork, int *info);
int zhesv_rook_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, int *lwork, int *info);
int zhesvx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex *work, int *lwork, double *rwork, int *info);
int zhesvxx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int * nparams, double *params, dcomplex *work, double *rwork, int *info);
int zheswapr_check(char *uplo, int *n, dcomplex *a, int *lda, int *i1, int *i2);
int zhetd2_check(char *uplo, int *n, dcomplex *a, int *lda, double *d__, double *e, dcomplex *tau, int *info);
int zhetf2_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, int *info);
int zhetf2_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, int *info);
int zhetrd_check(char *uplo, int *n, dcomplex *a, int *lda, double *d__, double *e, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zhetrf_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *lwork, int *info);
int zhetrf_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *lwork, int *info);
int zhetri_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *info);
int zhetri2_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *lwork, int *info);
int zhetri2x_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *nb, int *info);
int zhetri_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *info);
int zhetrs_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *info);
int zhetrs2_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, int *info);
int zhetrs_rook_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *info);
int zhfrk_check(char *transr, char *uplo, char *trans, int *n, int *k, double *alpha, dcomplex *a, int *lda, double *beta, dcomplex *c__);
int zhgeqz_check(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *t, int *ldt, dcomplex *alpha, dcomplex * beta, dcomplex *q, int *ldq, dcomplex *z__, int * ldz, dcomplex *work, int *lwork, double *rwork, int * info);
int zhpcon_check(char *uplo, int *n, dcomplex *ap, int *ipiv, double *anorm, double *rcond, dcomplex * work, int *info);
int zhpev_check(char *jobz, char *uplo, int *n, dcomplex *ap, double *w, dcomplex *z__, int *ldz, dcomplex * work, double *rwork, int *info);
int zhpevd_check(char *jobz, char *uplo, int *n, dcomplex *ap, double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double *rwork, int * lrwork, int *iwork, int *liwork, int *info);
int zhpevx_check(char *jobz, char *range, char *uplo, int *n, dcomplex *ap, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, dcomplex *z__, int *ldz, dcomplex *work, double * rwork, int *iwork, int *ifail, int *info);
int zhpgst_check(int *itype, char *uplo, int *n, dcomplex *ap, dcomplex *bp, int *info);
int zhpgv_check(int *itype, char *jobz, char *uplo, int * n, dcomplex *ap, dcomplex *bp, double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, int * info);
int zhpgvd_check(int *itype, char *jobz, char *uplo, int * n, dcomplex *ap, dcomplex *bp, double *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double * rwork, int *lrwork, int *iwork, int *liwork, int * info);
int zhpgvx_check(int *itype, char *jobz, char *range, char * uplo, int *n, dcomplex *ap, dcomplex *bp, double * vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, dcomplex *z__, int *ldz, dcomplex *work, double *rwork, int *iwork, int * ifail, int *info);
int zhprfs_check(char *uplo, int *n, int *nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex * b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int * info);
int zhpsv_check(char *uplo, int *n, int *nrhs, dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, int *info);
int zhpsvx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex * work, double *rwork, int *info);
int zhptrd_check(char *uplo, int *n, dcomplex *ap, double *d__, double *e, dcomplex *tau, int *info);
int zhptrf_check(char *uplo, int *n, dcomplex *ap, int *ipiv, int *info);
int zhptri_check(char *uplo, int *n, dcomplex *ap, int *ipiv, dcomplex *work, int *info);
int zhptrs_check(char *uplo, int *n, int *nrhs, dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, int *info);
int zhsein_check(char *side, char *eigsrc, char *initv, logical * select, int *n, dcomplex *h__, int *ldh, dcomplex * w, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, int *mm, int *m, dcomplex *work, double *rwork, int *ifaill, int *ifailr, int *info);
int zhseqr_check(char *job, char *compz, int *n, int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *w, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, int *info);
int zla_gbamv_check(int *trans, int *m, int *n, int *kl, int *ku, double *alpha, dcomplex *ab, int *ldab, dcomplex *x, int *incx, double *beta, double *y, int *incy);
double zla_gbrcond_c_check(char *trans, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *ipiv, double *c__, logical *capply, int *info, dcomplex *work, double *rwork);
double zla_gbrcond_x_check(char *trans, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *ipiv, dcomplex *x, int *info, dcomplex *work, double *rwork);
int zla_gbrfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *kl, int *ku, int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, int *ipiv, logical *colequ, double *c__, dcomplex *b, int *ldb, dcomplex *y, int *ldy, double *berr_out__, int *n_norms__, double *err_bnds_norm__, double * err_bnds_comp__, dcomplex *res, double *ayb, dcomplex * dy, dcomplex *y_tail__, double *rcond, int *ithresh, double *rthresh, double *dz_ub__, logical *ignore_cwise__, int *info);
double zla_gbrpvgrw_check(int *n, int *kl, int *ku, int * ncols, dcomplex *ab, int *ldab, dcomplex *afb, int * ldafb);
int zla_geamv_check(int *trans, int *m, int *n, double *alpha, dcomplex *a, int *lda, dcomplex *x, int *incx, double *beta, double *y, int *incy);
double zla_gercond_c_check(char *trans, int *n, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, double * c__, logical *capply, int *info, dcomplex *work, double * rwork);
double zla_gercond_x_check(char *trans, int *n, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex * x, int *info, dcomplex *work, double *rwork);
int zla_gerfsx_extended_check(int *prec_type__, int * trans_type__, int *n, int *nrhs, dcomplex *a, int * lda, dcomplex *af, int *ldaf, int *ipiv, logical *colequ, double *c__, dcomplex *b, int *ldb, dcomplex *y, int *ldy, double *berr_out__, int *n_norms__, double * errs_n__, double *errs_c__, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond, int *ithresh, double *rthresh, double *dz_ub__, logical * ignore_cwise__, int *info);
double zla_gerpvgrw_check(int *n, int *ncols, dcomplex *a, int *lda, dcomplex *af, int *ldaf);
int zla_heamv_check(int *uplo, int *n, double *alpha, dcomplex *a, int *lda, dcomplex *x, int *incx, double *beta, double *y, int *incy);
double zla_hercond_c_check(char *uplo, int *n, dcomplex *a, int * lda, dcomplex *af, int *ldaf, int *ipiv, double *c__, logical *capply, int *info, dcomplex *work, double * rwork);
double zla_hercond_x_check(char *uplo, int *n, dcomplex *a, int * lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex * x, int *info, dcomplex *work, double *rwork);
int zla_herfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, logical *colequ, double *c__, dcomplex *b, int *ldb, dcomplex *y, int *ldy, double *berr_out__, int *n_norms__, double * err_bnds_norm__, double *err_bnds_comp__, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond, int *ithresh, double *rthresh, double * dz_ub__, logical *ignore_cwise__, int *info);
double zla_herpvgrw_check(char *uplo, int *n, int *info, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, double *work);
int zla_lin_berr_check(int *n, int *nz, int *nrhs, dcomplex *res, double *ayb, double *berr);
double zla_porcond_c_check(char *uplo, int *n, dcomplex *a, int * lda, dcomplex *af, int *ldaf, double *c__, logical * capply, int *info, dcomplex *work, double *rwork);
double zla_porcond_x_check(char *uplo, int *n, dcomplex *a, int * lda, dcomplex *af, int *ldaf, dcomplex *x, int * info, dcomplex *work, double *rwork);
int zla_porfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, logical *colequ, double *c__, dcomplex *b, int *ldb, dcomplex *y, int *ldy, double *berr_out__, int *n_norms__, double * err_bnds_norm__, double *err_bnds_comp__, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond, int *ithresh, double *rthresh, double * dz_ub__, logical *ignore_cwise__, int *info);
double zla_porpvgrw_check(char *uplo, int *ncols, dcomplex *a, int *lda, dcomplex *af, int *ldaf, double *work);
int zla_syamv_check(int *uplo, int *n, double *alpha, dcomplex *a, int *lda, dcomplex *x, int *incx, double *beta, double *y, int *incy);
double zla_syrcond_c_check(char *uplo, int *n, dcomplex *a, int * lda, dcomplex *af, int *ldaf, int *ipiv, double *c__, logical *capply, int *info, dcomplex *work, double * rwork);
double zla_syrcond_x_check(char *uplo, int *n, dcomplex *a, int * lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex * x, int *info, dcomplex *work, double *rwork);
int zla_syrfsx_extended_check(int *prec_type__, char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, logical *colequ, double *c__, dcomplex *b, int *ldb, dcomplex *y, int *ldy, double *berr_out__, int *n_norms__, double * err_bnds_norm__, double *err_bnds_comp__, dcomplex *res, double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond, int *ithresh, double *rthresh, double * dz_ub__, logical *ignore_cwise__, int *info);
double zla_syrpvgrw_check(char *uplo, int *n, int *info, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, double *work);
int zla_wwaddw_check(int *n, dcomplex *x, dcomplex *y, dcomplex *w);
int zlabrd_check(int *m, int *n, int *nb, dcomplex *a, int *lda, double *d__, double *e, dcomplex *tauq, dcomplex *taup, dcomplex *x, int * ldx, dcomplex *y, int *ldy);
int zlacgv_check(int *n, dcomplex *x, int *incx);
int zlacn2_check(int *n, dcomplex *v, dcomplex *x, double *est, int *kase, int *isave);
int zlacon_check(int *n, dcomplex *v, dcomplex *x, double *est, int *kase);
int zlacp2_check(char *uplo, int *m, int *n, double * a, int *lda, dcomplex *b, int *ldb);
int zlacpy_check(char *uplo, int *m, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb);
int zlacrm_check(int *m, int *n, dcomplex *a, int *lda, double *b, int *ldb, dcomplex *c__, int *ldc, double *rwork);
int zlacrt_check(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy, dcomplex *c__, dcomplex * s);
VOID zladiv_check(dcomplex * ret_val, dcomplex *x, dcomplex *y);
int zlaed0_check(int *qsiz, int *n, double *d__, double *e, dcomplex *q, int *ldq, dcomplex *qstore, int *ldqs, double *rwork, int *iwork, int *info);
int zlaed7_check(int *n, int *cutpnt, int *qsiz, int *tlvls, int *curlvl, int *curpbm, double *d__, dcomplex *q, int *ldq, double *rho, int *indxq, double *qstore, int *qptr, int *prmptr, int *perm, int *givptr, int *givcol, double *givnum, dcomplex * work, double *rwork, int *iwork, int *info);
int zlaed8_check(int *k, int *n, int *qsiz, dcomplex *q, int *ldq, double *d__, double *rho, int *cutpnt, double *z__, double *dlamda, dcomplex * q2, int *ldq2, double *w, int *indxp, int *indx, int *indxq, int *perm, int *givptr, int *givcol, double *givnum, int *info);
int zlaein_check(logical *rightv, logical *noinit, int *n, dcomplex *h__, int *ldh, dcomplex *w, dcomplex *v, dcomplex *b, int *ldb, double *rwork, double *eps3, double *smlnum, int *info);
int zlaesy_check(dcomplex *a, dcomplex *b, dcomplex *c__, dcomplex *rt1, dcomplex *rt2, dcomplex *evscal, dcomplex *cs1, dcomplex *sn1);
int zlaev2_check(dcomplex *a, dcomplex *b, dcomplex *c__, double *rt1, double *rt2, double *cs1, dcomplex *sn1);
int zlag2c_check(int *m, int *n, dcomplex *a, int *lda, scomplex *sa, int *ldsa, int *info);
int zlags2_check(logical *upper, double *a1, dcomplex * a2, double *a3, double *b1, dcomplex *b2, double *b3, double *csu, dcomplex *snu, double *csv, dcomplex * snv, double *csq, dcomplex *snq);
int zlagtm_check(char *trans, int *n, int *nrhs, double *alpha, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *x, int *ldx, double *beta, dcomplex *b, int *ldb);
int zlahef_check(char *uplo, int *n, int *nb, int *kb, dcomplex *a, int *lda, int *ipiv, dcomplex *w, int *ldw, int *info);
int zlahef_rook_check(char *uplo, int *n, int *nb, int *kb, dcomplex *a, int *lda, int *ipiv, dcomplex *w, int *ldw, int *info);
int zlahqr_check(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *w, int *iloz, int *ihiz, dcomplex *z__, int *ldz, int *info);
int zlahr2_check(int *n, int *k, int *nb, dcomplex *a, int *lda, dcomplex *tau, dcomplex *t, int *ldt, dcomplex *y, int *ldy);
int zlahrd_check(int *n, int *k, int *nb, dcomplex *a, int *lda, dcomplex *tau, dcomplex *t, int *ldt, dcomplex *y, int *ldy);
int zlaic1_check(int *job, int *j, dcomplex *x, double *sest, dcomplex *w, dcomplex *gamma, double * sestpr, dcomplex *s, dcomplex *c__);
int zlals0_check(int *icompq, int *nl, int *nr, int *sqre, int *nrhs, dcomplex *b, int *ldb, dcomplex *bx, int *ldbx, int *perm, int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, double *poles, double *difl, double *difr, double * z__, int *k, double *c__, double *s, double *rwork, int *info);
int zlalsa_check(int *icompq, int *smlsiz, int *n, int *nrhs, dcomplex *b, int *ldb, dcomplex *bx, int *ldbx, double *u, int *ldu, double *vt, int * k, double *difl, double *difr, double *z__, double * poles, int *givptr, int *givcol, int *ldgcol, int * perm, double *givnum, double *c__, double *s, double * rwork, int *iwork, int *info);
int zlalsd_check(char *uplo, int *smlsiz, int *n, int *nrhs, double *d__, double *e, dcomplex *b, int *ldb, double *rcond, int *rank, dcomplex *work, double * rwork, int *iwork, int *info);
double zlangb_check(char *norm, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, double *work);
double zlange_check(char *norm, int *m, int *n, dcomplex *a, int *lda, double *work);
double zlangt_check(char *norm, int *n, dcomplex *dl, dcomplex * d__, dcomplex *du);
double zlanhb_check(char *norm, char *uplo, int *n, int *k, dcomplex *ab, int *ldab, double *work);
double zlanhe_check(char *norm, char *uplo, int *n, dcomplex *a, int *lda, double *work);
double zlanhf_check(char *norm, char *transr, char *uplo, int *n, dcomplex *a, double *work);
double zlanhp_check(char *norm, char *uplo, int *n, dcomplex *ap, double *work);
double zlanhs_check(char *norm, int *n, dcomplex *a, int *lda, double *work);
double zlanht_check(char *norm, int *n, double *d__, dcomplex *e);
double zlansb_check(char *norm, char *uplo, int *n, int *k, dcomplex *ab, int *ldab, double *work);
double zlansp_check(char *norm, char *uplo, int *n, dcomplex *ap, double *work);
double zlansy_check(char *norm, char *uplo, int *n, dcomplex *a, int *lda, double *work);
double zlantb_check(char *norm, char *uplo, char *diag, int *n, int *k, dcomplex *ab, int *ldab, double *work);
double zlantp_check(char *norm, char *uplo, char *diag, int *n, dcomplex *ap, double *work);
double zlantr_check(char *norm, char *uplo, char *diag, int *m, int *n, dcomplex *a, int *lda, double *work);
int zlapll_check(int *n, dcomplex *x, int *incx, dcomplex *y, int *incy, double *ssmin);
int zlapmr_check(logical *forwrd, int *m, int *n, dcomplex *x, int *ldx, int *k);
int zlapmt_check(logical *forwrd, int *m, int *n, dcomplex *x, int *ldx, int *k);
int zlaqgb_check(int *m, int *n, int *kl, int *ku, dcomplex *ab, int *ldab, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, char *equed);
int zlaqge_check(int *m, int *n, dcomplex *a, int *lda, double *r__, double *c__, double *rowcnd, double *colcnd, double *amax, char *equed);
int zlaqhb_check(char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, double *s, double *scond, double *amax, char *equed);
int zlaqhe_check(char *uplo, int *n, dcomplex *a, int *lda, double *s, double *scond, double *amax, char *equed);
int zlaqhp_check(char *uplo, int *n, dcomplex *ap, double *s, double *scond, double *amax, char *equed);
int zlaqp2_check(int *m, int *n, int *offset, dcomplex *a, int *lda, int *jpvt, dcomplex *tau, double *vn1, double *vn2, dcomplex *work);
int zlaqps_check(int *m, int *n, int *offset, int *nb, int *kb, dcomplex *a, int *lda, int *jpvt, dcomplex *tau, double *vn1, double *vn2, dcomplex * auxv, dcomplex *f, int *ldf);
int zlaqr0_check(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *w, int *iloz, int *ihiz, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, int *info);
int zlaqr1_check(int *n, dcomplex *h__, int *ldh, dcomplex *s1, dcomplex *s2, dcomplex *v);
int zlaqr2_check(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, dcomplex *h__, int *ldh, int *iloz, int *ihiz, dcomplex *z__, int *ldz, int *ns, int *nd, dcomplex *sh, dcomplex *v, int *ldv, int *nh, dcomplex *t, int *ldt, int *nv, dcomplex *wv, int *ldwv, dcomplex *work, int *lwork);
int zlaqr3_check(logical *wantt, logical *wantz, int *n, int *ktop, int *kbot, int *nw, dcomplex *h__, int *ldh, int *iloz, int *ihiz, dcomplex *z__, int *ldz, int *ns, int *nd, dcomplex *sh, dcomplex *v, int *ldv, int *nh, dcomplex *t, int *ldt, int *nv, dcomplex *wv, int *ldwv, dcomplex *work, int *lwork);
int zlaqr4_check(logical *wantt, logical *wantz, int *n, int *ilo, int *ihi, dcomplex *h__, int *ldh, dcomplex *w, int *iloz, int *ihiz, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, int *info);
int zlaqr5_check(logical *wantt, logical *wantz, int *kacc22, int *n, int *ktop, int *kbot, int *nshfts, dcomplex *s, dcomplex *h__, int *ldh, int *iloz, int *ihiz, dcomplex *z__, int *ldz, dcomplex *v, int *ldv, dcomplex *u, int *ldu, int *nv, dcomplex *wv, int *ldwv, int *nh, dcomplex *wh, int *ldwh);
int zlaqsb_check(char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, double *s, double *scond, double *amax, char *equed);
int zlaqsp_check(char *uplo, int *n, dcomplex *ap, double *s, double *scond, double *amax, char *equed);
int zlaqsy_check(char *uplo, int *n, dcomplex *a, int *lda, double *s, double *scond, double *amax, char *equed);
int zlar1v_check(int *n, int *b1, int *bn, double *lambda, double *d__, double *l, double *ld, double * lld, double *pivmin, double *gaptol, dcomplex *z__, logical *wantnc, int *negcnt, double *ztz, double *mingma, int *r__, int *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);
int zlar2v_check(int *n, dcomplex *x, dcomplex *y, dcomplex *z__, int *incx, double *c__, dcomplex *s, int *incc);
int zlarcm_check(int *m, int *n, double *a, int * lda, dcomplex *b, int *ldb, dcomplex *c__, int *ldc, double *rwork);
int zlarf_check(char *side, int *m, int *n, dcomplex *v, int *incv, dcomplex *tau, dcomplex *c__, int * ldc, dcomplex *work);
int zlarfb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, dcomplex *v, int *ldv, dcomplex *t, int *ldt, dcomplex *c__, int * ldc, dcomplex *work, int *ldwork);
int zlarfg_check(int *n, dcomplex *alpha, dcomplex * x, int *incx, dcomplex *tau);
int zlarfgp_check(int *n, dcomplex *alpha, dcomplex *x, int *incx, dcomplex *tau);
int zlarft_check(char *direct, char *storev, int *n, int * k, dcomplex *v, int *ldv, dcomplex *tau, dcomplex * t, int *ldt);
int zlarfx_check(char *side, int *m, int *n, dcomplex *v, dcomplex *tau, dcomplex *c__, int * ldc, dcomplex *work);
int zlargv_check(int *n, dcomplex *x, int *incx, dcomplex *y, int *incy, double *c__, int *incc);
int zlarnv_check(int *idist, int *iseed, int *n, dcomplex *x);
int zlarrv_check(int *n, double *vl, double *vu, double *d__, double *l, double *pivmin, int *isplit, int *m, int *dol, int *dou, double *minrgp, double *rtol1, double *rtol2, double *w, double *werr, double *wgap, int *iblock, int *indexw, double *gers, dcomplex *z__, int *ldz, int *isuppz, double *work, int *iwork, int *info);
int zlarscl2_check(int *m, int *n, double *d__, dcomplex *x, int *ldx);
int zlartg_check(dcomplex *f, dcomplex *g, double * cs, dcomplex *sn, dcomplex *r__);
int zlartv_check(int *n, dcomplex *x, int *incx, dcomplex *y, int *incy, double *c__, dcomplex *s, int *incc);
int zlarz_check(char *side, int *m, int *n, int *l, dcomplex *v, int *incv, dcomplex *tau, dcomplex * c__, int *ldc, dcomplex *work);
int zlarzb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, int *l, dcomplex *v, int *ldv, dcomplex *t, int *ldt, dcomplex *c__, int *ldc, dcomplex *work, int *ldwork);
int zlarzt_check(char *direct, char *storev, int *n, int * k, dcomplex *v, int *ldv, dcomplex *tau, dcomplex * t, int *ldt);
int zlascl_check(char *type__, int *kl, int *ku, double *cfrom, double *cto, int *m, int *n, dcomplex *a, int *lda, int *info);
int zlascl2_check(int *m, int *n, double *d__, dcomplex *x, int *ldx);
int zlaset_check(char *uplo, int *m, int *n, dcomplex *alpha, dcomplex *beta, dcomplex *a, int * lda);
int zlasr_check(char *side, char *pivot, char *direct, int *m, int *n, double *c__, double *s, dcomplex *a, int *lda);
int zlassq_check(int *n, dcomplex *x, int *incx, double *scale, double *sumsq);
int zlaswp_check(int *n, dcomplex *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);
int zlasyf_check(char *uplo, int *n, int *nb, int *kb, dcomplex *a, int *lda, int *ipiv, dcomplex *w, int *ldw, int *info);
int zlasyf_rook_check(char *uplo, int *n, int *nb, int *kb, dcomplex *a, int *lda, int *ipiv, dcomplex *w, int *ldw, int *info);
int zlat2c_check(char *uplo, int *n, dcomplex *a, int *lda, scomplex *sa, int *ldsa, int *info);
int zlatbs_check(char *uplo, char *trans, char *diag, char * normin, int *n, int *kd, dcomplex *ab, int *ldab, dcomplex *x, double *scale, double *cnorm, int *info);
int zlatdf_check(int *ijob, int *n, dcomplex *z__, int *ldz, dcomplex *rhs, double *rdsum, double * rdscal, int *ipiv, int *jpiv);
int zlatps_check(char *uplo, char *trans, char *diag, char * normin, int *n, dcomplex *ap, dcomplex *x, double * scale, double *cnorm, int *info);
int zlatrd_check(char *uplo, int *n, int *nb, dcomplex *a, int *lda, double *e, dcomplex *tau, dcomplex *w, int *ldw);
int zlatrs_check(char *uplo, char *trans, char *diag, char * normin, int *n, dcomplex *a, int *lda, dcomplex *x, double *scale, double *cnorm, int *info);
int zlatrz_check(int *m, int *n, int *l, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work);
int zlatzm_check(char *side, int *m, int *n, dcomplex *v, int *incv, dcomplex *tau, dcomplex * c1, dcomplex *c2, int *ldc, dcomplex *work);
int zlauu2_check(char *uplo, int *n, dcomplex *a, int *lda, int *info);
int zlauum_check(char *uplo, int *n, dcomplex *a, int *lda, int *info);
int zpbcon_check(char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, double *anorm, double * rcond, dcomplex *work, double *rwork, int *info);
int zpbequ_check(char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, double *s, double *scond, double *amax, int *info);
int zpbrfs_check(char *uplo, int *n, int *kd, int * nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int * ldafb, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double * rwork, int *info);
int zpbstf_check(char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, int *info);
int zpbsv_check(char *uplo, int *n, int *kd, int * nrhs, dcomplex *ab, int *ldab, dcomplex *b, int * ldb, int *info);
int zpbsvx_check(char *fact, char *uplo, int *n, int *kd, int *nrhs, dcomplex *ab, int *ldab, dcomplex *afb, int *ldafb, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double * ferr, double *berr, dcomplex *work, double *rwork, int *info);
int zpbtf2_check(char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, int *info);
int zpbtrf_check(char *uplo, int *n, int *kd, dcomplex *ab, int *ldab, int *info);
int zpbtrs_check(char *uplo, int *n, int *kd, int * nrhs, dcomplex *ab, int *ldab, dcomplex *b, int * ldb, int *info);
int zpftrf_check(char *transr, char *uplo, int *n, dcomplex *a, int *info);
int zpftri_check(char *transr, char *uplo, int *n, dcomplex *a, int *info);
int zpftrs_check(char *transr, char *uplo, int *n, int * nrhs, dcomplex *a, dcomplex *b, int *ldb, int *info);
int zpocon_check(char *uplo, int *n, dcomplex *a, int *lda, double *anorm, double *rcond, dcomplex * work, double *rwork, int *info);
int zpoequ_check(int *n, dcomplex *a, int *lda, double *s, double *scond, double *amax, int *info);
int zpoequb_check(int *n, dcomplex *a, int *lda, double *s, double *scond, double *amax, int *info);
int zporfs_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double * rwork, int *info);
int zporfsx_check(char *uplo, char *equed, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *berr, int * n_err_bnds__, double *err_bnds_norm__, double * err_bnds_comp__, int *nparams, double *params, dcomplex * work, double *rwork, int *info);
int zposv_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *info);
int zposvx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex *work, double *rwork, int * info);
int zposvxx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int *nparams, double *params, dcomplex *work, double *rwork, int *info);
int zpotf2_check(char *uplo, int *n, dcomplex *a, int *lda, int *info);
int zpotrf_check(char *uplo, int *n, dcomplex *a, int *lda, int *info);
int zpotri_check(char *uplo, int *n, dcomplex *a, int *lda, int *info);
int zpotrs_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *info);
int zppcon_check(char *uplo, int *n, dcomplex *ap, double *anorm, double *rcond, dcomplex *work, double *rwork, int *info);
int zppequ_check(char *uplo, int *n, dcomplex *ap, double *s, double *scond, double *amax, int *info);
int zpprfs_check(char *uplo, int *n, int *nrhs, dcomplex *ap, dcomplex *afp, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *info);
int zppsv_check(char *uplo, int *n, int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, int *info);
int zppsvx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *ap, dcomplex *afp, char *equed, double * s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex * work, double *rwork, int *info);
int zpptrf_check(char *uplo, int *n, dcomplex *ap, int *info);
int zpptri_check(char *uplo, int *n, dcomplex *ap, int *info);
int zpptrs_check(char *uplo, int *n, int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, int *info);
int zpstf2_check(char *uplo, int *n, dcomplex *a, int *lda, int *piv, int *rank, double *tol, double *work, int *info);
int zpstrf_check(char *uplo, int *n, dcomplex *a, int *lda, int *piv, int *rank, double *tol, double *work, int *info);
int zptcon_check(int *n, double *d__, dcomplex *e, double *anorm, double *rcond, double *rwork, int * info);
int zpteqr_check(char *compz, int *n, double *d__, double *e, dcomplex *z__, int *ldz, double *work, int *info);
int zptrfs_check(char *uplo, int *n, int *nrhs, double *d__, dcomplex *e, double *df, dcomplex *ef, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double * rwork, int *info);
int zptsv_check(int *n, int *nrhs, double *d__, dcomplex *e, dcomplex *b, int *ldb, int *info);
int zptsvx_check(char *fact, int *n, int *nrhs, double *d__, dcomplex *e, double *df, dcomplex *ef, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex * work, double *rwork, int *info);
int zpttrf_check(int *n, double *d__, dcomplex *e, int *info);
int zpttrs_check(char *uplo, int *n, int *nrhs, double *d__, dcomplex *e, dcomplex *b, int *ldb, int *info);
int zptts2_check(int *iuplo, int *n, int *nrhs, double *d__, dcomplex *e, dcomplex *b, int *ldb);
int zrot_check(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy, double *c__, dcomplex *s);
int zspcon_check(char *uplo, int *n, dcomplex *ap, int *ipiv, double *anorm, double *rcond, dcomplex * work, int *info);
int zspmv_check(char *uplo, int *n, dcomplex *alpha, dcomplex *ap, dcomplex *x, int *incx, dcomplex * beta, dcomplex *y, int *incy);
int zspr_check(char *uplo, int *n, dcomplex *alpha, dcomplex *x, int *incx, dcomplex *ap);
int zsprfs_check(char *uplo, int *n, int *nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex * b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int * info);
int zspsv_check(char *uplo, int *n, int *nrhs, dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, int *info);
int zspsvx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *ap, dcomplex *afp, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex * work, double *rwork, int *info);
int zsptrf_check(char *uplo, int *n, dcomplex *ap, int *ipiv, int *info);
int zsptri_check(char *uplo, int *n, dcomplex *ap, int *ipiv, dcomplex *work, int *info);
int zsptrs_check(char *uplo, int *n, int *nrhs, dcomplex *ap, int *ipiv, dcomplex *b, int *ldb, int *info);
int zstedc_check(char *compz, int *n, double *d__, double *e, dcomplex *z__, int *ldz, dcomplex *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);
int zstegr_check(char *jobz, char *range, int *n, double * d__, double *e, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, dcomplex *z__, int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);
int zstein_check(int *n, double *d__, double *e, int *m, double *w, int *iblock, int *isplit, dcomplex *z__, int *ldz, double *work, int *iwork, int *ifail, int *info);
int zstemr_check(char *jobz, char *range, int *n, double * d__, double *e, double *vl, double *vu, int *il, int *iu, int *m, double *w, dcomplex *z__, int * ldz, int *nzc, int *isuppz, logical *tryrac, double *work, int *lwork, int *iwork, int *liwork, int *info);
int zsteqr_check(char *compz, int *n, double *d__, double *e, dcomplex *z__, int *ldz, double *work, int *info);
int zsycon_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, double *anorm, double *rcond, dcomplex *work, int *info);
int zsycon_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, double *anorm, double *rcond, dcomplex *work, int *info);
int zsyconv_check(char *uplo, char *way, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *info);
int zsyequb_check(char *uplo, int *n, dcomplex *a, int *lda, double *s, double *scond, double *amax, dcomplex *work, int *info);
int zsymv_check(char *uplo, int *n, dcomplex *alpha, dcomplex *a, int *lda, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y, int *incy);
int zsyr_check(char *uplo, int *n, dcomplex *alpha, dcomplex *x, int *incx, dcomplex *a, int *lda);
int zsyrfs_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *af, int *ldaf, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *info);
int zsyrfsx_check(char *uplo, char *equed, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double * err_bnds_comp__, int *nparams, double *params, dcomplex * work, double *rwork, int *info);
int zsysv_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, int *lwork, int *info);
int zsysv_rook_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, int *lwork, int *info);
int zsysvx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *ferr, double *berr, dcomplex *work, int *lwork, double *rwork, int *info);
int zsysvxx_check(char *fact, char *uplo, int *n, int * nrhs, dcomplex *a, int *lda, dcomplex *af, int * ldaf, int *ipiv, char *equed, double *s, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *rcond, double *rpvgrw, double *berr, int *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__, int * nparams, double *params, dcomplex *work, double *rwork, int *info);
int zsyswapr_check(char *uplo, int *n, dcomplex *a, int *lda, int *i1, int *i2);
int zsytf2_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, int *info);
int zsytf2_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, int *info);
int zsytrf_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *lwork, int *info);
int zsytrf_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *lwork, int *info);
int zsytri_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *info);
int zsytri2_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *lwork, int *info);
int zsytri2x_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *nb, int *info);
int zsytri_rook_check(char *uplo, int *n, dcomplex *a, int *lda, int *ipiv, dcomplex *work, int *info);
int zsytrs_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *info);
int zsytrs2_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, dcomplex *work, int *info);
int zsytrs_rook_check(char *uplo, int *n, int *nrhs, dcomplex *a, int *lda, int *ipiv, dcomplex *b, int *ldb, int *info);
int ztbcon_check(char *norm, char *uplo, char *diag, int *n, int *kd, dcomplex *ab, int *ldab, double *rcond, dcomplex *work, double *rwork, int *info);
int ztbrfs_check(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, dcomplex *ab, int *ldab, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double * rwork, int *info);
int ztbtrs_check(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, dcomplex *ab, int *ldab, dcomplex *b, int *ldb, int *info);
int ztfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, int *m, int *n, dcomplex *alpha, dcomplex *a, dcomplex *b, int *ldb);
int ztftri_check(char *transr, char *uplo, char *diag, int *n, dcomplex *a, int *info);
int ztfttp_check(char *transr, char *uplo, int *n, dcomplex *arf, dcomplex *ap, int *info);
int ztfttr_check(char *transr, char *uplo, int *n, dcomplex *arf, dcomplex *a, int *lda, int *info);
int ztgevc_check(char *side, char *howmny, logical *select, int *n, dcomplex *s, int *lds, dcomplex *p, int *ldp, dcomplex *vl, int *ldvl, dcomplex *vr, int * ldvr, int *mm, int *m, dcomplex *work, double *rwork, int *info);
int ztgex2_check(logical *wantq, logical *wantz, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *q, int *ldq, dcomplex *z__, int *ldz, int *j1, int *info);
int ztgexc_check(logical *wantq, logical *wantz, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *q, int *ldq, dcomplex *z__, int *ldz, int *ifst, int *ilst, int *info);
int ztgsen_check(int *ijob, logical *wantq, logical *wantz, logical *select, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *alpha, dcomplex * beta, dcomplex *q, int *ldq, dcomplex *z__, int * ldz, int *m, double *pl, double *pr, double *dif, dcomplex *work, int *lwork, int *iwork, int *liwork, int *info);
int ztgsja_check(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, dcomplex *a, int *lda, dcomplex *b, int *ldb, double *tola, double *tolb, double *alpha, double *beta, dcomplex * u, int *ldu, dcomplex *v, int *ldv, dcomplex *q, int *ldq, dcomplex *work, int *ncycle, int *info);
int ztgsna_check(char *job, char *howmny, logical *select, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *vl, int *ldvl, dcomplex *vr, int * ldvr, double *s, double *dif, int *mm, int *m, dcomplex *work, int *lwork, int *iwork, int *info);
int ztgsy2_check(char *trans, int *ijob, int *m, int * n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *c__, int *ldc, dcomplex *d__, int *ldd, dcomplex *e, int *lde, dcomplex *f, int *ldf, double *scale, double *rdsum, double *rdscal, int * info);
int ztgsyl_check(char *trans, int *ijob, int *m, int * n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *c__, int *ldc, dcomplex *d__, int *ldd, dcomplex *e, int *lde, dcomplex *f, int *ldf, double *scale, double *dif, dcomplex *work, int * lwork, int *iwork, int *info);
int ztpcon_check(char *norm, char *uplo, char *diag, int *n, dcomplex *ap, double *rcond, dcomplex *work, double *rwork, int *info);
int ztpmqrt_check(char *side, char *trans, int *m, int *n, int *k, int *l, int *nb, dcomplex *v, int *ldv, dcomplex *t, int *ldt, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *work, int *info);
int ztpqrt_check(int *m, int *n, int *l, int *nb, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *t, int *ldt, dcomplex *work, int *info);
int ztpqrt2_check(int *m, int *n, int *l, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *t, int *ldt, int *info);
int ztprfb_check(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, int *l, dcomplex *v, int *ldv, dcomplex *t, int *ldt, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *work, int *ldwork);
int ztprfs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int *info);
int ztptri_check(char *uplo, char *diag, int *n, dcomplex *ap, int *info);
int ztptrs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, dcomplex *ap, dcomplex *b, int *ldb, int *info);
int ztpttf_check(char *transr, char *uplo, int *n, dcomplex *ap, dcomplex *arf, int *info);
int ztpttr_check(char *uplo, int *n, dcomplex *ap, dcomplex *a, int *lda, int *info);
int ztrcon_check(char *norm, char *uplo, char *diag, int *n, dcomplex *a, int *lda, double *rcond, dcomplex * work, double *rwork, int *info);
int ztrevc_check(char *side, char *howmny, logical *select, int *n, dcomplex *t, int *ldt, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, int *mm, int *m, dcomplex *work, double *rwork, int *info);
int ztrexc_check(char *compq, int *n, dcomplex *t, int *ldt, dcomplex *q, int *ldq, int *ifst, int * ilst, int *info);
int ztrrfs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *x, int *ldx, double *ferr, double *berr, dcomplex *work, double *rwork, int * info);
int ztrsen_check(char *job, char *compq, logical *select, int *n, dcomplex *t, int *ldt, dcomplex *q, int *ldq, dcomplex *w, int *m, double *s, double *sep, dcomplex *work, int *lwork, int *info);
int ztrsna_check(char *job, char *howmny, logical *select, int *n, dcomplex *t, int *ldt, dcomplex *vl, int *ldvl, dcomplex *vr, int *ldvr, double *s, double *sep, int *mm, int *m, dcomplex *work, int *ldwork, double *rwork, int *info);
int ztrsyl_check(char *trana, char *tranb, int *isgn, int *m, int *n, dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *c__, int *ldc, double *scale, int *info);
int ztrti2_check(char *uplo, char *diag, int *n, dcomplex *a, int *lda, int *info);
int ztrtri_check(char *uplo, char *diag, int *n, dcomplex *a, int *lda, int *info);
int ztrtrs_check(char *uplo, char *trans, char *diag, int *n, int *nrhs, dcomplex *a, int *lda, dcomplex *b, int *ldb, int *info);
int ztrttf_check(char *transr, char *uplo, int *n, dcomplex *a, int *lda, dcomplex *arf, int *info);
int ztrttp_check(char *uplo, int *n, dcomplex *a, int *lda, dcomplex *ap, int *info);
int ztzrqf_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, int *info);
int ztzrzf_check(int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zunbdb_check(char *trans, char *signs, int *m, int *p, int *q, dcomplex *x11, int *ldx11, dcomplex *x12, int *ldx12, dcomplex *x21, int *ldx21, dcomplex * x22, int *ldx22, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *tauq2, dcomplex *work, int *lwork, int * info);
int zunbdb1_check(int *m, int *p, int *q, dcomplex *x11, int *ldx11, dcomplex *x21, int * ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, int *lwork, int *info);
int zunbdb2_check(int *m, int *p, int *q, dcomplex *x11, int *ldx11, dcomplex *x21, int * ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, int *lwork, int *info);
int zunbdb3_check(int *m, int *p, int *q, dcomplex *x11, int *ldx11, dcomplex *x21, int * ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *work, int *lwork, int *info);
int zunbdb4_check(int *m, int *p, int *q, dcomplex *x11, int *ldx11, dcomplex *x21, int * ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2, dcomplex *tauq1, dcomplex *phantom, dcomplex *work, int *lwork, int *info);
int zunbdb5_check(int *m1, int *m2, int *n, dcomplex *x1, int *incx1, dcomplex *x2, int *incx2, dcomplex *q1, int *ldq1, dcomplex *q2, int *ldq2, dcomplex *work, int *lwork, int *info);
int zunbdb6_check(int *m1, int *m2, int *n, dcomplex *x1, int *incx1, dcomplex *x2, int *incx2, dcomplex *q1, int *ldq1, dcomplex *q2, int *ldq2, dcomplex *work, int *lwork, int *info);
int zuncsd_check(char *jobu1, char *jobu2, char *jobv1t, char * jobv2t, char *trans, char *signs, int *m, int *p, int *q, dcomplex *x11, int *ldx11, dcomplex *x12, int * ldx12, dcomplex *x21, int *ldx21, dcomplex *x22, int *ldx22, double *theta, dcomplex *u1, int *ldu1, dcomplex *u2, int *ldu2, dcomplex *v1t, int *ldv1t, dcomplex *v2t, int *ldv2t, dcomplex *work, int * lwork, double *rwork, int *lrwork, int *iwork, int * info);
int zuncsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, int *m, int *p, int *q, dcomplex *x11, int * ldx11, dcomplex *x21, int *ldx21, double *theta, dcomplex *u1, int *ldu1, dcomplex *u2, int *ldu2, dcomplex *v1t, int *ldv1t, dcomplex *work, int * lwork, double *rwork, int *lrwork, int *iwork, int * info);
int zung2l_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *info);
int zung2r_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *info);
int zungbr_check(char *vect, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info);
int zunghr_check(int *n, int *ilo, int *ihi, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info);
int zungl2_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *info);
int zunglq_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info);
int zungql_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info);
int zungqr_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info);
int zungr2_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *info);
int zungrq_check(int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex * work, int *lwork, int *info);
int zungtr_check(char *uplo, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *work, int *lwork, int *info);
int zunm2l_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *info);
int zunm2r_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *info);
int zunmbr_check(char *vect, char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int * lwork, int *info);
int zunmhr_check(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex * work, int *lwork, int *info);
int zunml2_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *info);
int zunmlq_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *info);
int zunmql_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *info);
int zunmqr_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *info);
int zunmr2_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *info);
int zunmr3_check(char *side, char *trans, int *m, int *n, int *k, int *l, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int * info);
int zunmrq_check(char *side, char *trans, int *m, int *n, int *k, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *info);
int zunmrz_check(char *side, char *trans, int *m, int *n, int *k, int *l, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int * lwork, int *info);
int zunmtr_check(char *side, char *uplo, char *trans, int *m, int *n, dcomplex *a, int *lda, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *lwork, int *info);
int zupgtr_check(char *uplo, int *n, dcomplex *ap, dcomplex *tau, dcomplex *q, int *ldq, dcomplex * work, int *info);
int zupmtr_check(char *side, char *uplo, char *trans, int *m, int *n, dcomplex *ap, dcomplex *tau, dcomplex *c__, int *ldc, dcomplex *work, int *info);
