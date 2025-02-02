/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
int cbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m,
                 integer *p, integer *q, float *theta, float *phi, scomplex *u1, integer *ldu1,
                 scomplex *u2, integer *ldu2, scomplex *v1t, integer *ldv1t, scomplex *v2t,
                 integer *ldv2t, float *b11d, float *b11e, float *b12d, float *b12e, float *b21d,
                 float *b21e, float *b22d, float *b22e, float *rwork, integer *lrwork,
                 integer *info);
int cbdsqr_check(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, float *d__,
                 float *e, scomplex *vt, integer *ldvt, scomplex *u, integer *ldu, scomplex *c__,
                 integer *ldc, float *rwork, integer *info);
int cgbbrd_check(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku,
                 scomplex *ab, integer *ldab, float *d__, float *e, scomplex *q, integer *ldq,
                 scomplex *pt, integer *ldpt, scomplex *c__, integer *ldc, scomplex *work,
                 float *rwork, integer *info);
int cgbcon_check(char *norm, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab,
                 integer *ipiv, float *anorm, float *rcond, scomplex *work, float *rwork,
                 integer *info);
int cgbequ_check(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, integer *info);
int cgbequb_check(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab,
                  float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, integer *info);
int cgbrfs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, scomplex *ab,
                 integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv, scomplex *b,
                 integer *ldb, scomplex *x, integer *ldx, float *ferr, float *berr, scomplex *work,
                 float *rwork, integer *info);
int cgbrfsx_check(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs,
                  scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv,
                  float *r__, float *c__, scomplex *b, integer *ldb, scomplex *x, integer *ldx,
                  float *rcond, float *berr, integer *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, integer *nparams, float *params, scomplex *work,
                  float *rwork, integer *info);
int cgbsv_check(integer *n, integer *kl, integer *ku, integer *nrhs, scomplex *ab, integer *ldab,
                integer *ipiv, scomplex *b, integer *ldb, integer *info);
int cgbsvx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                 scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv,
                 char *equed, float *r__, float *c__, scomplex *b, integer *ldb, scomplex *x,
                 integer *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork,
                 integer *info);
int cgbsvxx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                  scomplex *ab, integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv,
                  char *equed, float *r__, float *c__, scomplex *b, integer *ldb, scomplex *x,
                  integer *ldx, float *rcond, float *rpvgrw, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  scomplex *work, float *rwork, integer *info);
int cgbtf2_check(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab,
                 integer *ipiv, integer *info);
int cgbtrf_check(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab,
                 integer *ipiv, integer *info);
int cgbtrs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, scomplex *ab,
                 integer *ldab, integer *ipiv, scomplex *b, integer *ldb, integer *info);
int cgebak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, float *scale,
                 integer *m, scomplex *v, integer *ldv, integer *info);
int cgebal_check(char *job, integer *n, scomplex *a, integer *lda, integer *ilo, integer *ihi,
                 float *scale, integer *info);
int cgebd2_check(integer *m, integer *n, scomplex *a, integer *lda, float *d__, float *e,
                 scomplex *tauq, scomplex *taup, scomplex *work, integer *info);
int cgebrd_check(integer *m, integer *n, scomplex *a, integer *lda, float *d__, float *e,
                 scomplex *tauq, scomplex *taup, scomplex *work, integer *lwork, integer *info);
int cgecon_check(char *norm, integer *n, scomplex *a, integer *lda, float *anorm, float *rcond,
                 scomplex *work, float *rwork, integer *info);
int cgeequ_check(integer *m, integer *n, scomplex *a, integer *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, integer *info);
int cgeequb_check(integer *m, integer *n, scomplex *a, integer *lda, float *r__, float *c__,
                  float *rowcnd, float *colcnd, float *amax, integer *info);
int cgees_check(char *jobvs, char *sort, L_fp select, integer *n, scomplex *a, integer *lda,
                integer *sdim, scomplex *w, scomplex *vs, integer *ldvs, scomplex *work,
                integer *lwork, float *rwork, logical *bwork, integer *info);
int cgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, integer *n, scomplex *a,
                 integer *lda, integer *sdim, scomplex *w, scomplex *vs, integer *ldvs,
                 float *rconde, float *rcondv, scomplex *work, integer *lwork, float *rwork,
                 logical *bwork, integer *info);
int cgeev_check(char *jobvl, char *jobvr, integer *n, scomplex *a, integer *lda, scomplex *w,
                scomplex *vl, integer *ldvl, scomplex *vr, integer *ldvr, scomplex *work,
                integer *lwork, float *rwork, integer *info);
int cgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, scomplex *a,
                 integer *lda, scomplex *w, scomplex *vl, integer *ldvl, scomplex *vr,
                 integer *ldvr, integer *ilo, integer *ihi, float *scale, float *abnrm,
                 float *rconde, float *rcondv, scomplex *work, integer *lwork, float *rwork,
                 integer *info);
int cgegs_check(char *jobvsl, char *jobvsr, integer *n, scomplex *a, integer *lda, scomplex *b,
                integer *ldb, scomplex *alpha, scomplex *beta, scomplex *vsl, integer *ldvsl,
                scomplex *vsr, integer *ldvsr, scomplex *work, integer *lwork, float *rwork,
                integer *info);
int cgegv_check(char *jobvl, char *jobvr, integer *n, scomplex *a, integer *lda, scomplex *b,
                integer *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, integer *ldvl,
                scomplex *vr, integer *ldvr, scomplex *work, integer *lwork, float *rwork,
                integer *info);
int cgehd2_check(integer *n, integer *ilo, integer *ihi, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *info);
int cgehrd_check(integer *n, integer *ilo, integer *ihi, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *lwork, integer *info);
int cgelq2_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *info);
int cgelqf_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *lwork, integer *info);
int cgels_check(char *trans, integer *m, integer *n, integer *nrhs, scomplex *a, integer *lda,
                scomplex *b, integer *ldb, scomplex *work, integer *lwork, integer *info);
int cgelsd_check(integer *m, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, float *s, float *rcond, integer *rank, scomplex *work,
                 integer *lwork, float *rwork, integer *iwork, integer *info);
int cgelss_check(integer *m, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, float *s, float *rcond, integer *rank, scomplex *work,
                 integer *lwork, float *rwork, integer *info);
int cgelsx_check(integer *m, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, integer *jpvt, float *rcond, integer *rank, scomplex *work,
                 float *rwork, integer *info);
int cgelsy_check(integer *m, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, integer *jpvt, float *rcond, integer *rank, scomplex *work,
                 integer *lwork, float *rwork, integer *info);
int cgemqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *nb,
                  scomplex *v, integer *ldv, scomplex *t, integer *ldt, scomplex *c__, integer *ldc,
                  scomplex *work, integer *info);
int cgeql2_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *info);
int cgeqlf_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *lwork, integer *info);
int cgeqp3_check(integer *m, integer *n, scomplex *a, integer *lda, integer *jpvt, scomplex *tau,
                 scomplex *work, integer *lwork, float *rwork, integer *info);
int cgeqpf_check(integer *m, integer *n, scomplex *a, integer *lda, integer *jpvt, scomplex *tau,
                 scomplex *work, float *rwork, integer *info);
int cgeqr2_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *info);
int cgeqr2p_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                  integer *info);
int cgeqrf_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *lwork, integer *info);
int cgeqrfp_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                  integer *lwork, integer *info);
int cgeqrt_check(integer *m, integer *n, integer *nb, scomplex *a, integer *lda, scomplex *t,
                 integer *ldt, scomplex *work, integer *info);
int cgeqrt2_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *t, integer *ldt,
                  integer *info);
int cgeqrt3_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *t, integer *ldt,
                  integer *info);
int cgerfs_check(char *trans, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af,
                 integer *ldaf, integer *ipiv, scomplex *b, integer *ldb, scomplex *x, integer *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int cgerfsx_check(char *trans, char *equed, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, integer *ipiv, float *r__, float *c__, scomplex *b,
                  integer *ldb, scomplex *x, integer *ldx, float *rcond, float *berr,
                  integer *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  integer *nparams, float *params, scomplex *work, float *rwork, integer *info);
int cgerq2_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *info);
int cgerqf_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *lwork, integer *info);
int cgesc2_check(integer *n, scomplex *a, integer *lda, scomplex *rhs, integer *ipiv, integer *jpiv,
                 float *scale);
int cgesdd_check(char *jobz, integer *m, integer *n, scomplex *a, integer *lda, float *s,
                 scomplex *u, integer *ldu, scomplex *vt, integer *ldvt, scomplex *work,
                 integer *lwork, float *rwork, integer *iwork, integer *info);
int cgesv_check(integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv, scomplex *b,
                integer *ldb, integer *info);
int cgesvd_check(char *jobu, char *jobvt, integer *m, integer *n, scomplex *a, integer *lda,
                 float *s, scomplex *u, integer *ldu, scomplex *vt, integer *ldvt, scomplex *work,
                 integer *lwork, float *rwork, integer *info);
int cgesvx_check(char *fact, char *trans, integer *n, integer *nrhs, scomplex *a, integer *lda,
                 scomplex *af, integer *ldaf, integer *ipiv, char *equed, float *r__, float *c__,
                 scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *rcond, float *ferr,
                 float *berr, scomplex *work, float *rwork, integer *info);
int cgesvxx_check(char *fact, char *trans, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, integer *ipiv, char *equed, float *r__, float *c__,
                  scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *rcond, float *rpvgrw,
                  float *berr, integer *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, integer *nparams, float *params, scomplex *work,
                  float *rwork, integer *info);
int cgetc2_check(integer *n, scomplex *a, integer *lda, integer *ipiv, integer *jpiv,
                 integer *info);
int cgetf2_check(integer *m, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info);
int cgetrf_check(integer *m, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info);
int cgetrfnp_check(integer *m, integer *n, scomplex *a, integer *lda, integer *info);
int cgetrfnpi_check(integer *m, integer *n, integer *nfact, scomplex *a, integer *lda,
                    integer *info);
int cgetri_check(integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                 integer *lwork, integer *info);
int cgetrs_check(char *trans, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv,
                 scomplex *b, integer *ldb, integer *info);
int cggbak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, float *lscale,
                 float *rscale, integer *m, scomplex *v, integer *ldv, integer *info);
int cggbal_check(char *job, integer *n, scomplex *a, integer *lda, scomplex *b, integer *ldb,
                 integer *ilo, integer *ihi, float *lscale, float *rscale, float *work,
                 integer *info);
int cgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, scomplex *a,
                integer *lda, scomplex *b, integer *ldb, integer *sdim, scomplex *alpha,
                scomplex *beta, scomplex *vsl, integer *ldvsl, scomplex *vsr, integer *ldvsr,
                scomplex *work, integer *lwork, float *rwork, logical *bwork, integer *info);
int cggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n,
                 scomplex *a, integer *lda, scomplex *b, integer *ldb, integer *sdim,
                 scomplex *alpha, scomplex *beta, scomplex *vsl, integer *ldvsl, scomplex *vsr,
                 integer *ldvsr, float *rconde, float *rcondv, scomplex *work, integer *lwork,
                 float *rwork, integer *iwork, integer *liwork, logical *bwork, integer *info);
int cggev_check(char *jobvl, char *jobvr, integer *n, scomplex *a, integer *lda, scomplex *b,
                integer *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, integer *ldvl,
                scomplex *vr, integer *ldvr, scomplex *work, integer *lwork, float *rwork,
                integer *info);
int cggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, scomplex *a,
                 integer *lda, scomplex *b, integer *ldb, scomplex *alpha, scomplex *beta,
                 scomplex *vl, integer *ldvl, scomplex *vr, integer *ldvr, integer *ilo,
                 integer *ihi, float *lscale, float *rscale, float *abnrm, float *bbnrm,
                 float *rconde, float *rcondv, scomplex *work, integer *lwork, float *rwork,
                 integer *iwork, logical *bwork, integer *info);
int cggglm_check(integer *n, integer *m, integer *p, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, scomplex *d__, scomplex *x, scomplex *y, scomplex *work,
                 integer *lwork, integer *info);
int cgghrd_check(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, scomplex *a,
                 integer *lda, scomplex *b, integer *ldb, scomplex *q, integer *ldq, scomplex *z__,
                 integer *ldz, integer *info);
int cgglse_check(integer *m, integer *n, integer *p, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, scomplex *c__, scomplex *d__, scomplex *x, scomplex *work,
                 integer *lwork, integer *info);
int cggqrf_check(integer *n, integer *m, integer *p, scomplex *a, integer *lda, scomplex *taua,
                 scomplex *b, integer *ldb, scomplex *taub, scomplex *work, integer *lwork,
                 integer *info);
int cggrqf_check(integer *m, integer *p, integer *n, scomplex *a, integer *lda, scomplex *taua,
                 scomplex *b, integer *ldb, scomplex *taub, scomplex *work, integer *lwork,
                 integer *info);
int cggsvd_check(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k,
                 integer *l, scomplex *a, integer *lda, scomplex *b, integer *ldb, float *alpha,
                 float *beta, scomplex *u, integer *ldu, scomplex *v, integer *ldv, scomplex *q,
                 integer *ldq, scomplex *work, float *rwork, integer *iwork, integer *info);
int cggsvp_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n,
                 scomplex *a, integer *lda, scomplex *b, integer *ldb, float *tola, float *tolb,
                 integer *k, integer *l, scomplex *u, integer *ldu, scomplex *v, integer *ldv,
                 scomplex *q, integer *ldq, integer *iwork, float *rwork, scomplex *tau,
                 scomplex *work, integer *info);
int cgtcon_check(char *norm, integer *n, scomplex *dl, scomplex *d__, scomplex *du, scomplex *du2,
                 integer *ipiv, float *anorm, float *rcond, scomplex *work, integer *info);
int cgtrfs_check(char *trans, integer *n, integer *nrhs, scomplex *dl, scomplex *d__, scomplex *du,
                 scomplex *dlf, scomplex *df, scomplex *duf, scomplex *du2, integer *ipiv,
                 scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, integer *info);
int cgtsv_check(integer *n, integer *nrhs, scomplex *dl, scomplex *d__, scomplex *du, scomplex *b,
                integer *ldb, integer *info);
int cgtsvx_check(char *fact, char *trans, integer *n, integer *nrhs, scomplex *dl, scomplex *d__,
                 scomplex *du, scomplex *dlf, scomplex *df, scomplex *duf, scomplex *du2,
                 integer *ipiv, scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int cgttrf_check(integer *n, scomplex *dl, scomplex *d__, scomplex *du, scomplex *du2,
                 integer *ipiv, integer *info);
int cgttrs_check(char *trans, integer *n, integer *nrhs, scomplex *dl, scomplex *d__, scomplex *du,
                 scomplex *du2, integer *ipiv, scomplex *b, integer *ldb, integer *info);
int cgtts2_check(integer *itrans, integer *n, integer *nrhs, scomplex *dl, scomplex *d__,
                 scomplex *du, scomplex *du2, integer *ipiv, scomplex *b, integer *ldb);
int chbev_check(char *jobz, char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab,
                float *w, scomplex *z__, integer *ldz, scomplex *work, float *rwork, integer *info);
int chbevd_check(char *jobz, char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab,
                 float *w, scomplex *z__, integer *ldz, scomplex *work, integer *lwork,
                 float *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);
int chbevx_check(char *jobz, char *range, char *uplo, integer *n, integer *kd, scomplex *ab,
                 integer *ldab, scomplex *q, integer *ldq, float *vl, float *vu, integer *il,
                 integer *iu, float *abstol, integer *m, float *w, scomplex *z__, integer *ldz,
                 scomplex *work, float *rwork, integer *iwork, integer *ifail, integer *info);
int chbgst_check(char *vect, char *uplo, integer *n, integer *ka, integer *kb, scomplex *ab,
                 integer *ldab, scomplex *bb, integer *ldbb, scomplex *x, integer *ldx,
                 scomplex *work, float *rwork, integer *info);
int chbgv_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, scomplex *ab,
                integer *ldab, scomplex *bb, integer *ldbb, float *w, scomplex *z__, integer *ldz,
                scomplex *work, float *rwork, integer *info);
int chbgvd_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, scomplex *ab,
                 integer *ldab, scomplex *bb, integer *ldbb, float *w, scomplex *z__, integer *ldz,
                 scomplex *work, integer *lwork, float *rwork, integer *lrwork, integer *iwork,
                 integer *liwork, integer *info);
int chbgvx_check(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb,
                 scomplex *ab, integer *ldab, scomplex *bb, integer *ldbb, scomplex *q,
                 integer *ldq, float *vl, float *vu, integer *il, integer *iu, float *abstol,
                 integer *m, float *w, scomplex *z__, integer *ldz, scomplex *work, float *rwork,
                 integer *iwork, integer *ifail, integer *info);
int chbtrd_check(char *vect, char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab,
                 float *d__, float *e, scomplex *q, integer *ldq, scomplex *work, integer *info);
int checon_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, float *anorm,
                 float *rcond, scomplex *work, integer *info);
int checon_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      float *anorm, float *rcond, scomplex *work, integer *info);
int cheequb_check(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond,
                  float *amax, scomplex *work, integer *info);
int cheev_check(char *jobz, char *uplo, integer *n, scomplex *a, integer *lda, float *w,
                scomplex *work, integer *lwork, float *rwork, integer *info);
int cheevd_check(char *jobz, char *uplo, integer *n, scomplex *a, integer *lda, float *w,
                 scomplex *work, integer *lwork, float *rwork, integer *lrwork, integer *iwork,
                 integer *liwork, integer *info);
int cheevr_check(char *jobz, char *range, char *uplo, integer *n, scomplex *a, integer *lda,
                 float *vl, float *vu, integer *il, integer *iu, float *abstol, integer *m,
                 float *w, scomplex *z__, integer *ldz, integer *isuppz, scomplex *work,
                 integer *lwork, float *rwork, integer *lrwork, integer *iwork, integer *liwork,
                 integer *info);
int cheevx_check(char *jobz, char *range, char *uplo, integer *n, scomplex *a, integer *lda,
                 float *vl, float *vu, integer *il, integer *iu, float *abstol, integer *m,
                 float *w, scomplex *z__, integer *ldz, scomplex *work, integer *lwork,
                 float *rwork, integer *iwork, integer *ifail, integer *info);
int chegs2_check(integer *itype, char *uplo, integer *n, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, integer *info);
int chegst_check(integer *itype, char *uplo, integer *n, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, integer *info);
int chegv_check(integer *itype, char *jobz, char *uplo, integer *n, scomplex *a, integer *lda,
                scomplex *b, integer *ldb, float *w, scomplex *work, integer *lwork, float *rwork,
                integer *info);
int chegvd_check(integer *itype, char *jobz, char *uplo, integer *n, scomplex *a, integer *lda,
                 scomplex *b, integer *ldb, float *w, scomplex *work, integer *lwork, float *rwork,
                 integer *lrwork, integer *iwork, integer *liwork, integer *info);
int chegvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, scomplex *a,
                 integer *lda, scomplex *b, integer *ldb, float *vl, float *vu, integer *il,
                 integer *iu, float *abstol, integer *m, float *w, scomplex *z__, integer *ldz,
                 scomplex *work, integer *lwork, float *rwork, integer *iwork, integer *ifail,
                 integer *info);
int cherfs_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af,
                 integer *ldaf, integer *ipiv, scomplex *b, integer *ldb, scomplex *x, integer *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int cherfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, integer *ipiv, float *s, scomplex *b, integer *ldb,
                  scomplex *x, integer *ldx, float *rcond, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  scomplex *work, float *rwork, integer *info);
int chesv_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv,
                scomplex *b, integer *ldb, scomplex *work, integer *lwork, integer *info);
int chesv_rook_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                     integer *ipiv, scomplex *b, integer *ldb, scomplex *work, integer *lwork,
                     integer *info);
int chesvx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                 scomplex *af, integer *ldaf, integer *ipiv, scomplex *b, integer *ldb, scomplex *x,
                 integer *ldx, float *rcond, float *ferr, float *berr, scomplex *work,
                 integer *lwork, float *rwork, integer *info);
int chesvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, integer *ipiv, char *equed, float *s, scomplex *b,
                  integer *ldb, scomplex *x, integer *ldx, float *rcond, float *rpvgrw, float *berr,
                  integer *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  integer *nparams, float *params, scomplex *work, float *rwork, integer *info);
int cheswapr_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *i1, integer *i2);
int chetd2_check(char *uplo, integer *n, scomplex *a, integer *lda, float *d__, float *e,
                 scomplex *tau, integer *info);
int chetf2_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info);
int chetf2_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      integer *info);
int chetrd_check(char *uplo, integer *n, scomplex *a, integer *lda, float *d__, float *e,
                 scomplex *tau, scomplex *work, integer *lwork, integer *info);
int chetrf_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                 integer *lwork, integer *info);
int chetrf_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      scomplex *work, integer *lwork, integer *info);
int chetri_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                 integer *info);
int chetri2_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                  integer *lwork, integer *info);
int chetri2x_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                   integer *nb, integer *info);
int chetri_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      scomplex *work, integer *info);
int chetrs_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv,
                 scomplex *b, integer *ldb, integer *info);
int chetrs2_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv,
                  scomplex *b, integer *ldb, scomplex *work, integer *info);
int chetrs_rook_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                      integer *ipiv, scomplex *b, integer *ldb, integer *info);
int chfrk_check(char *transr, char *uplo, char *trans, integer *n, integer *k, float *alpha,
                scomplex *a, integer *lda, float *beta, scomplex *c__);
int chgeqz_check(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi,
                 scomplex *h__, integer *ldh, scomplex *t, integer *ldt, scomplex *alpha,
                 scomplex *beta, scomplex *q, integer *ldq, scomplex *z__, integer *ldz,
                 scomplex *work, integer *lwork, float *rwork, integer *info);
VOID chla_transtype_check(char *ret_val, integer *trans);
int chpcon_check(char *uplo, integer *n, scomplex *ap, integer *ipiv, float *anorm, float *rcond,
                 scomplex *work, integer *info);
int chpev_check(char *jobz, char *uplo, integer *n, scomplex *ap, float *w, scomplex *z__,
                integer *ldz, scomplex *work, float *rwork, integer *info);
int chpevd_check(char *jobz, char *uplo, integer *n, scomplex *ap, float *w, scomplex *z__,
                 integer *ldz, scomplex *work, integer *lwork, float *rwork, integer *lrwork,
                 integer *iwork, integer *liwork, integer *info);
int chpevx_check(char *jobz, char *range, char *uplo, integer *n, scomplex *ap, float *vl,
                 float *vu, integer *il, integer *iu, float *abstol, integer *m, float *w,
                 scomplex *z__, integer *ldz, scomplex *work, float *rwork, integer *iwork,
                 integer *ifail, integer *info);
int chpgst_check(integer *itype, char *uplo, integer *n, scomplex *ap, scomplex *bp, integer *info);
int chpgv_check(integer *itype, char *jobz, char *uplo, integer *n, scomplex *ap, scomplex *bp,
                float *w, scomplex *z__, integer *ldz, scomplex *work, float *rwork, integer *info);
int chpgvd_check(integer *itype, char *jobz, char *uplo, integer *n, scomplex *ap, scomplex *bp,
                 float *w, scomplex *z__, integer *ldz, scomplex *work, integer *lwork,
                 float *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);
int chpgvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, scomplex *ap,
                 scomplex *bp, float *vl, float *vu, integer *il, integer *iu, float *abstol,
                 integer *m, float *w, scomplex *z__, integer *ldz, scomplex *work, float *rwork,
                 integer *iwork, integer *ifail, integer *info);
int chprfs_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *afp, integer *ipiv,
                 scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, integer *info);
int chpsv_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, integer *ipiv, scomplex *b,
                integer *ldb, integer *info);
int chpsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *afp,
                 integer *ipiv, scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int chptrd_check(char *uplo, integer *n, scomplex *ap, float *d__, float *e, scomplex *tau,
                 integer *info);
int chptrf_check(char *uplo, integer *n, scomplex *ap, integer *ipiv, integer *info);
int chptri_check(char *uplo, integer *n, scomplex *ap, integer *ipiv, scomplex *work,
                 integer *info);
int chptrs_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, integer *ipiv, scomplex *b,
                 integer *ldb, integer *info);
int chsein_check(char *side, char *eigsrc, char *initv, logical *select, integer *n, scomplex *h__,
                 integer *ldh, scomplex *w, scomplex *vl, integer *ldvl, scomplex *vr,
                 integer *ldvr, integer *mm, integer *m, scomplex *work, float *rwork,
                 integer *ifaill, integer *ifailr, integer *info);
int chseqr_check(char *job, char *compz, integer *n, integer *ilo, integer *ihi, scomplex *h__,
                 integer *ldh, scomplex *w, scomplex *z__, integer *ldz, scomplex *work,
                 integer *lwork, integer *info);
int cla_gbamv_check(integer *trans, integer *m, integer *n, integer *kl, integer *ku, float *alpha,
                    scomplex *ab, integer *ldab, scomplex *x, integer *incx, float *beta, float *y,
                    integer *incy);
float cla_gbrcond_c_check(char *trans, integer *n, integer *kl, integer *ku, scomplex *ab,
                          integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv, float *c__,
                          logical *capply, integer *info, scomplex *work, float *rwork);
float cla_gbrcond_x_check(char *trans, integer *n, integer *kl, integer *ku, scomplex *ab,
                          integer *ldab, scomplex *afb, integer *ldafb, integer *ipiv, scomplex *x,
                          integer *info, scomplex *work, float *rwork);
int cla_gbrfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n, integer *kl,
                              integer *ku, integer *nrhs, scomplex *ab, integer *ldab,
                              scomplex *afb, integer *ldafb, integer *ipiv, logical *colequ,
                              float *c__, scomplex *b, integer *ldb, scomplex *y, integer *ldy,
                              float *berr_out__, integer *n_norms__, float *err_bnds_norm__,
                              float *err_bnds_comp__, scomplex *res, float *ayb, scomplex *dy,
                              scomplex *y_tail__, float *rcond, integer *ithresh, float *rthresh,
                              float *dz_ub__, logical *ignore_cwise__, integer *info);
float cla_gbrpvgrw_check(integer *n, integer *kl, integer *ku, integer *ncols, scomplex *ab,
                         integer *ldab, scomplex *afb, integer *ldafb);
int cla_geamv_check(integer *trans, integer *m, integer *n, float *alpha, scomplex *a, integer *lda,
                    scomplex *x, integer *incx, float *beta, float *y, integer *incy);
float cla_gercond_c_check(char *trans, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, integer *ipiv, float *c__, logical *capply, integer *info,
                          scomplex *work, float *rwork);
float cla_gercond_x_check(char *trans, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work,
                          float *rwork);
int cla_gerfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n,
                              integer *nrhs, scomplex *a, integer *lda, scomplex *af, integer *ldaf,
                              integer *ipiv, logical *colequ, float *c__, scomplex *b, integer *ldb,
                              scomplex *y, integer *ldy, float *berr_out__, integer *n_norms__,
                              float *errs_n__, float *errs_c__, scomplex *res, float *ayb,
                              scomplex *dy, scomplex *y_tail__, float *rcond, integer *ithresh,
                              float *rthresh, float *dz_ub__, logical *ignore_cwise__,
                              integer *info);
float cla_gerpvgrw_check(integer *n, integer *ncols, scomplex *a, integer *lda, scomplex *af,
                         integer *ldaf);
int cla_heamv_check(integer *uplo, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x,
                    integer *incx, float *beta, float *y, integer *incy);
float cla_hercond_c_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, integer *ipiv, float *c__, logical *capply, integer *info,
                          scomplex *work, float *rwork);
float cla_hercond_x_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work,
                          float *rwork);
int cla_herfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv,
                              logical *colequ, float *c__, scomplex *b, integer *ldb, scomplex *y,
                              integer *ldy, float *berr_out__, integer *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res,
                              float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond,
                              integer *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, integer *info);
float cla_herpvgrw_check(char *uplo, integer *n, integer *info, scomplex *a, integer *lda,
                         scomplex *af, integer *ldaf, integer *ipiv, float *work);
int cla_lin_berr_check(integer *n, integer *nz, integer *nrhs, scomplex *res, float *ayb,
                       float *berr);
float cla_porcond_c_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, float *c__, logical *capply, integer *info, scomplex *work,
                          float *rwork);
float cla_porcond_x_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, scomplex *x, integer *info, scomplex *work, float *rwork);
int cla_porfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              scomplex *a, integer *lda, scomplex *af, integer *ldaf,
                              logical *colequ, float *c__, scomplex *b, integer *ldb, scomplex *y,
                              integer *ldy, float *berr_out__, integer *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res,
                              float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond,
                              integer *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, integer *info);
float cla_porpvgrw_check(char *uplo, integer *ncols, scomplex *a, integer *lda, scomplex *af,
                         integer *ldaf, float *work);
int cla_syamv_check(integer *uplo, integer *n, float *alpha, scomplex *a, integer *lda, scomplex *x,
                    integer *incx, float *beta, float *y, integer *incy);
float cla_syrcond_c_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, integer *ipiv, float *c__, logical *capply, integer *info,
                          scomplex *work, float *rwork);
float cla_syrcond_x_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *af,
                          integer *ldaf, integer *ipiv, scomplex *x, integer *info, scomplex *work,
                          float *rwork);
int cla_syrfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              scomplex *a, integer *lda, scomplex *af, integer *ldaf, integer *ipiv,
                              logical *colequ, float *c__, scomplex *b, integer *ldb, scomplex *y,
                              integer *ldy, float *berr_out__, integer *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res,
                              float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond,
                              integer *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, integer *info);
float cla_syrpvgrw_check(char *uplo, integer *n, integer *info, scomplex *a, integer *lda,
                         scomplex *af, integer *ldaf, integer *ipiv, float *work);
int cla_wwaddw_check(integer *n, scomplex *x, scomplex *y, scomplex *w);
int clabrd_check(integer *m, integer *n, integer *nb, scomplex *a, integer *lda, float *d__,
                 float *e, scomplex *tauq, scomplex *taup, scomplex *x, integer *ldx, scomplex *y,
                 integer *ldy);
int clacgv_check(integer *n, scomplex *x, integer *incx);
int clacn2_check(integer *n, scomplex *v, scomplex *x, float *est, integer *kase, integer *isave);
int clacon_check(integer *n, scomplex *v, scomplex *x, float *est, integer *kase);
int clacp2_check(char *uplo, integer *m, integer *n, float *a, integer *lda, scomplex *b,
                 integer *ldb);
int clacpy_check(char *uplo, integer *m, integer *n, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb);
int clacrm_check(integer *m, integer *n, scomplex *a, integer *lda, float *b, integer *ldb,
                 scomplex *c__, integer *ldc, float *rwork);
int clacrt_check(integer *n, scomplex *cx, integer *incx, scomplex *cy, integer *incy,
                 scomplex *c__, scomplex *s);
VOID cladiv_check(scomplex *ret_val, scomplex *x, scomplex *y);
int claed0_check(integer *qsiz, integer *n, float *d__, float *e, scomplex *q, integer *ldq,
                 scomplex *qstore, integer *ldqs, float *rwork, integer *iwork, integer *info);
int claed7_check(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl,
                 integer *curpbm, float *d__, scomplex *q, integer *ldq, float *rho, integer *indxq,
                 float *qstore, integer *qptr, integer *prmptr, integer *perm, integer *givptr,
                 integer *givcol, float *givnum, scomplex *work, float *rwork, integer *iwork,
                 integer *info);
int claed8_check(integer *k, integer *n, integer *qsiz, scomplex *q, integer *ldq, float *d__,
                 float *rho, integer *cutpnt, float *z__, float *dlamda, scomplex *q2,
                 integer *ldq2, float *w, integer *indxp, integer *indx, integer *indxq,
                 integer *perm, integer *givptr, integer *givcol, float *givnum, integer *info);
int claein_check(logical *rightv, logical *noinit, integer *n, scomplex *h__, integer *ldh,
                 scomplex *w, scomplex *v, scomplex *b, integer *ldb, float *rwork, float *eps3,
                 float *smlnum, integer *info);
int claesy_check(scomplex *a, scomplex *b, scomplex *c__, scomplex *rt1, scomplex *rt2,
                 scomplex *evscal, scomplex *cs1, scomplex *sn1);
int claev2_check(scomplex *a, scomplex *b, scomplex *c__, float *rt1, float *rt2, float *cs1,
                 scomplex *sn1);
int clag2z_check(integer *m, integer *n, scomplex *sa, integer *ldsa, dcomplex *a, integer *lda,
                 integer *info);
int clags2_check(logical *upper, float *a1, scomplex *a2, float *a3, float *b1, scomplex *b2,
                 float *b3, float *csu, scomplex *snu, float *csv, scomplex *snv, float *csq,
                 scomplex *snq);
int clagtm_check(char *trans, integer *n, integer *nrhs, float *alpha, scomplex *dl, scomplex *d__,
                 scomplex *du, scomplex *x, integer *ldx, float *beta, scomplex *b, integer *ldb);
int clahef_check(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda,
                 integer *ipiv, scomplex *w, integer *ldw, integer *info);
int clahef_rook_check(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda,
                      integer *ipiv, scomplex *w, integer *ldw, integer *info);
int clahqr_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 scomplex *h__, integer *ldh, scomplex *w, integer *iloz, integer *ihiz,
                 scomplex *z__, integer *ldz, integer *info);
int clahr2_check(integer *n, integer *k, integer *nb, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *t, integer *ldt, scomplex *y, integer *ldy);
int clahrd_check(integer *n, integer *k, integer *nb, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *t, integer *ldt, scomplex *y, integer *ldy);
int claic1_check(integer *job, integer *j, scomplex *x, float *sest, scomplex *w, scomplex *gamma,
                 float *sestpr, scomplex *s, scomplex *c__);
int clals0_check(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs,
                 scomplex *b, integer *ldb, scomplex *bx, integer *ldbx, integer *perm,
                 integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum,
                 float *poles, float *difl, float *difr, float *z__, integer *k, float *c__,
                 float *s, float *rwork, integer *info);
int clalsa_check(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, scomplex *b,
                 integer *ldb, scomplex *bx, integer *ldbx, float *u, integer *ldu, float *vt,
                 integer *k, float *difl, float *difr, float *z__, float *poles, integer *givptr,
                 integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c__,
                 float *s, float *rwork, integer *iwork, integer *info);
int clalsd_check(char *uplo, integer *smlsiz, integer *n, integer *nrhs, float *d__, float *e,
                 scomplex *b, integer *ldb, float *rcond, integer *rank, scomplex *work,
                 float *rwork, integer *iwork, integer *info);
float clangb_check(char *norm, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab,
                   float *work);
float clange_check(char *norm, integer *m, integer *n, scomplex *a, integer *lda, float *work);
float clangt_check(char *norm, integer *n, scomplex *dl, scomplex *d__, scomplex *du);
float clanhb_check(char *norm, char *uplo, integer *n, integer *k, scomplex *ab, integer *ldab,
                   float *work);
float clanhe_check(char *norm, char *uplo, integer *n, scomplex *a, integer *lda, float *work);
float clanhf_check(char *norm, char *transr, char *uplo, integer *n, scomplex *a, float *work);
float clanhp_check(char *norm, char *uplo, integer *n, scomplex *ap, float *work);
float clanhs_check(char *norm, integer *n, scomplex *a, integer *lda, float *work);
float clanht_check(char *norm, integer *n, float *d__, scomplex *e);
float clansb_check(char *norm, char *uplo, integer *n, integer *k, scomplex *ab, integer *ldab,
                   float *work);
float clansp_check(char *norm, char *uplo, integer *n, scomplex *ap, float *work);
float clansy_check(char *norm, char *uplo, integer *n, scomplex *a, integer *lda, float *work);
float clantb_check(char *norm, char *uplo, char *diag, integer *n, integer *k, scomplex *ab,
                   integer *ldab, float *work);
float clantp_check(char *norm, char *uplo, char *diag, integer *n, scomplex *ap, float *work);
float clantr_check(char *norm, char *uplo, char *diag, integer *m, integer *n, scomplex *a,
                   integer *lda, float *work);
int clapll_check(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *ssmin);
int clapmr_check(logical *forwrd, integer *m, integer *n, scomplex *x, integer *ldx, integer *k);
int clapmt_check(logical *forwrd, integer *m, integer *n, scomplex *x, integer *ldx, integer *k);
int claqgb_check(integer *m, integer *n, integer *kl, integer *ku, scomplex *ab, integer *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *equed);
int claqge_check(integer *m, integer *n, scomplex *a, integer *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, char *equed);
int claqhb_check(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *s,
                 float *scond, float *amax, char *equed);
int claqhe_check(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond,
                 float *amax, char *equed);
int claqhp_check(char *uplo, integer *n, scomplex *ap, float *s, float *scond, float *amax,
                 char *equed);
int claqp2_check(integer *m, integer *n, integer *offset, scomplex *a, integer *lda, integer *jpvt,
                 scomplex *tau, float *vn1, float *vn2, scomplex *work);
int claqps_check(integer *m, integer *n, integer *offset, integer *nb, integer *kb, scomplex *a,
                 integer *lda, integer *jpvt, scomplex *tau, float *vn1, float *vn2, scomplex *auxv,
                 scomplex *f, integer *ldf);
int claqr0_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 scomplex *h__, integer *ldh, scomplex *w, integer *iloz, integer *ihiz,
                 scomplex *z__, integer *ldz, scomplex *work, integer *lwork, integer *info);
int claqr1_check(integer *n, scomplex *h__, integer *ldh, scomplex *s1, scomplex *s2, scomplex *v);
int claqr2_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, scomplex *h__, integer *ldh, integer *iloz, integer *ihiz,
                 scomplex *z__, integer *ldz, integer *ns, integer *nd, scomplex *sh, scomplex *v,
                 integer *ldv, integer *nh, scomplex *t, integer *ldt, integer *nv, scomplex *wv,
                 integer *ldwv, scomplex *work, integer *lwork);
int claqr3_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, scomplex *h__, integer *ldh, integer *iloz, integer *ihiz,
                 scomplex *z__, integer *ldz, integer *ns, integer *nd, scomplex *sh, scomplex *v,
                 integer *ldv, integer *nh, scomplex *t, integer *ldt, integer *nv, scomplex *wv,
                 integer *ldwv, scomplex *work, integer *lwork);
int claqr4_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 scomplex *h__, integer *ldh, scomplex *w, integer *iloz, integer *ihiz,
                 scomplex *z__, integer *ldz, scomplex *work, integer *lwork, integer *info);
int claqr5_check(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
                 integer *kbot, integer *nshfts, scomplex *s, scomplex *h__, integer *ldh,
                 integer *iloz, integer *ihiz, scomplex *z__, integer *ldz, scomplex *v,
                 integer *ldv, scomplex *u, integer *ldu, integer *nv, scomplex *wv, integer *ldwv,
                 integer *nh, scomplex *wh, integer *ldwh);
int claqsb_check(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *s,
                 float *scond, float *amax, char *equed);
int claqsp_check(char *uplo, integer *n, scomplex *ap, float *s, float *scond, float *amax,
                 char *equed);
int claqsy_check(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond,
                 float *amax, char *equed);
int clar1v_check(integer *n, integer *b1, integer *bn, float *lambda, float *d__, float *l,
                 float *ld, float *lld, float *pivmin, float *gaptol, scomplex *z__,
                 logical *wantnc, integer *negcnt, float *ztz, float *mingma, integer *r__,
                 integer *isuppz, float *nrminv, float *resid, float *rqcorr, float *work);
int clar2v_check(integer *n, scomplex *x, scomplex *y, scomplex *z__, integer *incx, float *c__,
                 scomplex *s, integer *incc);
int clarcm_check(integer *m, integer *n, float *a, integer *lda, scomplex *b, integer *ldb,
                 scomplex *c__, integer *ldc, float *rwork);
int clarf_check(char *side, integer *m, integer *n, scomplex *v, integer *incv, scomplex *tau,
                scomplex *c__, integer *ldc, scomplex *work);
int clarfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, scomplex *v, integer *ldv, scomplex *t, integer *ldt, scomplex *c__,
                 integer *ldc, scomplex *work, integer *ldwork);
int clarfg_check(integer *n, scomplex *alpha, scomplex *x, integer *incx, scomplex *tau);
int clarfgp_check(integer *n, scomplex *alpha, scomplex *x, integer *incx, scomplex *tau);
int clarft_check(char *direct, char *storev, integer *n, integer *k, scomplex *v, integer *ldv,
                 scomplex *tau, scomplex *t, integer *ldt);
int clarfx_check(char *side, integer *m, integer *n, scomplex *v, scomplex *tau, scomplex *c__,
                 integer *ldc, scomplex *work);
int clargv_check(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *c__,
                 integer *incc);
int clarnv_check(integer *idist, integer *iseed, integer *n, scomplex *x);
int clarrv_check(integer *n, float *vl, float *vu, float *d__, float *l, float *pivmin,
                 integer *isplit, integer *m, integer *dol, integer *dou, float *minrgp,
                 float *rtol1, float *rtol2, float *w, float *werr, float *wgap, integer *iblock,
                 integer *indexw, float *gers, scomplex *z__, integer *ldz, integer *isuppz,
                 float *work, integer *iwork, integer *info);
int clarscl2_check(integer *m, integer *n, float *d__, scomplex *x, integer *ldx);
int clartg_check(scomplex *f, scomplex *g, float *cs, scomplex *sn, scomplex *r__);
int clartv_check(integer *n, scomplex *x, integer *incx, scomplex *y, integer *incy, float *c__,
                 scomplex *s, integer *incc);
int clarz_check(char *side, integer *m, integer *n, integer *l, scomplex *v, integer *incv,
                scomplex *tau, scomplex *c__, integer *ldc, scomplex *work);
int clarzb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, scomplex *v, integer *ldv, scomplex *t, integer *ldt,
                 scomplex *c__, integer *ldc, scomplex *work, integer *ldwork);
int clarzt_check(char *direct, char *storev, integer *n, integer *k, scomplex *v, integer *ldv,
                 scomplex *tau, scomplex *t, integer *ldt);
int clascl_check(char *type__, integer *kl, integer *ku, float *cfrom, float *cto, integer *m,
                 integer *n, scomplex *a, integer *lda, integer *info);
int clascl2_check(integer *m, integer *n, float *d__, scomplex *x, integer *ldx);
int claset_check(char *uplo, integer *m, integer *n, scomplex *alpha, scomplex *beta, scomplex *a,
                 integer *lda);
int clasr_check(char *side, char *pivot, char *direct, integer *m, integer *n, float *c__, float *s,
                scomplex *a, integer *lda);
int classq_check(integer *n, scomplex *x, integer *incx, float *scale, float *sumsq);
int claswp_check(integer *n, scomplex *a, integer *lda, integer *k1, integer *k2, integer *ipiv,
                 integer *incx);
int clasyf_check(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda,
                 integer *ipiv, scomplex *w, integer *ldw, integer *info);
int clasyf_rook_check(char *uplo, integer *n, integer *nb, integer *kb, scomplex *a, integer *lda,
                      integer *ipiv, scomplex *w, integer *ldw, integer *info);
int clatbs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd,
                 scomplex *ab, integer *ldab, scomplex *x, float *scale, float *cnorm,
                 integer *info);
int clatdf_check(integer *ijob, integer *n, scomplex *z__, integer *ldz, scomplex *rhs,
                 float *rdsum, float *rdscal, integer *ipiv, integer *jpiv);
int clatps_check(char *uplo, char *trans, char *diag, char *normin, integer *n, scomplex *ap,
                 scomplex *x, float *scale, float *cnorm, integer *info);
int clatrd_check(char *uplo, integer *n, integer *nb, scomplex *a, integer *lda, float *e,
                 scomplex *tau, scomplex *w, integer *ldw);
int clatrs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, scomplex *a,
                 integer *lda, scomplex *x, float *scale, float *cnorm, integer *info);
int clatrz_check(integer *m, integer *n, integer *l, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work);
int clatzm_check(char *side, integer *m, integer *n, scomplex *v, integer *incv, scomplex *tau,
                 scomplex *c1, scomplex *c2, integer *ldc, scomplex *work);
int clauu2_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *info);
int clauum_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *info);
int cpbcon_check(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *anorm,
                 float *rcond, scomplex *work, float *rwork, integer *info);
int cpbequ_check(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, float *s,
                 float *scond, float *amax, integer *info);
int cpbrfs_check(char *uplo, integer *n, integer *kd, integer *nrhs, scomplex *ab, integer *ldab,
                 scomplex *afb, integer *ldafb, scomplex *b, integer *ldb, scomplex *x,
                 integer *ldx, float *ferr, float *berr, scomplex *work, float *rwork,
                 integer *info);
int cpbstf_check(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, integer *info);
int cpbsv_check(char *uplo, integer *n, integer *kd, integer *nrhs, scomplex *ab, integer *ldab,
                scomplex *b, integer *ldb, integer *info);
int cpbsvx_check(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, scomplex *ab,
                 integer *ldab, scomplex *afb, integer *ldafb, char *equed, float *s, scomplex *b,
                 integer *ldb, scomplex *x, integer *ldx, float *rcond, float *ferr, float *berr,
                 scomplex *work, float *rwork, integer *info);
int cpbtf2_check(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, integer *info);
int cpbtrf_check(char *uplo, integer *n, integer *kd, scomplex *ab, integer *ldab, integer *info);
int cpbtrs_check(char *uplo, integer *n, integer *kd, integer *nrhs, scomplex *ab, integer *ldab,
                 scomplex *b, integer *ldb, integer *info);
int cpftrf_check(char *transr, char *uplo, integer *n, scomplex *a, integer *info);
int cpftri_check(char *transr, char *uplo, integer *n, scomplex *a, integer *info);
int cpftrs_check(char *transr, char *uplo, integer *n, integer *nrhs, scomplex *a, scomplex *b,
                 integer *ldb, integer *info);
int cpocon_check(char *uplo, integer *n, scomplex *a, integer *lda, float *anorm, float *rcond,
                 scomplex *work, float *rwork, integer *info);
int cpoequ_check(integer *n, scomplex *a, integer *lda, float *s, float *scond, float *amax,
                 integer *info);
int cpoequb_check(integer *n, scomplex *a, integer *lda, float *s, float *scond, float *amax,
                  integer *info);
int cporfs_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af,
                 integer *ldaf, scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *ferr,
                 float *berr, scomplex *work, float *rwork, integer *info);
int cporfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, float *s, scomplex *b, integer *ldb, scomplex *x,
                  integer *ldx, float *rcond, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  scomplex *work, float *rwork, integer *info);
int cposv_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *b,
                integer *ldb, integer *info);
int cposvx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                 scomplex *af, integer *ldaf, char *equed, float *s, scomplex *b, integer *ldb,
                 scomplex *x, integer *ldx, float *rcond, float *ferr, float *berr, scomplex *work,
                 float *rwork, integer *info);
int cposvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, char *equed, float *s, scomplex *b, integer *ldb,
                  scomplex *x, integer *ldx, float *rcond, float *rpvgrw, float *berr,
                  integer *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  integer *nparams, float *params, scomplex *work, float *rwork, integer *info);
int cpotf2_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *info);
int cpotrf_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *info);
int cpotri_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *info);
int cpotrs_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, integer *info);
int cppcon_check(char *uplo, integer *n, scomplex *ap, float *anorm, float *rcond, scomplex *work,
                 float *rwork, integer *info);
int cppequ_check(char *uplo, integer *n, scomplex *ap, float *s, float *scond, float *amax,
                 integer *info);
int cpprfs_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *afp, scomplex *b,
                 integer *ldb, scomplex *x, integer *ldx, float *ferr, float *berr, scomplex *work,
                 float *rwork, integer *info);
int cppsv_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *b, integer *ldb,
                integer *info);
int cppsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *afp,
                 char *equed, float *s, scomplex *b, integer *ldb, scomplex *x, integer *ldx,
                 float *rcond, float *ferr, float *berr, scomplex *work, float *rwork,
                 integer *info);
int cpptrf_check(char *uplo, integer *n, scomplex *ap, integer *info);
int cpptri_check(char *uplo, integer *n, scomplex *ap, integer *info);
int cpptrs_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *b, integer *ldb,
                 integer *info);
int cpstf2_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *piv, integer *rank,
                 float *tol, float *work, integer *info);
int cpstrf_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *piv, integer *rank,
                 float *tol, float *work, integer *info);
int cptcon_check(integer *n, float *d__, scomplex *e, float *anorm, float *rcond, float *rwork,
                 integer *info);
int cpteqr_check(char *compz, integer *n, float *d__, float *e, scomplex *z__, integer *ldz,
                 float *work, integer *info);
int cptrfs_check(char *uplo, integer *n, integer *nrhs, float *d__, scomplex *e, float *df,
                 scomplex *ef, scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *ferr,
                 float *berr, scomplex *work, float *rwork, integer *info);
int cptsv_check(integer *n, integer *nrhs, float *d__, scomplex *e, scomplex *b, integer *ldb,
                integer *info);
int cptsvx_check(char *fact, integer *n, integer *nrhs, float *d__, scomplex *e, float *df,
                 scomplex *ef, scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int cpttrf_check(integer *n, float *d__, scomplex *e, integer *info);
int cpttrs_check(char *uplo, integer *n, integer *nrhs, float *d__, scomplex *e, scomplex *b,
                 integer *ldb, integer *info);
int cptts2_check(integer *iuplo, integer *n, integer *nrhs, float *d__, scomplex *e, scomplex *b,
                 integer *ldb);
int crot_check(integer *n, scomplex *cx, integer *incx, scomplex *cy, integer *incy, float *c__,
               scomplex *s);
int cspcon_check(char *uplo, integer *n, scomplex *ap, integer *ipiv, float *anorm, float *rcond,
                 scomplex *work, integer *info);
int cspmv_check(char *uplo, integer *n, scomplex *alpha, scomplex *ap, scomplex *x, integer *incx,
                scomplex *beta, scomplex *y, integer *incy);
int cspr_check(char *uplo, integer *n, scomplex *alpha, scomplex *x, integer *incx, scomplex *ap);
int csprfs_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *afp, integer *ipiv,
                 scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, integer *info);
int cspsv_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, integer *ipiv, scomplex *b,
                integer *ldb, integer *info);
int cspsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *ap, scomplex *afp,
                 integer *ipiv, scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int csptrf_check(char *uplo, integer *n, scomplex *ap, integer *ipiv, integer *info);
int csptri_check(char *uplo, integer *n, scomplex *ap, integer *ipiv, scomplex *work,
                 integer *info);
int csptrs_check(char *uplo, integer *n, integer *nrhs, scomplex *ap, integer *ipiv, scomplex *b,
                 integer *ldb, integer *info);
int csrscl_check(integer *n, float *sa, scomplex *sx, integer *incx);
int cstedc_check(char *compz, integer *n, float *d__, float *e, scomplex *z__, integer *ldz,
                 scomplex *work, integer *lwork, float *rwork, integer *lrwork, integer *iwork,
                 integer *liwork, integer *info);
int cstegr_check(char *jobz, char *range, integer *n, float *d__, float *e, float *vl, float *vu,
                 integer *il, integer *iu, float *abstol, integer *m, float *w, scomplex *z__,
                 integer *ldz, integer *isuppz, float *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int cstein_check(integer *n, float *d__, float *e, integer *m, float *w, integer *iblock,
                 integer *isplit, scomplex *z__, integer *ldz, float *work, integer *iwork,
                 integer *ifail, integer *info);
int cstemr_check(char *jobz, char *range, integer *n, float *d__, float *e, float *vl, float *vu,
                 integer *il, integer *iu, integer *m, float *w, scomplex *z__, integer *ldz,
                 integer *nzc, integer *isuppz, logical *tryrac, float *work, integer *lwork,
                 integer *iwork, integer *liwork, integer *info);
int csteqr_check(char *compz, integer *n, float *d__, float *e, scomplex *z__, integer *ldz,
                 float *work, integer *info);
int csycon_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, float *anorm,
                 float *rcond, scomplex *work, integer *info);
int csycon_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      float *anorm, float *rcond, scomplex *work, integer *info);
int csyconv_check(char *uplo, char *way, integer *n, scomplex *a, integer *lda, integer *ipiv,
                  scomplex *work, integer *info);
int csyequb_check(char *uplo, integer *n, scomplex *a, integer *lda, float *s, float *scond,
                  float *amax, scomplex *work, integer *info);
int csymv_check(char *uplo, integer *n, scomplex *alpha, scomplex *a, integer *lda, scomplex *x,
                integer *incx, scomplex *beta, scomplex *y, integer *incy);
int csyr_check(char *uplo, integer *n, scomplex *alpha, scomplex *x, integer *incx, scomplex *a,
               integer *lda);
int csyrfs_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, scomplex *af,
                 integer *ldaf, integer *ipiv, scomplex *b, integer *ldb, scomplex *x, integer *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int csyrfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, integer *ipiv, float *s, scomplex *b, integer *ldb,
                  scomplex *x, integer *ldx, float *rcond, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  scomplex *work, float *rwork, integer *info);
int csysv_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv,
                scomplex *b, integer *ldb, scomplex *work, integer *lwork, integer *info);
int csysv_rook_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                     integer *ipiv, scomplex *b, integer *ldb, scomplex *work, integer *lwork,
                     integer *info);
int csysvx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                 scomplex *af, integer *ldaf, integer *ipiv, scomplex *b, integer *ldb, scomplex *x,
                 integer *ldx, float *rcond, float *ferr, float *berr, scomplex *work,
                 integer *lwork, float *rwork, integer *info);
int csysvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                  scomplex *af, integer *ldaf, integer *ipiv, char *equed, float *s, scomplex *b,
                  integer *ldb, scomplex *x, integer *ldx, float *rcond, float *rpvgrw, float *berr,
                  integer *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  integer *nparams, float *params, scomplex *work, float *rwork, integer *info);
int csyswapr_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *i1, integer *i2);
int csytf2_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, integer *info);
int csytf2_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      integer *info);
int csytrf_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                 integer *lwork, integer *info);
int csytrf_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      scomplex *work, integer *lwork, integer *info);
int csytri_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                 integer *info);
int csytri2_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                  integer *lwork, integer *info);
int csytri2x_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv, scomplex *work,
                   integer *nb, integer *info);
int csytri_rook_check(char *uplo, integer *n, scomplex *a, integer *lda, integer *ipiv,
                      scomplex *work, integer *info);
int csytrs_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv,
                 scomplex *b, integer *ldb, integer *info);
int csytrs2_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda, integer *ipiv,
                  scomplex *b, integer *ldb, scomplex *work, integer *info);
int csytrs_rook_check(char *uplo, integer *n, integer *nrhs, scomplex *a, integer *lda,
                      integer *ipiv, scomplex *b, integer *ldb, integer *info);
int ctbcon_check(char *norm, char *uplo, char *diag, integer *n, integer *kd, scomplex *ab,
                 integer *ldab, float *rcond, scomplex *work, float *rwork, integer *info);
int ctbrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 scomplex *ab, integer *ldab, scomplex *b, integer *ldb, scomplex *x, integer *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, integer *info);
int ctbtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 scomplex *ab, integer *ldab, scomplex *b, integer *ldb, integer *info);
int ctfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m,
                integer *n, scomplex *alpha, scomplex *a, scomplex *b, integer *ldb);
int ctftri_check(char *transr, char *uplo, char *diag, integer *n, scomplex *a, integer *info);
int ctfttp_check(char *transr, char *uplo, integer *n, scomplex *arf, scomplex *ap, integer *info);
int ctfttr_check(char *transr, char *uplo, integer *n, scomplex *arf, scomplex *a, integer *lda,
                 integer *info);
int ctgevc_check(char *side, char *howmny, logical *select, integer *n, scomplex *s, integer *lds,
                 scomplex *p, integer *ldp, scomplex *vl, integer *ldvl, scomplex *vr,
                 integer *ldvr, integer *mm, integer *m, scomplex *work, float *rwork,
                 integer *info);
int ctgex2_check(logical *wantq, logical *wantz, integer *n, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, scomplex *q, integer *ldq, scomplex *z__, integer *ldz, integer *j1,
                 integer *info);
int ctgexc_check(logical *wantq, logical *wantz, integer *n, scomplex *a, integer *lda, scomplex *b,
                 integer *ldb, scomplex *q, integer *ldq, scomplex *z__, integer *ldz,
                 integer *ifst, integer *ilst, integer *info);
int ctgsen_check(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n,
                 scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *alpha,
                 scomplex *beta, scomplex *q, integer *ldq, scomplex *z__, integer *ldz, integer *m,
                 float *pl, float *pr, float *dif, scomplex *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int ctgsja_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k,
                 integer *l, scomplex *a, integer *lda, scomplex *b, integer *ldb, float *tola,
                 float *tolb, float *alpha, float *beta, scomplex *u, integer *ldu, scomplex *v,
                 integer *ldv, scomplex *q, integer *ldq, scomplex *work, integer *ncycle,
                 integer *info);
int ctgsna_check(char *job, char *howmny, logical *select, integer *n, scomplex *a, integer *lda,
                 scomplex *b, integer *ldb, scomplex *vl, integer *ldvl, scomplex *vr,
                 integer *ldvr, float *s, float *dif, integer *mm, integer *m, scomplex *work,
                 integer *lwork, integer *iwork, integer *info);
int ctgsy2_check(char *trans, integer *ijob, integer *m, integer *n, scomplex *a, integer *lda,
                 scomplex *b, integer *ldb, scomplex *c__, integer *ldc, scomplex *d__,
                 integer *ldd, scomplex *e, integer *lde, scomplex *f, integer *ldf, float *scale,
                 float *rdsum, float *rdscal, integer *info);
int ctgsyl_check(char *trans, integer *ijob, integer *m, integer *n, scomplex *a, integer *lda,
                 scomplex *b, integer *ldb, scomplex *c__, integer *ldc, scomplex *d__,
                 integer *ldd, scomplex *e, integer *lde, scomplex *f, integer *ldf, float *scale,
                 float *dif, scomplex *work, integer *lwork, integer *iwork, integer *info);
int ctpcon_check(char *norm, char *uplo, char *diag, integer *n, scomplex *ap, float *rcond,
                 scomplex *work, float *rwork, integer *info);
int ctpmqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                  integer *nb, scomplex *v, integer *ldv, scomplex *t, integer *ldt, scomplex *a,
                  integer *lda, scomplex *b, integer *ldb, scomplex *work, integer *info);
int ctpqrt_check(integer *m, integer *n, integer *l, integer *nb, scomplex *a, integer *lda,
                 scomplex *b, integer *ldb, scomplex *t, integer *ldt, scomplex *work,
                 integer *info);
int ctpqrt2_check(integer *m, integer *n, integer *l, scomplex *a, integer *lda, scomplex *b,
                  integer *ldb, scomplex *t, integer *ldt, integer *info);
int ctprfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, scomplex *v, integer *ldv, scomplex *t, integer *ldt,
                 scomplex *a, integer *lda, scomplex *b, integer *ldb, scomplex *work,
                 integer *ldwork);
int ctprfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, scomplex *ap,
                 scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, integer *info);
int ctptri_check(char *uplo, char *diag, integer *n, scomplex *ap, integer *info);
int ctptrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, scomplex *ap,
                 scomplex *b, integer *ldb, integer *info);
int ctpttf_check(char *transr, char *uplo, integer *n, scomplex *ap, scomplex *arf, integer *info);
int ctpttr_check(char *uplo, integer *n, scomplex *ap, scomplex *a, integer *lda, integer *info);
int ctrcon_check(char *norm, char *uplo, char *diag, integer *n, scomplex *a, integer *lda,
                 float *rcond, scomplex *work, float *rwork, integer *info);
int ctrevc_check(char *side, char *howmny, logical *select, integer *n, scomplex *t, integer *ldt,
                 scomplex *vl, integer *ldvl, scomplex *vr, integer *ldvr, integer *mm, integer *m,
                 scomplex *work, float *rwork, integer *info);
int ctrexc_check(char *compq, integer *n, scomplex *t, integer *ldt, scomplex *q, integer *ldq,
                 integer *ifst, integer *ilst, integer *info);
int ctrrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, scomplex *a,
                 integer *lda, scomplex *b, integer *ldb, scomplex *x, integer *ldx, float *ferr,
                 float *berr, scomplex *work, float *rwork, integer *info);
int ctrsen_check(char *job, char *compq, logical *select, integer *n, scomplex *t, integer *ldt,
                 scomplex *q, integer *ldq, scomplex *w, integer *m, float *s, float *sep,
                 scomplex *work, integer *lwork, integer *info);
int ctrsna_check(char *job, char *howmny, logical *select, integer *n, scomplex *t, integer *ldt,
                 scomplex *vl, integer *ldvl, scomplex *vr, integer *ldvr, float *s, float *sep,
                 integer *mm, integer *m, scomplex *work, integer *ldwork, float *rwork,
                 integer *info);
int ctrsyl_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, scomplex *a,
                 integer *lda, scomplex *b, integer *ldb, scomplex *c__, integer *ldc, float *scale,
                 integer *info);
int ctrsyl3_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, complex *a,
                  integer *lda, complex *b, integer *ldb, complex *c__, integer *ldc, real *scale,
                  real *swork, integer *ldswork, integer *info);
int ctrti2_check(char *uplo, char *diag, integer *n, scomplex *a, integer *lda, integer *info);
int ctrtri_check(char *uplo, char *diag, integer *n, scomplex *a, integer *lda, integer *info);
int ctrtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, scomplex *a,
                 integer *lda, scomplex *b, integer *ldb, integer *info);
int ctrttf_check(char *transr, char *uplo, integer *n, scomplex *a, integer *lda, scomplex *arf,
                 integer *info);
int ctrttp_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *ap, integer *info);
int ctzrqf_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, integer *info);
int ctzrzf_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *lwork, integer *info);
int cunbdb_check(char *trans, char *signs, integer *m, integer *p, integer *q, scomplex *x11,
                 integer *ldx11, scomplex *x12, integer *ldx12, scomplex *x21, integer *ldx21,
                 scomplex *x22, integer *ldx22, float *theta, float *phi, scomplex *taup1,
                 scomplex *taup2, scomplex *tauq1, scomplex *tauq2, scomplex *work, integer *lwork,
                 integer *info);
int cunbdb1_check(integer *m, integer *p, integer *q, scomplex *x11, integer *ldx11, scomplex *x21,
                  integer *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *work, integer *lwork, integer *info);
int cunbdb2_check(integer *m, integer *p, integer *q, scomplex *x11, integer *ldx11, scomplex *x21,
                  integer *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *work, integer *lwork, integer *info);
int cunbdb3_check(integer *m, integer *p, integer *q, scomplex *x11, integer *ldx11, scomplex *x21,
                  integer *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *work, integer *lwork, integer *info);
int cunbdb4_check(integer *m, integer *p, integer *q, scomplex *x11, integer *ldx11, scomplex *x21,
                  integer *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *phantom, scomplex *work, integer *lwork,
                  integer *info);
int cunbdb5_check(integer *m1, integer *m2, integer *n, scomplex *x1, integer *incx1, scomplex *x2,
                  integer *incx2, scomplex *q1, integer *ldq1, scomplex *q2, integer *ldq2,
                  scomplex *work, integer *lwork, integer *info);
int cunbdb6_check(integer *m1, integer *m2, integer *n, scomplex *x1, integer *incx1, scomplex *x2,
                  integer *incx2, scomplex *q1, integer *ldq1, scomplex *q2, integer *ldq2,
                  scomplex *work, integer *lwork, integer *info);
int cuncsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 integer *m, integer *p, integer *q, scomplex *x11, integer *ldx11, scomplex *x12,
                 integer *ldx12, scomplex *x21, integer *ldx21, scomplex *x22, integer *ldx22,
                 float *theta, scomplex *u1, integer *ldu1, scomplex *u2, integer *ldu2,
                 scomplex *v1t, integer *ldv1t, scomplex *v2t, integer *ldv2t, scomplex *work,
                 integer *lwork, float *rwork, integer *lrwork, integer *iwork, integer *info);
int cuncsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, integer *m, integer *p, integer *q,
                     scomplex *x11, integer *ldx11, scomplex *x21, integer *ldx21, float *theta,
                     scomplex *u1, integer *ldu1, scomplex *u2, integer *ldu2, scomplex *v1t,
                     integer *ldv1t, scomplex *work, integer *lwork, float *rwork, integer *lrwork,
                     integer *iwork, integer *info);
int cung2l_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *info);
int cung2r_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *info);
int cungbr_check(char *vect, integer *m, integer *n, integer *k, scomplex *a, integer *lda,
                 scomplex *tau, scomplex *work, integer *lwork, integer *info);
int cunghr_check(integer *n, integer *ilo, integer *ihi, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *lwork, integer *info);
int cungl2_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *info);
int cunglq_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *lwork, integer *info);
int cungql_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *lwork, integer *info);
int cungqr_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *lwork, integer *info);
int cungr2_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *info);
int cungrq_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau,
                 scomplex *work, integer *lwork, integer *info);
int cungtr_check(char *uplo, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                 integer *lwork, integer *info);
int cunm2l_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *info);
int cunm2r_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *info);
int cunmbr_check(char *vect, char *side, char *trans, integer *m, integer *n, integer *k,
                 scomplex *a, integer *lda, scomplex *tau, scomplex *c__, integer *ldc,
                 scomplex *work, integer *lwork, integer *info);
int cunmhr_check(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi,
                 scomplex *a, integer *lda, scomplex *tau, scomplex *c__, integer *ldc,
                 scomplex *work, integer *lwork, integer *info);
int cunml2_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *info);
int cunmlq_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *lwork, integer *info);
int cunmql_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *lwork, integer *info);
int cunmqr_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *lwork, integer *info);
int cunmr2_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *info);
int cunmr3_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                 scomplex *a, integer *lda, scomplex *tau, scomplex *c__, integer *ldc,
                 scomplex *work, integer *info);
int cunmrq_check(char *side, char *trans, integer *m, integer *n, integer *k, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *lwork, integer *info);
int cunmrz_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                 scomplex *a, integer *lda, scomplex *tau, scomplex *c__, integer *ldc,
                 scomplex *work, integer *lwork, integer *info);
int cunmtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, scomplex *a,
                 integer *lda, scomplex *tau, scomplex *c__, integer *ldc, scomplex *work,
                 integer *lwork, integer *info);
int cupgtr_check(char *uplo, integer *n, scomplex *ap, scomplex *tau, scomplex *q, integer *ldq,
                 scomplex *work, integer *info);
int cupmtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, scomplex *ap,
                 scomplex *tau, scomplex *c__, integer *ldc, scomplex *work, integer *info);
int dbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m,
                 integer *p, integer *q, double *theta, double *phi, double *u1, integer *ldu1,
                 double *u2, integer *ldu2, double *v1t, integer *ldv1t, double *v2t,
                 integer *ldv2t, double *b11d, double *b11e, double *b12d, double *b12e,
                 double *b21d, double *b21e, double *b22d, double *b22e, double *work,
                 integer *lwork, integer *info);
int dbdsdc_check(char *uplo, char *compq, integer *n, double *d__, double *e, double *u,
                 integer *ldu, double *vt, integer *ldvt, double *q, integer *iq, double *work,
                 integer *iwork, integer *info);
int dbdsqr_check(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, double *d__,
                 double *e, double *vt, integer *ldvt, double *u, integer *ldu, double *c__,
                 integer *ldc, double *work, integer *info);
int ddisna_check(char *job, integer *m, integer *n, double *d__, double *sep, integer *info);
int dgbbrd_check(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku,
                 double *ab, integer *ldab, double *d__, double *e, double *q, integer *ldq,
                 double *pt, integer *ldpt, double *c__, integer *ldc, double *work, integer *info);
int dgbcon_check(char *norm, integer *n, integer *kl, integer *ku, double *ab, integer *ldab,
                 integer *ipiv, double *anorm, double *rcond, double *work, integer *iwork,
                 integer *info);
int dgbequ_check(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 integer *info);
int dgbequb_check(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab,
                  double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                  integer *info);
int dgbrfs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, double *ab,
                 integer *ldab, double *afb, integer *ldafb, integer *ipiv, double *b, integer *ldb,
                 double *x, integer *ldx, double *ferr, double *berr, double *work, integer *iwork,
                 integer *info);
int dgbrfsx_check(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs,
                  double *ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv,
                  double *r__, double *c__, double *b, integer *ldb, double *x, integer *ldx,
                  double *rcond, double *berr, integer *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, integer *nparams, double *params, double *work,
                  integer *iwork, integer *info);
int dgbsv_check(integer *n, integer *kl, integer *ku, integer *nrhs, double *ab, integer *ldab,
                integer *ipiv, double *b, integer *ldb, integer *info);
int dgbsvx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                 double *ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv, char *equed,
                 double *r__, double *c__, double *b, integer *ldb, double *x, integer *ldx,
                 double *rcond, double *ferr, double *berr, double *work, integer *iwork,
                 integer *info);
int dgbsvxx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                  double *ab, integer *ldab, double *afb, integer *ldafb, integer *ipiv,
                  char *equed, double *r__, double *c__, double *b, integer *ldb, double *x,
                  integer *ldx, double *rcond, double *rpvgrw, double *berr, integer *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams,
                  double *params, double *work, integer *iwork, integer *info);
int dgbtf2_check(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab,
                 integer *ipiv, integer *info);
int dgbtrf_check(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab,
                 integer *ipiv, integer *info);
int dgbtrs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, double *ab,
                 integer *ldab, integer *ipiv, double *b, integer *ldb, integer *info);
int dgebak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, double *scale,
                 integer *m, double *v, integer *ldv, integer *info);
int dgebal_check(char *job, integer *n, double *a, integer *lda, integer *ilo, integer *ihi,
                 double *scale, integer *info);
int dgebd2_check(integer *m, integer *n, double *a, integer *lda, double *d__, double *e,
                 double *tauq, double *taup, double *work, integer *info);
int dgebrd_check(integer *m, integer *n, double *a, integer *lda, double *d__, double *e,
                 double *tauq, double *taup, double *work, integer *lwork, integer *info);
int dgecon_check(char *norm, integer *n, double *a, integer *lda, double *anorm, double *rcond,
                 double *work, integer *iwork, integer *info);
int dgeequ_check(integer *m, integer *n, double *a, integer *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, integer *info);
int dgeequb_check(integer *m, integer *n, double *a, integer *lda, double *r__, double *c__,
                  double *rowcnd, double *colcnd, double *amax, integer *info);
int dgees_check(char *jobvs, char *sort, L_fp select, integer *n, double *a, integer *lda,
                integer *sdim, double *wr, double *wi, double *vs, integer *ldvs, double *work,
                integer *lwork, logical *bwork, integer *info);
int dgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, integer *n, double *a,
                 integer *lda, integer *sdim, double *wr, double *wi, double *vs, integer *ldvs,
                 double *rconde, double *rcondv, double *work, integer *lwork, integer *iwork,
                 integer *liwork, logical *bwork, integer *info);
int dgeev_check(char *jobvl, char *jobvr, integer *n, double *a, integer *lda, double *wr,
                double *wi, double *vl, integer *ldvl, double *vr, integer *ldvr, double *work,
                integer *lwork, integer *info);
int dgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, double *a,
                 integer *lda, double *wr, double *wi, double *vl, integer *ldvl, double *vr,
                 integer *ldvr, integer *ilo, integer *ihi, double *scale, double *abnrm,
                 double *rconde, double *rcondv, double *work, integer *lwork, integer *iwork,
                 integer *info);
int dgegs_check(char *jobvsl, char *jobvsr, integer *n, double *a, integer *lda, double *b,
                integer *ldb, double *alphar, double *alphai, double *beta, double *vsl,
                integer *ldvsl, double *vsr, integer *ldvsr, double *work, integer *lwork,
                integer *info);
int dgegv_check(char *jobvl, char *jobvr, integer *n, double *a, integer *lda, double *b,
                integer *ldb, double *alphar, double *alphai, double *beta, double *vl,
                integer *ldvl, double *vr, integer *ldvr, double *work, integer *lwork,
                integer *info);
int dgehd2_check(integer *n, integer *ilo, integer *ihi, double *a, integer *lda, double *tau,
                 double *work, integer *info);
int dgehrd_check(integer *n, integer *ilo, integer *ihi, double *a, integer *lda, double *tau,
                 double *work, integer *lwork, integer *info);
int dgejsv_check(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, integer *m,
                 integer *n, double *a, integer *lda, double *sva, double *u, integer *ldu,
                 double *v, integer *ldv, double *work, integer *lwork, integer *iwork,
                 integer *info);
int dgelq2_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *info);
int dgelqf_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *lwork, integer *info);
int dgels_check(char *trans, integer *m, integer *n, integer *nrhs, double *a, integer *lda,
                double *b, integer *ldb, double *work, integer *lwork, integer *info);
int dgelsd_check(integer *m, integer *n, integer *nrhs, double *a, integer *lda, double *b,
                 integer *ldb, double *s, double *rcond, integer *rank, double *work,
                 integer *lwork, integer *iwork, integer *info);
int dgelss_check(integer *m, integer *n, integer *nrhs, double *a, integer *lda, double *b,
                 integer *ldb, double *s, double *rcond, integer *rank, double *work,
                 integer *lwork, integer *info);
int dgelsx_check(integer *m, integer *n, integer *nrhs, double *a, integer *lda, double *b,
                 integer *ldb, integer *jpvt, double *rcond, integer *rank, double *work,
                 integer *info);
int dgelsy_check(integer *m, integer *n, integer *nrhs, double *a, integer *lda, double *b,
                 integer *ldb, integer *jpvt, double *rcond, integer *rank, double *work,
                 integer *lwork, integer *info);
int dgemqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *nb,
                  double *v, integer *ldv, double *t, integer *ldt, double *c__, integer *ldc,
                  double *work, integer *info);
int dgeql2_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *info);
int dgeqlf_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *lwork, integer *info);
int dgeqp3_check(integer *m, integer *n, double *a, integer *lda, integer *jpvt, double *tau,
                 double *work, integer *lwork, integer *info);
int dgeqpf_check(integer *m, integer *n, double *a, integer *lda, integer *jpvt, double *tau,
                 double *work, integer *info);
int dgeqr2_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *info);
int dgeqr2p_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                  integer *info);
int dgeqrf_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *lwork, integer *info);
int dgeqrfp_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                  integer *lwork, integer *info);
int dgeqrt_check(integer *m, integer *n, integer *nb, double *a, integer *lda, double *t,
                 integer *ldt, double *work, integer *info);
int dgeqrt2_check(integer *m, integer *n, double *a, integer *lda, double *t, integer *ldt,
                  integer *info);
int dgeqrt3_check(integer *m, integer *n, double *a, integer *lda, double *t, integer *ldt,
                  integer *info);
int dgerfs_check(char *trans, integer *n, integer *nrhs, double *a, integer *lda, double *af,
                 integer *ldaf, integer *ipiv, double *b, integer *ldb, double *x, integer *ldx,
                 double *ferr, double *berr, double *work, integer *iwork, integer *info);
int dgerfsx_check(char *trans, char *equed, integer *n, integer *nrhs, double *a, integer *lda,
                  double *af, integer *ldaf, integer *ipiv, double *r__, double *c__, double *b,
                  integer *ldb, double *x, integer *ldx, double *rcond, double *berr,
                  integer *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  integer *nparams, double *params, double *work, integer *iwork, integer *info);
int dgerq2_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *info);
int dgerqf_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *lwork, integer *info);
int dgesc2_check(integer *n, double *a, integer *lda, double *rhs, integer *ipiv, integer *jpiv,
                 double *scale);
int dgesdd_check(char *jobz, integer *m, integer *n, double *a, integer *lda, double *s, double *u,
                 integer *ldu, double *vt, integer *ldvt, double *work, integer *lwork,
                 integer *iwork, integer *info);
int dgesdd_fla_check(char *jobu, char *jobvt, integer *m, integer *n, double *a, integer *lda,
                     double *s, double *u, integer *ldu, double *vt, integer *ldvt, double *work,
                     integer *lwork, integer *info);
int dgesv_check(integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv, double *b,
                integer *ldb, integer *info);
int dgesvd_check(char *jobu, char *jobvt, integer *m, integer *n, double *a, integer *lda,
                 double *s, double *u, integer *ldu, double *vt, integer *ldvt, double *work,
                 integer *lwork, integer *info);
int dgesvj_check(char *joba, char *jobu, char *jobv, integer *m, integer *n, double *a,
                 integer *lda, double *sva, integer *mv, double *v, integer *ldv, double *work,
                 integer *lwork, integer *info);
int dgesvx_check(char *fact, char *trans, integer *n, integer *nrhs, double *a, integer *lda,
                 double *af, integer *ldaf, integer *ipiv, char *equed, double *r__, double *c__,
                 double *b, integer *ldb, double *x, integer *ldx, double *rcond, double *ferr,
                 double *berr, double *work, integer *iwork, integer *info);
int dgesvxx_check(char *fact, char *trans, integer *n, integer *nrhs, double *a, integer *lda,
                  double *af, integer *ldaf, integer *ipiv, char *equed, double *r__, double *c__,
                  double *b, integer *ldb, double *x, integer *ldx, double *rcond, double *rpvgrw,
                  double *berr, integer *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, integer *nparams, double *params, double *work,
                  integer *iwork, integer *info);
int dgetc2_check(integer *n, double *a, integer *lda, integer *ipiv, integer *jpiv, integer *info);
int dgetf2_check(integer *m, integer *n, double *a, integer *lda, integer *ipiv, integer *info);
int dgetrf_check(integer *m, integer *n, double *a, integer *lda, integer *ipiv, integer *info);
int dgetrfnp_check(integer *m, integer *n, double *a, integer *lda, integer *info);
int dgetrfnpi_check(integer *m, integer *n, integer *nfact, double *a, integer *lda, integer *info);
int dgetri_check(integer *n, double *a, integer *lda, integer *ipiv, double *work, integer *lwork,
                 integer *info);
int dgetrs_check(char *trans, integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv,
                 double *b, integer *ldb, integer *info);
int dggbak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, double *lscale,
                 double *rscale, integer *m, double *v, integer *ldv, integer *info);
int dggbal_check(char *job, integer *n, double *a, integer *lda, double *b, integer *ldb,
                 integer *ilo, integer *ihi, double *lscale, double *rscale, double *work,
                 integer *info);
int dgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, double *a,
                integer *lda, double *b, integer *ldb, integer *sdim, double *alphar,
                double *alphai, double *beta, double *vsl, integer *ldvsl, double *vsr,
                integer *ldvsr, double *work, integer *lwork, logical *bwork, integer *info);
int dggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n,
                 double *a, integer *lda, double *b, integer *ldb, integer *sdim, double *alphar,
                 double *alphai, double *beta, double *vsl, integer *ldvsl, double *vsr,
                 integer *ldvsr, double *rconde, double *rcondv, double *work, integer *lwork,
                 integer *iwork, integer *liwork, logical *bwork, integer *info);
int dggev_check(char *jobvl, char *jobvr, integer *n, double *a, integer *lda, double *b,
                integer *ldb, double *alphar, double *alphai, double *beta, double *vl,
                integer *ldvl, double *vr, integer *ldvr, double *work, integer *lwork,
                integer *info);
int dggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, double *a,
                 integer *lda, double *b, integer *ldb, double *alphar, double *alphai,
                 double *beta, double *vl, integer *ldvl, double *vr, integer *ldvr, integer *ilo,
                 integer *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm,
                 double *rconde, double *rcondv, double *work, integer *lwork, integer *iwork,
                 logical *bwork, integer *info);
int dggglm_check(integer *n, integer *m, integer *p, double *a, integer *lda, double *b,
                 integer *ldb, double *d__, double *x, double *y, double *work, integer *lwork,
                 integer *info);
int dgghrd_check(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, double *a,
                 integer *lda, double *b, integer *ldb, double *q, integer *ldq, double *z__,
                 integer *ldz, integer *info);
int dgglse_check(integer *m, integer *n, integer *p, double *a, integer *lda, double *b,
                 integer *ldb, double *c__, double *d__, double *x, double *work, integer *lwork,
                 integer *info);
int dggqrf_check(integer *n, integer *m, integer *p, double *a, integer *lda, double *taua,
                 double *b, integer *ldb, double *taub, double *work, integer *lwork,
                 integer *info);
int dggrqf_check(integer *m, integer *p, integer *n, double *a, integer *lda, double *taua,
                 double *b, integer *ldb, double *taub, double *work, integer *lwork,
                 integer *info);
int dggsvd_check(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k,
                 integer *l, double *a, integer *lda, double *b, integer *ldb, double *alpha,
                 double *beta, double *u, integer *ldu, double *v, integer *ldv, double *q,
                 integer *ldq, double *work, integer *iwork, integer *info);
int dggsvp_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, double *a,
                 integer *lda, double *b, integer *ldb, double *tola, double *tolb, integer *k,
                 integer *l, double *u, integer *ldu, double *v, integer *ldv, double *q,
                 integer *ldq, integer *iwork, double *tau, double *work, integer *info);
int dgsvj0_check(char *jobv, integer *m, integer *n, double *a, integer *lda, double *d__,
                 double *sva, integer *mv, double *v, integer *ldv, double *eps, double *sfmin,
                 double *tol, integer *nsweep, double *work, integer *lwork, integer *info);
int dgsvj1_check(char *jobv, integer *m, integer *n, integer *n1, double *a, integer *lda,
                 double *d__, double *sva, integer *mv, double *v, integer *ldv, double *eps,
                 double *sfmin, double *tol, integer *nsweep, double *work, integer *lwork,
                 integer *info);
int dgtcon_check(char *norm, integer *n, double *dl, double *d__, double *du, double *du2,
                 integer *ipiv, double *anorm, double *rcond, double *work, integer *iwork,
                 integer *info);
int dgtrfs_check(char *trans, integer *n, integer *nrhs, double *dl, double *d__, double *du,
                 double *dlf, double *df, double *duf, double *du2, integer *ipiv, double *b,
                 integer *ldb, double *x, integer *ldx, double *ferr, double *berr, double *work,
                 integer *iwork, integer *info);
int dgtsv_check(integer *n, integer *nrhs, double *dl, double *d__, double *du, double *b,
                integer *ldb, integer *info);
int dgtsvx_check(char *fact, char *trans, integer *n, integer *nrhs, double *dl, double *d__,
                 double *du, double *dlf, double *df, double *duf, double *du2, integer *ipiv,
                 double *b, integer *ldb, double *x, integer *ldx, double *rcond, double *ferr,
                 double *berr, double *work, integer *iwork, integer *info);
int dgttrf_check(integer *n, double *dl, double *d__, double *du, double *du2, integer *ipiv,
                 integer *info);
int dgttrs_check(char *trans, integer *n, integer *nrhs, double *dl, double *d__, double *du,
                 double *du2, integer *ipiv, double *b, integer *ldb, integer *info);
int dgtts2_check(integer *itrans, integer *n, integer *nrhs, double *dl, double *d__, double *du,
                 double *du2, integer *ipiv, double *b, integer *ldb);
int dhgeqz_check(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi,
                 double *h__, integer *ldh, double *t, integer *ldt, double *alphar, double *alphai,
                 double *beta, double *q, integer *ldq, double *z__, integer *ldz, double *work,
                 integer *lwork, integer *info);
int dhsein_check(char *side, char *eigsrc, char *initv, logical *select, integer *n, double *h__,
                 integer *ldh, double *wr, double *wi, double *vl, integer *ldvl, double *vr,
                 integer *ldvr, integer *mm, integer *m, double *work, integer *ifaill,
                 integer *ifailr, integer *info);
int dhseqr_check(char *job, char *compz, integer *n, integer *ilo, integer *ihi, double *h__,
                 integer *ldh, double *wr, double *wi, double *z__, integer *ldz, double *work,
                 integer *lwork, integer *info);
logical disnan_check(double *din);
int dla_gbamv_check(integer *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha,
                    double *ab, integer *ldab, double *x, integer *incx, double *beta, double *y,
                    integer *incy);
double dla_gbrcond_check(char *trans, integer *n, integer *kl, integer *ku, double *ab,
                         integer *ldab, double *afb, integer *ldafb, integer *ipiv, integer *cmode,
                         double *c__, integer *info, double *work, integer *iwork);
int dla_gbrfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n, integer *kl,
                              integer *ku, integer *nrhs, double *ab, integer *ldab, double *afb,
                              integer *ldafb, integer *ipiv, logical *colequ, double *c__,
                              double *b, integer *ldb, double *y, integer *ldy, double *berr_out__,
                              integer *n_norms__, double *err_bnds_norm__, double *err_bnds_comp__,
                              double *res, double *ayb, double *dy, double *y_tail__, double *rcond,
                              integer *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, integer *info);
double dla_gbrpvgrw_check(integer *n, integer *kl, integer *ku, integer *ncols, double *ab,
                          integer *ldab, double *afb, integer *ldafb);
int dla_geamv_check(integer *trans, integer *m, integer *n, double *alpha, double *a, integer *lda,
                    double *x, integer *incx, double *beta, double *y, integer *incy);
double dla_gercond_check(char *trans, integer *n, double *a, integer *lda, double *af,
                         integer *ldaf, integer *ipiv, integer *cmode, double *c__, integer *info,
                         double *work, integer *iwork);
int dla_gerfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n,
                              integer *nrhs, double *a, integer *lda, double *af, integer *ldaf,
                              integer *ipiv, logical *colequ, double *c__, double *b, integer *ldb,
                              double *y, integer *ldy, double *berr_out__, integer *n_norms__,
                              double *errs_n__, double *errs_c__, double *res, double *ayb,
                              double *dy, double *y_tail__, double *rcond, integer *ithresh,
                              double *rthresh, double *dz_ub__, logical *ignore_cwise__,
                              integer *info);
double dla_gerpvgrw_check(integer *n, integer *ncols, double *a, integer *lda, double *af,
                          integer *ldaf);
int dla_lin_berr_check(integer *n, integer *nz, integer *nrhs, double *res, double *ayb,
                       double *berr);
double dla_porcond_check(char *uplo, integer *n, double *a, integer *lda, double *af, integer *ldaf,
                         integer *cmode, double *c__, integer *info, double *work, integer *iwork);
int dla_porfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              double *a, integer *lda, double *af, integer *ldaf, logical *colequ,
                              double *c__, double *b, integer *ldb, double *y, integer *ldy,
                              double *berr_out__, integer *n_norms__, double *err_bnds_norm__,
                              double *err_bnds_comp__, double *res, double *ayb, double *dy,
                              double *y_tail__, double *rcond, integer *ithresh, double *rthresh,
                              double *dz_ub__, logical *ignore_cwise__, integer *info);
double dla_porpvgrw_check(char *uplo, integer *ncols, double *a, integer *lda, double *af,
                          integer *ldaf, double *work);
int dla_syamv_check(integer *uplo, integer *n, double *alpha, double *a, integer *lda, double *x,
                    integer *incx, double *beta, double *y, integer *incy);
double dla_syrcond_check(char *uplo, integer *n, double *a, integer *lda, double *af, integer *ldaf,
                         integer *ipiv, integer *cmode, double *c__, integer *info, double *work,
                         integer *iwork);
int dla_syrfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              double *a, integer *lda, double *af, integer *ldaf, integer *ipiv,
                              logical *colequ, double *c__, double *b, integer *ldb, double *y,
                              integer *ldy, double *berr_out__, integer *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, double *res,
                              double *ayb, double *dy, double *y_tail__, double *rcond,
                              integer *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, integer *info);
double dla_syrpvgrw_check(char *uplo, integer *n, integer *info, double *a, integer *lda,
                          double *af, integer *ldaf, integer *ipiv, double *work);
int dla_wwaddw_check(integer *n, double *x, double *y, double *w);
int dlabad_check(double *small_, double *large);
int dlabrd_check(integer *m, integer *n, integer *nb, double *a, integer *lda, double *d__,
                 double *e, double *tauq, double *taup, double *x, integer *ldx, double *y,
                 integer *ldy);
int dlacn2_check(integer *n, double *v, double *x, integer *isgn, double *est, integer *kase,
                 integer *isave);
int dlacon_check(integer *n, double *v, double *x, integer *isgn, double *est, integer *kase);
int dlacpy_check(char *uplo, integer *m, integer *n, double *a, integer *lda, double *b,
                 integer *ldb);
int dladiv_check(double *a, double *b, double *c__, double *d__, double *p, double *q);
int dlae2_check(double *a, double *b, double *c__, double *rt1, double *rt2);
int dlaebz_check(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp,
                 integer *nbmin, double *abstol, double *reltol, double *pivmin, double *d__,
                 double *e, double *e2, integer *nval, double *ab, double *c__, integer *mout,
                 integer *nab, double *work, integer *iwork, integer *info);
int dlaed0_check(integer *icompq, integer *qsiz, integer *n, double *d__, double *e, double *q,
                 integer *ldq, double *qstore, integer *ldqs, double *work, integer *iwork,
                 integer *info);
int dlaed1_check(integer *n, double *d__, double *q, integer *ldq, integer *indxq, double *rho,
                 integer *cutpnt, double *work, integer *iwork, integer *info);
int dlaed2_check(integer *k, integer *n, integer *n1, double *d__, double *q, integer *ldq,
                 integer *indxq, double *rho, double *z__, double *dlamda, double *w, double *q2,
                 integer *indx, integer *indxc, integer *indxp, integer *coltyp, integer *info);
int dlaed3_check(integer *k, integer *n, integer *n1, double *d__, double *q, integer *ldq,
                 double *rho, double *dlamda, double *q2, integer *indx, integer *ctot, double *w,
                 double *s, integer *info);
int dlaed4_check(integer *n, integer *i__, double *d__, double *z__, double *delta, double *rho,
                 double *dlam, integer *info);
int dlaed5_check(integer *i__, double *d__, double *z__, double *delta, double *rho, double *dlam);
int dlaed6_check(integer *kniter, logical *orgati, double *rho, double *d__, double *z__,
                 double *finit, double *tau, integer *info);
int dlaed7_check(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl,
                 integer *curpbm, double *d__, double *q, integer *ldq, integer *indxq, double *rho,
                 integer *cutpnt, double *qstore, integer *qptr, integer *prmptr, integer *perm,
                 integer *givptr, integer *givcol, double *givnum, double *work, integer *iwork,
                 integer *info);
int dlaed8_check(integer *icompq, integer *k, integer *n, integer *qsiz, double *d__, double *q,
                 integer *ldq, integer *indxq, double *rho, integer *cutpnt, double *z__,
                 double *dlamda, double *q2, integer *ldq2, double *w, integer *perm,
                 integer *givptr, integer *givcol, double *givnum, integer *indxp, integer *indx,
                 integer *info);
int dlaed9_check(integer *k, integer *kstart, integer *kstop, integer *n, double *d__, double *q,
                 integer *ldq, double *rho, double *dlamda, double *w, double *s, integer *lds,
                 integer *info);
int dlaeda_check(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr,
                 integer *perm, integer *givptr, integer *givcol, double *givnum, double *q,
                 integer *qptr, double *z__, double *ztemp, integer *info);
int dlaein_check(logical *rightv, logical *noinit, integer *n, double *h__, integer *ldh,
                 double *wr, double *wi, double *vr, double *vi, double *b, integer *ldb,
                 double *work, double *eps3, double *smlnum, double *bignum, integer *info);
int dlaev2_check(double *a, double *b, double *c__, double *rt1, double *rt2, double *cs1,
                 double *sn1);
int dlaexc_check(logical *wantq, integer *n, double *t, integer *ldt, double *q, integer *ldq,
                 integer *j1, integer *n1, integer *n2, double *work, integer *info);
int dlag2_check(double *a, integer *lda, double *b, integer *ldb, double *safmin, double *scale1,
                double *scale2, double *wr1, double *wr2, double *wi);
int dlag2s_check(integer *m, integer *n, double *a, integer *lda, float *sa, integer *ldsa,
                 integer *info);
int dlags2_check(logical *upper, double *a1, double *a2, double *a3, double *b1, double *b2,
                 double *b3, double *csu, double *snu, double *csv, double *snv, double *csq,
                 double *snq);
int dlagtf_check(integer *n, double *a, double *lambda, double *b, double *c__, double *tol,
                 double *d__, integer *in, integer *info);
int dlagtm_check(char *trans, integer *n, integer *nrhs, double *alpha, double *dl, double *d__,
                 double *du, double *x, integer *ldx, double *beta, double *b, integer *ldb);
int dlagts_check(integer *job, integer *n, double *a, double *b, double *c__, double *d__,
                 integer *in, double *y, double *tol, integer *info);
int dlagv2_check(double *a, integer *lda, double *b, integer *ldb, double *alphar, double *alphai,
                 double *beta, double *csl, double *snl, double *csr, double *snr);
int dlahqr_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 double *h__, integer *ldh, double *wr, double *wi, integer *iloz, integer *ihiz,
                 double *z__, integer *ldz, integer *info);
int dlahr2_check(integer *n, integer *k, integer *nb, double *a, integer *lda, double *tau,
                 double *t, integer *ldt, double *y, integer *ldy);
int dlahrd_check(integer *n, integer *k, integer *nb, double *a, integer *lda, double *tau,
                 double *t, integer *ldt, double *y, integer *ldy);
int dlaic1_check(integer *job, integer *j, double *x, double *sest, double *w, double *gamma,
                 double *sestpr, double *s, double *c__);
logical dlaisnan_check(double *din1, double *din2);
int dlaln2_check(logical *ltrans, integer *na, integer *nw, double *smin, double *ca, double *a,
                 integer *lda, double *d1, double *d2, double *b, integer *ldb, double *wr,
                 double *wi, double *x, integer *ldx, double *scale, double *xnorm, integer *info);
int dlals0_check(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, double *b,
                 integer *ldb, double *bx, integer *ldbx, integer *perm, integer *givptr,
                 integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum, double *poles,
                 double *difl, double *difr, double *z__, integer *k, double *c__, double *s,
                 double *work, integer *info);
int dlalsa_check(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, double *b,
                 integer *ldb, double *bx, integer *ldbx, double *u, integer *ldu, double *vt,
                 integer *k, double *difl, double *difr, double *z__, double *poles,
                 integer *givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum,
                 double *c__, double *s, double *work, integer *iwork, integer *info);
int dlalsd_check(char *uplo, integer *smlsiz, integer *n, integer *nrhs, double *d__, double *e,
                 double *b, integer *ldb, double *rcond, integer *rank, double *work,
                 integer *iwork, integer *info);
int dlamrg_check(integer *n1, integer *n2, double *a, integer *dtrd1, integer *dtrd2,
                 integer *index);
int dlaneg_check(integer *n, double *d__, double *lld, double *sigma, double *pivmin, integer *r__);
double dlangb_check(char *norm, integer *n, integer *kl, integer *ku, double *ab, integer *ldab,
                    double *work);
double dlange_check(char *norm, integer *m, integer *n, double *a, integer *lda, double *work);
double dlangt_check(char *norm, integer *n, double *dl, double *d__, double *du);
double dlanhs_check(char *norm, integer *n, double *a, integer *lda, double *work);
double dlansb_check(char *norm, char *uplo, integer *n, integer *k, double *ab, integer *ldab,
                    double *work);
double dlansf_check(char *norm, char *transr, char *uplo, integer *n, double *a, double *work);
double dlansp_check(char *norm, char *uplo, integer *n, double *ap, double *work);
double dlanst_check(char *norm, integer *n, double *d__, double *e);
double dlansy_check(char *norm, char *uplo, integer *n, double *a, integer *lda, double *work);
double dlantb_check(char *norm, char *uplo, char *diag, integer *n, integer *k, double *ab,
                    integer *ldab, double *work);
double dlantp_check(char *norm, char *uplo, char *diag, integer *n, double *ap, double *work);
double dlantr_check(char *norm, char *uplo, char *diag, integer *m, integer *n, double *a,
                    integer *lda, double *work);
int dlanv2_check(double *a, double *b, double *c__, double *d__, double *rt1r, double *rt1i,
                 double *rt2r, double *rt2i, double *cs, double *sn);
int dlapll_check(integer *n, double *x, integer *incx, double *y, integer *incy, double *ssmin);
int dlapmr_check(logical *forwrd, integer *m, integer *n, double *x, integer *ldx, integer *k);
int dlapmt_check(logical *forwrd, integer *m, integer *n, double *x, integer *ldx, integer *k);
double dlapy2_check(double *x, double *y);
double dlapy3_check(double *x, double *y, double *z__);
int dlaqgb_check(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 char *equed);
int dlaqge_check(integer *m, integer *n, double *a, integer *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, char *equed);
int dlaqp2_check(integer *m, integer *n, integer *offset, double *a, integer *lda, integer *jpvt,
                 double *tau, double *vn1, double *vn2, double *work);
int dlaqps_check(integer *m, integer *n, integer *offset, integer *nb, integer *kb, double *a,
                 integer *lda, integer *jpvt, double *tau, double *vn1, double *vn2, double *auxv,
                 double *f, integer *ldf);
int dlaqr0_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 double *h__, integer *ldh, double *wr, double *wi, integer *iloz, integer *ihiz,
                 double *z__, integer *ldz, double *work, integer *lwork, integer *info);
int dlaqr1_check(integer *n, double *h__, integer *ldh, double *sr1, double *si1, double *sr2,
                 double *si2, double *v);
int dlaqr2_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, double *h__, integer *ldh, integer *iloz, integer *ihiz, double *z__,
                 integer *ldz, integer *ns, integer *nd, double *sr, double *si, double *v,
                 integer *ldv, integer *nh, double *t, integer *ldt, integer *nv, double *wv,
                 integer *ldwv, double *work, integer *lwork);
int dlaqr3_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, double *h__, integer *ldh, integer *iloz, integer *ihiz, double *z__,
                 integer *ldz, integer *ns, integer *nd, double *sr, double *si, double *v,
                 integer *ldv, integer *nh, double *t, integer *ldt, integer *nv, double *wv,
                 integer *ldwv, double *work, integer *lwork);
int dlaqr4_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 double *h__, integer *ldh, double *wr, double *wi, integer *iloz, integer *ihiz,
                 double *z__, integer *ldz, double *work, integer *lwork, integer *info);
int dlaqr5_check(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
                 integer *kbot, integer *nshfts, double *sr, double *si, double *h__, integer *ldh,
                 integer *iloz, integer *ihiz, double *z__, integer *ldz, double *v, integer *ldv,
                 double *u, integer *ldu, integer *nv, double *wv, integer *ldwv, integer *nh,
                 double *wh, integer *ldwh);
int dlaqsb_check(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, double *s,
                 double *scond, double *amax, char *equed);
int dlaqsp_check(char *uplo, integer *n, double *ap, double *s, double *scond, double *amax,
                 char *equed);
int dlaqsy_check(char *uplo, integer *n, double *a, integer *lda, double *s, double *scond,
                 double *amax, char *equed);
int dlaqtr_check(logical *ltran, logical *lfloat, integer *n, double *t, integer *ldt, double *b,
                 double *w, double *scale, double *x, double *work, integer *info);
int dlar1v_check(integer *n, integer *b1, integer *bn, double *lambda, double *d__, double *l,
                 double *ld, double *lld, double *pivmin, double *gaptol, double *z__,
                 logical *wantnc, integer *negcnt, double *ztz, double *mingma, integer *r__,
                 integer *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);
int dlar2v_check(integer *n, double *x, double *y, double *z__, integer *incx, double *c__,
                 double *s, integer *incc);
int dlarf_check(char *side, integer *m, integer *n, double *v, integer *incv, double *tau,
                double *c__, integer *ldc, double *work);
int dlarfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, double *v, integer *ldv, double *t, integer *ldt, double *c__,
                 integer *ldc, double *work, integer *ldwork);
int dlarfg_check(integer *n, double *alpha, double *x, integer *incx, double *tau);
int dlarfgp_check(integer *n, double *alpha, double *x, integer *incx, double *tau);
int dlarft_check(char *direct, char *storev, integer *n, integer *k, double *v, integer *ldv,
                 double *tau, double *t, integer *ldt);
int dlarfx_check(char *side, integer *m, integer *n, double *v, double *tau, double *c__,
                 integer *ldc, double *work);
int dlargv_check(integer *n, double *x, integer *incx, double *y, integer *incy, double *c__,
                 integer *incc);
int dlarnv_check(integer *idist, integer *iseed, integer *n, double *x);
int dlarra_check(integer *n, double *d__, double *e, double *e2, double *spltol, double *tnrm,
                 integer *nsplit, integer *isplit, integer *info);
int dlarrb_check(integer *n, double *d__, double *lld, integer *ifirst, integer *ilast,
                 double *rtol1, double *rtol2, integer *offset, double *w, double *wgap,
                 double *werr, double *work, integer *iwork, double *pivmin, double *spdiam,
                 integer *twist, integer *info);
int dlarrc_check(char *jobt, integer *n, double *vl, double *vu, double *d__, double *e,
                 double *pivmin, integer *eigcnt, integer *lcnt, integer *rcnt, integer *info);
int dlarrd_check(char *range, char *order, integer *n, double *vl, double *vu, integer *il,
                 integer *iu, double *gers, double *reltol, double *d__, double *e, double *e2,
                 double *pivmin, integer *nsplit, integer *isplit, integer *m, double *w,
                 double *werr, double *wl, double *wu, integer *iblock, integer *indexw,
                 double *work, integer *iwork, integer *info);
int dlarre_check(char *range, integer *n, double *vl, double *vu, integer *il, integer *iu,
                 double *d__, double *e, double *e2, double *rtol1, double *rtol2, double *spltol,
                 integer *nsplit, integer *isplit, integer *m, double *w, double *werr,
                 double *wgap, integer *iblock, integer *indexw, double *gers, double *pivmin,
                 double *work, integer *iwork, integer *info);
int dlarrf_check(integer *n, double *d__, double *l, double *ld, integer *clstrt, integer *clend,
                 double *w, double *wgap, double *werr, double *spdiam, double *clgapl,
                 double *clgapr, double *pivmin, double *sigma, double *dplus, double *lplus,
                 double *work, integer *info);
int dlarrj_check(integer *n, double *d__, double *e2, integer *ifirst, integer *ilast, double *rtol,
                 integer *offset, double *w, double *werr, double *work, integer *iwork,
                 double *pivmin, double *spdiam, integer *info);
int dlarrk_check(integer *n, integer *iw, double *gl, double *gu, double *d__, double *e2,
                 double *pivmin, double *reltol, double *w, double *werr, integer *info);
int dlarrr_check(integer *n, double *d__, double *e, integer *info);
int dlarrv_check(integer *n, double *vl, double *vu, double *d__, double *l, double *pivmin,
                 integer *isplit, integer *m, integer *dol, integer *dou, double *minrgp,
                 double *rtol1, double *rtol2, double *w, double *werr, double *wgap,
                 integer *iblock, integer *indexw, double *gers, double *z__, integer *ldz,
                 integer *isuppz, double *work, integer *iwork, integer *info);
int dlarscl2_check(integer *m, integer *n, double *d__, double *x, integer *ldx);
int dlartg_check(double *f, double *g, double *cs, double *sn, double *r__);
int dlartgp_check(double *f, double *g, double *cs, double *sn, double *r__);
int dlartgs_check(double *x, double *y, double *sigma, double *cs, double *sn);
int dlartv_check(integer *n, double *x, integer *incx, double *y, integer *incy, double *c__,
                 double *s, integer *incc);
int dlaruv_check(integer *iseed, integer *n, double *x);
int dlarz_check(char *side, integer *m, integer *n, integer *l, double *v, integer *incv,
                double *tau, double *c__, integer *ldc, double *work);
int dlarzb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, double *v, integer *ldv, double *t, integer *ldt,
                 double *c__, integer *ldc, double *work, integer *ldwork);
int dlarzt_check(char *direct, char *storev, integer *n, integer *k, double *v, integer *ldv,
                 double *tau, double *t, integer *ldt);
int dlas2_check(double *f, double *g, double *h__, double *ssmin, double *ssmax);
int dlascl_check(char *type__, integer *kl, integer *ku, double *cfrom, double *cto, integer *m,
                 integer *n, double *a, integer *lda, integer *info);
int dlascl2_check(integer *m, integer *n, double *d__, double *x, integer *ldx);
int dlasd0_check(integer *n, integer *sqre, double *d__, double *e, double *u, integer *ldu,
                 double *vt, integer *ldvt, integer *smlsiz, integer *iwork, double *work,
                 integer *info);
int dlasd1_check(integer *nl, integer *nr, integer *sqre, double *d__, double *alpha, double *beta,
                 double *u, integer *ldu, double *vt, integer *ldvt, integer *idxq, integer *iwork,
                 double *work, integer *info);
int dlasd2_check(integer *nl, integer *nr, integer *sqre, integer *k, double *d__, double *z__,
                 double *alpha, double *beta, double *u, integer *ldu, double *vt, integer *ldvt,
                 double *dsigma, double *u2, integer *ldu2, double *vt2, integer *ldvt2,
                 integer *idxp, integer *idx, integer *idxc, integer *idxq, integer *coltyp,
                 integer *info);
int dlasd3_check(integer *nl, integer *nr, integer *sqre, integer *k, double *d__, double *q,
                 integer *ldq, double *dsigma, double *u, integer *ldu, double *u2, integer *ldu2,
                 double *vt, integer *ldvt, double *vt2, integer *ldvt2, integer *idxc,
                 integer *ctot, double *z__, integer *info);
int dlasd4_check(integer *n, integer *i__, double *d__, double *z__, double *delta, double *rho,
                 double *sigma, double *work, integer *info);
int dlasd5_check(integer *i__, double *d__, double *z__, double *delta, double *rho, double *dsigma,
                 double *work);
int dlasd6_check(integer *icompq, integer *nl, integer *nr, integer *sqre, double *d__, double *vf,
                 double *vl, double *alpha, double *beta, integer *idxq, integer *perm,
                 integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum,
                 double *poles, double *difl, double *difr, double *z__, integer *k, double *c__,
                 double *s, double *work, integer *iwork, integer *info);
int dlasd7_check(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, double *d__,
                 double *z__, double *zw, double *vf, double *vfw, double *vl, double *vlw,
                 double *alpha, double *beta, double *dsigma, integer *idx, integer *idxp,
                 integer *idxq, integer *perm, integer *givptr, integer *givcol, integer *ldgcol,
                 double *givnum, integer *ldgnum, double *c__, double *s, integer *info);
int dlasd8_check(integer *icompq, integer *k, double *d__, double *z__, double *vf, double *vl,
                 double *difl, double *difr, integer *lddifr, double *dsigma, double *work,
                 integer *info);
int dlasda_check(integer *icompq, integer *smlsiz, integer *n, integer *sqre, double *d__,
                 double *e, double *u, integer *ldu, double *vt, integer *k, double *difl,
                 double *difr, double *z__, double *poles, integer *givptr, integer *givcol,
                 integer *ldgcol, integer *perm, double *givnum, double *c__, double *s,
                 double *work, integer *iwork, integer *info);
int dlasdq_check(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc,
                 double *d__, double *e, double *vt, integer *ldvt, double *u, integer *ldu,
                 double *c__, integer *ldc, double *work, integer *info);
int dlasdt_check(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml,
                 integer *ndimr, integer *msub);
int dlaset_check(char *uplo, integer *m, integer *n, double *alpha, double *beta, double *a,
                 integer *lda);
int dlasq1_check(integer *n, double *d__, double *e, double *work, integer *info);
int dlasq2_check(integer *n, double *z__, integer *info);
int dlasq3_check(integer *i0, integer *n0, double *z__, integer *pp, double *dmin__, double *sigma,
                 double *desig, double *qmax, integer *nfail, integer *iter, integer *ndiv,
                 logical *ieee, integer *ttype, double *dmin1, double *dmin2, double *dn,
                 double *dn1, double *dn2, double *g, double *tau);
int dlasq4_check(integer *i0, integer *n0, double *z__, integer *pp, integer *n0in, double *dmin__,
                 double *dmin1, double *dmin2, double *dn, double *dn1, double *dn2, double *tau,
                 integer *ttype, double *g);
int dlasq5_check(integer *i0, integer *n0, double *z__, integer *pp, double *tau, double *sigma,
                 double *dmin__, double *dmin1, double *dmin2, double *dn, double *dnm1,
                 double *dnm2, logical *ieee, double *eps);
int dlasq6_check(integer *i0, integer *n0, double *z__, integer *pp, double *dmin__, double *dmin1,
                 double *dmin2, double *dn, double *dnm1, double *dnm2);
int dlasr_check(char *side, char *pivot, char *direct, integer *m, integer *n, double *c__,
                double *s, double *a, integer *lda);
int dlasrt_check(char *id, integer *n, double *d__, integer *info);
int dlassq_check(integer *n, double *x, integer *incx, double *scale, double *sumsq);
int dlasv2_check(double *f, double *g, double *h__, double *ssmin, double *ssmax, double *snr,
                 double *csr, double *snl, double *csl);
int dlaswp_check(integer *n, double *a, integer *lda, integer *k1, integer *k2, integer *ipiv,
                 integer *incx);
int dlasy2_check(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2,
                 double *tl, integer *ldtl, double *tr, integer *ldtr, double *b, integer *ldb,
                 double *scale, double *x, integer *ldx, double *xnorm, integer *info);
int dlasyf_check(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda,
                 integer *ipiv, double *w, integer *ldw, integer *info);
int dlasyf_rook_check(char *uplo, integer *n, integer *nb, integer *kb, double *a, integer *lda,
                      integer *ipiv, double *w, integer *ldw, integer *info);
int dlat2s_check(char *uplo, integer *n, double *a, integer *lda, float *sa, integer *ldsa,
                 integer *info);
int dlatbs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd,
                 double *ab, integer *ldab, double *x, double *scale, double *cnorm, integer *info);
int dlatdf_check(integer *ijob, integer *n, double *z__, integer *ldz, double *rhs, double *rdsum,
                 double *rdscal, integer *ipiv, integer *jpiv);
int dlatps_check(char *uplo, char *trans, char *diag, char *normin, integer *n, double *ap,
                 double *x, double *scale, double *cnorm, integer *info);
int dlatrd_check(char *uplo, integer *n, integer *nb, double *a, integer *lda, double *e,
                 double *tau, double *w, integer *ldw);
int dlatrs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, double *a,
                 integer *lda, double *x, double *scale, double *cnorm, integer *info);
int dlatrz_check(integer *m, integer *n, integer *l, double *a, integer *lda, double *tau,
                 double *work);
int dlatzm_check(char *side, integer *m, integer *n, double *v, integer *incv, double *tau,
                 double *c1, double *c2, integer *ldc, double *work);
int dlauu2_check(char *uplo, integer *n, double *a, integer *lda, integer *info);
int dlauum_check(char *uplo, integer *n, double *a, integer *lda, integer *info);
int dopgtr_check(char *uplo, integer *n, double *ap, double *tau, double *q, integer *ldq,
                 double *work, integer *info);
int dopmtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, double *ap,
                 double *tau, double *c__, integer *ldc, double *work, integer *info);
int dorbdb_check(char *trans, char *signs, integer *m, integer *p, integer *q, double *x11,
                 integer *ldx11, double *x12, integer *ldx12, double *x21, integer *ldx21,
                 double *x22, integer *ldx22, double *theta, double *phi, double *taup1,
                 double *taup2, double *tauq1, double *tauq2, double *work, integer *lwork,
                 integer *info);
int dorbdb1_check(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21,
                  integer *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *work, integer *lwork, integer *info);
int dorbdb2_check(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21,
                  integer *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *work, integer *lwork, integer *info);
int dorbdb3_check(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21,
                  integer *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *work, integer *lwork, integer *info);
int dorbdb4_check(integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x21,
                  integer *ldx21, double *theta, double *phi, double *taup1, double *taup2,
                  double *tauq1, double *phantom, double *work, integer *lwork, integer *info);
int dorbdb5_check(integer *m1, integer *m2, integer *n, double *x1, integer *incx1, double *x2,
                  integer *incx2, double *q1, integer *ldq1, double *q2, integer *ldq2,
                  double *work, integer *lwork, integer *info);
int dorbdb6_check(integer *m1, integer *m2, integer *n, double *x1, integer *incx1, double *x2,
                  integer *incx2, double *q1, integer *ldq1, double *q2, integer *ldq2,
                  double *work, integer *lwork, integer *info);
int dorcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 integer *m, integer *p, integer *q, double *x11, integer *ldx11, double *x12,
                 integer *ldx12, double *x21, integer *ldx21, double *x22, integer *ldx22,
                 double *theta, double *u1, integer *ldu1, double *u2, integer *ldu2, double *v1t,
                 integer *ldv1t, double *v2t, integer *ldv2t, double *work, integer *lwork,
                 integer *iwork, integer *info);
int dorcsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, integer *m, integer *p, integer *q,
                     double *x11, integer *ldx11, double *x21, integer *ldx21, double *theta,
                     double *u1, integer *ldu1, double *u2, integer *ldu2, double *v1t,
                     integer *ldv1t, double *work, integer *lwork, integer *iwork, integer *info);
int dorg2l_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *info);
int dorg2r_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *info);
int dorgbr_check(char *vect, integer *m, integer *n, integer *k, double *a, integer *lda,
                 double *tau, double *work, integer *lwork, integer *info);
int dorghr_check(integer *n, integer *ilo, integer *ihi, double *a, integer *lda, double *tau,
                 double *work, integer *lwork, integer *info);
int dorgl2_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *info);
int dorglq_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *lwork, integer *info);
int dorgql_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *lwork, integer *info);
int dorgqr_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *lwork, integer *info);
int dorgr2_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *info);
int dorgrq_check(integer *m, integer *n, integer *k, double *a, integer *lda, double *tau,
                 double *work, integer *lwork, integer *info);
int dorgtr_check(char *uplo, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *lwork, integer *info);
int dorm2l_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *info);
int dorm2r_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *info);
int dormbr_check(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *lwork,
                 integer *info);
int dormhr_check(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi,
                 double *a, integer *lda, double *tau, double *c__, integer *ldc, double *work,
                 integer *lwork, integer *info);
int dorml2_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *info);
int dormlq_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *lwork,
                 integer *info);
int dormql_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *lwork,
                 integer *info);
int dormqr_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *lwork,
                 integer *info);
int dormr2_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *info);
int dormr3_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *info);
int dormrq_check(char *side, char *trans, integer *m, integer *n, integer *k, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *lwork,
                 integer *info);
int dormrz_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *lwork,
                 integer *info);
int dormtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, double *a,
                 integer *lda, double *tau, double *c__, integer *ldc, double *work, integer *lwork,
                 integer *info);
int dpbcon_check(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, double *anorm,
                 double *rcond, double *work, integer *iwork, integer *info);
int dpbequ_check(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, double *s,
                 double *scond, double *amax, integer *info);
int dpbrfs_check(char *uplo, integer *n, integer *kd, integer *nrhs, double *ab, integer *ldab,
                 double *afb, integer *ldafb, double *b, integer *ldb, double *x, integer *ldx,
                 double *ferr, double *berr, double *work, integer *iwork, integer *info);
int dpbstf_check(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, integer *info);
int dpbsv_check(char *uplo, integer *n, integer *kd, integer *nrhs, double *ab, integer *ldab,
                double *b, integer *ldb, integer *info);
int dpbsvx_check(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, double *ab,
                 integer *ldab, double *afb, integer *ldafb, char *equed, double *s, double *b,
                 integer *ldb, double *x, integer *ldx, double *rcond, double *ferr, double *berr,
                 double *work, integer *iwork, integer *info);
int dpbtf2_check(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, integer *info);
int dpbtrf_check(char *uplo, integer *n, integer *kd, double *ab, integer *ldab, integer *info);
int dpbtrs_check(char *uplo, integer *n, integer *kd, integer *nrhs, double *ab, integer *ldab,
                 double *b, integer *ldb, integer *info);
int dpftrf_check(char *transr, char *uplo, integer *n, double *a, integer *info);
int dpftri_check(char *transr, char *uplo, integer *n, double *a, integer *info);
int dpftrs_check(char *transr, char *uplo, integer *n, integer *nrhs, double *a, double *b,
                 integer *ldb, integer *info);
int dpocon_check(char *uplo, integer *n, double *a, integer *lda, double *anorm, double *rcond,
                 double *work, integer *iwork, integer *info);
int dpoequ_check(integer *n, double *a, integer *lda, double *s, double *scond, double *amax,
                 integer *info);
int dpoequb_check(integer *n, double *a, integer *lda, double *s, double *scond, double *amax,
                  integer *info);
int dporfs_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *af,
                 integer *ldaf, double *b, integer *ldb, double *x, integer *ldx, double *ferr,
                 double *berr, double *work, integer *iwork, integer *info);
int dporfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, double *a, integer *lda,
                  double *af, integer *ldaf, double *s, double *b, integer *ldb, double *x,
                  integer *ldx, double *rcond, double *berr, integer *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams,
                  double *params, double *work, integer *iwork, integer *info);
int dposv_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *b,
                integer *ldb, integer *info);
int dposvx_check(char *fact, char *uplo, integer *n, integer *nrhs, double *a, integer *lda,
                 double *af, integer *ldaf, char *equed, double *s, double *b, integer *ldb,
                 double *x, integer *ldx, double *rcond, double *ferr, double *berr, double *work,
                 integer *iwork, integer *info);
int dposvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, double *a, integer *lda,
                  double *af, integer *ldaf, char *equed, double *s, double *b, integer *ldb,
                  double *x, integer *ldx, double *rcond, double *rpvgrw, double *berr,
                  integer *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  integer *nparams, double *params, double *work, integer *iwork, integer *info);
int dpotf2_check(char *uplo, integer *n, double *a, integer *lda, integer *info);
int dpotrf_check(char *uplo, integer *n, double *a, integer *lda, integer *info);
int dpotri_check(char *uplo, integer *n, double *a, integer *lda, integer *info);
int dpotrs_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *b,
                 integer *ldb, integer *info);
int dppcon_check(char *uplo, integer *n, double *ap, double *anorm, double *rcond, double *work,
                 integer *iwork, integer *info);
int dppequ_check(char *uplo, integer *n, double *ap, double *s, double *scond, double *amax,
                 integer *info);
int dpprfs_check(char *uplo, integer *n, integer *nrhs, double *ap, double *afp, double *b,
                 integer *ldb, double *x, integer *ldx, double *ferr, double *berr, double *work,
                 integer *iwork, integer *info);
int dppsv_check(char *uplo, integer *n, integer *nrhs, double *ap, double *b, integer *ldb,
                integer *info);
int dppsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, double *ap, double *afp,
                 char *equed, double *s, double *b, integer *ldb, double *x, integer *ldx,
                 double *rcond, double *ferr, double *berr, double *work, integer *iwork,
                 integer *info);
int dpptrf_check(char *uplo, integer *n, double *ap, integer *info);
int dpptri_check(char *uplo, integer *n, double *ap, integer *info);
int dpptrs_check(char *uplo, integer *n, integer *nrhs, double *ap, double *b, integer *ldb,
                 integer *info);
int dpstf2_check(char *uplo, integer *n, double *a, integer *lda, integer *piv, integer *rank,
                 double *tol, double *work, integer *info);
int dpstrf_check(char *uplo, integer *n, double *a, integer *lda, integer *piv, integer *rank,
                 double *tol, double *work, integer *info);
int dptcon_check(integer *n, double *d__, double *e, double *anorm, double *rcond, double *work,
                 integer *info);
int dpteqr_check(char *compz, integer *n, double *d__, double *e, double *z__, integer *ldz,
                 double *work, integer *info);
int dptrfs_check(integer *n, integer *nrhs, double *d__, double *e, double *df, double *ef,
                 double *b, integer *ldb, double *x, integer *ldx, double *ferr, double *berr,
                 double *work, integer *info);
int dptsv_check(integer *n, integer *nrhs, double *d__, double *e, double *b, integer *ldb,
                integer *info);
int dptsvx_check(char *fact, integer *n, integer *nrhs, double *d__, double *e, double *df,
                 double *ef, double *b, integer *ldb, double *x, integer *ldx, double *rcond,
                 double *ferr, double *berr, double *work, integer *info);
int dpttrf_check(integer *n, double *d__, double *e, integer *info);
int dpttrs_check(integer *n, integer *nrhs, double *d__, double *e, double *b, integer *ldb,
                 integer *info);
int dptts2_check(integer *n, integer *nrhs, double *d__, double *e, double *b, integer *ldb);
int drscl_check(integer *n, double *sa, double *sx, integer *incx);
int dsbev_check(char *jobz, char *uplo, integer *n, integer *kd, double *ab, integer *ldab,
                double *w, double *z__, integer *ldz, double *work, integer *info);
int dsbevd_check(char *jobz, char *uplo, integer *n, integer *kd, double *ab, integer *ldab,
                 double *w, double *z__, integer *ldz, double *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int dsbevx_check(char *jobz, char *range, char *uplo, integer *n, integer *kd, double *ab,
                 integer *ldab, double *q, integer *ldq, double *vl, double *vu, integer *il,
                 integer *iu, double *abstol, integer *m, double *w, double *z__, integer *ldz,
                 double *work, integer *iwork, integer *ifail, integer *info);
int dsbgst_check(char *vect, char *uplo, integer *n, integer *ka, integer *kb, double *ab,
                 integer *ldab, double *bb, integer *ldbb, double *x, integer *ldx, double *work,
                 integer *info);
int dsbgv_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, double *ab,
                integer *ldab, double *bb, integer *ldbb, double *w, double *z__, integer *ldz,
                double *work, integer *info);
int dsbgvd_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, double *ab,
                 integer *ldab, double *bb, integer *ldbb, double *w, double *z__, integer *ldz,
                 double *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int dsbgvx_check(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb,
                 double *ab, integer *ldab, double *bb, integer *ldbb, double *q, integer *ldq,
                 double *vl, double *vu, integer *il, integer *iu, double *abstol, integer *m,
                 double *w, double *z__, integer *ldz, double *work, integer *iwork, integer *ifail,
                 integer *info);
int dsbtrd_check(char *vect, char *uplo, integer *n, integer *kd, double *ab, integer *ldab,
                 double *d__, double *e, double *q, integer *ldq, double *work, integer *info);
int dsfrk_check(char *transr, char *uplo, char *trans, integer *n, integer *k, double *alpha,
                double *a, integer *lda, double *beta, double *c__);
int dsgesv_check(integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv, double *b,
                 integer *ldb, double *x, integer *ldx, double *work, float *swork, integer *iter,
                 integer *info);
int dspcon_check(char *uplo, integer *n, double *ap, integer *ipiv, double *anorm, double *rcond,
                 double *work, integer *iwork, integer *info);
int dspev_check(char *jobz, char *uplo, integer *n, double *ap, double *w, double *z__,
                integer *ldz, double *work, integer *info);
int dspevd_check(char *jobz, char *uplo, integer *n, double *ap, double *w, double *z__,
                 integer *ldz, double *work, integer *lwork, integer *iwork, integer *liwork,
                 integer *info);
int dspevx_check(char *jobz, char *range, char *uplo, integer *n, double *ap, double *vl,
                 double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w,
                 double *z__, integer *ldz, double *work, integer *iwork, integer *ifail,
                 integer *info);
int dspgst_check(integer *itype, char *uplo, integer *n, double *ap, double *bp, integer *info);
int dspgv_check(integer *itype, char *jobz, char *uplo, integer *n, double *ap, double *bp,
                double *w, double *z__, integer *ldz, double *work, integer *info);
int dspgvd_check(integer *itype, char *jobz, char *uplo, integer *n, double *ap, double *bp,
                 double *w, double *z__, integer *ldz, double *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int dspgvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, double *ap,
                 double *bp, double *vl, double *vu, integer *il, integer *iu, double *abstol,
                 integer *m, double *w, double *z__, integer *ldz, double *work, integer *iwork,
                 integer *ifail, integer *info);
int dsposv_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *b,
                 integer *ldb, double *x, integer *ldx, double *work, float *swork, integer *iter,
                 integer *info);
int dsprfs_check(char *uplo, integer *n, integer *nrhs, double *ap, double *afp, integer *ipiv,
                 double *b, integer *ldb, double *x, integer *ldx, double *ferr, double *berr,
                 double *work, integer *iwork, integer *info);
int dspsv_check(char *uplo, integer *n, integer *nrhs, double *ap, integer *ipiv, double *b,
                integer *ldb, integer *info);
int dspsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, double *ap, double *afp,
                 integer *ipiv, double *b, integer *ldb, double *x, integer *ldx, double *rcond,
                 double *ferr, double *berr, double *work, integer *iwork, integer *info);
int dsptrd_check(char *uplo, integer *n, double *ap, double *d__, double *e, double *tau,
                 integer *info);
int dsptrf_check(char *uplo, integer *n, double *ap, integer *ipiv, integer *info);
int dsptri_check(char *uplo, integer *n, double *ap, integer *ipiv, double *work, integer *info);
int dsptrs_check(char *uplo, integer *n, integer *nrhs, double *ap, integer *ipiv, double *b,
                 integer *ldb, integer *info);
int dstebz_check(char *range, char *order, integer *n, double *vl, double *vu, integer *il,
                 integer *iu, double *abstol, double *d__, double *e, integer *m, integer *nsplit,
                 double *w, integer *iblock, integer *isplit, double *work, integer *iwork,
                 integer *info);
int dstedc_check(char *compz, integer *n, double *d__, double *e, double *z__, integer *ldz,
                 double *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int dstegr_check(char *jobz, char *range, integer *n, double *d__, double *e, double *vl,
                 double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w,
                 double *z__, integer *ldz, integer *isuppz, double *work, integer *lwork,
                 integer *iwork, integer *liwork, integer *info);
int dstein_check(integer *n, double *d__, double *e, integer *m, double *w, integer *iblock,
                 integer *isplit, double *z__, integer *ldz, double *work, integer *iwork,
                 integer *ifail, integer *info);
int dstemr_check(char *jobz, char *range, integer *n, double *d__, double *e, double *vl,
                 double *vu, integer *il, integer *iu, integer *m, double *w, double *z__,
                 integer *ldz, integer *nzc, integer *isuppz, logical *tryrac, double *work,
                 integer *lwork, integer *iwork, integer *liwork, integer *info);
int dsteqr_check(char *compz, integer *n, double *d__, double *e, double *z__, integer *ldz,
                 double *work, integer *info);
int dsterf_check(integer *n, double *d__, double *e, integer *info);
int dstev_check(char *jobz, integer *n, double *d__, double *e, double *z__, integer *ldz,
                double *work, integer *info);
int dstevd_check(char *jobz, integer *n, double *d__, double *e, double *z__, integer *ldz,
                 double *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int dstevr_check(char *jobz, char *range, integer *n, double *d__, double *e, double *vl,
                 double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w,
                 double *z__, integer *ldz, integer *isuppz, double *work, integer *lwork,
                 integer *iwork, integer *liwork, integer *info);
int dstevx_check(char *jobz, char *range, integer *n, double *d__, double *e, double *vl,
                 double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w,
                 double *z__, integer *ldz, double *work, integer *iwork, integer *ifail,
                 integer *info);
int dsycon_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *anorm,
                 double *rcond, double *work, integer *iwork, integer *info);
int dsycon_rook_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *anorm,
                      double *rcond, double *work, integer *iwork, integer *info);
int dsyconv_check(char *uplo, char *way, integer *n, double *a, integer *lda, integer *ipiv,
                  double *work, integer *info);
int dsyequb_check(char *uplo, integer *n, double *a, integer *lda, double *s, double *scond,
                  double *amax, double *work, integer *info);
int dsyev_check(char *jobz, char *uplo, integer *n, double *a, integer *lda, double *w,
                double *work, integer *lwork, integer *info);
int dsyevd_check(char *jobz, char *uplo, integer *n, double *a, integer *lda, double *w,
                 double *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int dsyevr_check(char *jobz, char *range, char *uplo, integer *n, double *a, integer *lda,
                 double *vl, double *vu, integer *il, integer *iu, double *abstol, integer *m,
                 double *w, double *z__, integer *ldz, integer *isuppz, double *work,
                 integer *lwork, integer *iwork, integer *liwork, integer *info);
int dsyevx_check(char *jobz, char *range, char *uplo, integer *n, double *a, integer *lda,
                 double *vl, double *vu, integer *il, integer *iu, double *abstol, integer *m,
                 double *w, double *z__, integer *ldz, double *work, integer *lwork, integer *iwork,
                 integer *ifail, integer *info);
int dsygs2_check(integer *itype, char *uplo, integer *n, double *a, integer *lda, double *b,
                 integer *ldb, integer *info);
int dsygst_check(integer *itype, char *uplo, integer *n, double *a, integer *lda, double *b,
                 integer *ldb, integer *info);
int dsygv_check(integer *itype, char *jobz, char *uplo, integer *n, double *a, integer *lda,
                double *b, integer *ldb, double *w, double *work, integer *lwork, integer *info);
int dsygvd_check(integer *itype, char *jobz, char *uplo, integer *n, double *a, integer *lda,
                 double *b, integer *ldb, double *w, double *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int dsygvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, double *a,
                 integer *lda, double *b, integer *ldb, double *vl, double *vu, integer *il,
                 integer *iu, double *abstol, integer *m, double *w, double *z__, integer *ldz,
                 double *work, integer *lwork, integer *iwork, integer *ifail, integer *info);
int dsyrfs_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, double *af,
                 integer *ldaf, integer *ipiv, double *b, integer *ldb, double *x, integer *ldx,
                 double *ferr, double *berr, double *work, integer *iwork, integer *info);
int dsyrfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, double *a, integer *lda,
                  double *af, integer *ldaf, integer *ipiv, double *s, double *b, integer *ldb,
                  double *x, integer *ldx, double *rcond, double *berr, integer *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams,
                  double *params, double *work, integer *iwork, integer *info);
int dsysv_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv,
                double *b, integer *ldb, double *work, integer *lwork, integer *info);
int dsysv_rook_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv,
                     double *b, integer *ldb, double *work, integer *lwork, integer *info);
int dsysvx_check(char *fact, char *uplo, integer *n, integer *nrhs, double *a, integer *lda,
                 double *af, integer *ldaf, integer *ipiv, double *b, integer *ldb, double *x,
                 integer *ldx, double *rcond, double *ferr, double *berr, double *work,
                 integer *lwork, integer *iwork, integer *info);
int dsysvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, double *a, integer *lda,
                  double *af, integer *ldaf, integer *ipiv, char *equed, double *s, double *b,
                  integer *ldb, double *x, integer *ldx, double *rcond, double *rpvgrw,
                  double *berr, integer *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, integer *nparams, double *params, double *work,
                  integer *iwork, integer *info);
int dsyswapr_check(char *uplo, integer *n, double *a, integer *lda, integer *i1, integer *i2);
int dsytd2_check(char *uplo, integer *n, double *a, integer *lda, double *d__, double *e,
                 double *tau, integer *info);
int dsytf2_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, integer *info);
int dsytf2_rook_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv,
                      integer *info);
int dsytrd_check(char *uplo, integer *n, double *a, integer *lda, double *d__, double *e,
                 double *tau, double *work, integer *lwork, integer *info);
int dsytrf_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *work,
                 integer *lwork, integer *info);
int dsytrf_rook_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *work,
                      integer *lwork, integer *info);
int dsytri_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *work,
                 integer *info);
int dsytri2_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *work,
                  integer *lwork, integer *info);
int dsytri2x_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *work,
                   integer *nb, integer *info);
int dsytri_rook_check(char *uplo, integer *n, double *a, integer *lda, integer *ipiv, double *work,
                      integer *info);
int dsytrs_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv,
                 double *b, integer *ldb, integer *info);
int dsytrs2_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv,
                  double *b, integer *ldb, double *work, integer *info);
int dsytrs_rook_check(char *uplo, integer *n, integer *nrhs, double *a, integer *lda, integer *ipiv,
                      double *b, integer *ldb, integer *info);
int dtbcon_check(char *norm, char *uplo, char *diag, integer *n, integer *kd, double *ab,
                 integer *ldab, double *rcond, double *work, integer *iwork, integer *info);
int dtbrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 double *ab, integer *ldab, double *b, integer *ldb, double *x, integer *ldx,
                 double *ferr, double *berr, double *work, integer *iwork, integer *info);
int dtbtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 double *ab, integer *ldab, double *b, integer *ldb, integer *info);
int dtfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m,
                integer *n, double *alpha, double *a, double *b, integer *ldb);
int dtftri_check(char *transr, char *uplo, char *diag, integer *n, double *a, integer *info);
int dtfttp_check(char *transr, char *uplo, integer *n, double *arf, double *ap, integer *info);
int dtfttr_check(char *transr, char *uplo, integer *n, double *arf, double *a, integer *lda,
                 integer *info);
int dtgevc_check(char *side, char *howmny, logical *select, integer *n, double *s, integer *lds,
                 double *p, integer *ldp, double *vl, integer *ldvl, double *vr, integer *ldvr,
                 integer *mm, integer *m, double *work, integer *info);
int dtgex2_check(logical *wantq, logical *wantz, integer *n, double *a, integer *lda, double *b,
                 integer *ldb, double *q, integer *ldq, double *z__, integer *ldz, integer *j1,
                 integer *n1, integer *n2, double *work, integer *lwork, integer *info);
int dtgexc_check(logical *wantq, logical *wantz, integer *n, double *a, integer *lda, double *b,
                 integer *ldb, double *q, integer *ldq, double *z__, integer *ldz, integer *ifst,
                 integer *ilst, double *work, integer *lwork, integer *info);
int dtgsen_check(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n,
                 double *a, integer *lda, double *b, integer *ldb, double *alphar, double *alphai,
                 double *beta, double *q, integer *ldq, double *z__, integer *ldz, integer *m,
                 double *pl, double *pr, double *dif, double *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int dtgsja_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k,
                 integer *l, double *a, integer *lda, double *b, integer *ldb, double *tola,
                 double *tolb, double *alpha, double *beta, double *u, integer *ldu, double *v,
                 integer *ldv, double *q, integer *ldq, double *work, integer *ncycle,
                 integer *info);
int dtgsna_check(char *job, char *howmny, logical *select, integer *n, double *a, integer *lda,
                 double *b, integer *ldb, double *vl, integer *ldvl, double *vr, integer *ldvr,
                 double *s, double *dif, integer *mm, integer *m, double *work, integer *lwork,
                 integer *iwork, integer *info);
int dtgsy2_check(char *trans, integer *ijob, integer *m, integer *n, double *a, integer *lda,
                 double *b, integer *ldb, double *c__, integer *ldc, double *d__, integer *ldd,
                 double *e, integer *lde, double *f, integer *ldf, double *scale, double *rdsum,
                 double *rdscal, integer *iwork, integer *pq, integer *info);
int dtgsyl_check(char *trans, integer *ijob, integer *m, integer *n, double *a, integer *lda,
                 double *b, integer *ldb, double *c__, integer *ldc, double *d__, integer *ldd,
                 double *e, integer *lde, double *f, integer *ldf, double *scale, double *dif,
                 double *work, integer *lwork, integer *iwork, integer *info);
int dtpcon_check(char *norm, char *uplo, char *diag, integer *n, double *ap, double *rcond,
                 double *work, integer *iwork, integer *info);
int dtpmqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                  integer *nb, double *v, integer *ldv, double *t, integer *ldt, double *a,
                  integer *lda, double *b, integer *ldb, double *work, integer *info);
int dtpqrt_check(integer *m, integer *n, integer *l, integer *nb, double *a, integer *lda,
                 double *b, integer *ldb, double *t, integer *ldt, double *work, integer *info);
int dtpqrt2_check(integer *m, integer *n, integer *l, double *a, integer *lda, double *b,
                  integer *ldb, double *t, integer *ldt, integer *info);
int dtprfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, double *v, integer *ldv, double *t, integer *ldt,
                 double *a, integer *lda, double *b, integer *ldb, double *work, integer *ldwork);
int dtprfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, double *ap,
                 double *b, integer *ldb, double *x, integer *ldx, double *ferr, double *berr,
                 double *work, integer *iwork, integer *info);
int dtptri_check(char *uplo, char *diag, integer *n, double *ap, integer *info);
int dtptrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, double *ap,
                 double *b, integer *ldb, integer *info);
int dtpttf_check(char *transr, char *uplo, integer *n, double *ap, double *arf, integer *info);
int dtpttr_check(char *uplo, integer *n, double *ap, double *a, integer *lda, integer *info);
int dtrcon_check(char *norm, char *uplo, char *diag, integer *n, double *a, integer *lda,
                 double *rcond, double *work, integer *iwork, integer *info);
int dtrevc_check(char *side, char *howmny, logical *select, integer *n, double *t, integer *ldt,
                 double *vl, integer *ldvl, double *vr, integer *ldvr, integer *mm, integer *m,
                 double *work, integer *info);
int dtrexc_check(char *compq, integer *n, double *t, integer *ldt, double *q, integer *ldq,
                 integer *ifst, integer *ilst, double *work, integer *info);
int dtrrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, double *a,
                 integer *lda, double *b, integer *ldb, double *x, integer *ldx, double *ferr,
                 double *berr, double *work, integer *iwork, integer *info);
int dtrsen_check(char *job, char *compq, logical *select, integer *n, double *t, integer *ldt,
                 double *q, integer *ldq, double *wr, double *wi, integer *m, double *s,
                 double *sep, double *work, integer *lwork, integer *iwork, integer *liwork,
                 integer *info);
int dtrsna_check(char *job, char *howmny, logical *select, integer *n, double *t, integer *ldt,
                 double *vl, integer *ldvl, double *vr, integer *ldvr, double *s, double *sep,
                 integer *mm, integer *m, double *work, integer *ldwork, integer *iwork,
                 integer *info);
int dtrsyl_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, double *a,
                 integer *lda, double *b, integer *ldb, double *c__, integer *ldc, double *scale,
                 integer *info);
int dtrsyl3_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, doublereal *a,
                  integer *lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc,
                  doublereal *scale, integer *iwork, integer *liwork, doublereal *swork,
                  integer *ldswork, integer *info);
int dtrti2_check(char *uplo, char *diag, integer *n, double *a, integer *lda, integer *info);
int dtrtri_check(char *uplo, char *diag, integer *n, double *a, integer *lda, integer *info);
int dtrtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, double *a,
                 integer *lda, double *b, integer *ldb, integer *info);
int dtrttf_check(char *transr, char *uplo, integer *n, double *a, integer *lda, double *arf,
                 integer *info);
int dtrttp_check(char *uplo, integer *n, double *a, integer *lda, double *ap, integer *info);
int dtzrqf_check(integer *m, integer *n, double *a, integer *lda, double *tau, integer *info);
int dtzrzf_check(integer *m, integer *n, double *a, integer *lda, double *tau, double *work,
                 integer *lwork, integer *info);
double dzsum1_check(integer *n, dcomplex *cx, integer *incx);
int icmax1_check(integer *n, scomplex *cx, integer *incx);
int ilaclc_check(integer *m, integer *n, scomplex *a, integer *lda);
int ilaclr_check(integer *m, integer *n, scomplex *a, integer *lda);
int iladiag_check(char *diag);
int iladlc_check(integer *m, integer *n, double *a, integer *lda);
int iladlr_check(integer *m, integer *n, double *a, integer *lda);
int ilaprec_check(char *prec);
int ilaslc_check(integer *m, integer *n, float *a, integer *lda);
int ilaslr_check(integer *m, integer *n, float *a, integer *lda);
int ilatrans_check(char *trans);
int ilauplo_check(char *uplo);
int ilaver_check(integer *vers_major__, integer *vers_minor__, integer *vers_patch__);
int ilazlc_check(integer *m, integer *n, dcomplex *a, integer *lda);
int ilazlr_check(integer *m, integer *n, dcomplex *a, integer *lda);
int izmax1_check(integer *n, dcomplex *cx, integer *incx);
int sbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m,
                 integer *p, integer *q, float *theta, float *phi, float *u1, integer *ldu1,
                 float *u2, integer *ldu2, float *v1t, integer *ldv1t, float *v2t, integer *ldv2t,
                 float *b11d, float *b11e, float *b12d, float *b12e, float *b21d, float *b21e,
                 float *b22d, float *b22e, float *work, integer *lwork, integer *info);
int sbdsdc_check(char *uplo, char *compq, integer *n, float *d__, float *e, float *u, integer *ldu,
                 float *vt, integer *ldvt, float *q, integer *iq, float *work, integer *iwork,
                 integer *info);
int sbdsqr_check(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, float *d__,
                 float *e, float *vt, integer *ldvt, float *u, integer *ldu, float *c__,
                 integer *ldc, float *work, integer *info);
float scsum1_check(integer *n, scomplex *cx, integer *incx);
int sdisna_check(char *job, integer *m, integer *n, float *d__, float *sep, integer *info);
int sgbbrd_check(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku,
                 float *ab, integer *ldab, float *d__, float *e, float *q, integer *ldq, float *pt,
                 integer *ldpt, float *c__, integer *ldc, float *work, integer *info);
int sgbcon_check(char *norm, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                 integer *ipiv, float *anorm, float *rcond, float *work, integer *iwork,
                 integer *info);
int sgbequ_check(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, integer *info);
int sgbequb_check(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                  float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, integer *info);
int sgbrfs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, float *ab,
                 integer *ldab, float *afb, integer *ldafb, integer *ipiv, float *b, integer *ldb,
                 float *x, integer *ldx, float *ferr, float *berr, float *work, integer *iwork,
                 integer *info);
int sgbrfsx_check(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs,
                  float *ab, integer *ldab, float *afb, integer *ldafb, integer *ipiv, float *r__,
                  float *c__, float *b, integer *ldb, float *x, integer *ldx, float *rcond,
                  float *berr, integer *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, integer *nparams, float *params, float *work,
                  integer *iwork, integer *info);
int sgbsv_check(integer *n, integer *kl, integer *ku, integer *nrhs, float *ab, integer *ldab,
                integer *ipiv, float *b, integer *ldb, integer *info);
int sgbsvx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                 float *ab, integer *ldab, float *afb, integer *ldafb, integer *ipiv, char *equed,
                 float *r__, float *c__, float *b, integer *ldb, float *x, integer *ldx,
                 float *rcond, float *ferr, float *berr, float *work, integer *iwork,
                 integer *info);
int sgbsvxx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                  float *ab, integer *ldab, float *afb, integer *ldafb, integer *ipiv, char *equed,
                  float *r__, float *c__, float *b, integer *ldb, float *x, integer *ldx,
                  float *rcond, float *rpvgrw, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  float *work, integer *iwork, integer *info);
int sgbtf2_check(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                 integer *ipiv, integer *info);
int sgbtrf_check(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                 integer *ipiv, integer *info);
int sgbtrs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, float *ab,
                 integer *ldab, integer *ipiv, float *b, integer *ldb, integer *info);
int sgebak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, float *scale,
                 integer *m, float *v, integer *ldv, integer *info);
int sgebal_check(char *job, integer *n, float *a, integer *lda, integer *ilo, integer *ihi,
                 float *scale, integer *info);
int sgebd2_check(integer *m, integer *n, float *a, integer *lda, float *d__, float *e, float *tauq,
                 float *taup, float *work, integer *info);
int sgebrd_check(integer *m, integer *n, float *a, integer *lda, float *d__, float *e, float *tauq,
                 float *taup, float *work, integer *lwork, integer *info);
int sgecon_check(char *norm, integer *n, float *a, integer *lda, float *anorm, float *rcond,
                 float *work, integer *iwork, integer *info);
int sgeequ_check(integer *m, integer *n, float *a, integer *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, integer *info);
int sgeequb_check(integer *m, integer *n, float *a, integer *lda, float *r__, float *c__,
                  float *rowcnd, float *colcnd, float *amax, integer *info);
int sgees_check(char *jobvs, char *sort, L_fp select, integer *n, float *a, integer *lda,
                integer *sdim, float *wr, float *wi, float *vs, integer *ldvs, float *work,
                integer *lwork, logical *bwork, integer *info);
int sgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, integer *n, float *a,
                 integer *lda, integer *sdim, float *wr, float *wi, float *vs, integer *ldvs,
                 float *rconde, float *rcondv, float *work, integer *lwork, integer *iwork,
                 integer *liwork, logical *bwork, integer *info);
int sgeev_check(char *jobvl, char *jobvr, integer *n, float *a, integer *lda, float *wr, float *wi,
                float *vl, integer *ldvl, float *vr, integer *ldvr, float *work, integer *lwork,
                integer *info);
int sgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, float *a,
                 integer *lda, float *wr, float *wi, float *vl, integer *ldvl, float *vr,
                 integer *ldvr, integer *ilo, integer *ihi, float *scale, float *abnrm,
                 float *rconde, float *rcondv, float *work, integer *lwork, integer *iwork,
                 integer *info);
int sgegs_check(char *jobvsl, char *jobvsr, integer *n, float *a, integer *lda, float *b,
                integer *ldb, float *alphar, float *alphai, float *beta, float *vsl, integer *ldvsl,
                float *vsr, integer *ldvsr, float *work, integer *lwork, integer *info);
int sgegv_check(char *jobvl, char *jobvr, integer *n, float *a, integer *lda, float *b,
                integer *ldb, float *alphar, float *alphai, float *beta, float *vl, integer *ldvl,
                float *vr, integer *ldvr, float *work, integer *lwork, integer *info);
int sgehd2_check(integer *n, integer *ilo, integer *ihi, float *a, integer *lda, float *tau,
                 float *work, integer *info);
int sgehrd_check(integer *n, integer *ilo, integer *ihi, float *a, integer *lda, float *tau,
                 float *work, integer *lwork, integer *info);
int sgejsv_check(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, integer *m,
                 integer *n, float *a, integer *lda, float *sva, float *u, integer *ldu, float *v,
                 integer *ldv, float *work, integer *lwork, integer *iwork, integer *info);
int sgelq2_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *info);
int sgelqf_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *lwork, integer *info);
int sgels_check(char *trans, integer *m, integer *n, integer *nrhs, float *a, integer *lda,
                float *b, integer *ldb, float *work, integer *lwork, integer *info);
int sgelsd_check(integer *m, integer *n, integer *nrhs, float *a, integer *lda, float *b,
                 integer *ldb, float *s, float *rcond, integer *rank, float *work, integer *lwork,
                 integer *iwork, integer *info);
int sgelss_check(integer *m, integer *n, integer *nrhs, float *a, integer *lda, float *b,
                 integer *ldb, float *s, float *rcond, integer *rank, float *work, integer *lwork,
                 integer *info);
int sgelsx_check(integer *m, integer *n, integer *nrhs, float *a, integer *lda, float *b,
                 integer *ldb, integer *jpvt, float *rcond, integer *rank, float *work,
                 integer *info);
int sgelsy_check(integer *m, integer *n, integer *nrhs, float *a, integer *lda, float *b,
                 integer *ldb, integer *jpvt, float *rcond, integer *rank, float *work,
                 integer *lwork, integer *info);
int sgemqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *nb,
                  float *v, integer *ldv, float *t, integer *ldt, float *c__, integer *ldc,
                  float *work, integer *info);
int sgeql2_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *info);
int sgeqlf_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *lwork, integer *info);
int sgeqp3_check(integer *m, integer *n, float *a, integer *lda, integer *jpvt, float *tau,
                 float *work, integer *lwork, integer *info);
int sgeqpf_check(integer *m, integer *n, float *a, integer *lda, integer *jpvt, float *tau,
                 float *work, integer *info);
int sgeqr2_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *info);
int sgeqr2p_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                  integer *info);
int sgeqrf_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *lwork, integer *info);
int sgeqrfp_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                  integer *lwork, integer *info);
int sgeqrt_check(integer *m, integer *n, integer *nb, float *a, integer *lda, float *t,
                 integer *ldt, float *work, integer *info);
int sgeqrt2_check(integer *m, integer *n, float *a, integer *lda, float *t, integer *ldt,
                  integer *info);
int sgeqrt3_check(integer *m, integer *n, float *a, integer *lda, float *t, integer *ldt,
                  integer *info);
int sgerfs_check(char *trans, integer *n, integer *nrhs, float *a, integer *lda, float *af,
                 integer *ldaf, integer *ipiv, float *b, integer *ldb, float *x, integer *ldx,
                 float *ferr, float *berr, float *work, integer *iwork, integer *info);
int sgerfsx_check(char *trans, char *equed, integer *n, integer *nrhs, float *a, integer *lda,
                  float *af, integer *ldaf, integer *ipiv, float *r__, float *c__, float *b,
                  integer *ldb, float *x, integer *ldx, float *rcond, float *berr,
                  integer *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  integer *nparams, float *params, float *work, integer *iwork, integer *info);
int sgerq2_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *info);
int sgerqf_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *lwork, integer *info);
int sgesc2_check(integer *n, float *a, integer *lda, float *rhs, integer *ipiv, integer *jpiv,
                 float *scale);
int sgesdd_check(char *jobz, integer *m, integer *n, float *a, integer *lda, float *s, float *u,
                 integer *ldu, float *vt, integer *ldvt, float *work, integer *lwork,
                 integer *iwork, integer *info);
int sgesv_check(integer *n, integer *nrhs, float *a, integer *lda, integer *ipiv, float *b,
                integer *ldb, integer *info);
int sgesvd_check(char *jobu, char *jobvt, integer *m, integer *n, float *a, integer *lda, float *s,
                 float *u, integer *ldu, float *vt, integer *ldvt, float *work, integer *lwork,
                 integer *info);
int sgesvj_check(char *joba, char *jobu, char *jobv, integer *m, integer *n, float *a, integer *lda,
                 float *sva, integer *mv, float *v, integer *ldv, float *work, integer *lwork,
                 integer *info);
int sgesvx_check(char *fact, char *trans, integer *n, integer *nrhs, float *a, integer *lda,
                 float *af, integer *ldaf, integer *ipiv, char *equed, float *r__, float *c__,
                 float *b, integer *ldb, float *x, integer *ldx, float *rcond, float *ferr,
                 float *berr, float *work, integer *iwork, integer *info);
int sgesvxx_check(char *fact, char *trans, integer *n, integer *nrhs, float *a, integer *lda,
                  float *af, integer *ldaf, integer *ipiv, char *equed, float *r__, float *c__,
                  float *b, integer *ldb, float *x, integer *ldx, float *rcond, float *rpvgrw,
                  float *berr, integer *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, integer *nparams, float *params, float *work,
                  integer *iwork, integer *info);
int sgetc2_check(integer *n, float *a, integer *lda, integer *ipiv, integer *jpiv, integer *info);
int sgetf2_check(integer *m, integer *n, float *a, integer *lda, integer *ipiv, integer *info);
int sgetrf_check(integer *m, integer *n, float *a, integer *lda, integer *ipiv, integer *info);
int sgetrfnp_check(integer *m, integer *n, float *a, integer *lda, integer *info);
int sgetrfnpi_check(integer *m, integer *n, integer *nfact, float *a, integer *lda, integer *info);
int sgetri_check(integer *n, float *a, integer *lda, integer *ipiv, float *work, integer *lwork,
                 integer *info);
int sgetrs_check(char *trans, integer *n, integer *nrhs, float *a, integer *lda, integer *ipiv,
                 float *b, integer *ldb, integer *info);
int sggbak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, float *lscale,
                 float *rscale, integer *m, float *v, integer *ldv, integer *info);
int sggbal_check(char *job, integer *n, float *a, integer *lda, float *b, integer *ldb,
                 integer *ilo, integer *ihi, float *lscale, float *rscale, float *work,
                 integer *info);
int sgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, float *a,
                integer *lda, float *b, integer *ldb, integer *sdim, float *alphar, float *alphai,
                float *beta, float *vsl, integer *ldvsl, float *vsr, integer *ldvsr, float *work,
                integer *lwork, logical *bwork, integer *info);
int sggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n,
                 float *a, integer *lda, float *b, integer *ldb, integer *sdim, float *alphar,
                 float *alphai, float *beta, float *vsl, integer *ldvsl, float *vsr, integer *ldvsr,
                 float *rconde, float *rcondv, float *work, integer *lwork, integer *iwork,
                 integer *liwork, logical *bwork, integer *info);
int sggev_check(char *jobvl, char *jobvr, integer *n, float *a, integer *lda, float *b,
                integer *ldb, float *alphar, float *alphai, float *beta, float *vl, integer *ldvl,
                float *vr, integer *ldvr, float *work, integer *lwork, integer *info);
int sggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, float *a,
                 integer *lda, float *b, integer *ldb, float *alphar, float *alphai, float *beta,
                 float *vl, integer *ldvl, float *vr, integer *ldvr, integer *ilo, integer *ihi,
                 float *lscale, float *rscale, float *abnrm, float *bbnrm, float *rconde,
                 float *rcondv, float *work, integer *lwork, integer *iwork, logical *bwork,
                 integer *info);
int sggglm_check(integer *n, integer *m, integer *p, float *a, integer *lda, float *b, integer *ldb,
                 float *d__, float *x, float *y, float *work, integer *lwork, integer *info);
int sgghrd_check(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, float *a,
                 integer *lda, float *b, integer *ldb, float *q, integer *ldq, float *z__,
                 integer *ldz, integer *info);
int sgglse_check(integer *m, integer *n, integer *p, float *a, integer *lda, float *b, integer *ldb,
                 float *c__, float *d__, float *x, float *work, integer *lwork, integer *info);
int sggqrf_check(integer *n, integer *m, integer *p, float *a, integer *lda, float *taua, float *b,
                 integer *ldb, float *taub, float *work, integer *lwork, integer *info);
int sggrqf_check(integer *m, integer *p, integer *n, float *a, integer *lda, float *taua, float *b,
                 integer *ldb, float *taub, float *work, integer *lwork, integer *info);
int sggsvd_check(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k,
                 integer *l, float *a, integer *lda, float *b, integer *ldb, float *alpha,
                 float *beta, float *u, integer *ldu, float *v, integer *ldv, float *q,
                 integer *ldq, float *work, integer *iwork, integer *info);
int sggsvp_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, float *a,
                 integer *lda, float *b, integer *ldb, float *tola, float *tolb, integer *k,
                 integer *l, float *u, integer *ldu, float *v, integer *ldv, float *q, integer *ldq,
                 integer *iwork, float *tau, float *work, integer *info);
int sgsvj0_check(char *jobv, integer *m, integer *n, float *a, integer *lda, float *d__, float *sva,
                 integer *mv, float *v, integer *ldv, float *eps, float *sfmin, float *tol,
                 integer *nsweep, float *work, integer *lwork, integer *info);
int sgsvj1_check(char *jobv, integer *m, integer *n, integer *n1, float *a, integer *lda,
                 float *d__, float *sva, integer *mv, float *v, integer *ldv, float *eps,
                 float *sfmin, float *tol, integer *nsweep, float *work, integer *lwork,
                 integer *info);
int sgtcon_check(char *norm, integer *n, float *dl, float *d__, float *du, float *du2,
                 integer *ipiv, float *anorm, float *rcond, float *work, integer *iwork,
                 integer *info);
int sgtrfs_check(char *trans, integer *n, integer *nrhs, float *dl, float *d__, float *du,
                 float *dlf, float *df, float *duf, float *du2, integer *ipiv, float *b,
                 integer *ldb, float *x, integer *ldx, float *ferr, float *berr, float *work,
                 integer *iwork, integer *info);
int sgtsv_check(integer *n, integer *nrhs, float *dl, float *d__, float *du, float *b, integer *ldb,
                integer *info);
int sgtsvx_check(char *fact, char *trans, integer *n, integer *nrhs, float *dl, float *d__,
                 float *du, float *dlf, float *df, float *duf, float *du2, integer *ipiv, float *b,
                 integer *ldb, float *x, integer *ldx, float *rcond, float *ferr, float *berr,
                 float *work, integer *iwork, integer *info);
int sgttrf_check(integer *n, float *dl, float *d__, float *du, float *du2, integer *ipiv,
                 integer *info);
int sgttrs_check(char *trans, integer *n, integer *nrhs, float *dl, float *d__, float *du,
                 float *du2, integer *ipiv, float *b, integer *ldb, integer *info);
int sgtts2_check(integer *itrans, integer *n, integer *nrhs, float *dl, float *d__, float *du,
                 float *du2, integer *ipiv, float *b, integer *ldb);
int shgeqz_check(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi,
                 float *h__, integer *ldh, float *t, integer *ldt, float *alphar, float *alphai,
                 float *beta, float *q, integer *ldq, float *z__, integer *ldz, float *work,
                 integer *lwork, integer *info);
int shsein_check(char *side, char *eigsrc, char *initv, logical *select, integer *n, float *h__,
                 integer *ldh, float *wr, float *wi, float *vl, integer *ldvl, float *vr,
                 integer *ldvr, integer *mm, integer *m, float *work, integer *ifaill,
                 integer *ifailr, integer *info);
int shseqr_check(char *job, char *compz, integer *n, integer *ilo, integer *ihi, float *h__,
                 integer *ldh, float *wr, float *wi, float *z__, integer *ldz, float *work,
                 integer *lwork, integer *info);
logical sisnan_check(float *sin__);
int sla_gbamv_check(integer *trans, integer *m, integer *n, integer *kl, integer *ku, float *alpha,
                    float *ab, integer *ldab, float *x, integer *incx, float *beta, float *y,
                    integer *incy);
float sla_gbrcond_check(char *trans, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                        float *afb, integer *ldafb, integer *ipiv, integer *cmode, float *c__,
                        integer *info, float *work, integer *iwork);
int sla_gbrfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n, integer *kl,
                              integer *ku, integer *nrhs, float *ab, integer *ldab, float *afb,
                              integer *ldafb, integer *ipiv, logical *colequ, float *c__, float *b,
                              integer *ldb, float *y, integer *ldy, float *berr_out__,
                              integer *n_norms__, float *err_bnds_norm__, float *err_bnds_comp__,
                              float *res, float *ayb, float *dy, float *y_tail__, float *rcond,
                              integer *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, integer *info);
float sla_gbrpvgrw_check(integer *n, integer *kl, integer *ku, integer *ncols, float *ab,
                         integer *ldab, float *afb, integer *ldafb);
int sla_geamv_check(integer *trans, integer *m, integer *n, float *alpha, float *a, integer *lda,
                    float *x, integer *incx, float *beta, float *y, integer *incy);
float sla_gercond_check(char *trans, integer *n, float *a, integer *lda, float *af, integer *ldaf,
                        integer *ipiv, integer *cmode, float *c__, integer *info, float *work,
                        integer *iwork);
int sla_gerfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n,
                              integer *nrhs, float *a, integer *lda, float *af, integer *ldaf,
                              integer *ipiv, logical *colequ, float *c__, float *b, integer *ldb,
                              float *y, integer *ldy, float *berr_out__, integer *n_norms__,
                              float *errs_n__, float *errs_c__, float *res, float *ayb, float *dy,
                              float *y_tail__, float *rcond, integer *ithresh, float *rthresh,
                              float *dz_ub__, logical *ignore_cwise__, integer *info);
float sla_gerpvgrw_check(integer *n, integer *ncols, float *a, integer *lda, float *af,
                         integer *ldaf);
int sla_lin_berr_check(integer *n, integer *nz, integer *nrhs, float *res, float *ayb, float *berr);
float sla_porcond_check(char *uplo, integer *n, float *a, integer *lda, float *af, integer *ldaf,
                        integer *cmode, float *c__, integer *info, float *work, integer *iwork);
int sla_porfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs, float *a,
                              integer *lda, float *af, integer *ldaf, logical *colequ, float *c__,
                              float *b, integer *ldb, float *y, integer *ldy, float *berr_out__,
                              integer *n_norms__, float *err_bnds_norm__, float *err_bnds_comp__,
                              float *res, float *ayb, float *dy, float *y_tail__, float *rcond,
                              integer *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, integer *info);
float sla_porpvgrw_check(char *uplo, integer *ncols, float *a, integer *lda, float *af,
                         integer *ldaf, float *work);
int sla_syamv_check(integer *uplo, integer *n, float *alpha, float *a, integer *lda, float *x,
                    integer *incx, float *beta, float *y, integer *incy);
float sla_syrcond_check(char *uplo, integer *n, float *a, integer *lda, float *af, integer *ldaf,
                        integer *ipiv, integer *cmode, float *c__, integer *info, float *work,
                        integer *iwork);
int sla_syrfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs, float *a,
                              integer *lda, float *af, integer *ldaf, integer *ipiv,
                              logical *colequ, float *c__, float *b, integer *ldb, float *y,
                              integer *ldy, float *berr_out__, integer *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, float *res,
                              float *ayb, float *dy, float *y_tail__, float *rcond,
                              integer *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, integer *info);
float sla_syrpvgrw_check(char *uplo, integer *n, integer *info, float *a, integer *lda, float *af,
                         integer *ldaf, integer *ipiv, float *work);
int sla_wwaddw_check(integer *n, float *x, float *y, float *w);
int slabad_check(float *small_, float *large);
int slabrd_check(integer *m, integer *n, integer *nb, float *a, integer *lda, float *d__, float *e,
                 float *tauq, float *taup, float *x, integer *ldx, float *y, integer *ldy);
int slacn2_check(integer *n, float *v, float *x, integer *isgn, float *est, integer *kase,
                 integer *isave);
int slacon_check(integer *n, float *v, float *x, integer *isgn, float *est, integer *kase);
int slacpy_check(char *uplo, integer *m, integer *n, float *a, integer *lda, float *b,
                 integer *ldb);
int sladiv_check(float *a, float *b, float *c__, float *d__, float *p, float *q);
int slae2_check(float *a, float *b, float *c__, float *rt1, float *rt2);
int slaebz_check(integer *ijob, integer *nitmax, integer *n, integer *mmax, integer *minp,
                 integer *nbmin, float *abstol, float *reltol, float *pivmin, float *d__, float *e,
                 float *e2, integer *nval, float *ab, float *c__, integer *mout, integer *nab,
                 float *work, integer *iwork, integer *info);
int slaed0_check(integer *icompq, integer *qsiz, integer *n, float *d__, float *e, float *q,
                 integer *ldq, float *qstore, integer *ldqs, float *work, integer *iwork,
                 integer *info);
int slaed1_check(integer *n, float *d__, float *q, integer *ldq, integer *indxq, float *rho,
                 integer *cutpnt, float *work, integer *iwork, integer *info);
int slaed2_check(integer *k, integer *n, integer *n1, float *d__, float *q, integer *ldq,
                 integer *indxq, float *rho, float *z__, float *dlamda, float *w, float *q2,
                 integer *indx, integer *indxc, integer *indxp, integer *coltyp, integer *info);
int slaed3_check(integer *k, integer *n, integer *n1, float *d__, float *q, integer *ldq,
                 float *rho, float *dlamda, float *q2, integer *indx, integer *ctot, float *w,
                 float *s, integer *info);
int slaed4_check(integer *n, integer *i__, float *d__, float *z__, float *delta, float *rho,
                 float *dlam, integer *info);
int slaed5_check(integer *i__, float *d__, float *z__, float *delta, float *rho, float *dlam);
int slaed6_check(integer *kniter, logical *orgati, float *rho, float *d__, float *z__, float *finit,
                 float *tau, integer *info);
int slaed7_check(integer *icompq, integer *n, integer *qsiz, integer *tlvls, integer *curlvl,
                 integer *curpbm, float *d__, float *q, integer *ldq, integer *indxq, float *rho,
                 integer *cutpnt, float *qstore, integer *qptr, integer *prmptr, integer *perm,
                 integer *givptr, integer *givcol, float *givnum, float *work, integer *iwork,
                 integer *info);
int slaed8_check(integer *icompq, integer *k, integer *n, integer *qsiz, float *d__, float *q,
                 integer *ldq, integer *indxq, float *rho, integer *cutpnt, float *z__,
                 float *dlamda, float *q2, integer *ldq2, float *w, integer *perm, integer *givptr,
                 integer *givcol, float *givnum, integer *indxp, integer *indx, integer *info);
int slaed9_check(integer *k, integer *kstart, integer *kstop, integer *n, float *d__, float *q,
                 integer *ldq, float *rho, float *dlamda, float *w, float *s, integer *lds,
                 integer *info);
int slaeda_check(integer *n, integer *tlvls, integer *curlvl, integer *curpbm, integer *prmptr,
                 integer *perm, integer *givptr, integer *givcol, float *givnum, float *q,
                 integer *qptr, float *z__, float *ztemp, integer *info);
int slaein_check(logical *rightv, logical *noinit, integer *n, float *h__, integer *ldh, float *wr,
                 float *wi, float *vr, float *vi, float *b, integer *ldb, float *work, float *eps3,
                 float *smlnum, float *bignum, integer *info);
int slaev2_check(float *a, float *b, float *c__, float *rt1, float *rt2, float *cs1, float *sn1);
int slaexc_check(logical *wantq, integer *n, float *t, integer *ldt, float *q, integer *ldq,
                 integer *j1, integer *n1, integer *n2, float *work, integer *info);
int slag2_check(float *a, integer *lda, float *b, integer *ldb, float *safmin, float *scale1,
                float *scale2, float *wr1, float *wr2, float *wi);
int slag2d_check(integer *m, integer *n, float *sa, integer *ldsa, double *a, integer *lda,
                 integer *info);
int slags2_check(logical *upper, float *a1, float *a2, float *a3, float *b1, float *b2, float *b3,
                 float *csu, float *snu, float *csv, float *snv, float *csq, float *snq);
int slagtf_check(integer *n, float *a, float *lambda, float *b, float *c__, float *tol, float *d__,
                 integer *in, integer *info);
int slagtm_check(char *trans, integer *n, integer *nrhs, float *alpha, float *dl, float *d__,
                 float *du, float *x, integer *ldx, float *beta, float *b, integer *ldb);
int slagts_check(integer *job, integer *n, float *a, float *b, float *c__, float *d__, integer *in,
                 float *y, float *tol, integer *info);
int slagv2_check(float *a, integer *lda, float *b, integer *ldb, float *alphar, float *alphai,
                 float *beta, float *csl, float *snl, float *csr, float *snr);
int slahqr_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h__,
                 integer *ldh, float *wr, float *wi, integer *iloz, integer *ihiz, float *z__,
                 integer *ldz, integer *info);
int slahr2_check(integer *n, integer *k, integer *nb, float *a, integer *lda, float *tau, float *t,
                 integer *ldt, float *y, integer *ldy);
int slahrd_check(integer *n, integer *k, integer *nb, float *a, integer *lda, float *tau, float *t,
                 integer *ldt, float *y, integer *ldy);
int slaic1_check(integer *job, integer *j, float *x, float *sest, float *w, float *gamma,
                 float *sestpr, float *s, float *c__);
logical slaisnan_check(float *sin1, float *sin2);
int slaln2_check(logical *ltrans, integer *na, integer *nw, float *smin, float *ca, float *a,
                 integer *lda, float *d1, float *d2, float *b, integer *ldb, float *wr, float *wi,
                 float *x, integer *ldx, float *scale, float *xnorm, integer *info);
int slals0_check(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, float *b,
                 integer *ldb, float *bx, integer *ldbx, integer *perm, integer *givptr,
                 integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum, float *poles,
                 float *difl, float *difr, float *z__, integer *k, float *c__, float *s,
                 float *work, integer *info);
int slalsa_check(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, float *b,
                 integer *ldb, float *bx, integer *ldbx, float *u, integer *ldu, float *vt,
                 integer *k, float *difl, float *difr, float *z__, float *poles, integer *givptr,
                 integer *givcol, integer *ldgcol, integer *perm, float *givnum, float *c__,
                 float *s, float *work, integer *iwork, integer *info);
int slalsd_check(char *uplo, integer *smlsiz, integer *n, integer *nrhs, float *d__, float *e,
                 float *b, integer *ldb, float *rcond, integer *rank, float *work, integer *iwork,
                 integer *info);
int slamrg_check(integer *n1, integer *n2, float *a, integer *strd1, integer *strd2,
                 integer *index);
int slaneg_check(integer *n, float *d__, float *lld, float *sigma, float *pivmin, integer *r__);
float slangb_check(char *norm, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                   float *work);
float slange_check(char *norm, integer *m, integer *n, float *a, integer *lda, float *work);
float slangt_check(char *norm, integer *n, float *dl, float *d__, float *du);
float slanhs_check(char *norm, integer *n, float *a, integer *lda, float *work);
float slansb_check(char *norm, char *uplo, integer *n, integer *k, float *ab, integer *ldab,
                   float *work);
float slansf_check(char *norm, char *transr, char *uplo, integer *n, float *a, float *work);
float slansp_check(char *norm, char *uplo, integer *n, float *ap, float *work);
float slanst_check(char *norm, integer *n, float *d__, float *e);
float slansy_check(char *norm, char *uplo, integer *n, float *a, integer *lda, float *work);
float slantb_check(char *norm, char *uplo, char *diag, integer *n, integer *k, float *ab,
                   integer *ldab, float *work);
float slantp_check(char *norm, char *uplo, char *diag, integer *n, float *ap, float *work);
float slantr_check(char *norm, char *uplo, char *diag, integer *m, integer *n, float *a,
                   integer *lda, float *work);
int slanv2_check(float *a, float *b, float *c__, float *d__, float *rt1r, float *rt1i, float *rt2r,
                 float *rt2i, float *cs, float *sn);
int slapll_check(integer *n, float *x, integer *incx, float *y, integer *incy, float *ssmin);
int slapmr_check(logical *forwrd, integer *m, integer *n, float *x, integer *ldx, integer *k);
int slapmt_check(logical *forwrd, integer *m, integer *n, float *x, integer *ldx, integer *k);
float slapy2_check(float *x, float *y);
float slapy3_check(float *x, float *y, float *z__);
int slaqgb_check(integer *m, integer *n, integer *kl, integer *ku, float *ab, integer *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *equed);
int slaqge_check(integer *m, integer *n, float *a, integer *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, char *equed);
int slaqp2_check(integer *m, integer *n, integer *offset, float *a, integer *lda, integer *jpvt,
                 float *tau, float *vn1, float *vn2, float *work);
int slaqps_check(integer *m, integer *n, integer *offset, integer *nb, integer *kb, float *a,
                 integer *lda, integer *jpvt, float *tau, float *vn1, float *vn2, float *auxv,
                 float *f, integer *ldf);
int slaqr0_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h__,
                 integer *ldh, float *wr, float *wi, integer *iloz, integer *ihiz, float *z__,
                 integer *ldz, float *work, integer *lwork, integer *info);
int slaqr1_check(integer *n, float *h__, integer *ldh, float *sr1, float *si1, float *sr2,
                 float *si2, float *v);
int slaqr2_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, float *h__, integer *ldh, integer *iloz, integer *ihiz, float *z__,
                 integer *ldz, integer *ns, integer *nd, float *sr, float *si, float *v,
                 integer *ldv, integer *nh, float *t, integer *ldt, integer *nv, float *wv,
                 integer *ldwv, float *work, integer *lwork);
int slaqr3_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, float *h__, integer *ldh, integer *iloz, integer *ihiz, float *z__,
                 integer *ldz, integer *ns, integer *nd, float *sr, float *si, float *v,
                 integer *ldv, integer *nh, float *t, integer *ldt, integer *nv, float *wv,
                 integer *ldwv, float *work, integer *lwork);
int slaqr4_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi, float *h__,
                 integer *ldh, float *wr, float *wi, integer *iloz, integer *ihiz, float *z__,
                 integer *ldz, float *work, integer *lwork, integer *info);
int slaqr5_check(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
                 integer *kbot, integer *nshfts, float *sr, float *si, float *h__, integer *ldh,
                 integer *iloz, integer *ihiz, float *z__, integer *ldz, float *v, integer *ldv,
                 float *u, integer *ldu, integer *nv, float *wv, integer *ldwv, integer *nh,
                 float *wh, integer *ldwh);
int slaqsb_check(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *s,
                 float *scond, float *amax, char *equed);
int slaqsp_check(char *uplo, integer *n, float *ap, float *s, float *scond, float *amax,
                 char *equed);
int slaqsy_check(char *uplo, integer *n, float *a, integer *lda, float *s, float *scond,
                 float *amax, char *equed);
int slaqtr_check(logical *ltran, logical *lfloat, integer *n, float *t, integer *ldt, float *b,
                 float *w, float *scale, float *x, float *work, integer *info);
int slar1v_check(integer *n, integer *b1, integer *bn, float *lambda, float *d__, float *l,
                 float *ld, float *lld, float *pivmin, float *gaptol, float *z__, logical *wantnc,
                 integer *negcnt, float *ztz, float *mingma, integer *r__, integer *isuppz,
                 float *nrminv, float *resid, float *rqcorr, float *work);
int slar2v_check(integer *n, float *x, float *y, float *z__, integer *incx, float *c__, float *s,
                 integer *incc);
int slarf_check(char *side, integer *m, integer *n, float *v, integer *incv, float *tau, float *c__,
                integer *ldc, float *work);
int slarfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, float *v, integer *ldv, float *t, integer *ldt, float *c__,
                 integer *ldc, float *work, integer *ldwork);
int slarfg_check(integer *n, float *alpha, float *x, integer *incx, float *tau);
int slarfgp_check(integer *n, float *alpha, float *x, integer *incx, float *tau);
int slarft_check(char *direct, char *storev, integer *n, integer *k, float *v, integer *ldv,
                 float *tau, float *t, integer *ldt);
int slarfx_check(char *side, integer *m, integer *n, float *v, float *tau, float *c__, integer *ldc,
                 float *work);
int slargv_check(integer *n, float *x, integer *incx, float *y, integer *incy, float *c__,
                 integer *incc);
int slarnv_check(integer *idist, integer *iseed, integer *n, float *x);
int slarra_check(integer *n, float *d__, float *e, float *e2, float *spltol, float *tnrm,
                 integer *nsplit, integer *isplit, integer *info);
int slarrb_check(integer *n, float *d__, float *lld, integer *ifirst, integer *ilast, float *rtol1,
                 float *rtol2, integer *offset, float *w, float *wgap, float *werr, float *work,
                 integer *iwork, float *pivmin, float *spdiam, integer *twist, integer *info);
int slarrc_check(char *jobt, integer *n, float *vl, float *vu, float *d__, float *e, float *pivmin,
                 integer *eigcnt, integer *lcnt, integer *rcnt, integer *info);
int slarrd_check(char *range, char *order, integer *n, float *vl, float *vu, integer *il,
                 integer *iu, float *gers, float *reltol, float *d__, float *e, float *e2,
                 float *pivmin, integer *nsplit, integer *isplit, integer *m, float *w, float *werr,
                 float *wl, float *wu, integer *iblock, integer *indexw, float *work,
                 integer *iwork, integer *info);
int slarre_check(char *range, integer *n, float *vl, float *vu, integer *il, integer *iu,
                 float *d__, float *e, float *e2, float *rtol1, float *rtol2, float *spltol,
                 integer *nsplit, integer *isplit, integer *m, float *w, float *werr, float *wgap,
                 integer *iblock, integer *indexw, float *gers, float *pivmin, float *work,
                 integer *iwork, integer *info);
int slarrf_check(integer *n, float *d__, float *l, float *ld, integer *clstrt, integer *clend,
                 float *w, float *wgap, float *werr, float *spdiam, float *clgapl, float *clgapr,
                 float *pivmin, float *sigma, float *dplus, float *lplus, float *work,
                 integer *info);
int slarrj_check(integer *n, float *d__, float *e2, integer *ifirst, integer *ilast, float *rtol,
                 integer *offset, float *w, float *werr, float *work, integer *iwork, float *pivmin,
                 float *spdiam, integer *info);
int slarrk_check(integer *n, integer *iw, float *gl, float *gu, float *d__, float *e2,
                 float *pivmin, float *reltol, float *w, float *werr, integer *info);
int slarrr_check(integer *n, float *d__, float *e, integer *info);
int slarrv_check(integer *n, float *vl, float *vu, float *d__, float *l, float *pivmin,
                 integer *isplit, integer *m, integer *dol, integer *dou, float *minrgp,
                 float *rtol1, float *rtol2, float *w, float *werr, float *wgap, integer *iblock,
                 integer *indexw, float *gers, float *z__, integer *ldz, integer *isuppz,
                 float *work, integer *iwork, integer *info);
int slarscl2_check(integer *m, integer *n, float *d__, float *x, integer *ldx);
int slartg_check(float *f, float *g, float *cs, float *sn, float *r__);
int slartgp_check(float *f, float *g, float *cs, float *sn, float *r__);
int slartgs_check(float *x, float *y, float *sigma, float *cs, float *sn);
int slartv_check(integer *n, float *x, integer *incx, float *y, integer *incy, float *c__, float *s,
                 integer *incc);
int slaruv_check(integer *iseed, integer *n, float *x);
int slarz_check(char *side, integer *m, integer *n, integer *l, float *v, integer *incv, float *tau,
                float *c__, integer *ldc, float *work);
int slarzb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, float *v, integer *ldv, float *t, integer *ldt, float *c__,
                 integer *ldc, float *work, integer *ldwork);
int slarzt_check(char *direct, char *storev, integer *n, integer *k, float *v, integer *ldv,
                 float *tau, float *t, integer *ldt);
int slas2_check(float *f, float *g, float *h__, float *ssmin, float *ssmax);
int slascl_check(char *type__, integer *kl, integer *ku, float *cfrom, float *cto, integer *m,
                 integer *n, float *a, integer *lda, integer *info);
int slascl2_check(integer *m, integer *n, float *d__, float *x, integer *ldx);
int slasd0_check(integer *n, integer *sqre, float *d__, float *e, float *u, integer *ldu, float *vt,
                 integer *ldvt, integer *smlsiz, integer *iwork, float *work, integer *info);
int slasd1_check(integer *nl, integer *nr, integer *sqre, float *d__, float *alpha, float *beta,
                 float *u, integer *ldu, float *vt, integer *ldvt, integer *idxq, integer *iwork,
                 float *work, integer *info);
int slasd2_check(integer *nl, integer *nr, integer *sqre, integer *k, float *d__, float *z__,
                 float *alpha, float *beta, float *u, integer *ldu, float *vt, integer *ldvt,
                 float *dsigma, float *u2, integer *ldu2, float *vt2, integer *ldvt2, integer *idxp,
                 integer *idx, integer *idxc, integer *idxq, integer *coltyp, integer *info);
int slasd3_check(integer *nl, integer *nr, integer *sqre, integer *k, float *d__, float *q,
                 integer *ldq, float *dsigma, float *u, integer *ldu, float *u2, integer *ldu2,
                 float *vt, integer *ldvt, float *vt2, integer *ldvt2, integer *idxc, integer *ctot,
                 float *z__, integer *info);
int slasd4_check(integer *n, integer *i__, float *d__, float *z__, float *delta, float *rho,
                 float *sigma, float *work, integer *info);
int slasd5_check(integer *i__, float *d__, float *z__, float *delta, float *rho, float *dsigma,
                 float *work);
int slasd6_check(integer *icompq, integer *nl, integer *nr, integer *sqre, float *d__, float *vf,
                 float *vl, float *alpha, float *beta, integer *idxq, integer *perm,
                 integer *givptr, integer *givcol, integer *ldgcol, float *givnum, integer *ldgnum,
                 float *poles, float *difl, float *difr, float *z__, integer *k, float *c__,
                 float *s, float *work, integer *iwork, integer *info);
int slasd7_check(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, float *d__,
                 float *z__, float *zw, float *vf, float *vfw, float *vl, float *vlw, float *alpha,
                 float *beta, float *dsigma, integer *idx, integer *idxp, integer *idxq,
                 integer *perm, integer *givptr, integer *givcol, integer *ldgcol, float *givnum,
                 integer *ldgnum, float *c__, float *s, integer *info);
int slasd8_check(integer *icompq, integer *k, float *d__, float *z__, float *vf, float *vl,
                 float *difl, float *difr, integer *lddifr, float *dsigma, float *work,
                 integer *info);
int slasda_check(integer *icompq, integer *smlsiz, integer *n, integer *sqre, float *d__, float *e,
                 float *u, integer *ldu, float *vt, integer *k, float *difl, float *difr,
                 float *z__, float *poles, integer *givptr, integer *givcol, integer *ldgcol,
                 integer *perm, float *givnum, float *c__, float *s, float *work, integer *iwork,
                 integer *info);
int slasdq_check(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc,
                 float *d__, float *e, float *vt, integer *ldvt, float *u, integer *ldu, float *c__,
                 integer *ldc, float *work, integer *info);
int slasdt_check(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml,
                 integer *ndimr, integer *msub);
int slaset_check(char *uplo, integer *m, integer *n, float *alpha, float *beta, float *a,
                 integer *lda);
int slasq1_check(integer *n, float *d__, float *e, float *work, integer *info);
int slasq2_check(integer *n, float *z__, integer *info);
int slasq3_check(integer *i0, integer *n0, float *z__, integer *pp, float *dmin__, float *sigma,
                 float *desig, float *qmax, integer *nfail, integer *iter, integer *ndiv,
                 logical *ieee, integer *ttype, float *dmin1, float *dmin2, float *dn, float *dn1,
                 float *dn2, float *g, float *tau);
int slasq4_check(integer *i0, integer *n0, float *z__, integer *pp, integer *n0in, float *dmin__,
                 float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau,
                 integer *ttype, float *g);
int slasq5_check(integer *i0, integer *n0, float *z__, integer *pp, float *tau, float *sigma,
                 float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float *dnm2,
                 logical *ieee, float *eps);
int slasq6_check(integer *i0, integer *n0, float *z__, integer *pp, float *dmin__, float *dmin1,
                 float *dmin2, float *dn, float *dnm1, float *dnm2);
int slasr_check(char *side, char *pivot, char *direct, integer *m, integer *n, float *c__, float *s,
                float *a, integer *lda);
int slasrt_check(char *id, integer *n, float *d__, integer *info);
int slassq_check(integer *n, float *x, integer *incx, float *scale, float *sumsq);
int slasv2_check(float *f, float *g, float *h__, float *ssmin, float *ssmax, float *snr, float *csr,
                 float *snl, float *csl);
int slaswp_check(integer *n, float *a, integer *lda, integer *k1, integer *k2, integer *ipiv,
                 integer *incx);
int slasy2_check(logical *ltranl, logical *ltranr, integer *isgn, integer *n1, integer *n2,
                 float *tl, integer *ldtl, float *tr, integer *ldtr, float *b, integer *ldb,
                 float *scale, float *x, integer *ldx, float *xnorm, integer *info);
int slasyf_check(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda,
                 integer *ipiv, float *w, integer *ldw, integer *info);
int slasyf_rook_check(char *uplo, integer *n, integer *nb, integer *kb, float *a, integer *lda,
                      integer *ipiv, float *w, integer *ldw, integer *info);
int slatbs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd,
                 float *ab, integer *ldab, float *x, float *scale, float *cnorm, integer *info);
int slatdf_check(integer *ijob, integer *n, float *z__, integer *ldz, float *rhs, float *rdsum,
                 float *rdscal, integer *ipiv, integer *jpiv);
int slatps_check(char *uplo, char *trans, char *diag, char *normin, integer *n, float *ap, float *x,
                 float *scale, float *cnorm, integer *info);
int slatrd_check(char *uplo, integer *n, integer *nb, float *a, integer *lda, float *e, float *tau,
                 float *w, integer *ldw);
int slatrs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, float *a,
                 integer *lda, float *x, float *scale, float *cnorm, integer *info);
int slatrz_check(integer *m, integer *n, integer *l, float *a, integer *lda, float *tau,
                 float *work);
int slatzm_check(char *side, integer *m, integer *n, float *v, integer *incv, float *tau, float *c1,
                 float *c2, integer *ldc, float *work);
int slauu2_check(char *uplo, integer *n, float *a, integer *lda, integer *info);
int slauum_check(char *uplo, integer *n, float *a, integer *lda, integer *info);
int sopgtr_check(char *uplo, integer *n, float *ap, float *tau, float *q, integer *ldq, float *work,
                 integer *info);
int sopmtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, float *ap, float *tau,
                 float *c__, integer *ldc, float *work, integer *info);
int sorbdb_check(char *trans, char *signs, integer *m, integer *p, integer *q, float *x11,
                 integer *ldx11, float *x12, integer *ldx12, float *x21, integer *ldx21, float *x22,
                 integer *ldx22, float *theta, float *phi, float *taup1, float *taup2, float *tauq1,
                 float *tauq2, float *work, integer *lwork, integer *info);
int sorbdb1_check(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21,
                  integer *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *work, integer *lwork, integer *info);
int sorbdb2_check(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21,
                  integer *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *work, integer *lwork, integer *info);
int sorbdb3_check(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21,
                  integer *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *work, integer *lwork, integer *info);
int sorbdb4_check(integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x21,
                  integer *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *phantom, float *work, integer *lwork, integer *info);
int sorbdb5_check(integer *m1, integer *m2, integer *n, float *x1, integer *incx1, float *x2,
                  integer *incx2, float *q1, integer *ldq1, float *q2, integer *ldq2, float *work,
                  integer *lwork, integer *info);
int sorbdb6_check(integer *m1, integer *m2, integer *n, float *x1, integer *incx1, float *x2,
                  integer *incx2, float *q1, integer *ldq1, float *q2, integer *ldq2, float *work,
                  integer *lwork, integer *info);
int sorcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 integer *m, integer *p, integer *q, float *x11, integer *ldx11, float *x12,
                 integer *ldx12, float *x21, integer *ldx21, float *x22, integer *ldx22,
                 float *theta, float *u1, integer *ldu1, float *u2, integer *ldu2, float *v1t,
                 integer *ldv1t, float *v2t, integer *ldv2t, float *work, integer *lwork,
                 integer *iwork, integer *info);
int sorcsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, integer *m, integer *p, integer *q,
                     float *x11, integer *ldx11, float *x21, integer *ldx21, float *theta,
                     float *u1, integer *ldu1, float *u2, integer *ldu2, float *v1t, integer *ldv1t,
                     float *work, integer *lwork, integer *iwork, integer *info);
int sorg2l_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *info);
int sorg2r_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *info);
int sorgbr_check(char *vect, integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *lwork, integer *info);
int sorghr_check(integer *n, integer *ilo, integer *ihi, float *a, integer *lda, float *tau,
                 float *work, integer *lwork, integer *info);
int sorgl2_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *info);
int sorglq_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *lwork, integer *info);
int sorgql_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *lwork, integer *info);
int sorgqr_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *lwork, integer *info);
int sorgr2_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *info);
int sorgrq_check(integer *m, integer *n, integer *k, float *a, integer *lda, float *tau,
                 float *work, integer *lwork, integer *info);
int sorgtr_check(char *uplo, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *lwork, integer *info);
int sorm2l_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *info);
int sorm2r_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *info);
int sormbr_check(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork,
                 integer *info);
int sormhr_check(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi,
                 float *a, integer *lda, float *tau, float *c__, integer *ldc, float *work,
                 integer *lwork, integer *info);
int sorml2_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *info);
int sormlq_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork,
                 integer *info);
int sormql_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork,
                 integer *info);
int sormqr_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork,
                 integer *info);
int sormr2_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *info);
int sormr3_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *info);
int sormrq_check(char *side, char *trans, integer *m, integer *n, integer *k, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork,
                 integer *info);
int sormrz_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork,
                 integer *info);
int sormtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, float *a,
                 integer *lda, float *tau, float *c__, integer *ldc, float *work, integer *lwork,
                 integer *info);
int spbcon_check(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *anorm,
                 float *rcond, float *work, integer *iwork, integer *info);
int spbequ_check(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *s,
                 float *scond, float *amax, integer *info);
int spbrfs_check(char *uplo, integer *n, integer *kd, integer *nrhs, float *ab, integer *ldab,
                 float *afb, integer *ldafb, float *b, integer *ldb, float *x, integer *ldx,
                 float *ferr, float *berr, float *work, integer *iwork, integer *info);
int spbstf_check(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, integer *info);
int spbsv_check(char *uplo, integer *n, integer *kd, integer *nrhs, float *ab, integer *ldab,
                float *b, integer *ldb, integer *info);
int spbsvx_check(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, float *ab,
                 integer *ldab, float *afb, integer *ldafb, char *equed, float *s, float *b,
                 integer *ldb, float *x, integer *ldx, float *rcond, float *ferr, float *berr,
                 float *work, integer *iwork, integer *info);
int spbtf2_check(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, integer *info);
int spbtrf_check(char *uplo, integer *n, integer *kd, float *ab, integer *ldab, integer *info);
int spbtrs_check(char *uplo, integer *n, integer *kd, integer *nrhs, float *ab, integer *ldab,
                 float *b, integer *ldb, integer *info);
int spftrf_check(char *transr, char *uplo, integer *n, float *a, integer *info);
int spftri_check(char *transr, char *uplo, integer *n, float *a, integer *info);
int spftrs_check(char *transr, char *uplo, integer *n, integer *nrhs, float *a, float *b,
                 integer *ldb, integer *info);
int spocon_check(char *uplo, integer *n, float *a, integer *lda, float *anorm, float *rcond,
                 float *work, integer *iwork, integer *info);
int spoequ_check(integer *n, float *a, integer *lda, float *s, float *scond, float *amax,
                 integer *info);
int spoequb_check(integer *n, float *a, integer *lda, float *s, float *scond, float *amax,
                  integer *info);
int sporfs_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *af,
                 integer *ldaf, float *b, integer *ldb, float *x, integer *ldx, float *ferr,
                 float *berr, float *work, integer *iwork, integer *info);
int sporfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, float *a, integer *lda,
                  float *af, integer *ldaf, float *s, float *b, integer *ldb, float *x,
                  integer *ldx, float *rcond, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  float *work, integer *iwork, integer *info);
int sposv_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *b,
                integer *ldb, integer *info);
int sposvx_check(char *fact, char *uplo, integer *n, integer *nrhs, float *a, integer *lda,
                 float *af, integer *ldaf, char *equed, float *s, float *b, integer *ldb, float *x,
                 integer *ldx, float *rcond, float *ferr, float *berr, float *work, integer *iwork,
                 integer *info);
int sposvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, float *a, integer *lda,
                  float *af, integer *ldaf, char *equed, float *s, float *b, integer *ldb, float *x,
                  integer *ldx, float *rcond, float *rpvgrw, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  float *work, integer *iwork, integer *info);
int spotf2_check(char *uplo, integer *n, float *a, integer *lda, integer *info);
int spotrf_check(char *uplo, integer *n, float *a, integer *lda, integer *info);
int spotri_check(char *uplo, integer *n, float *a, integer *lda, integer *info);
int spotrs_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *b,
                 integer *ldb, integer *info);
int sppcon_check(char *uplo, integer *n, float *ap, float *anorm, float *rcond, float *work,
                 integer *iwork, integer *info);
int sppequ_check(char *uplo, integer *n, float *ap, float *s, float *scond, float *amax,
                 integer *info);
int spprfs_check(char *uplo, integer *n, integer *nrhs, float *ap, float *afp, float *b,
                 integer *ldb, float *x, integer *ldx, float *ferr, float *berr, float *work,
                 integer *iwork, integer *info);
int sppsv_check(char *uplo, integer *n, integer *nrhs, float *ap, float *b, integer *ldb,
                integer *info);
int sppsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, float *ap, float *afp,
                 char *equed, float *s, float *b, integer *ldb, float *x, integer *ldx,
                 float *rcond, float *ferr, float *berr, float *work, integer *iwork,
                 integer *info);
int spptrf_check(char *uplo, integer *n, float *ap, integer *info);
int spptri_check(char *uplo, integer *n, float *ap, integer *info);
int spptrs_check(char *uplo, integer *n, integer *nrhs, float *ap, float *b, integer *ldb,
                 integer *info);
int spstf2_check(char *uplo, integer *n, float *a, integer *lda, integer *piv, integer *rank,
                 float *tol, float *work, integer *info);
int spstrf_check(char *uplo, integer *n, float *a, integer *lda, integer *piv, integer *rank,
                 float *tol, float *work, integer *info);
int sptcon_check(integer *n, float *d__, float *e, float *anorm, float *rcond, float *work,
                 integer *info);
int spteqr_check(char *compz, integer *n, float *d__, float *e, float *z__, integer *ldz,
                 float *work, integer *info);
int sptrfs_check(integer *n, integer *nrhs, float *d__, float *e, float *df, float *ef, float *b,
                 integer *ldb, float *x, integer *ldx, float *ferr, float *berr, float *work,
                 integer *info);
int sptsv_check(integer *n, integer *nrhs, float *d__, float *e, float *b, integer *ldb,
                integer *info);
int sptsvx_check(char *fact, integer *n, integer *nrhs, float *d__, float *e, float *df, float *ef,
                 float *b, integer *ldb, float *x, integer *ldx, float *rcond, float *ferr,
                 float *berr, float *work, integer *info);
int spttrf_check(integer *n, float *d__, float *e, integer *info);
int spttrs_check(integer *n, integer *nrhs, float *d__, float *e, float *b, integer *ldb,
                 integer *info);
int sptts2_check(integer *n, integer *nrhs, float *d__, float *e, float *b, integer *ldb);
int srscl_check(integer *n, float *sa, float *sx, integer *incx);
int ssbev_check(char *jobz, char *uplo, integer *n, integer *kd, float *ab, integer *ldab, float *w,
                float *z__, integer *ldz, float *work, integer *info);
int ssbevd_check(char *jobz, char *uplo, integer *n, integer *kd, float *ab, integer *ldab,
                 float *w, float *z__, integer *ldz, float *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int ssbevx_check(char *jobz, char *range, char *uplo, integer *n, integer *kd, float *ab,
                 integer *ldab, float *q, integer *ldq, float *vl, float *vu, integer *il,
                 integer *iu, float *abstol, integer *m, float *w, float *z__, integer *ldz,
                 float *work, integer *iwork, integer *ifail, integer *info);
int ssbgst_check(char *vect, char *uplo, integer *n, integer *ka, integer *kb, float *ab,
                 integer *ldab, float *bb, integer *ldbb, float *x, integer *ldx, float *work,
                 integer *info);
int ssbgv_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, float *ab,
                integer *ldab, float *bb, integer *ldbb, float *w, float *z__, integer *ldz,
                float *work, integer *info);
int ssbgvd_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, float *ab,
                 integer *ldab, float *bb, integer *ldbb, float *w, float *z__, integer *ldz,
                 float *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int ssbgvx_check(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb,
                 float *ab, integer *ldab, float *bb, integer *ldbb, float *q, integer *ldq,
                 float *vl, float *vu, integer *il, integer *iu, float *abstol, integer *m,
                 float *w, float *z__, integer *ldz, float *work, integer *iwork, integer *ifail,
                 integer *info);
int ssbtrd_check(char *vect, char *uplo, integer *n, integer *kd, float *ab, integer *ldab,
                 float *d__, float *e, float *q, integer *ldq, float *work, integer *info);
int ssfrk_check(char *transr, char *uplo, char *trans, integer *n, integer *k, float *alpha,
                float *a, integer *lda, float *beta, float *c__);
int sspcon_check(char *uplo, integer *n, float *ap, integer *ipiv, float *anorm, float *rcond,
                 float *work, integer *iwork, integer *info);
int sspev_check(char *jobz, char *uplo, integer *n, float *ap, float *w, float *z__, integer *ldz,
                float *work, integer *info);
int sspevd_check(char *jobz, char *uplo, integer *n, float *ap, float *w, float *z__, integer *ldz,
                 float *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int sspevx_check(char *jobz, char *range, char *uplo, integer *n, float *ap, float *vl, float *vu,
                 integer *il, integer *iu, float *abstol, integer *m, float *w, float *z__,
                 integer *ldz, float *work, integer *iwork, integer *ifail, integer *info);
int sspgst_check(integer *itype, char *uplo, integer *n, float *ap, float *bp, integer *info);
int sspgv_check(integer *itype, char *jobz, char *uplo, integer *n, float *ap, float *bp, float *w,
                float *z__, integer *ldz, float *work, integer *info);
int sspgvd_check(integer *itype, char *jobz, char *uplo, integer *n, float *ap, float *bp, float *w,
                 float *z__, integer *ldz, float *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int sspgvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, float *ap,
                 float *bp, float *vl, float *vu, integer *il, integer *iu, float *abstol,
                 integer *m, float *w, float *z__, integer *ldz, float *work, integer *iwork,
                 integer *ifail, integer *info);
int ssprfs_check(char *uplo, integer *n, integer *nrhs, float *ap, float *afp, integer *ipiv,
                 float *b, integer *ldb, float *x, integer *ldx, float *ferr, float *berr,
                 float *work, integer *iwork, integer *info);
int sspsv_check(char *uplo, integer *n, integer *nrhs, float *ap, integer *ipiv, float *b,
                integer *ldb, integer *info);
int sspsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, float *ap, float *afp,
                 integer *ipiv, float *b, integer *ldb, float *x, integer *ldx, float *rcond,
                 float *ferr, float *berr, float *work, integer *iwork, integer *info);
int ssptrd_check(char *uplo, integer *n, float *ap, float *d__, float *e, float *tau,
                 integer *info);
int ssptrf_check(char *uplo, integer *n, float *ap, integer *ipiv, integer *info);
int ssptri_check(char *uplo, integer *n, float *ap, integer *ipiv, float *work, integer *info);
int ssptrs_check(char *uplo, integer *n, integer *nrhs, float *ap, integer *ipiv, float *b,
                 integer *ldb, integer *info);
int sstebz_check(char *range, char *order, integer *n, float *vl, float *vu, integer *il,
                 integer *iu, float *abstol, float *d__, float *e, integer *m, integer *nsplit,
                 float *w, integer *iblock, integer *isplit, float *work, integer *iwork,
                 integer *info);
int sstedc_check(char *compz, integer *n, float *d__, float *e, float *z__, integer *ldz,
                 float *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int sstegr_check(char *jobz, char *range, integer *n, float *d__, float *e, float *vl, float *vu,
                 integer *il, integer *iu, float *abstol, integer *m, float *w, float *z__,
                 integer *ldz, integer *isuppz, float *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int sstein_check(integer *n, float *d__, float *e, integer *m, float *w, integer *iblock,
                 integer *isplit, float *z__, integer *ldz, float *work, integer *iwork,
                 integer *ifail, integer *info);
int sstemr_check(char *jobz, char *range, integer *n, float *d__, float *e, float *vl, float *vu,
                 integer *il, integer *iu, integer *m, float *w, float *z__, integer *ldz,
                 integer *nzc, integer *isuppz, logical *tryrac, float *work, integer *lwork,
                 integer *iwork, integer *liwork, integer *info);
int ssteqr_check(char *compz, integer *n, float *d__, float *e, float *z__, integer *ldz,
                 float *work, integer *info);
int ssterf_check(integer *n, float *d__, float *e, integer *info);
int sstev_check(char *jobz, integer *n, float *d__, float *e, float *z__, integer *ldz, float *work,
                integer *info);
int sstevd_check(char *jobz, integer *n, float *d__, float *e, float *z__, integer *ldz,
                 float *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int sstevr_check(char *jobz, char *range, integer *n, float *d__, float *e, float *vl, float *vu,
                 integer *il, integer *iu, float *abstol, integer *m, float *w, float *z__,
                 integer *ldz, integer *isuppz, float *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int sstevx_check(char *jobz, char *range, integer *n, float *d__, float *e, float *vl, float *vu,
                 integer *il, integer *iu, float *abstol, integer *m, float *w, float *z__,
                 integer *ldz, float *work, integer *iwork, integer *ifail, integer *info);
int ssycon_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *anorm,
                 float *rcond, float *work, integer *iwork, integer *info);
int ssycon_rook_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *anorm,
                      float *rcond, float *work, integer *iwork, integer *info);
int ssyconv_check(char *uplo, char *way, integer *n, float *a, integer *lda, integer *ipiv,
                  float *work, integer *info);
int ssyequb_check(char *uplo, integer *n, float *a, integer *lda, float *s, float *scond,
                  float *amax, float *work, integer *info);
int ssyev_check(char *jobz, char *uplo, integer *n, float *a, integer *lda, float *w, float *work,
                integer *lwork, integer *info);
int ssyevd_check(char *jobz, char *uplo, integer *n, float *a, integer *lda, float *w, float *work,
                 integer *lwork, integer *iwork, integer *liwork, integer *info);
int ssyevr_check(char *jobz, char *range, char *uplo, integer *n, float *a, integer *lda, float *vl,
                 float *vu, integer *il, integer *iu, float *abstol, integer *m, float *w,
                 float *z__, integer *ldz, integer *isuppz, float *work, integer *lwork,
                 integer *iwork, integer *liwork, integer *info);
int ssyevx_check(char *jobz, char *range, char *uplo, integer *n, float *a, integer *lda, float *vl,
                 float *vu, integer *il, integer *iu, float *abstol, integer *m, float *w,
                 float *z__, integer *ldz, float *work, integer *lwork, integer *iwork,
                 integer *ifail, integer *info);
int ssygs2_check(integer *itype, char *uplo, integer *n, float *a, integer *lda, float *b,
                 integer *ldb, integer *info);
int ssygst_check(integer *itype, char *uplo, integer *n, float *a, integer *lda, float *b,
                 integer *ldb, integer *info);
int ssygv_check(integer *itype, char *jobz, char *uplo, integer *n, float *a, integer *lda,
                float *b, integer *ldb, float *w, float *work, integer *lwork, integer *info);
int ssygvd_check(integer *itype, char *jobz, char *uplo, integer *n, float *a, integer *lda,
                 float *b, integer *ldb, float *w, float *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int ssygvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, float *a,
                 integer *lda, float *b, integer *ldb, float *vl, float *vu, integer *il,
                 integer *iu, float *abstol, integer *m, float *w, float *z__, integer *ldz,
                 float *work, integer *lwork, integer *iwork, integer *ifail, integer *info);
int ssyrfs_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, float *af,
                 integer *ldaf, integer *ipiv, float *b, integer *ldb, float *x, integer *ldx,
                 float *ferr, float *berr, float *work, integer *iwork, integer *info);
int ssyrfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, float *a, integer *lda,
                  float *af, integer *ldaf, integer *ipiv, float *s, float *b, integer *ldb,
                  float *x, integer *ldx, float *rcond, float *berr, integer *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, integer *nparams, float *params,
                  float *work, integer *iwork, integer *info);
int ssysv_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, integer *ipiv,
                float *b, integer *ldb, float *work, integer *lwork, integer *info);
int ssysv_rook_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, integer *ipiv,
                     float *b, integer *ldb, float *work, integer *lwork, integer *info);
int ssysvx_check(char *fact, char *uplo, integer *n, integer *nrhs, float *a, integer *lda,
                 float *af, integer *ldaf, integer *ipiv, float *b, integer *ldb, float *x,
                 integer *ldx, float *rcond, float *ferr, float *berr, float *work, integer *lwork,
                 integer *iwork, integer *info);
int ssysvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, float *a, integer *lda,
                  float *af, integer *ldaf, integer *ipiv, char *equed, float *s, float *b,
                  integer *ldb, float *x, integer *ldx, float *rcond, float *rpvgrw, float *berr,
                  integer *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  integer *nparams, float *params, float *work, integer *iwork, integer *info);
int ssyswapr_check(char *uplo, integer *n, float *a, integer *lda, integer *i1, integer *i2);
int ssytd2_check(char *uplo, integer *n, float *a, integer *lda, float *d__, float *e, float *tau,
                 integer *info);
int ssytf2_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, integer *info);
int ssytf2_rook_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, integer *info);
int ssytrd_check(char *uplo, integer *n, float *a, integer *lda, float *d__, float *e, float *tau,
                 float *work, integer *lwork, integer *info);
int ssytrf_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *work,
                 integer *lwork, integer *info);
int ssytrf_rook_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *work,
                      integer *lwork, integer *info);
int ssytri_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *work,
                 integer *info);
int ssytri2_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *work,
                  integer *lwork, integer *info);
int ssytri2x_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *work,
                   integer *nb, integer *info);
int ssytri_rook_check(char *uplo, integer *n, float *a, integer *lda, integer *ipiv, float *work,
                      integer *info);
int ssytrs_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, integer *ipiv,
                 float *b, integer *ldb, integer *info);
int ssytrs2_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, integer *ipiv,
                  float *b, integer *ldb, float *work, integer *info);
int ssytrs_rook_check(char *uplo, integer *n, integer *nrhs, float *a, integer *lda, integer *ipiv,
                      float *b, integer *ldb, integer *info);
int stbcon_check(char *norm, char *uplo, char *diag, integer *n, integer *kd, float *ab,
                 integer *ldab, float *rcond, float *work, integer *iwork, integer *info);
int stbrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 float *ab, integer *ldab, float *b, integer *ldb, float *x, integer *ldx,
                 float *ferr, float *berr, float *work, integer *iwork, integer *info);
int stbtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 float *ab, integer *ldab, float *b, integer *ldb, integer *info);
int stfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m,
                integer *n, float *alpha, float *a, float *b, integer *ldb);
int stftri_check(char *transr, char *uplo, char *diag, integer *n, float *a, integer *info);
int stfttp_check(char *transr, char *uplo, integer *n, float *arf, float *ap, integer *info);
int stfttr_check(char *transr, char *uplo, integer *n, float *arf, float *a, integer *lda,
                 integer *info);
int stgevc_check(char *side, char *howmny, logical *select, integer *n, float *s, integer *lds,
                 float *p, integer *ldp, float *vl, integer *ldvl, float *vr, integer *ldvr,
                 integer *mm, integer *m, float *work, integer *info);
int stgex2_check(logical *wantq, logical *wantz, integer *n, float *a, integer *lda, float *b,
                 integer *ldb, float *q, integer *ldq, float *z__, integer *ldz, integer *j1,
                 integer *n1, integer *n2, float *work, integer *lwork, integer *info);
int stgexc_check(logical *wantq, logical *wantz, integer *n, float *a, integer *lda, float *b,
                 integer *ldb, float *q, integer *ldq, float *z__, integer *ldz, integer *ifst,
                 integer *ilst, float *work, integer *lwork, integer *info);
int stgsen_check(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n,
                 float *a, integer *lda, float *b, integer *ldb, float *alphar, float *alphai,
                 float *beta, float *q, integer *ldq, float *z__, integer *ldz, integer *m,
                 float *pl, float *pr, float *dif, float *work, integer *lwork, integer *iwork,
                 integer *liwork, integer *info);
int stgsja_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k,
                 integer *l, float *a, integer *lda, float *b, integer *ldb, float *tola,
                 float *tolb, float *alpha, float *beta, float *u, integer *ldu, float *v,
                 integer *ldv, float *q, integer *ldq, float *work, integer *ncycle, integer *info);
int stgsna_check(char *job, char *howmny, logical *select, integer *n, float *a, integer *lda,
                 float *b, integer *ldb, float *vl, integer *ldvl, float *vr, integer *ldvr,
                 float *s, float *dif, integer *mm, integer *m, float *work, integer *lwork,
                 integer *iwork, integer *info);
int stgsy2_check(char *trans, integer *ijob, integer *m, integer *n, float *a, integer *lda,
                 float *b, integer *ldb, float *c__, integer *ldc, float *d__, integer *ldd,
                 float *e, integer *lde, float *f, integer *ldf, float *scale, float *rdsum,
                 float *rdscal, integer *iwork, integer *pq, integer *info);
int stgsyl_check(char *trans, integer *ijob, integer *m, integer *n, float *a, integer *lda,
                 float *b, integer *ldb, float *c__, integer *ldc, float *d__, integer *ldd,
                 float *e, integer *lde, float *f, integer *ldf, float *scale, float *dif,
                 float *work, integer *lwork, integer *iwork, integer *info);
int stpcon_check(char *norm, char *uplo, char *diag, integer *n, float *ap, float *rcond,
                 float *work, integer *iwork, integer *info);
int stpmqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                  integer *nb, float *v, integer *ldv, float *t, integer *ldt, float *a,
                  integer *lda, float *b, integer *ldb, float *work, integer *info);
int stpqrt_check(integer *m, integer *n, integer *l, integer *nb, float *a, integer *lda, float *b,
                 integer *ldb, float *t, integer *ldt, float *work, integer *info);
int stpqrt2_check(integer *m, integer *n, integer *l, float *a, integer *lda, float *b,
                  integer *ldb, float *t, integer *ldt, integer *info);
int stprfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, float *v, integer *ldv, float *t, integer *ldt, float *a,
                 integer *lda, float *b, integer *ldb, float *work, integer *ldwork);
int stprfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, float *ap,
                 float *b, integer *ldb, float *x, integer *ldx, float *ferr, float *berr,
                 float *work, integer *iwork, integer *info);
int stptri_check(char *uplo, char *diag, integer *n, float *ap, integer *info);
int stptrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, float *ap,
                 float *b, integer *ldb, integer *info);
int stpttf_check(char *transr, char *uplo, integer *n, float *ap, float *arf, integer *info);
int stpttr_check(char *uplo, integer *n, float *ap, float *a, integer *lda, integer *info);
int strcon_check(char *norm, char *uplo, char *diag, integer *n, float *a, integer *lda,
                 float *rcond, float *work, integer *iwork, integer *info);
int strevc_check(char *side, char *howmny, logical *select, integer *n, float *t, integer *ldt,
                 float *vl, integer *ldvl, float *vr, integer *ldvr, integer *mm, integer *m,
                 float *work, integer *info);
int strexc_check(char *compq, integer *n, float *t, integer *ldt, float *q, integer *ldq,
                 integer *ifst, integer *ilst, float *work, integer *info);
int strrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, float *a,
                 integer *lda, float *b, integer *ldb, float *x, integer *ldx, float *ferr,
                 float *berr, float *work, integer *iwork, integer *info);
int strsen_check(char *job, char *compq, logical *select, integer *n, float *t, integer *ldt,
                 float *q, integer *ldq, float *wr, float *wi, integer *m, float *s, float *sep,
                 float *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
int strsna_check(char *job, char *howmny, logical *select, integer *n, float *t, integer *ldt,
                 float *vl, integer *ldvl, float *vr, integer *ldvr, float *s, float *sep,
                 integer *mm, integer *m, float *work, integer *ldwork, integer *iwork,
                 integer *info);
int strsyl_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, float *a,
                 integer *lda, float *b, integer *ldb, float *c__, integer *ldc, float *scale,
                 integer *info);
int strsyl3_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, real *a,
                  integer *lda, real *b, integer *ldb, real *c__, integer *ldc, real *scale,
                  integer *iwork, integer *liwork, real *swork, integer *ldswork, integer *info);
int strti2_check(char *uplo, char *diag, integer *n, float *a, integer *lda, integer *info);
int strtri_check(char *uplo, char *diag, integer *n, float *a, integer *lda, integer *info);
int strtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, float *a,
                 integer *lda, float *b, integer *ldb, integer *info);
int strttf_check(char *transr, char *uplo, integer *n, float *a, integer *lda, float *arf,
                 integer *info);
int strttp_check(char *uplo, integer *n, float *a, integer *lda, float *ap, integer *info);
int stzrqf_check(integer *m, integer *n, float *a, integer *lda, float *tau, integer *info);
int stzrzf_check(integer *m, integer *n, float *a, integer *lda, float *tau, float *work,
                 integer *lwork, integer *info);
int zbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, integer *m,
                 integer *p, integer *q, double *theta, double *phi, dcomplex *u1, integer *ldu1,
                 dcomplex *u2, integer *ldu2, dcomplex *v1t, integer *ldv1t, dcomplex *v2t,
                 integer *ldv2t, double *b11d, double *b11e, double *b12d, double *b12e,
                 double *b21d, double *b21e, double *b22d, double *b22e, double *rwork,
                 integer *lrwork, integer *info);
int zbdsqr_check(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, double *d__,
                 double *e, dcomplex *vt, integer *ldvt, dcomplex *u, integer *ldu, dcomplex *c__,
                 integer *ldc, double *rwork, integer *info);
int zcgesv_check(integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv, dcomplex *b,
                 integer *ldb, dcomplex *x, integer *ldx, dcomplex *work, scomplex *swork,
                 double *rwork, integer *iter, integer *info);
int zcposv_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, dcomplex *x, integer *ldx, dcomplex *work, scomplex *swork,
                 double *rwork, integer *iter, integer *info);
int zdrscl_check(integer *n, double *sa, dcomplex *sx, integer *incx);
int zgbbrd_check(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku,
                 dcomplex *ab, integer *ldab, double *d__, double *e, dcomplex *q, integer *ldq,
                 dcomplex *pt, integer *ldpt, dcomplex *c__, integer *ldc, dcomplex *work,
                 double *rwork, integer *info);
int zgbcon_check(char *norm, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab,
                 integer *ipiv, double *anorm, double *rcond, dcomplex *work, double *rwork,
                 integer *info);
int zgbequ_check(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 integer *info);
int zgbequb_check(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab,
                  double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                  integer *info);
int zgbrfs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, dcomplex *ab,
                 integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv, dcomplex *b,
                 integer *ldb, dcomplex *x, integer *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int zgbrfsx_check(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs,
                  dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv,
                  double *r__, double *c__, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx,
                  double *rcond, double *berr, integer *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, integer *nparams, double *params, dcomplex *work,
                  double *rwork, integer *info);
int zgbsv_check(integer *n, integer *kl, integer *ku, integer *nrhs, dcomplex *ab, integer *ldab,
                integer *ipiv, dcomplex *b, integer *ldb, integer *info);
int zgbsvx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                 dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv,
                 char *equed, double *r__, double *c__, dcomplex *b, integer *ldb, dcomplex *x,
                 integer *ldx, double *rcond, double *ferr, double *berr, dcomplex *work,
                 double *rwork, integer *info);
int zgbsvxx_check(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs,
                  dcomplex *ab, integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv,
                  char *equed, double *r__, double *c__, dcomplex *b, integer *ldb, dcomplex *x,
                  integer *ldx, double *rcond, double *rpvgrw, double *berr, integer *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams,
                  double *params, dcomplex *work, double *rwork, integer *info);
int zgbtf2_check(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab,
                 integer *ipiv, integer *info);
int zgbtrf_check(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab,
                 integer *ipiv, integer *info);
int zgbtrs_check(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, dcomplex *ab,
                 integer *ldab, integer *ipiv, dcomplex *b, integer *ldb, integer *info);
int zgebak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, double *scale,
                 integer *m, dcomplex *v, integer *ldv, integer *info);
int zgebal_check(char *job, integer *n, dcomplex *a, integer *lda, integer *ilo, integer *ihi,
                 double *scale, integer *info);
int zgebd2_check(integer *m, integer *n, dcomplex *a, integer *lda, double *d__, double *e,
                 dcomplex *tauq, dcomplex *taup, dcomplex *work, integer *info);
int zgebrd_check(integer *m, integer *n, dcomplex *a, integer *lda, double *d__, double *e,
                 dcomplex *tauq, dcomplex *taup, dcomplex *work, integer *lwork, integer *info);
int zgecon_check(char *norm, integer *n, dcomplex *a, integer *lda, double *anorm, double *rcond,
                 dcomplex *work, double *rwork, integer *info);
int zgeequ_check(integer *m, integer *n, dcomplex *a, integer *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, integer *info);
int zgeequb_check(integer *m, integer *n, dcomplex *a, integer *lda, double *r__, double *c__,
                  double *rowcnd, double *colcnd, double *amax, integer *info);
int zgees_check(char *jobvs, char *sort, L_fp select, integer *n, dcomplex *a, integer *lda,
                integer *sdim, dcomplex *w, dcomplex *vs, integer *ldvs, dcomplex *work,
                integer *lwork, double *rwork, logical *bwork, integer *info);
int zgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, integer *n, dcomplex *a,
                 integer *lda, integer *sdim, dcomplex *w, dcomplex *vs, integer *ldvs,
                 double *rconde, double *rcondv, dcomplex *work, integer *lwork, double *rwork,
                 logical *bwork, integer *info);
int zgeev_check(char *jobvl, char *jobvr, integer *n, dcomplex *a, integer *lda, dcomplex *w,
                dcomplex *vl, integer *ldvl, dcomplex *vr, integer *ldvr, dcomplex *work,
                integer *lwork, double *rwork, integer *info);
int zgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, dcomplex *a,
                 integer *lda, dcomplex *w, dcomplex *vl, integer *ldvl, dcomplex *vr,
                 integer *ldvr, integer *ilo, integer *ihi, double *scale, double *abnrm,
                 double *rconde, double *rcondv, dcomplex *work, integer *lwork, double *rwork,
                 integer *info);
int zgegs_check(char *jobvsl, char *jobvsr, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                integer *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vsl, integer *ldvsl,
                dcomplex *vsr, integer *ldvsr, dcomplex *work, integer *lwork, double *rwork,
                integer *info);
int zgegv_check(char *jobvl, char *jobvr, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                integer *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vl, integer *ldvl,
                dcomplex *vr, integer *ldvr, dcomplex *work, integer *lwork, double *rwork,
                integer *info);
int zgehd2_check(integer *n, integer *ilo, integer *ihi, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *info);
int zgehrd_check(integer *n, integer *ilo, integer *ihi, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *lwork, integer *info);
int zgelq2_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *info);
int zgelqf_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *lwork, integer *info);
int zgels_check(char *trans, integer *m, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                dcomplex *b, integer *ldb, dcomplex *work, integer *lwork, integer *info);
int zgelsd_check(integer *m, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, double *s, double *rcond, integer *rank, dcomplex *work,
                 integer *lwork, double *rwork, integer *iwork, integer *info);
int zgelss_check(integer *m, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, double *s, double *rcond, integer *rank, dcomplex *work,
                 integer *lwork, double *rwork, integer *info);
int zgelsx_check(integer *m, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, integer *jpvt, double *rcond, integer *rank, dcomplex *work,
                 double *rwork, integer *info);
int zgelsy_check(integer *m, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, integer *jpvt, double *rcond, integer *rank, dcomplex *work,
                 integer *lwork, double *rwork, integer *info);
int zgemqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *nb,
                  dcomplex *v, integer *ldv, dcomplex *t, integer *ldt, dcomplex *c__, integer *ldc,
                  dcomplex *work, integer *info);
int zgeql2_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *info);
int zgeqlf_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *lwork, integer *info);
int zgeqp3_check(integer *m, integer *n, dcomplex *a, integer *lda, integer *jpvt, dcomplex *tau,
                 dcomplex *work, integer *lwork, double *rwork, integer *info);
int zgeqpf_check(integer *m, integer *n, dcomplex *a, integer *lda, integer *jpvt, dcomplex *tau,
                 dcomplex *work, double *rwork, integer *info);
int zgeqr2_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *info);
int zgeqr2p_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                  integer *info);
int zgeqrf_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *lwork, integer *info);
int zgeqrfp_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                  integer *lwork, integer *info);
int zgeqrt_check(integer *m, integer *n, integer *nb, dcomplex *a, integer *lda, dcomplex *t,
                 integer *ldt, dcomplex *work, integer *info);
int zgeqrt2_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *t, integer *ldt,
                  integer *info);
int zgeqrt3_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *t, integer *ldt,
                  integer *info);
int zgerfs_check(char *trans, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af,
                 integer *ldaf, integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int zgerfsx_check(char *trans, char *equed, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, integer *ipiv, double *r__, double *c__, dcomplex *b,
                  integer *ldb, dcomplex *x, integer *ldx, double *rcond, double *berr,
                  integer *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  integer *nparams, double *params, dcomplex *work, double *rwork, integer *info);
int zgerq2_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *info);
int zgerqf_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *lwork, integer *info);
int zgesc2_check(integer *n, dcomplex *a, integer *lda, dcomplex *rhs, integer *ipiv, integer *jpiv,
                 double *scale);
int zgesdd_check(char *jobz, integer *m, integer *n, dcomplex *a, integer *lda, double *s,
                 dcomplex *u, integer *ldu, dcomplex *vt, integer *ldvt, dcomplex *work,
                 integer *lwork, double *rwork, integer *iwork, integer *info);
int zgesv_check(integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv, dcomplex *b,
                integer *ldb, integer *info);
int zgesvd_check(char *jobu, char *jobvt, integer *m, integer *n, dcomplex *a, integer *lda,
                 double *s, dcomplex *u, integer *ldu, dcomplex *vt, integer *ldvt, dcomplex *work,
                 integer *lwork, double *rwork, integer *info);
int zgesvx_check(char *fact, char *trans, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                 dcomplex *af, integer *ldaf, integer *ipiv, char *equed, double *r__, double *c__,
                 dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *rcond, double *ferr,
                 double *berr, dcomplex *work, double *rwork, integer *info);
int zgesvxx_check(char *fact, char *trans, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, integer *ipiv, char *equed, double *r__, double *c__,
                  dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *rcond,
                  double *rpvgrw, double *berr, integer *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, integer *nparams, double *params, dcomplex *work,
                  double *rwork, integer *info);
int zgetc2_check(integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *jpiv,
                 integer *info);
int zgetf2_check(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info);
int zgetrf_check(integer *m, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info);
int zgetrfnp_check(integer *m, integer *n, dcomplex *a, integer *lda, integer *info);
int zgetrfnpi_check(integer *m, integer *n, integer *nfact, dcomplex *a, integer *lda,
                    integer *info);
int zgetri_check(integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                 integer *lwork, integer *info);
int zgetrs_check(char *trans, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv,
                 dcomplex *b, integer *ldb, integer *info);
int zggbak_check(char *job, char *side, integer *n, integer *ilo, integer *ihi, double *lscale,
                 double *rscale, integer *m, dcomplex *v, integer *ldv, integer *info);
int zggbal_check(char *job, integer *n, dcomplex *a, integer *lda, dcomplex *b, integer *ldb,
                 integer *ilo, integer *ihi, double *lscale, double *rscale, double *work,
                 integer *info);
int zgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, dcomplex *a,
                integer *lda, dcomplex *b, integer *ldb, integer *sdim, dcomplex *alpha,
                dcomplex *beta, dcomplex *vsl, integer *ldvsl, dcomplex *vsr, integer *ldvsr,
                dcomplex *work, integer *lwork, double *rwork, logical *bwork, integer *info);
int zggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n,
                 dcomplex *a, integer *lda, dcomplex *b, integer *ldb, integer *sdim,
                 dcomplex *alpha, dcomplex *beta, dcomplex *vsl, integer *ldvsl, dcomplex *vsr,
                 integer *ldvsr, double *rconde, double *rcondv, dcomplex *work, integer *lwork,
                 double *rwork, integer *iwork, integer *liwork, logical *bwork, integer *info);
int zggev_check(char *jobvl, char *jobvr, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                integer *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vl, integer *ldvl,
                dcomplex *vr, integer *ldvr, dcomplex *work, integer *lwork, double *rwork,
                integer *info);
int zggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, dcomplex *a,
                 integer *lda, dcomplex *b, integer *ldb, dcomplex *alpha, dcomplex *beta,
                 dcomplex *vl, integer *ldvl, dcomplex *vr, integer *ldvr, integer *ilo,
                 integer *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm,
                 double *rconde, double *rcondv, dcomplex *work, integer *lwork, double *rwork,
                 integer *iwork, logical *bwork, integer *info);
int zggglm_check(integer *n, integer *m, integer *p, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, dcomplex *d__, dcomplex *x, dcomplex *y, dcomplex *work,
                 integer *lwork, integer *info);
int zgghrd_check(char *compq, char *compz, integer *n, integer *ilo, integer *ihi, dcomplex *a,
                 integer *lda, dcomplex *b, integer *ldb, dcomplex *q, integer *ldq, dcomplex *z__,
                 integer *ldz, integer *info);
int zgglse_check(integer *m, integer *n, integer *p, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, dcomplex *c__, dcomplex *d__, dcomplex *x, dcomplex *work,
                 integer *lwork, integer *info);
int zggqrf_check(integer *n, integer *m, integer *p, dcomplex *a, integer *lda, dcomplex *taua,
                 dcomplex *b, integer *ldb, dcomplex *taub, dcomplex *work, integer *lwork,
                 integer *info);
int zggrqf_check(integer *m, integer *p, integer *n, dcomplex *a, integer *lda, dcomplex *taua,
                 dcomplex *b, integer *ldb, dcomplex *taub, dcomplex *work, integer *lwork,
                 integer *info);
int zggsvd_check(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k,
                 integer *l, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, double *alpha,
                 double *beta, dcomplex *u, integer *ldu, dcomplex *v, integer *ldv, dcomplex *q,
                 integer *ldq, dcomplex *work, double *rwork, integer *iwork, integer *info);
int zggsvp_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n,
                 dcomplex *a, integer *lda, dcomplex *b, integer *ldb, double *tola, double *tolb,
                 integer *k, integer *l, dcomplex *u, integer *ldu, dcomplex *v, integer *ldv,
                 dcomplex *q, integer *ldq, integer *iwork, double *rwork, dcomplex *tau,
                 dcomplex *work, integer *info);
int zgtcon_check(char *norm, integer *n, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2,
                 integer *ipiv, double *anorm, double *rcond, dcomplex *work, integer *info);
int zgtrfs_check(char *trans, integer *n, integer *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du,
                 dcomplex *dlf, dcomplex *df, dcomplex *duf, dcomplex *du2, integer *ipiv,
                 dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int zgtsv_check(integer *n, integer *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *b,
                integer *ldb, integer *info);
int zgtsvx_check(char *fact, char *trans, integer *n, integer *nrhs, dcomplex *dl, dcomplex *d__,
                 dcomplex *du, dcomplex *dlf, dcomplex *df, dcomplex *duf, dcomplex *du2,
                 integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int zgttrf_check(integer *n, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2,
                 integer *ipiv, integer *info);
int zgttrs_check(char *trans, integer *n, integer *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du,
                 dcomplex *du2, integer *ipiv, dcomplex *b, integer *ldb, integer *info);
int zgtts2_check(integer *itrans, integer *n, integer *nrhs, dcomplex *dl, dcomplex *d__,
                 dcomplex *du, dcomplex *du2, integer *ipiv, dcomplex *b, integer *ldb);
int zhbev_check(char *jobz, char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab,
                double *w, dcomplex *z__, integer *ldz, dcomplex *work, double *rwork,
                integer *info);
int zhbevd_check(char *jobz, char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab,
                 double *w, dcomplex *z__, integer *ldz, dcomplex *work, integer *lwork,
                 double *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);
int zhbevx_check(char *jobz, char *range, char *uplo, integer *n, integer *kd, dcomplex *ab,
                 integer *ldab, dcomplex *q, integer *ldq, double *vl, double *vu, integer *il,
                 integer *iu, double *abstol, integer *m, double *w, dcomplex *z__, integer *ldz,
                 dcomplex *work, double *rwork, integer *iwork, integer *ifail, integer *info);
int zhbgst_check(char *vect, char *uplo, integer *n, integer *ka, integer *kb, dcomplex *ab,
                 integer *ldab, dcomplex *bb, integer *ldbb, dcomplex *x, integer *ldx,
                 dcomplex *work, double *rwork, integer *info);
int zhbgv_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, dcomplex *ab,
                integer *ldab, dcomplex *bb, integer *ldbb, double *w, dcomplex *z__, integer *ldz,
                dcomplex *work, double *rwork, integer *info);
int zhbgvd_check(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, dcomplex *ab,
                 integer *ldab, dcomplex *bb, integer *ldbb, double *w, dcomplex *z__, integer *ldz,
                 dcomplex *work, integer *lwork, double *rwork, integer *lrwork, integer *iwork,
                 integer *liwork, integer *info);
int zhbgvx_check(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb,
                 dcomplex *ab, integer *ldab, dcomplex *bb, integer *ldbb, dcomplex *q,
                 integer *ldq, double *vl, double *vu, integer *il, integer *iu, double *abstol,
                 integer *m, double *w, dcomplex *z__, integer *ldz, dcomplex *work, double *rwork,
                 integer *iwork, integer *ifail, integer *info);
int zhbtrd_check(char *vect, char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab,
                 double *d__, double *e, dcomplex *q, integer *ldq, dcomplex *work, integer *info);
int zhecon_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, double *anorm,
                 double *rcond, dcomplex *work, integer *info);
int zhecon_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      double *anorm, double *rcond, dcomplex *work, integer *info);
int zheequb_check(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond,
                  double *amax, dcomplex *work, integer *info);
int zheev_check(char *jobz, char *uplo, integer *n, dcomplex *a, integer *lda, double *w,
                dcomplex *work, integer *lwork, double *rwork, integer *info);
int zheevd_check(char *jobz, char *uplo, integer *n, dcomplex *a, integer *lda, double *w,
                 dcomplex *work, integer *lwork, double *rwork, integer *lrwork, integer *iwork,
                 integer *liwork, integer *info);
int zheevr_check(char *jobz, char *range, char *uplo, integer *n, dcomplex *a, integer *lda,
                 double *vl, double *vu, integer *il, integer *iu, double *abstol, integer *m,
                 double *w, dcomplex *z__, integer *ldz, integer *isuppz, dcomplex *work,
                 integer *lwork, double *rwork, integer *lrwork, integer *iwork, integer *liwork,
                 integer *info);
int zheevx_check(char *jobz, char *range, char *uplo, integer *n, dcomplex *a, integer *lda,
                 double *vl, double *vu, integer *il, integer *iu, double *abstol, integer *m,
                 double *w, dcomplex *z__, integer *ldz, dcomplex *work, integer *lwork,
                 double *rwork, integer *iwork, integer *ifail, integer *info);
int zhegs2_check(integer *itype, char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, integer *info);
int zhegst_check(integer *itype, char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, integer *info);
int zhegv_check(integer *itype, char *jobz, char *uplo, integer *n, dcomplex *a, integer *lda,
                dcomplex *b, integer *ldb, double *w, dcomplex *work, integer *lwork, double *rwork,
                integer *info);
int zhegvd_check(integer *itype, char *jobz, char *uplo, integer *n, dcomplex *a, integer *lda,
                 dcomplex *b, integer *ldb, double *w, dcomplex *work, integer *lwork,
                 double *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);
int zhegvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, dcomplex *a,
                 integer *lda, dcomplex *b, integer *ldb, double *vl, double *vu, integer *il,
                 integer *iu, double *abstol, integer *m, double *w, dcomplex *z__, integer *ldz,
                 dcomplex *work, integer *lwork, double *rwork, integer *iwork, integer *ifail,
                 integer *info);
int zherfs_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af,
                 integer *ldaf, integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int zherfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, integer *ipiv, double *s, dcomplex *b, integer *ldb,
                  dcomplex *x, integer *ldx, double *rcond, double *berr, integer *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams,
                  double *params, dcomplex *work, double *rwork, integer *info);
int zhesv_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv,
                dcomplex *b, integer *ldb, dcomplex *work, integer *lwork, integer *info);
int zhesv_rook_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                     integer *ipiv, dcomplex *b, integer *ldb, dcomplex *work, integer *lwork,
                     integer *info);
int zhesvx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                 dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x,
                 integer *ldx, double *rcond, double *ferr, double *berr, dcomplex *work,
                 integer *lwork, double *rwork, integer *info);
int zhesvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, integer *ipiv, char *equed, double *s, dcomplex *b,
                  integer *ldb, dcomplex *x, integer *ldx, double *rcond, double *rpvgrw,
                  double *berr, integer *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, integer *nparams, double *params, dcomplex *work,
                  double *rwork, integer *info);
int zheswapr_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *i1, integer *i2);
int zhetd2_check(char *uplo, integer *n, dcomplex *a, integer *lda, double *d__, double *e,
                 dcomplex *tau, integer *info);
int zhetf2_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info);
int zhetf2_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      integer *info);
int zhetrd_check(char *uplo, integer *n, dcomplex *a, integer *lda, double *d__, double *e,
                 dcomplex *tau, dcomplex *work, integer *lwork, integer *info);
int zhetrf_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                 integer *lwork, integer *info);
int zhetrf_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      dcomplex *work, integer *lwork, integer *info);
int zhetri_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                 integer *info);
int zhetri2_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                  integer *lwork, integer *info);
int zhetri2x_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                   integer *nb, integer *info);
int zhetri_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      dcomplex *work, integer *info);
int zhetrs_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv,
                 dcomplex *b, integer *ldb, integer *info);
int zhetrs2_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv,
                  dcomplex *b, integer *ldb, dcomplex *work, integer *info);
int zhetrs_rook_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                      integer *ipiv, dcomplex *b, integer *ldb, integer *info);
int zhfrk_check(char *transr, char *uplo, char *trans, integer *n, integer *k, double *alpha,
                dcomplex *a, integer *lda, double *beta, dcomplex *c__);
int zhgeqz_check(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi,
                 dcomplex *h__, integer *ldh, dcomplex *t, integer *ldt, dcomplex *alpha,
                 dcomplex *beta, dcomplex *q, integer *ldq, dcomplex *z__, integer *ldz,
                 dcomplex *work, integer *lwork, double *rwork, integer *info);
int zhpcon_check(char *uplo, integer *n, dcomplex *ap, integer *ipiv, double *anorm, double *rcond,
                 dcomplex *work, integer *info);
int zhpev_check(char *jobz, char *uplo, integer *n, dcomplex *ap, double *w, dcomplex *z__,
                integer *ldz, dcomplex *work, double *rwork, integer *info);
int zhpevd_check(char *jobz, char *uplo, integer *n, dcomplex *ap, double *w, dcomplex *z__,
                 integer *ldz, dcomplex *work, integer *lwork, double *rwork, integer *lrwork,
                 integer *iwork, integer *liwork, integer *info);
int zhpevx_check(char *jobz, char *range, char *uplo, integer *n, dcomplex *ap, double *vl,
                 double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w,
                 dcomplex *z__, integer *ldz, dcomplex *work, double *rwork, integer *iwork,
                 integer *ifail, integer *info);
int zhpgst_check(integer *itype, char *uplo, integer *n, dcomplex *ap, dcomplex *bp, integer *info);
int zhpgv_check(integer *itype, char *jobz, char *uplo, integer *n, dcomplex *ap, dcomplex *bp,
                double *w, dcomplex *z__, integer *ldz, dcomplex *work, double *rwork,
                integer *info);
int zhpgvd_check(integer *itype, char *jobz, char *uplo, integer *n, dcomplex *ap, dcomplex *bp,
                 double *w, dcomplex *z__, integer *ldz, dcomplex *work, integer *lwork,
                 double *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);
int zhpgvx_check(integer *itype, char *jobz, char *range, char *uplo, integer *n, dcomplex *ap,
                 dcomplex *bp, double *vl, double *vu, integer *il, integer *iu, double *abstol,
                 integer *m, double *w, dcomplex *z__, integer *ldz, dcomplex *work, double *rwork,
                 integer *iwork, integer *ifail, integer *info);
int zhprfs_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *afp, integer *ipiv,
                 dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int zhpsv_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, integer *ipiv, dcomplex *b,
                integer *ldb, integer *info);
int zhpsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *afp,
                 integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int zhptrd_check(char *uplo, integer *n, dcomplex *ap, double *d__, double *e, dcomplex *tau,
                 integer *info);
int zhptrf_check(char *uplo, integer *n, dcomplex *ap, integer *ipiv, integer *info);
int zhptri_check(char *uplo, integer *n, dcomplex *ap, integer *ipiv, dcomplex *work,
                 integer *info);
int zhptrs_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, integer *ipiv, dcomplex *b,
                 integer *ldb, integer *info);
int zhsein_check(char *side, char *eigsrc, char *initv, logical *select, integer *n, dcomplex *h__,
                 integer *ldh, dcomplex *w, dcomplex *vl, integer *ldvl, dcomplex *vr,
                 integer *ldvr, integer *mm, integer *m, dcomplex *work, double *rwork,
                 integer *ifaill, integer *ifailr, integer *info);
int zhseqr_check(char *job, char *compz, integer *n, integer *ilo, integer *ihi, dcomplex *h__,
                 integer *ldh, dcomplex *w, dcomplex *z__, integer *ldz, dcomplex *work,
                 integer *lwork, integer *info);
int zla_gbamv_check(integer *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha,
                    dcomplex *ab, integer *ldab, dcomplex *x, integer *incx, double *beta,
                    double *y, integer *incy);
double zla_gbrcond_c_check(char *trans, integer *n, integer *kl, integer *ku, dcomplex *ab,
                           integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv, double *c__,
                           logical *capply, integer *info, dcomplex *work, double *rwork);
double zla_gbrcond_x_check(char *trans, integer *n, integer *kl, integer *ku, dcomplex *ab,
                           integer *ldab, dcomplex *afb, integer *ldafb, integer *ipiv, dcomplex *x,
                           integer *info, dcomplex *work, double *rwork);
int zla_gbrfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n, integer *kl,
                              integer *ku, integer *nrhs, dcomplex *ab, integer *ldab,
                              dcomplex *afb, integer *ldafb, integer *ipiv, logical *colequ,
                              double *c__, dcomplex *b, integer *ldb, dcomplex *y, integer *ldy,
                              double *berr_out__, integer *n_norms__, double *err_bnds_norm__,
                              double *err_bnds_comp__, dcomplex *res, double *ayb, dcomplex *dy,
                              dcomplex *y_tail__, double *rcond, integer *ithresh, double *rthresh,
                              double *dz_ub__, logical *ignore_cwise__, integer *info);
double zla_gbrpvgrw_check(integer *n, integer *kl, integer *ku, integer *ncols, dcomplex *ab,
                          integer *ldab, dcomplex *afb, integer *ldafb);
int zla_geamv_check(integer *trans, integer *m, integer *n, double *alpha, dcomplex *a,
                    integer *lda, dcomplex *x, integer *incx, double *beta, double *y,
                    integer *incy);
double zla_gercond_c_check(char *trans, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, integer *ipiv, double *c__, logical *capply,
                           integer *info, dcomplex *work, double *rwork);
double zla_gercond_x_check(char *trans, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work,
                           double *rwork);
int zla_gerfsx_extended_check(integer *prec_type__, integer *trans_type__, integer *n,
                              integer *nrhs, dcomplex *a, integer *lda, dcomplex *af, integer *ldaf,
                              integer *ipiv, logical *colequ, double *c__, dcomplex *b,
                              integer *ldb, dcomplex *y, integer *ldy, double *berr_out__,
                              integer *n_norms__, double *errs_n__, double *errs_c__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              integer *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, integer *info);
double zla_gerpvgrw_check(integer *n, integer *ncols, dcomplex *a, integer *lda, dcomplex *af,
                          integer *ldaf);
int zla_heamv_check(integer *uplo, integer *n, double *alpha, dcomplex *a, integer *lda,
                    dcomplex *x, integer *incx, double *beta, double *y, integer *incy);
double zla_hercond_c_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, integer *ipiv, double *c__, logical *capply,
                           integer *info, dcomplex *work, double *rwork);
double zla_hercond_x_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work,
                           double *rwork);
int zla_herfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv,
                              logical *colequ, double *c__, dcomplex *b, integer *ldb, dcomplex *y,
                              integer *ldy, double *berr_out__, integer *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              integer *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, integer *info);
double zla_herpvgrw_check(char *uplo, integer *n, integer *info, dcomplex *a, integer *lda,
                          dcomplex *af, integer *ldaf, integer *ipiv, double *work);
int zla_lin_berr_check(integer *n, integer *nz, integer *nrhs, dcomplex *res, double *ayb,
                       double *berr);
double zla_porcond_c_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, double *c__, logical *capply, integer *info,
                           dcomplex *work, double *rwork);
double zla_porcond_x_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, dcomplex *x, integer *info, dcomplex *work,
                           double *rwork);
int zla_porfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              dcomplex *a, integer *lda, dcomplex *af, integer *ldaf,
                              logical *colequ, double *c__, dcomplex *b, integer *ldb, dcomplex *y,
                              integer *ldy, double *berr_out__, integer *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              integer *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, integer *info);
double zla_porpvgrw_check(char *uplo, integer *ncols, dcomplex *a, integer *lda, dcomplex *af,
                          integer *ldaf, double *work);
int zla_syamv_check(integer *uplo, integer *n, double *alpha, dcomplex *a, integer *lda,
                    dcomplex *x, integer *incx, double *beta, double *y, integer *incy);
double zla_syrcond_c_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, integer *ipiv, double *c__, logical *capply,
                           integer *info, dcomplex *work, double *rwork);
double zla_syrcond_x_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *af,
                           integer *ldaf, integer *ipiv, dcomplex *x, integer *info, dcomplex *work,
                           double *rwork);
int zla_syrfsx_extended_check(integer *prec_type__, char *uplo, integer *n, integer *nrhs,
                              dcomplex *a, integer *lda, dcomplex *af, integer *ldaf, integer *ipiv,
                              logical *colequ, double *c__, dcomplex *b, integer *ldb, dcomplex *y,
                              integer *ldy, double *berr_out__, integer *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              integer *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, integer *info);
double zla_syrpvgrw_check(char *uplo, integer *n, integer *info, dcomplex *a, integer *lda,
                          dcomplex *af, integer *ldaf, integer *ipiv, double *work);
int zla_wwaddw_check(integer *n, dcomplex *x, dcomplex *y, dcomplex *w);
int zlabrd_check(integer *m, integer *n, integer *nb, dcomplex *a, integer *lda, double *d__,
                 double *e, dcomplex *tauq, dcomplex *taup, dcomplex *x, integer *ldx, dcomplex *y,
                 integer *ldy);
int zlacgv_check(integer *n, dcomplex *x, integer *incx);
int zlacn2_check(integer *n, dcomplex *v, dcomplex *x, double *est, integer *kase, integer *isave);
int zlacon_check(integer *n, dcomplex *v, dcomplex *x, double *est, integer *kase);
int zlacp2_check(char *uplo, integer *m, integer *n, double *a, integer *lda, dcomplex *b,
                 integer *ldb);
int zlacpy_check(char *uplo, integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb);
int zlacrm_check(integer *m, integer *n, dcomplex *a, integer *lda, double *b, integer *ldb,
                 dcomplex *c__, integer *ldc, double *rwork);
int zlacrt_check(integer *n, dcomplex *cx, integer *incx, dcomplex *cy, integer *incy,
                 dcomplex *c__, dcomplex *s);
VOID zladiv_check(dcomplex *ret_val, dcomplex *x, dcomplex *y);
int zlaed0_check(integer *qsiz, integer *n, double *d__, double *e, dcomplex *q, integer *ldq,
                 dcomplex *qstore, integer *ldqs, double *rwork, integer *iwork, integer *info);
int zlaed7_check(integer *n, integer *cutpnt, integer *qsiz, integer *tlvls, integer *curlvl,
                 integer *curpbm, double *d__, dcomplex *q, integer *ldq, double *rho,
                 integer *indxq, double *qstore, integer *qptr, integer *prmptr, integer *perm,
                 integer *givptr, integer *givcol, double *givnum, dcomplex *work, double *rwork,
                 integer *iwork, integer *info);
int zlaed8_check(integer *k, integer *n, integer *qsiz, dcomplex *q, integer *ldq, double *d__,
                 double *rho, integer *cutpnt, double *z__, double *dlamda, dcomplex *q2,
                 integer *ldq2, double *w, integer *indxp, integer *indx, integer *indxq,
                 integer *perm, integer *givptr, integer *givcol, double *givnum, integer *info);
int zlaein_check(logical *rightv, logical *noinit, integer *n, dcomplex *h__, integer *ldh,
                 dcomplex *w, dcomplex *v, dcomplex *b, integer *ldb, double *rwork, double *eps3,
                 double *smlnum, integer *info);
int zlaesy_check(dcomplex *a, dcomplex *b, dcomplex *c__, dcomplex *rt1, dcomplex *rt2,
                 dcomplex *evscal, dcomplex *cs1, dcomplex *sn1);
int zlaev2_check(dcomplex *a, dcomplex *b, dcomplex *c__, double *rt1, double *rt2, double *cs1,
                 dcomplex *sn1);
int zlag2c_check(integer *m, integer *n, dcomplex *a, integer *lda, scomplex *sa, integer *ldsa,
                 integer *info);
int zlags2_check(logical *upper, double *a1, dcomplex *a2, double *a3, double *b1, dcomplex *b2,
                 double *b3, double *csu, dcomplex *snu, double *csv, dcomplex *snv, double *csq,
                 dcomplex *snq);
int zlagtm_check(char *trans, integer *n, integer *nrhs, double *alpha, dcomplex *dl, dcomplex *d__,
                 dcomplex *du, dcomplex *x, integer *ldx, double *beta, dcomplex *b, integer *ldb);
int zlahef_check(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda,
                 integer *ipiv, dcomplex *w, integer *ldw, integer *info);
int zlahef_rook_check(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda,
                      integer *ipiv, dcomplex *w, integer *ldw, integer *info);
int zlahqr_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 dcomplex *h__, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz,
                 dcomplex *z__, integer *ldz, integer *info);
int zlahr2_check(integer *n, integer *k, integer *nb, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *t, integer *ldt, dcomplex *y, integer *ldy);
int zlahrd_check(integer *n, integer *k, integer *nb, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *t, integer *ldt, dcomplex *y, integer *ldy);
int zlaic1_check(integer *job, integer *j, dcomplex *x, double *sest, dcomplex *w, dcomplex *gamma,
                 double *sestpr, dcomplex *s, dcomplex *c__);
int zlals0_check(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs,
                 dcomplex *b, integer *ldb, dcomplex *bx, integer *ldbx, integer *perm,
                 integer *givptr, integer *givcol, integer *ldgcol, double *givnum, integer *ldgnum,
                 double *poles, double *difl, double *difr, double *z__, integer *k, double *c__,
                 double *s, double *rwork, integer *info);
int zlalsa_check(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, dcomplex *b,
                 integer *ldb, dcomplex *bx, integer *ldbx, double *u, integer *ldu, double *vt,
                 integer *k, double *difl, double *difr, double *z__, double *poles,
                 integer *givptr, integer *givcol, integer *ldgcol, integer *perm, double *givnum,
                 double *c__, double *s, double *rwork, integer *iwork, integer *info);
int zlalsd_check(char *uplo, integer *smlsiz, integer *n, integer *nrhs, double *d__, double *e,
                 dcomplex *b, integer *ldb, double *rcond, integer *rank, dcomplex *work,
                 double *rwork, integer *iwork, integer *info);
double zlangb_check(char *norm, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab,
                    double *work);
double zlange_check(char *norm, integer *m, integer *n, dcomplex *a, integer *lda, double *work);
double zlangt_check(char *norm, integer *n, dcomplex *dl, dcomplex *d__, dcomplex *du);
double zlanhb_check(char *norm, char *uplo, integer *n, integer *k, dcomplex *ab, integer *ldab,
                    double *work);
double zlanhe_check(char *norm, char *uplo, integer *n, dcomplex *a, integer *lda, double *work);
double zlanhf_check(char *norm, char *transr, char *uplo, integer *n, dcomplex *a, double *work);
double zlanhp_check(char *norm, char *uplo, integer *n, dcomplex *ap, double *work);
double zlanhs_check(char *norm, integer *n, dcomplex *a, integer *lda, double *work);
double zlanht_check(char *norm, integer *n, double *d__, dcomplex *e);
double zlansb_check(char *norm, char *uplo, integer *n, integer *k, dcomplex *ab, integer *ldab,
                    double *work);
double zlansp_check(char *norm, char *uplo, integer *n, dcomplex *ap, double *work);
double zlansy_check(char *norm, char *uplo, integer *n, dcomplex *a, integer *lda, double *work);
double zlantb_check(char *norm, char *uplo, char *diag, integer *n, integer *k, dcomplex *ab,
                    integer *ldab, double *work);
double zlantp_check(char *norm, char *uplo, char *diag, integer *n, dcomplex *ap, double *work);
double zlantr_check(char *norm, char *uplo, char *diag, integer *m, integer *n, dcomplex *a,
                    integer *lda, double *work);
int zlapll_check(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *ssmin);
int zlapmr_check(logical *forwrd, integer *m, integer *n, dcomplex *x, integer *ldx, integer *k);
int zlapmt_check(logical *forwrd, integer *m, integer *n, dcomplex *x, integer *ldx, integer *k);
int zlaqgb_check(integer *m, integer *n, integer *kl, integer *ku, dcomplex *ab, integer *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 char *equed);
int zlaqge_check(integer *m, integer *n, dcomplex *a, integer *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, char *equed);
int zlaqhb_check(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *s,
                 double *scond, double *amax, char *equed);
int zlaqhe_check(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond,
                 double *amax, char *equed);
int zlaqhp_check(char *uplo, integer *n, dcomplex *ap, double *s, double *scond, double *amax,
                 char *equed);
int zlaqp2_check(integer *m, integer *n, integer *offset, dcomplex *a, integer *lda, integer *jpvt,
                 dcomplex *tau, double *vn1, double *vn2, dcomplex *work);
int zlaqps_check(integer *m, integer *n, integer *offset, integer *nb, integer *kb, dcomplex *a,
                 integer *lda, integer *jpvt, dcomplex *tau, double *vn1, double *vn2,
                 dcomplex *auxv, dcomplex *f, integer *ldf);
int zlaqr0_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 dcomplex *h__, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz,
                 dcomplex *z__, integer *ldz, dcomplex *work, integer *lwork, integer *info);
int zlaqr1_check(integer *n, dcomplex *h__, integer *ldh, dcomplex *s1, dcomplex *s2, dcomplex *v);
int zlaqr2_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, dcomplex *h__, integer *ldh, integer *iloz, integer *ihiz,
                 dcomplex *z__, integer *ldz, integer *ns, integer *nd, dcomplex *sh, dcomplex *v,
                 integer *ldv, integer *nh, dcomplex *t, integer *ldt, integer *nv, dcomplex *wv,
                 integer *ldwv, dcomplex *work, integer *lwork);
int zlaqr3_check(logical *wantt, logical *wantz, integer *n, integer *ktop, integer *kbot,
                 integer *nw, dcomplex *h__, integer *ldh, integer *iloz, integer *ihiz,
                 dcomplex *z__, integer *ldz, integer *ns, integer *nd, dcomplex *sh, dcomplex *v,
                 integer *ldv, integer *nh, dcomplex *t, integer *ldt, integer *nv, dcomplex *wv,
                 integer *ldwv, dcomplex *work, integer *lwork);
int zlaqr4_check(logical *wantt, logical *wantz, integer *n, integer *ilo, integer *ihi,
                 dcomplex *h__, integer *ldh, dcomplex *w, integer *iloz, integer *ihiz,
                 dcomplex *z__, integer *ldz, dcomplex *work, integer *lwork, integer *info);
int zlaqr5_check(logical *wantt, logical *wantz, integer *kacc22, integer *n, integer *ktop,
                 integer *kbot, integer *nshfts, dcomplex *s, dcomplex *h__, integer *ldh,
                 integer *iloz, integer *ihiz, dcomplex *z__, integer *ldz, dcomplex *v,
                 integer *ldv, dcomplex *u, integer *ldu, integer *nv, dcomplex *wv, integer *ldwv,
                 integer *nh, dcomplex *wh, integer *ldwh);
int zlaqsb_check(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *s,
                 double *scond, double *amax, char *equed);
int zlaqsp_check(char *uplo, integer *n, dcomplex *ap, double *s, double *scond, double *amax,
                 char *equed);
int zlaqsy_check(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond,
                 double *amax, char *equed);
int zlar1v_check(integer *n, integer *b1, integer *bn, double *lambda, double *d__, double *l,
                 double *ld, double *lld, double *pivmin, double *gaptol, dcomplex *z__,
                 logical *wantnc, integer *negcnt, double *ztz, double *mingma, integer *r__,
                 integer *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);
int zlar2v_check(integer *n, dcomplex *x, dcomplex *y, dcomplex *z__, integer *incx, double *c__,
                 dcomplex *s, integer *incc);
int zlarcm_check(integer *m, integer *n, double *a, integer *lda, dcomplex *b, integer *ldb,
                 dcomplex *c__, integer *ldc, double *rwork);
int zlarf_check(char *side, integer *m, integer *n, dcomplex *v, integer *incv, dcomplex *tau,
                dcomplex *c__, integer *ldc, dcomplex *work);
int zlarfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt, dcomplex *c__,
                 integer *ldc, dcomplex *work, integer *ldwork);
int zlarfg_check(integer *n, dcomplex *alpha, dcomplex *x, integer *incx, dcomplex *tau);
int zlarfgp_check(integer *n, dcomplex *alpha, dcomplex *x, integer *incx, dcomplex *tau);
int zlarft_check(char *direct, char *storev, integer *n, integer *k, dcomplex *v, integer *ldv,
                 dcomplex *tau, dcomplex *t, integer *ldt);
int zlarfx_check(char *side, integer *m, integer *n, dcomplex *v, dcomplex *tau, dcomplex *c__,
                 integer *ldc, dcomplex *work);
int zlargv_check(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *c__,
                 integer *incc);
int zlarnv_check(integer *idist, integer *iseed, integer *n, dcomplex *x);
int zlarrv_check(integer *n, double *vl, double *vu, double *d__, double *l, double *pivmin,
                 integer *isplit, integer *m, integer *dol, integer *dou, double *minrgp,
                 double *rtol1, double *rtol2, double *w, double *werr, double *wgap,
                 integer *iblock, integer *indexw, double *gers, dcomplex *z__, integer *ldz,
                 integer *isuppz, double *work, integer *iwork, integer *info);
int zlarscl2_check(integer *m, integer *n, double *d__, dcomplex *x, integer *ldx);
int zlartg_check(dcomplex *f, dcomplex *g, double *cs, dcomplex *sn, dcomplex *r__);
int zlartv_check(integer *n, dcomplex *x, integer *incx, dcomplex *y, integer *incy, double *c__,
                 dcomplex *s, integer *incc);
int zlarz_check(char *side, integer *m, integer *n, integer *l, dcomplex *v, integer *incv,
                dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work);
int zlarzb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt,
                 dcomplex *c__, integer *ldc, dcomplex *work, integer *ldwork);
int zlarzt_check(char *direct, char *storev, integer *n, integer *k, dcomplex *v, integer *ldv,
                 dcomplex *tau, dcomplex *t, integer *ldt);
int zlascl_check(char *type__, integer *kl, integer *ku, double *cfrom, double *cto, integer *m,
                 integer *n, dcomplex *a, integer *lda, integer *info);
int zlascl2_check(integer *m, integer *n, double *d__, dcomplex *x, integer *ldx);
int zlaset_check(char *uplo, integer *m, integer *n, dcomplex *alpha, dcomplex *beta, dcomplex *a,
                 integer *lda);
int zlasr_check(char *side, char *pivot, char *direct, integer *m, integer *n, double *c__,
                double *s, dcomplex *a, integer *lda);
int zlassq_check(integer *n, dcomplex *x, integer *incx, double *scale, double *sumsq);
int zlaswp_check(integer *n, dcomplex *a, integer *lda, integer *k1, integer *k2, integer *ipiv,
                 integer *incx);
int zlasyf_check(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda,
                 integer *ipiv, dcomplex *w, integer *ldw, integer *info);
int zlasyf_rook_check(char *uplo, integer *n, integer *nb, integer *kb, dcomplex *a, integer *lda,
                      integer *ipiv, dcomplex *w, integer *ldw, integer *info);
int zlat2c_check(char *uplo, integer *n, dcomplex *a, integer *lda, scomplex *sa, integer *ldsa,
                 integer *info);
int zlatbs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd,
                 dcomplex *ab, integer *ldab, dcomplex *x, double *scale, double *cnorm,
                 integer *info);
int zlatdf_check(integer *ijob, integer *n, dcomplex *z__, integer *ldz, dcomplex *rhs,
                 double *rdsum, double *rdscal, integer *ipiv, integer *jpiv);
int zlatps_check(char *uplo, char *trans, char *diag, char *normin, integer *n, dcomplex *ap,
                 dcomplex *x, double *scale, double *cnorm, integer *info);
int zlatrd_check(char *uplo, integer *n, integer *nb, dcomplex *a, integer *lda, double *e,
                 dcomplex *tau, dcomplex *w, integer *ldw);
int zlatrs_check(char *uplo, char *trans, char *diag, char *normin, integer *n, dcomplex *a,
                 integer *lda, dcomplex *x, double *scale, double *cnorm, integer *info);
int zlatrz_check(integer *m, integer *n, integer *l, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work);
int zlatzm_check(char *side, integer *m, integer *n, dcomplex *v, integer *incv, dcomplex *tau,
                 dcomplex *c1, dcomplex *c2, integer *ldc, dcomplex *work);
int zlauu2_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *info);
int zlauum_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *info);
int zpbcon_check(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *anorm,
                 double *rcond, dcomplex *work, double *rwork, integer *info);
int zpbequ_check(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, double *s,
                 double *scond, double *amax, integer *info);
int zpbrfs_check(char *uplo, integer *n, integer *kd, integer *nrhs, dcomplex *ab, integer *ldab,
                 dcomplex *afb, integer *ldafb, dcomplex *b, integer *ldb, dcomplex *x,
                 integer *ldx, double *ferr, double *berr, dcomplex *work, double *rwork,
                 integer *info);
int zpbstf_check(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, integer *info);
int zpbsv_check(char *uplo, integer *n, integer *kd, integer *nrhs, dcomplex *ab, integer *ldab,
                dcomplex *b, integer *ldb, integer *info);
int zpbsvx_check(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, dcomplex *ab,
                 integer *ldab, dcomplex *afb, integer *ldafb, char *equed, double *s, dcomplex *b,
                 integer *ldb, dcomplex *x, integer *ldx, double *rcond, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int zpbtf2_check(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, integer *info);
int zpbtrf_check(char *uplo, integer *n, integer *kd, dcomplex *ab, integer *ldab, integer *info);
int zpbtrs_check(char *uplo, integer *n, integer *kd, integer *nrhs, dcomplex *ab, integer *ldab,
                 dcomplex *b, integer *ldb, integer *info);
int zpftrf_check(char *transr, char *uplo, integer *n, dcomplex *a, integer *info);
int zpftri_check(char *transr, char *uplo, integer *n, dcomplex *a, integer *info);
int zpftrs_check(char *transr, char *uplo, integer *n, integer *nrhs, dcomplex *a, dcomplex *b,
                 integer *ldb, integer *info);
int zpocon_check(char *uplo, integer *n, dcomplex *a, integer *lda, double *anorm, double *rcond,
                 dcomplex *work, double *rwork, integer *info);
int zpoequ_check(integer *n, dcomplex *a, integer *lda, double *s, double *scond, double *amax,
                 integer *info);
int zpoequb_check(integer *n, dcomplex *a, integer *lda, double *s, double *scond, double *amax,
                  integer *info);
int zporfs_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af,
                 integer *ldaf, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *ferr,
                 double *berr, dcomplex *work, double *rwork, integer *info);
int zporfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, double *s, dcomplex *b, integer *ldb, dcomplex *x,
                  integer *ldx, double *rcond, double *berr, integer *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams,
                  double *params, dcomplex *work, double *rwork, integer *info);
int zposv_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *b,
                integer *ldb, integer *info);
int zposvx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                 dcomplex *af, integer *ldaf, char *equed, double *s, dcomplex *b, integer *ldb,
                 dcomplex *x, integer *ldx, double *rcond, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int zposvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, char *equed, double *s, dcomplex *b, integer *ldb,
                  dcomplex *x, integer *ldx, double *rcond, double *rpvgrw, double *berr,
                  integer *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  integer *nparams, double *params, dcomplex *work, double *rwork, integer *info);
int zpotf2_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *info);
int zpotrf_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *info);
int zpotri_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *info);
int zpotrs_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, integer *info);
int zppcon_check(char *uplo, integer *n, dcomplex *ap, double *anorm, double *rcond, dcomplex *work,
                 double *rwork, integer *info);
int zppequ_check(char *uplo, integer *n, dcomplex *ap, double *s, double *scond, double *amax,
                 integer *info);
int zpprfs_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *afp, dcomplex *b,
                 integer *ldb, dcomplex *x, integer *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int zppsv_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *b, integer *ldb,
                integer *info);
int zppsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *afp,
                 char *equed, double *s, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx,
                 double *rcond, double *ferr, double *berr, dcomplex *work, double *rwork,
                 integer *info);
int zpptrf_check(char *uplo, integer *n, dcomplex *ap, integer *info);
int zpptri_check(char *uplo, integer *n, dcomplex *ap, integer *info);
int zpptrs_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *b, integer *ldb,
                 integer *info);
int zpstf2_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *piv, integer *rank,
                 double *tol, double *work, integer *info);
int zpstrf_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *piv, integer *rank,
                 double *tol, double *work, integer *info);
int zptcon_check(integer *n, double *d__, dcomplex *e, double *anorm, double *rcond, double *rwork,
                 integer *info);
int zpteqr_check(char *compz, integer *n, double *d__, double *e, dcomplex *z__, integer *ldz,
                 double *work, integer *info);
int zptrfs_check(char *uplo, integer *n, integer *nrhs, double *d__, dcomplex *e, double *df,
                 dcomplex *ef, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *ferr,
                 double *berr, dcomplex *work, double *rwork, integer *info);
int zptsv_check(integer *n, integer *nrhs, double *d__, dcomplex *e, dcomplex *b, integer *ldb,
                integer *info);
int zptsvx_check(char *fact, integer *n, integer *nrhs, double *d__, dcomplex *e, double *df,
                 dcomplex *ef, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int zpttrf_check(integer *n, double *d__, dcomplex *e, integer *info);
int zpttrs_check(char *uplo, integer *n, integer *nrhs, double *d__, dcomplex *e, dcomplex *b,
                 integer *ldb, integer *info);
int zptts2_check(integer *iuplo, integer *n, integer *nrhs, double *d__, dcomplex *e, dcomplex *b,
                 integer *ldb);
int zrot_check(integer *n, dcomplex *cx, integer *incx, dcomplex *cy, integer *incy, double *c__,
               dcomplex *s);
int zspcon_check(char *uplo, integer *n, dcomplex *ap, integer *ipiv, double *anorm, double *rcond,
                 dcomplex *work, integer *info);
int zspmv_check(char *uplo, integer *n, dcomplex *alpha, dcomplex *ap, dcomplex *x, integer *incx,
                dcomplex *beta, dcomplex *y, integer *incy);
int zspr_check(char *uplo, integer *n, dcomplex *alpha, dcomplex *x, integer *incx, dcomplex *ap);
int zsprfs_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *afp, integer *ipiv,
                 dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int zspsv_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, integer *ipiv, dcomplex *b,
                integer *ldb, integer *info);
int zspsvx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *ap, dcomplex *afp,
                 integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int zsptrf_check(char *uplo, integer *n, dcomplex *ap, integer *ipiv, integer *info);
int zsptri_check(char *uplo, integer *n, dcomplex *ap, integer *ipiv, dcomplex *work,
                 integer *info);
int zsptrs_check(char *uplo, integer *n, integer *nrhs, dcomplex *ap, integer *ipiv, dcomplex *b,
                 integer *ldb, integer *info);
int zstedc_check(char *compz, integer *n, double *d__, double *e, dcomplex *z__, integer *ldz,
                 dcomplex *work, integer *lwork, double *rwork, integer *lrwork, integer *iwork,
                 integer *liwork, integer *info);
int zstegr_check(char *jobz, char *range, integer *n, double *d__, double *e, double *vl,
                 double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w,
                 dcomplex *z__, integer *ldz, integer *isuppz, double *work, integer *lwork,
                 integer *iwork, integer *liwork, integer *info);
int zstein_check(integer *n, double *d__, double *e, integer *m, double *w, integer *iblock,
                 integer *isplit, dcomplex *z__, integer *ldz, double *work, integer *iwork,
                 integer *ifail, integer *info);
int zstemr_check(char *jobz, char *range, integer *n, double *d__, double *e, double *vl,
                 double *vu, integer *il, integer *iu, integer *m, double *w, dcomplex *z__,
                 integer *ldz, integer *nzc, integer *isuppz, logical *tryrac, double *work,
                 integer *lwork, integer *iwork, integer *liwork, integer *info);
int zsteqr_check(char *compz, integer *n, double *d__, double *e, dcomplex *z__, integer *ldz,
                 double *work, integer *info);
int zsycon_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, double *anorm,
                 double *rcond, dcomplex *work, integer *info);
int zsycon_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      double *anorm, double *rcond, dcomplex *work, integer *info);
int zsyconv_check(char *uplo, char *way, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                  dcomplex *work, integer *info);
int zsyequb_check(char *uplo, integer *n, dcomplex *a, integer *lda, double *s, double *scond,
                  double *amax, dcomplex *work, integer *info);
int zsymv_check(char *uplo, integer *n, dcomplex *alpha, dcomplex *a, integer *lda, dcomplex *x,
                integer *incx, dcomplex *beta, dcomplex *y, integer *incy);
int zsyr_check(char *uplo, integer *n, dcomplex *alpha, dcomplex *x, integer *incx, dcomplex *a,
               integer *lda);
int zsyrfs_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, dcomplex *af,
                 integer *ldaf, integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int zsyrfsx_check(char *uplo, char *equed, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, integer *ipiv, double *s, dcomplex *b, integer *ldb,
                  dcomplex *x, integer *ldx, double *rcond, double *berr, integer *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, integer *nparams,
                  double *params, dcomplex *work, double *rwork, integer *info);
int zsysv_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv,
                dcomplex *b, integer *ldb, dcomplex *work, integer *lwork, integer *info);
int zsysv_rook_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                     integer *ipiv, dcomplex *b, integer *ldb, dcomplex *work, integer *lwork,
                     integer *info);
int zsysvx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                 dcomplex *af, integer *ldaf, integer *ipiv, dcomplex *b, integer *ldb, dcomplex *x,
                 integer *ldx, double *rcond, double *ferr, double *berr, dcomplex *work,
                 integer *lwork, double *rwork, integer *info);
int zsysvxx_check(char *fact, char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                  dcomplex *af, integer *ldaf, integer *ipiv, char *equed, double *s, dcomplex *b,
                  integer *ldb, dcomplex *x, integer *ldx, double *rcond, double *rpvgrw,
                  double *berr, integer *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, integer *nparams, double *params, dcomplex *work,
                  double *rwork, integer *info);
int zsyswapr_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *i1, integer *i2);
int zsytf2_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, integer *info);
int zsytf2_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      integer *info);
int zsytrf_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                 integer *lwork, integer *info);
int zsytrf_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      dcomplex *work, integer *lwork, integer *info);
int zsytri_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                 integer *info);
int zsytri2_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                  integer *lwork, integer *info);
int zsytri2x_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv, dcomplex *work,
                   integer *nb, integer *info);
int zsytri_rook_check(char *uplo, integer *n, dcomplex *a, integer *lda, integer *ipiv,
                      dcomplex *work, integer *info);
int zsytrs_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv,
                 dcomplex *b, integer *ldb, integer *info);
int zsytrs2_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda, integer *ipiv,
                  dcomplex *b, integer *ldb, dcomplex *work, integer *info);
int zsytrs_rook_check(char *uplo, integer *n, integer *nrhs, dcomplex *a, integer *lda,
                      integer *ipiv, dcomplex *b, integer *ldb, integer *info);
int ztbcon_check(char *norm, char *uplo, char *diag, integer *n, integer *kd, dcomplex *ab,
                 integer *ldab, double *rcond, dcomplex *work, double *rwork, integer *info);
int ztbrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 dcomplex *ab, integer *ldab, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, integer *info);
int ztbtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs,
                 dcomplex *ab, integer *ldab, dcomplex *b, integer *ldb, integer *info);
int ztfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m,
                integer *n, dcomplex *alpha, dcomplex *a, dcomplex *b, integer *ldb);
int ztftri_check(char *transr, char *uplo, char *diag, integer *n, dcomplex *a, integer *info);
int ztfttp_check(char *transr, char *uplo, integer *n, dcomplex *arf, dcomplex *ap, integer *info);
int ztfttr_check(char *transr, char *uplo, integer *n, dcomplex *arf, dcomplex *a, integer *lda,
                 integer *info);
int ztgevc_check(char *side, char *howmny, logical *select, integer *n, dcomplex *s, integer *lds,
                 dcomplex *p, integer *ldp, dcomplex *vl, integer *ldvl, dcomplex *vr,
                 integer *ldvr, integer *mm, integer *m, dcomplex *work, double *rwork,
                 integer *info);
int ztgex2_check(logical *wantq, logical *wantz, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, dcomplex *q, integer *ldq, dcomplex *z__, integer *ldz, integer *j1,
                 integer *info);
int ztgexc_check(logical *wantq, logical *wantz, integer *n, dcomplex *a, integer *lda, dcomplex *b,
                 integer *ldb, dcomplex *q, integer *ldq, dcomplex *z__, integer *ldz,
                 integer *ifst, integer *ilst, integer *info);
int ztgsen_check(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n,
                 dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *alpha,
                 dcomplex *beta, dcomplex *q, integer *ldq, dcomplex *z__, integer *ldz, integer *m,
                 double *pl, double *pr, double *dif, dcomplex *work, integer *lwork,
                 integer *iwork, integer *liwork, integer *info);
int ztgsja_check(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k,
                 integer *l, dcomplex *a, integer *lda, dcomplex *b, integer *ldb, double *tola,
                 double *tolb, double *alpha, double *beta, dcomplex *u, integer *ldu, dcomplex *v,
                 integer *ldv, dcomplex *q, integer *ldq, dcomplex *work, integer *ncycle,
                 integer *info);
int ztgsna_check(char *job, char *howmny, logical *select, integer *n, dcomplex *a, integer *lda,
                 dcomplex *b, integer *ldb, dcomplex *vl, integer *ldvl, dcomplex *vr,
                 integer *ldvr, double *s, double *dif, integer *mm, integer *m, dcomplex *work,
                 integer *lwork, integer *iwork, integer *info);
int ztgsy2_check(char *trans, integer *ijob, integer *m, integer *n, dcomplex *a, integer *lda,
                 dcomplex *b, integer *ldb, dcomplex *c__, integer *ldc, dcomplex *d__,
                 integer *ldd, dcomplex *e, integer *lde, dcomplex *f, integer *ldf, double *scale,
                 double *rdsum, double *rdscal, integer *info);
int ztgsyl_check(char *trans, integer *ijob, integer *m, integer *n, dcomplex *a, integer *lda,
                 dcomplex *b, integer *ldb, dcomplex *c__, integer *ldc, dcomplex *d__,
                 integer *ldd, dcomplex *e, integer *lde, dcomplex *f, integer *ldf, double *scale,
                 double *dif, dcomplex *work, integer *lwork, integer *iwork, integer *info);
int ztpcon_check(char *norm, char *uplo, char *diag, integer *n, dcomplex *ap, double *rcond,
                 dcomplex *work, double *rwork, integer *info);
int ztpmqrt_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                  integer *nb, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt, dcomplex *a,
                  integer *lda, dcomplex *b, integer *ldb, dcomplex *work, integer *info);
int ztpqrt_check(integer *m, integer *n, integer *l, integer *nb, dcomplex *a, integer *lda,
                 dcomplex *b, integer *ldb, dcomplex *t, integer *ldt, dcomplex *work,
                 integer *info);
int ztpqrt2_check(integer *m, integer *n, integer *l, dcomplex *a, integer *lda, dcomplex *b,
                  integer *ldb, dcomplex *t, integer *ldt, integer *info);
int ztprfb_check(char *side, char *trans, char *direct, char *storev, integer *m, integer *n,
                 integer *k, integer *l, dcomplex *v, integer *ldv, dcomplex *t, integer *ldt,
                 dcomplex *a, integer *lda, dcomplex *b, integer *ldb, dcomplex *work,
                 integer *ldwork);
int ztprfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, dcomplex *ap,
                 dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, integer *info);
int ztptri_check(char *uplo, char *diag, integer *n, dcomplex *ap, integer *info);
int ztptrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, dcomplex *ap,
                 dcomplex *b, integer *ldb, integer *info);
int ztpttf_check(char *transr, char *uplo, integer *n, dcomplex *ap, dcomplex *arf, integer *info);
int ztpttr_check(char *uplo, integer *n, dcomplex *ap, dcomplex *a, integer *lda, integer *info);
int ztrcon_check(char *norm, char *uplo, char *diag, integer *n, dcomplex *a, integer *lda,
                 double *rcond, dcomplex *work, double *rwork, integer *info);
int ztrevc_check(char *side, char *howmny, logical *select, integer *n, dcomplex *t, integer *ldt,
                 dcomplex *vl, integer *ldvl, dcomplex *vr, integer *ldvr, integer *mm, integer *m,
                 dcomplex *work, double *rwork, integer *info);
int ztrexc_check(char *compq, integer *n, dcomplex *t, integer *ldt, dcomplex *q, integer *ldq,
                 integer *ifst, integer *ilst, integer *info);
int ztrrfs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, dcomplex *a,
                 integer *lda, dcomplex *b, integer *ldb, dcomplex *x, integer *ldx, double *ferr,
                 double *berr, dcomplex *work, double *rwork, integer *info);
int ztrsen_check(char *job, char *compq, logical *select, integer *n, dcomplex *t, integer *ldt,
                 dcomplex *q, integer *ldq, dcomplex *w, integer *m, double *s, double *sep,
                 dcomplex *work, integer *lwork, integer *info);
int ztrsna_check(char *job, char *howmny, logical *select, integer *n, dcomplex *t, integer *ldt,
                 dcomplex *vl, integer *ldvl, dcomplex *vr, integer *ldvr, double *s, double *sep,
                 integer *mm, integer *m, dcomplex *work, integer *ldwork, double *rwork,
                 integer *info);
int ztrsyl_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, dcomplex *a,
                 integer *lda, dcomplex *b, integer *ldb, dcomplex *c__, integer *ldc,
                 double *scale, integer *info);
int ztrsyl3_check(char *trana, char *tranb, integer *isgn, integer *m, integer *n, doublecomplex *a,
                  integer *lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc,
                  doublereal *scale, doublereal *swork, integer *ldswork, integer *info);
int ztrti2_check(char *uplo, char *diag, integer *n, dcomplex *a, integer *lda, integer *info);
int ztrtri_check(char *uplo, char *diag, integer *n, dcomplex *a, integer *lda, integer *info);
int ztrtrs_check(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, dcomplex *a,
                 integer *lda, dcomplex *b, integer *ldb, integer *info);
int ztrttf_check(char *transr, char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *arf,
                 integer *info);
int ztrttp_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *ap, integer *info);
int ztzrqf_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, integer *info);
int ztzrzf_check(integer *m, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *lwork, integer *info);
int zunbdb_check(char *trans, char *signs, integer *m, integer *p, integer *q, dcomplex *x11,
                 integer *ldx11, dcomplex *x12, integer *ldx12, dcomplex *x21, integer *ldx21,
                 dcomplex *x22, integer *ldx22, double *theta, double *phi, dcomplex *taup1,
                 dcomplex *taup2, dcomplex *tauq1, dcomplex *tauq2, dcomplex *work, integer *lwork,
                 integer *info);
int zunbdb1_check(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21,
                  integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info);
int zunbdb2_check(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21,
                  integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info);
int zunbdb3_check(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21,
                  integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *work, integer *lwork, integer *info);
int zunbdb4_check(integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x21,
                  integer *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *phantom, dcomplex *work, integer *lwork,
                  integer *info);
int zunbdb5_check(integer *m1, integer *m2, integer *n, dcomplex *x1, integer *incx1, dcomplex *x2,
                  integer *incx2, dcomplex *q1, integer *ldq1, dcomplex *q2, integer *ldq2,
                  dcomplex *work, integer *lwork, integer *info);
int zunbdb6_check(integer *m1, integer *m2, integer *n, dcomplex *x1, integer *incx1, dcomplex *x2,
                  integer *incx2, dcomplex *q1, integer *ldq1, dcomplex *q2, integer *ldq2,
                  dcomplex *work, integer *lwork, integer *info);
int zuncsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 integer *m, integer *p, integer *q, dcomplex *x11, integer *ldx11, dcomplex *x12,
                 integer *ldx12, dcomplex *x21, integer *ldx21, dcomplex *x22, integer *ldx22,
                 double *theta, dcomplex *u1, integer *ldu1, dcomplex *u2, integer *ldu2,
                 dcomplex *v1t, integer *ldv1t, dcomplex *v2t, integer *ldv2t, dcomplex *work,
                 integer *lwork, double *rwork, integer *lrwork, integer *iwork, integer *info);
int zuncsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, integer *m, integer *p, integer *q,
                     dcomplex *x11, integer *ldx11, dcomplex *x21, integer *ldx21, double *theta,
                     dcomplex *u1, integer *ldu1, dcomplex *u2, integer *ldu2, dcomplex *v1t,
                     integer *ldv1t, dcomplex *work, integer *lwork, double *rwork, integer *lrwork,
                     integer *iwork, integer *info);
int zung2l_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *info);
int zung2r_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *info);
int zungbr_check(char *vect, integer *m, integer *n, integer *k, dcomplex *a, integer *lda,
                 dcomplex *tau, dcomplex *work, integer *lwork, integer *info);
int zunghr_check(integer *n, integer *ilo, integer *ihi, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *lwork, integer *info);
int zungl2_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *info);
int zunglq_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *lwork, integer *info);
int zungql_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *lwork, integer *info);
int zungqr_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *lwork, integer *info);
int zungr2_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *info);
int zungrq_check(integer *m, integer *n, integer *k, dcomplex *a, integer *lda, dcomplex *tau,
                 dcomplex *work, integer *lwork, integer *info);
int zungtr_check(char *uplo, integer *n, dcomplex *a, integer *lda, dcomplex *tau, dcomplex *work,
                 integer *lwork, integer *info);
int zunm2l_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *info);
int zunm2r_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *info);
int zunmbr_check(char *vect, char *side, char *trans, integer *m, integer *n, integer *k,
                 dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc,
                 dcomplex *work, integer *lwork, integer *info);
int zunmhr_check(char *side, char *trans, integer *m, integer *n, integer *ilo, integer *ihi,
                 dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc,
                 dcomplex *work, integer *lwork, integer *info);
int zunml2_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *info);
int zunmlq_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *lwork, integer *info);
int zunmql_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *lwork, integer *info);
int zunmqr_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *lwork, integer *info);
int zunmr2_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *info);
int zunmr3_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                 dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc,
                 dcomplex *work, integer *info);
int zunmrq_check(char *side, char *trans, integer *m, integer *n, integer *k, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *lwork, integer *info);
int zunmrz_check(char *side, char *trans, integer *m, integer *n, integer *k, integer *l,
                 dcomplex *a, integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc,
                 dcomplex *work, integer *lwork, integer *info);
int zunmtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, dcomplex *a,
                 integer *lda, dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work,
                 integer *lwork, integer *info);
int zupgtr_check(char *uplo, integer *n, dcomplex *ap, dcomplex *tau, dcomplex *q, integer *ldq,
                 dcomplex *work, integer *info);
int zupmtr_check(char *side, char *uplo, char *trans, integer *m, integer *n, dcomplex *ap,
                 dcomplex *tau, dcomplex *c__, integer *ldc, dcomplex *work, integer *info);
