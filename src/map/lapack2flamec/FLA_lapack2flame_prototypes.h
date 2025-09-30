/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

    Copyright (c) 2025 Advanced Micro Devices, Inc. All rights reserved.
*/
#ifndef LAPACK_2_FLAME_H
#define LAPACK_2_FLAME_H

int cbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, aocl_int64_t *m,
                 aocl_int64_t *p, aocl_int64_t *q, float *theta, float *phi, scomplex *u1, aocl_int64_t *ldu1,
                 scomplex *u2, aocl_int64_t *ldu2, scomplex *v1t, aocl_int64_t *ldv1t, scomplex *v2t,
                 aocl_int64_t *ldv2t, float *b11d, float *b11e, float *b12d, float *b12e, float *b21d,
                 float *b21e, float *b22d, float *b22e, float *rwork, aocl_int64_t *lrwork,
                 aocl_int64_t *info);
int cbdsqr_check(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru, aocl_int64_t *ncc, float *d__,
                 float *e, scomplex *vt, aocl_int64_t *ldvt, scomplex *u, aocl_int64_t *ldu, scomplex *c__,
                 aocl_int64_t *ldc, float *rwork, aocl_int64_t *info);
int cgbbrd_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ncc, aocl_int64_t *kl, aocl_int64_t *ku,
                 scomplex *ab, aocl_int64_t *ldab, float *d__, float *e, scomplex *q, aocl_int64_t *ldq,
                 scomplex *pt, aocl_int64_t *ldpt, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 float *rwork, aocl_int64_t *info);
int cgbcon_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, float *anorm, float *rcond, scomplex *work, float *rwork,
                 aocl_int64_t *info);
int cgbequ_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab, aocl_int64_t *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int cgbequb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab, aocl_int64_t *ldab,
                  float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int cgbrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, scomplex *ab,
                 aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, scomplex *b,
                 aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr, float *berr, scomplex *work,
                 float *rwork, aocl_int64_t *info);
int cgbrfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  scomplex *ab, aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                  float *r__, float *c__, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx,
                  float *rcond, float *berr, aocl_int64_t *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, aocl_int64_t *nparams, float *params, scomplex *work,
                  float *rwork, aocl_int64_t *info);
int cgbsv_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, scomplex *ab, aocl_int64_t *ldab,
                aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int cgbsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                 scomplex *ab, aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                 char *equed, float *r__, float *c__, scomplex *b, aocl_int64_t *ldb, scomplex *x,
                 aocl_int64_t *ldx, float *rcond, float *ferr, float *berr, scomplex *work, float *rwork,
                 aocl_int64_t *info);
int cgbsvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  scomplex *ab, aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                  char *equed, float *r__, float *c__, scomplex *b, aocl_int64_t *ldb, scomplex *x,
                  aocl_int64_t *ldx, float *rcond, float *rpvgrw, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  scomplex *work, float *rwork, aocl_int64_t *info);
int cgbtf2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int cgbtrf_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int cgbtrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, scomplex *ab,
                 aocl_int64_t *ldab, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int cgebak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *scale,
                 aocl_int64_t *m, scomplex *v, aocl_int64_t *ldv, aocl_int64_t *info);
int cgebal_check(char *job, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 float *scale, aocl_int64_t *info);
int cgebd2_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *d__, float *e,
                 scomplex *tauq, scomplex *taup, scomplex *work, aocl_int64_t *info);
int cgebrd_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *d__, float *e,
                 scomplex *tauq, scomplex *taup, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cgecon_check(char *norm, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *anorm, float *rcond,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int cgeequ_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int cgeequb_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *r__, float *c__,
                  float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int cgees_check(char *jobvs, char *sort, L_fp select, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                aocl_int64_t *sdim, scomplex *w, scomplex *vs, aocl_int64_t *ldvs, scomplex *work,
                aocl_int64_t *lwork, float *rwork, logical *bwork, aocl_int64_t *info);
int cgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, aocl_int64_t *n, scomplex *a,
                 aocl_int64_t *lda, aocl_int64_t *sdim, scomplex *w, scomplex *vs, aocl_int64_t *ldvs,
                 float *rconde, float *rcondv, scomplex *work, aocl_int64_t *lwork, float *rwork,
                 logical *bwork, aocl_int64_t *info);
int cgeev_check(char *jobvl, char *jobvr, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *w,
                scomplex *vl, aocl_int64_t *ldvl, scomplex *vr, aocl_int64_t *ldvr, scomplex *work,
                aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int cgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, scomplex *a,
                 aocl_int64_t *lda, scomplex *w, scomplex *vl, aocl_int64_t *ldvl, scomplex *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *ilo, aocl_int64_t *ihi, float *scale, float *abnrm,
                 float *rconde, float *rcondv, scomplex *work, aocl_int64_t *lwork, float *rwork,
                 aocl_int64_t *info);
int cgegs_check(char *jobvsl, char *jobvsr, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                aocl_int64_t *ldb, scomplex *alpha, scomplex *beta, scomplex *vsl, aocl_int64_t *ldvsl,
                scomplex *vsr, aocl_int64_t *ldvsr, scomplex *work, aocl_int64_t *lwork, float *rwork,
                aocl_int64_t *info);
int cgegv_check(char *jobvl, char *jobvr, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                aocl_int64_t *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, aocl_int64_t *ldvl,
                scomplex *vr, aocl_int64_t *ldvr, scomplex *work, aocl_int64_t *lwork, float *rwork,
                aocl_int64_t *info);
int cgehd2_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *info);
int cgehrd_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cgelq2_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *info);
int cgelqf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cgels_check(char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cgelsd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, float *s, float *rcond, aocl_int64_t *rank, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int cgelss_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, float *s, float *rcond, aocl_int64_t *rank, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int cgelsx_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, float *rcond, aocl_int64_t *rank, scomplex *work,
                 float *rwork, aocl_int64_t *info);
int cgelsy_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, float *rcond, aocl_int64_t *rank, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int cgemqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb,
                  scomplex *v, aocl_int64_t *ldv, scomplex *t, aocl_int64_t *ldt, scomplex *c__, aocl_int64_t *ldc,
                  scomplex *work, aocl_int64_t *info);
int cgeql2_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *info);
int cgeqlf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cgeqp3_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int cgeqpf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt, scomplex *tau,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int cgeqr2_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *info);
int cgeqr2p_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                  aocl_int64_t *info);
int cgeqrf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cgeqrfp_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int cgeqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, scomplex *a, aocl_int64_t *lda, scomplex *t,
                 aocl_int64_t *ldt, scomplex *work, aocl_int64_t *info);
int cgeqrt2_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int cgeqrt3_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int cgerfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int cgerfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *r__, float *c__, scomplex *b,
                  aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond, float *berr,
                  aocl_int64_t *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  aocl_int64_t *nparams, float *params, scomplex *work, float *rwork, aocl_int64_t *info);
int cgerq2_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *info);
int cgerqf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cgesc2_check(aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *rhs, aocl_int_t *ipiv, aocl_int64_t *jpiv,
                 float *scale);
int cgesdd_check(char *jobz, aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *s,
                 scomplex *u, aocl_int64_t *ldu, scomplex *vt, aocl_int64_t *ldvt, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int cgesv_check(aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int cgesvd_check(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 float *s, scomplex *u, aocl_int64_t *ldu, scomplex *vt, aocl_int64_t *ldvt, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int cgesvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                 scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, float *r__, float *c__,
                 scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond, float *ferr,
                 float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int cgesvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, float *r__, float *c__,
                  scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond, float *rpvgrw,
                  float *berr, aocl_int64_t *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, aocl_int64_t *nparams, float *params, scomplex *work,
                  float *rwork, aocl_int64_t *info);
int cgetc2_check(aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *jpiv,
                 aocl_int64_t *info);
int cgetf2_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int cgetrf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int cgetrfnp_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int cgetrfnpi_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nfact, scomplex *a, aocl_int64_t *lda,
                    aocl_int64_t *info);
int cgetri_check(aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cgetrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int cggbak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *lscale,
                 float *rscale, aocl_int64_t *m, scomplex *v, aocl_int64_t *ldv, aocl_int64_t *info);
int cggbal_check(char *job, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb,
                 aocl_int64_t *ilo, aocl_int64_t *ihi, float *lscale, float *rscale, float *work,
                 aocl_int64_t *info);
int cgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, aocl_int64_t *n, scomplex *a,
                aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *sdim, scomplex *alpha,
                scomplex *beta, scomplex *vsl, aocl_int64_t *ldvsl, scomplex *vsr, aocl_int64_t *ldvsr,
                scomplex *work, aocl_int64_t *lwork, float *rwork, logical *bwork, aocl_int64_t *info);
int cggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, aocl_int64_t *n,
                 scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *sdim,
                 scomplex *alpha, scomplex *beta, scomplex *vsl, aocl_int64_t *ldvsl, scomplex *vsr,
                 aocl_int64_t *ldvsr, float *rconde, float *rcondv, scomplex *work, aocl_int64_t *lwork,
                 float *rwork, aocl_int64_t *iwork, aocl_int64_t *liwork, logical *bwork, aocl_int64_t *info);
int cggev_check(char *jobvl, char *jobvr, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                aocl_int64_t *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, aocl_int64_t *ldvl,
                scomplex *vr, aocl_int64_t *ldvr, scomplex *work, aocl_int64_t *lwork, float *rwork,
                aocl_int64_t *info);
int cggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, scomplex *a,
                 aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *alpha, scomplex *beta,
                 scomplex *vl, aocl_int64_t *ldvl, scomplex *vr, aocl_int64_t *ldvr, aocl_int64_t *ilo,
                 aocl_int64_t *ihi, float *lscale, float *rscale, float *abnrm, float *bbnrm,
                 float *rconde, float *rcondv, scomplex *work, aocl_int64_t *lwork, float *rwork,
                 aocl_int64_t *iwork, logical *bwork, aocl_int64_t *info);
int cggglm_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, scomplex *d__, scomplex *x, scomplex *y, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cgghrd_check(char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, scomplex *a,
                 aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *q, aocl_int64_t *ldq, scomplex *z__,
                 aocl_int64_t *ldz, aocl_int64_t *info);
int cgglse_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, scomplex *c__, scomplex *d__, scomplex *x, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cggqrf_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, scomplex *a, aocl_int64_t *lda, scomplex *taua,
                 scomplex *b, aocl_int64_t *ldb, scomplex *taub, scomplex *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int cggrqf_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *taua,
                 scomplex *b, aocl_int64_t *ldb, scomplex *taub, scomplex *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int cggsvd_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, aocl_int64_t *k,
                 aocl_int64_t *l, scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, float *alpha,
                 float *beta, scomplex *u, aocl_int64_t *ldu, scomplex *v, aocl_int64_t *ldv, scomplex *q,
                 aocl_int64_t *ldq, scomplex *work, float *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int cggsvp_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n,
                 scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, float *tola, float *tolb,
                 aocl_int64_t *k, aocl_int64_t *l, scomplex *u, aocl_int64_t *ldu, scomplex *v, aocl_int64_t *ldv,
                 scomplex *q, aocl_int64_t *ldq, aocl_int64_t *iwork, float *rwork, scomplex *tau,
                 scomplex *work, aocl_int64_t *info);
int cgtcon_check(char *norm, aocl_int64_t *n, scomplex *dl, scomplex *d__, scomplex *du, scomplex *du2,
                 aocl_int_t *ipiv, float *anorm, float *rcond, scomplex *work, aocl_int64_t *info);
int cgtrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *dl, scomplex *d__, scomplex *du,
                 scomplex *dlf, scomplex *df, scomplex *duf, scomplex *du2, aocl_int_t *ipiv,
                 scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int cgtsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *dl, scomplex *d__, scomplex *du, scomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int cgtsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *dl, scomplex *d__,
                 scomplex *du, scomplex *dlf, scomplex *df, scomplex *duf, scomplex *du2,
                 aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int cgttrf_check(aocl_int64_t *n, scomplex *dl, scomplex *d__, scomplex *du, scomplex *du2,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int cgttrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *dl, scomplex *d__, scomplex *du,
                 scomplex *du2, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int cgtts2_check(aocl_int64_t *itrans, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *dl, scomplex *d__,
                 scomplex *du, scomplex *du2, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb);
int chbev_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab,
                float *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work, float *rwork, aocl_int64_t *info);
int chbevd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab,
                 float *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work, aocl_int64_t *lwork,
                 float *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int chbevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab,
                 aocl_int64_t *ldab, scomplex *q, aocl_int64_t *ldq, float *vl, float *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, scomplex *z__, aocl_int64_t *ldz,
                 scomplex *work, float *rwork, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int chbgst_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, scomplex *ab,
                 aocl_int64_t *ldab, scomplex *bb, aocl_int64_t *ldbb, scomplex *x, aocl_int64_t *ldx,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int chbgv_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, scomplex *ab,
                aocl_int64_t *ldab, scomplex *bb, aocl_int64_t *ldbb, float *w, scomplex *z__, aocl_int64_t *ldz,
                scomplex *work, float *rwork, aocl_int64_t *info);
int chbgvd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, scomplex *ab,
                 aocl_int64_t *ldab, scomplex *bb, aocl_int64_t *ldbb, float *w, scomplex *z__, aocl_int64_t *ldz,
                 scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int chbgvx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb,
                 scomplex *ab, aocl_int64_t *ldab, scomplex *bb, aocl_int64_t *ldbb, scomplex *q,
                 aocl_int64_t *ldq, float *vl, float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol,
                 aocl_int64_t *m, float *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work, float *rwork,
                 aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int chbtrd_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab,
                 float *d__, float *e, scomplex *q, aocl_int64_t *ldq, scomplex *work, aocl_int64_t *info);
int checon_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *anorm,
                 float *rcond, scomplex *work, aocl_int64_t *info);
int checon_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      float *anorm, float *rcond, scomplex *work, aocl_int64_t *info);
int cheequb_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *s, float *scond,
                  float *amax, scomplex *work, aocl_int64_t *info);
int cheev_check(char *jobz, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *w,
                scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int cheevd_check(char *jobz, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *w,
                 scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int cheevr_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 float *vl, float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m,
                 float *w, scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork,
                 aocl_int64_t *info);
int cheevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 float *vl, float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m,
                 float *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work, aocl_int64_t *lwork,
                 float *rwork, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int chegs2_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int chegst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int chegv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                scomplex *b, aocl_int64_t *ldb, float *w, scomplex *work, aocl_int64_t *lwork, float *rwork,
                aocl_int64_t *info);
int chegvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 scomplex *b, aocl_int64_t *ldb, float *w, scomplex *work, aocl_int64_t *lwork, float *rwork,
                 aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int chegvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, scomplex *a,
                 aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, float *vl, float *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, scomplex *z__, aocl_int64_t *ldz,
                 scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *iwork, aocl_int64_t *ifail,
                 aocl_int64_t *info);
int cherfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int cherfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *s, scomplex *b, aocl_int64_t *ldb,
                  scomplex *x, aocl_int64_t *ldx, float *rcond, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  scomplex *work, float *rwork, aocl_int64_t *info);
int chesv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int chesv_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                     aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *lwork,
                     aocl_int64_t *info);
int chesvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                 scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x,
                 aocl_int64_t *ldx, float *rcond, float *ferr, float *berr, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int chesvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, float *s, scomplex *b,
                  aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond, float *rpvgrw, float *berr,
                  aocl_int64_t *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  aocl_int64_t *nparams, float *params, scomplex *work, float *rwork, aocl_int64_t *info);
int cheswapr_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *i1, aocl_int64_t *i2);
int chetd2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *d__, float *e,
                 scomplex *tau, aocl_int64_t *info);
int chetf2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int chetf2_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      aocl_int64_t *info);
int chetrd_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *d__, float *e,
                 scomplex *tau, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int chetrf_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int chetrf_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int chetri_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                 aocl_int64_t *info);
int chetri2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int chetri2x_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                   aocl_int64_t *nb, aocl_int64_t *info);
int chetri_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      scomplex *work, aocl_int64_t *info);
int chetrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int chetrs2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *info);
int chetrs_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int chfrk_check(char *transr, char *uplo, char *trans, aocl_int64_t *n, aocl_int64_t *k, float *alpha,
                scomplex *a, aocl_int64_t *lda, float *beta, scomplex *c__);
int chgeqz_check(char *job, char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 scomplex *h__, aocl_int64_t *ldh, scomplex *t, aocl_int64_t *ldt, scomplex *alpha,
                 scomplex *beta, scomplex *q, aocl_int64_t *ldq, scomplex *z__, aocl_int64_t *ldz,
                 scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
void chla_transtype_check(char *ret_val, aocl_int64_t *trans);
int chpcon_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int_t *ipiv, float *anorm, float *rcond,
                 scomplex *work, aocl_int64_t *info);
int chpev_check(char *jobz, char *uplo, aocl_int64_t *n, scomplex *ap, float *w, scomplex *z__,
                aocl_int64_t *ldz, scomplex *work, float *rwork, aocl_int64_t *info);
int chpevd_check(char *jobz, char *uplo, aocl_int64_t *n, scomplex *ap, float *w, scomplex *z__,
                 aocl_int64_t *ldz, scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *lrwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int chpevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, scomplex *ap, float *vl,
                 float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w,
                 scomplex *z__, aocl_int64_t *ldz, scomplex *work, float *rwork, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int chpgst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, scomplex *ap, scomplex *bp, aocl_int64_t *info);
int chpgv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, scomplex *ap, scomplex *bp,
                float *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work, float *rwork, aocl_int64_t *info);
int chpgvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, scomplex *ap, scomplex *bp,
                 float *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work, aocl_int64_t *lwork,
                 float *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int chpgvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, scomplex *ap,
                 scomplex *bp, float *vl, float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol,
                 aocl_int64_t *m, float *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work, float *rwork,
                 aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int chprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *afp, aocl_int_t *ipiv,
                 scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int chpsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, aocl_int_t *ipiv, scomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int chpsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *afp,
                 aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int chptrd_check(char *uplo, aocl_int64_t *n, scomplex *ap, float *d__, float *e, scomplex *tau,
                 aocl_int64_t *info);
int chptrf_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int_t *ipiv, aocl_int64_t *info);
int chptri_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int_t *ipiv, scomplex *work,
                 aocl_int64_t *info);
int chptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, aocl_int_t *ipiv, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int chsein_check(char *side, char *eigsrc, char *initv, logical *select, aocl_int64_t *n, scomplex *h__,
                 aocl_int64_t *ldh, scomplex *w, scomplex *vl, aocl_int64_t *ldvl, scomplex *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m, scomplex *work, float *rwork,
                 aocl_int64_t *ifaill, aocl_int64_t *ifailr, aocl_int64_t *info);
int chseqr_check(char *job, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, scomplex *h__,
                 aocl_int64_t *ldh, scomplex *w, scomplex *z__, aocl_int64_t *ldz, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cla_gbamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *alpha,
                    scomplex *ab, aocl_int64_t *ldab, scomplex *x, aocl_int64_t *incx, float *beta, float *y,
                    aocl_int64_t *incy);
float cla_gbrcond_c_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab,
                          aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, float *c__,
                          logical *capply, aocl_int64_t *info, scomplex *work, float *rwork);
float cla_gbrcond_x_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab,
                          aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, scomplex *x,
                          aocl_int64_t *info, scomplex *work, float *rwork);
int cla_gbrfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n, aocl_int64_t *kl,
                              aocl_int64_t *ku, aocl_int64_t *nrhs, scomplex *ab, aocl_int64_t *ldab,
                              scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, logical *colequ,
                              float *c__, scomplex *b, aocl_int64_t *ldb, scomplex *y, aocl_int64_t *ldy,
                              float *berr_out__, aocl_int64_t *n_norms__, float *err_bnds_norm__,
                              float *err_bnds_comp__, scomplex *res, float *ayb, scomplex *dy,
                              scomplex *y_tail__, float *rcond, aocl_int64_t *ithresh, float *rthresh,
                              float *dz_ub__, logical *ignore_cwise__, aocl_int64_t *info);
float cla_gbrpvgrw_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *ncols, scomplex *ab,
                         aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb);
int cla_geamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, float *alpha, scomplex *a, aocl_int64_t *lda,
                    scomplex *x, aocl_int64_t *incx, float *beta, float *y, aocl_int64_t *incy);
float cla_gercond_c_check(char *trans, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, aocl_int_t *ipiv, float *c__, logical *capply, aocl_int64_t *info,
                          scomplex *work, float *rwork);
float cla_gercond_x_check(char *trans, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *x, aocl_int64_t *info, scomplex *work,
                          float *rwork);
int cla_gerfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n,
                              aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *af, aocl_int64_t *ldaf,
                              aocl_int_t *ipiv, logical *colequ, float *c__, scomplex *b, aocl_int64_t *ldb,
                              scomplex *y, aocl_int64_t *ldy, float *berr_out__, aocl_int64_t *n_norms__,
                              float *errs_n__, float *errs_c__, scomplex *res, float *ayb,
                              scomplex *dy, scomplex *y_tail__, float *rcond, aocl_int64_t *ithresh,
                              float *rthresh, float *dz_ub__, logical *ignore_cwise__,
                              aocl_int64_t *info);
float cla_gerpvgrw_check(aocl_int64_t *n, aocl_int64_t *ncols, scomplex *a, aocl_int64_t *lda, scomplex *af,
                         aocl_int64_t *ldaf);
int cla_heamv_check(aocl_int64_t *uplo, aocl_int64_t *n, float *alpha, scomplex *a, aocl_int64_t *lda, scomplex *x,
                    aocl_int64_t *incx, float *beta, float *y, aocl_int64_t *incy);
float cla_hercond_c_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, aocl_int_t *ipiv, float *c__, logical *capply, aocl_int64_t *info,
                          scomplex *work, float *rwork);
float cla_hercond_x_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *x, aocl_int64_t *info, scomplex *work,
                          float *rwork);
int cla_herfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              scomplex *a, aocl_int64_t *lda, scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv,
                              logical *colequ, float *c__, scomplex *b, aocl_int64_t *ldb, scomplex *y,
                              aocl_int64_t *ldy, float *berr_out__, aocl_int64_t *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res,
                              float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond,
                              aocl_int64_t *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
float cla_herpvgrw_check(char *uplo, aocl_int64_t *n, aocl_int64_t *info, scomplex *a, aocl_int64_t *lda,
                         scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *work);
int cla_lin_berr_check(aocl_int64_t *n, aocl_int64_t *nz, aocl_int64_t *nrhs, scomplex *res, float *ayb,
                       float *berr);
float cla_porcond_c_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, float *c__, logical *capply, aocl_int64_t *info, scomplex *work,
                          float *rwork);
float cla_porcond_x_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, scomplex *x, aocl_int64_t *info, scomplex *work, float *rwork);
int cla_porfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              scomplex *a, aocl_int64_t *lda, scomplex *af, aocl_int64_t *ldaf,
                              logical *colequ, float *c__, scomplex *b, aocl_int64_t *ldb, scomplex *y,
                              aocl_int64_t *ldy, float *berr_out__, aocl_int64_t *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res,
                              float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond,
                              aocl_int64_t *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
float cla_porpvgrw_check(char *uplo, aocl_int64_t *ncols, scomplex *a, aocl_int64_t *lda, scomplex *af,
                         aocl_int64_t *ldaf, float *work);
int cla_syamv_check(aocl_int64_t *uplo, aocl_int64_t *n, float *alpha, scomplex *a, aocl_int64_t *lda, scomplex *x,
                    aocl_int64_t *incx, float *beta, float *y, aocl_int64_t *incy);
float cla_syrcond_c_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, aocl_int_t *ipiv, float *c__, logical *capply, aocl_int64_t *info,
                          scomplex *work, float *rwork);
float cla_syrcond_x_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *af,
                          aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *x, aocl_int64_t *info, scomplex *work,
                          float *rwork);
int cla_syrfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              scomplex *a, aocl_int64_t *lda, scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv,
                              logical *colequ, float *c__, scomplex *b, aocl_int64_t *ldb, scomplex *y,
                              aocl_int64_t *ldy, float *berr_out__, aocl_int64_t *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, scomplex *res,
                              float *ayb, scomplex *dy, scomplex *y_tail__, float *rcond,
                              aocl_int64_t *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
float cla_syrpvgrw_check(char *uplo, aocl_int64_t *n, aocl_int64_t *info, scomplex *a, aocl_int64_t *lda,
                         scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *work);
int cla_wwaddw_check(aocl_int64_t *n, scomplex *x, scomplex *y, scomplex *w);
int clabrd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, scomplex *a, aocl_int64_t *lda, float *d__,
                 float *e, scomplex *tauq, scomplex *taup, scomplex *x, aocl_int64_t *ldx, scomplex *y,
                 aocl_int64_t *ldy);
int clacgv_check(aocl_int64_t *n, scomplex *x, aocl_int64_t *incx);
int clacn2_check(aocl_int64_t *n, scomplex *v, scomplex *x, float *est, aocl_int64_t *kase, aocl_int64_t *isave);
int clacon_check(aocl_int64_t *n, scomplex *v, scomplex *x, float *est, aocl_int64_t *kase);
int clacp2_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb);
int clacpy_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb);
int clacrm_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb,
                 scomplex *c__, aocl_int64_t *ldc, float *rwork);
int clacrt_check(aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx, scomplex *cy, aocl_int64_t *incy,
                 scomplex *c__, scomplex *s);
void cladiv_check(scomplex *ret_val, scomplex *x, scomplex *y);
int claed0_check(aocl_int64_t *qsiz, aocl_int64_t *n, float *d__, float *e, scomplex *q, aocl_int64_t *ldq,
                 scomplex *qstore, aocl_int64_t *ldqs, float *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int claed7_check(aocl_int64_t *n, aocl_int64_t *cutpnt, aocl_int64_t *qsiz, aocl_int64_t *tlvls, aocl_int64_t *curlvl,
                 aocl_int64_t *curpbm, float *d__, scomplex *q, aocl_int64_t *ldq, float *rho, aocl_int64_t *indxq,
                 float *qstore, aocl_int64_t *qptr, aocl_int64_t *prmptr, aocl_int64_t *perm, aocl_int64_t *givptr,
                 aocl_int64_t *givcol, float *givnum, scomplex *work, float *rwork, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int claed8_check(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *qsiz, scomplex *q, aocl_int64_t *ldq, float *d__,
                 float *rho, aocl_int64_t *cutpnt, float *z__, float *dlamda, scomplex *q2,
                 aocl_int64_t *ldq2, float *w, aocl_int64_t *indxp, aocl_int64_t *indx, aocl_int64_t *indxq,
                 aocl_int64_t *perm, aocl_int64_t *givptr, aocl_int64_t *givcol, float *givnum, aocl_int64_t *info);
int claein_check(logical *rightv, logical *noinit, aocl_int64_t *n, scomplex *h__, aocl_int64_t *ldh,
                 scomplex *w, scomplex *v, scomplex *b, aocl_int64_t *ldb, float *rwork, float *eps3,
                 float *smlnum, aocl_int64_t *info);
int claesy_check(scomplex *a, scomplex *b, scomplex *c__, scomplex *rt1, scomplex *rt2,
                 scomplex *evscal, scomplex *cs1, scomplex *sn1);
int claev2_check(scomplex *a, scomplex *b, scomplex *c__, float *rt1, float *rt2, float *cs1,
                 scomplex *sn1);
int clag2z_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *sa, aocl_int64_t *ldsa, dcomplex *a, aocl_int64_t *lda,
                 aocl_int64_t *info);
int clags2_check(logical *upper, float *a1, scomplex *a2, float *a3, float *b1, scomplex *b2,
                 float *b3, float *csu, scomplex *snu, float *csv, scomplex *snv, float *csq,
                 scomplex *snq);
int clagtm_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *alpha, scomplex *dl, scomplex *d__,
                 scomplex *du, scomplex *x, aocl_int64_t *ldx, float *beta, scomplex *b, aocl_int64_t *ldb);
int clahef_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, scomplex *a, aocl_int64_t *lda,
                 aocl_int_t *ipiv, scomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int clahef_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, scomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, scomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int clahqr_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 scomplex *h__, aocl_int64_t *ldh, scomplex *w, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *info);
int clahr2_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *t, aocl_int64_t *ldt, scomplex *y, aocl_int64_t *ldy);
int clahrd_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *t, aocl_int64_t *ldt, scomplex *y, aocl_int64_t *ldy);
int claic1_check(aocl_int64_t *job, aocl_int64_t *j, scomplex *x, float *sest, scomplex *w, scomplex *gamma,
                 float *sestpr, scomplex *s, scomplex *c__);
int clals0_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *nrhs,
                 scomplex *b, aocl_int64_t *ldb, scomplex *bx, aocl_int64_t *ldbx, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol, float *givnum, aocl_int64_t *ldgnum,
                 float *poles, float *difl, float *difr, float *z__, aocl_int64_t *k, float *c__,
                 float *s, float *rwork, aocl_int64_t *info);
int clalsa_check(aocl_int64_t *icompq, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *b,
                 aocl_int64_t *ldb, scomplex *bx, aocl_int64_t *ldbx, float *u, aocl_int64_t *ldu, float *vt,
                 aocl_int64_t *k, float *difl, float *difr, float *z__, float *poles, aocl_int64_t *givptr,
                 aocl_int64_t *givcol, aocl_int64_t *ldgcol, aocl_int64_t *perm, float *givnum, float *c__,
                 float *s, float *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int clalsd_check(char *uplo, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, float *e,
                 scomplex *b, aocl_int64_t *ldb, float *rcond, aocl_int64_t *rank, scomplex *work,
                 float *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
float clangb_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab, aocl_int64_t *ldab,
                   float *work);
float clange_check(char *norm, aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *work);
float clangt_check(char *norm, aocl_int64_t *n, scomplex *dl, scomplex *d__, scomplex *du);
float clanhb_check(char *norm, char *uplo, aocl_int64_t *n, aocl_int64_t *k, scomplex *ab, aocl_int64_t *ldab,
                   float *work);
float clanhe_check(char *norm, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *work);
float clanhf_check(char *norm, char *transr, char *uplo, aocl_int64_t *n, scomplex *a, float *work);
float clanhp_check(char *norm, char *uplo, aocl_int64_t *n, scomplex *ap, float *work);
float clanhs_check(char *norm, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *work);
float clanht_check(char *norm, aocl_int64_t *n, float *d__, scomplex *e);
float clansb_check(char *norm, char *uplo, aocl_int64_t *n, aocl_int64_t *k, scomplex *ab, aocl_int64_t *ldab,
                   float *work);
float clansp_check(char *norm, char *uplo, aocl_int64_t *n, scomplex *ap, float *work);
float clansy_check(char *norm, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *work);
float clantb_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *k, scomplex *ab,
                   aocl_int64_t *ldab, float *work);
float clantp_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, scomplex *ap, float *work);
float clantr_check(char *norm, char *uplo, char *diag, aocl_int64_t *m, aocl_int64_t *n, scomplex *a,
                   aocl_int64_t *lda, float *work);
int clapll_check(aocl_int64_t *n, scomplex *x, aocl_int64_t *incx, scomplex *y, aocl_int64_t *incy, float *ssmin);
int clapmr_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, scomplex *x, aocl_int64_t *ldx, aocl_int64_t *k);
int clapmt_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, scomplex *x, aocl_int64_t *ldx, aocl_int64_t *k);
int claqgb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, scomplex *ab, aocl_int64_t *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *equed);
int claqge_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, char *equed);
int claqhb_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab, float *s,
                 float *scond, float *amax, char *equed);
int claqhe_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *s, float *scond,
                 float *amax, char *equed);
int claqhp_check(char *uplo, aocl_int64_t *n, scomplex *ap, float *s, float *scond, float *amax,
                 char *equed);
int claqp2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, scomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt,
                 scomplex *tau, float *vn1, float *vn2, scomplex *work);
int claqps_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, aocl_int64_t *nb, aocl_int64_t *kb, scomplex *a,
                 aocl_int64_t *lda, aocl_int64_t *jpvt, scomplex *tau, float *vn1, float *vn2, scomplex *auxv,
                 scomplex *f, aocl_int64_t *ldf);
int claqr0_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 scomplex *h__, aocl_int64_t *ldh, scomplex *w, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 scomplex *z__, aocl_int64_t *ldz, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int claqr1_check(aocl_int64_t *n, scomplex *h__, aocl_int64_t *ldh, scomplex *s1, scomplex *s2, scomplex *v);
int claqr2_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, scomplex *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, scomplex *sh, scomplex *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, scomplex *t, aocl_int64_t *ldt, aocl_int64_t *nv, scomplex *wv,
                 aocl_int64_t *ldwv, scomplex *work, aocl_int64_t *lwork);
int claqr3_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, scomplex *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, scomplex *sh, scomplex *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, scomplex *t, aocl_int64_t *ldt, aocl_int64_t *nv, scomplex *wv,
                 aocl_int64_t *ldwv, scomplex *work, aocl_int64_t *lwork);
int claqr4_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 scomplex *h__, aocl_int64_t *ldh, scomplex *w, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 scomplex *z__, aocl_int64_t *ldz, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int claqr5_check(logical *wantt, logical *wantz, aocl_int64_t *kacc22, aocl_int64_t *n, aocl_int64_t *ktop,
                 aocl_int64_t *kbot, aocl_int64_t *nshfts, scomplex *s, scomplex *h__, aocl_int64_t *ldh,
                 aocl_int64_t *iloz, aocl_int64_t *ihiz, scomplex *z__, aocl_int64_t *ldz, scomplex *v,
                 aocl_int64_t *ldv, scomplex *u, aocl_int64_t *ldu, aocl_int64_t *nv, scomplex *wv, aocl_int64_t *ldwv,
                 aocl_int64_t *nh, scomplex *wh, aocl_int64_t *ldwh);
int claqsb_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab, float *s,
                 float *scond, float *amax, char *equed);
int claqsp_check(char *uplo, aocl_int64_t *n, scomplex *ap, float *s, float *scond, float *amax,
                 char *equed);
int claqsy_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *s, float *scond,
                 float *amax, char *equed);
int clar1v_check(aocl_int64_t *n, aocl_int64_t *b1, aocl_int64_t *bn, float *lambda, float *d__, float *l,
                 float *ld, float *lld, float *pivmin, float *gaptol, scomplex *z__,
                 logical *wantnc, aocl_int64_t *negcnt, float *ztz, float *mingma, aocl_int64_t *r__,
                 aocl_int64_t *isuppz, float *nrminv, float *resid, float *rqcorr, float *work);
int clar2v_check(aocl_int64_t *n, scomplex *x, scomplex *y, scomplex *z__, aocl_int64_t *incx, float *c__,
                 scomplex *s, aocl_int64_t *incc);
int clarcm_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb,
                 scomplex *c__, aocl_int64_t *ldc, float *rwork);
int clarf_check(char *side, aocl_int64_t *m, aocl_int64_t *n, scomplex *v, aocl_int64_t *incv, scomplex *tau,
                scomplex *c__, aocl_int64_t *ldc, scomplex *work);
int clarfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, scomplex *v, aocl_int64_t *ldv, scomplex *t, aocl_int64_t *ldt, scomplex *c__,
                 aocl_int64_t *ldc, scomplex *work, aocl_int64_t *ldwork);
int clarfg_check(aocl_int64_t *n, scomplex *alpha, scomplex *x, aocl_int64_t *incx, scomplex *tau);
int clarfgp_check(aocl_int64_t *n, scomplex *alpha, scomplex *x, aocl_int64_t *incx, scomplex *tau);
int clarft_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, scomplex *v, aocl_int64_t *ldv,
                 scomplex *tau, scomplex *t, aocl_int64_t *ldt);
int clarfx_check(char *side, aocl_int64_t *m, aocl_int64_t *n, scomplex *v, scomplex *tau, scomplex *c__,
                 aocl_int64_t *ldc, scomplex *work);
int clargv_check(aocl_int64_t *n, scomplex *x, aocl_int64_t *incx, scomplex *y, aocl_int64_t *incy, float *c__,
                 aocl_int64_t *incc);
int clarnv_check(aocl_int64_t *idist, aocl_int64_t *iseed, aocl_int64_t *n, scomplex *x);
int clarrv_check(aocl_int64_t *n, float *vl, float *vu, float *d__, float *l, float *pivmin,
                 aocl_int64_t *isplit, aocl_int64_t *m, aocl_int64_t *dol, aocl_int64_t *dou, float *minrgp,
                 float *rtol1, float *rtol2, float *w, float *werr, float *wgap, aocl_int64_t *iblock,
                 aocl_int64_t *indexw, float *gers, scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int clarscl2_check(aocl_int64_t *m, aocl_int64_t *n, float *d__, scomplex *x, aocl_int64_t *ldx);
int clartg_check(scomplex *f, scomplex *g, float *cs, scomplex *sn, scomplex *r__);
int clartv_check(aocl_int64_t *n, scomplex *x, aocl_int64_t *incx, scomplex *y, aocl_int64_t *incy, float *c__,
                 scomplex *s, aocl_int64_t *incc);
int clarz_check(char *side, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, scomplex *v, aocl_int64_t *incv,
                scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work);
int clarzb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, scomplex *v, aocl_int64_t *ldv, scomplex *t, aocl_int64_t *ldt,
                 scomplex *c__, aocl_int64_t *ldc, scomplex *work, aocl_int64_t *ldwork);
int clarzt_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, scomplex *v, aocl_int64_t *ldv,
                 scomplex *tau, scomplex *t, aocl_int64_t *ldt);
int clascl_check(char *type__, aocl_int64_t *kl, aocl_int64_t *ku, float *cfrom, float *cto, aocl_int64_t *m,
                 aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int clascl2_check(aocl_int64_t *m, aocl_int64_t *n, float *d__, scomplex *x, aocl_int64_t *ldx);
int claset_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, scomplex *alpha, scomplex *beta, scomplex *a,
                 aocl_int64_t *lda);
int clasr_check(char *side, char *pivot, char *direct, aocl_int64_t *m, aocl_int64_t *n, float *c__, float *s,
                scomplex *a, aocl_int64_t *lda);
int classq_check(aocl_int64_t *n, scomplex *x, aocl_int64_t *incx, float *scale, float *sumsq);
int claswp_check(aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *k1, aocl_int64_t *k2, aocl_int_t *ipiv,
                 aocl_int64_t *incx);
int clasyf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, scomplex *a, aocl_int64_t *lda,
                 aocl_int_t *ipiv, scomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int clasyf_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, scomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, scomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int clatbs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, aocl_int64_t *kd,
                 scomplex *ab, aocl_int64_t *ldab, scomplex *x, float *scale, float *cnorm,
                 aocl_int64_t *info);
int clatdf_check(aocl_int64_t *ijob, aocl_int64_t *n, scomplex *z__, aocl_int64_t *ldz, scomplex *rhs,
                 float *rdsum, float *rdscal, aocl_int_t *ipiv, aocl_int64_t *jpiv);
int clatps_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, scomplex *ap,
                 scomplex *x, float *scale, float *cnorm, aocl_int64_t *info);
int clatrd_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, scomplex *a, aocl_int64_t *lda, float *e,
                 scomplex *tau, scomplex *w, aocl_int64_t *ldw);
int clatrs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, scomplex *a,
                 aocl_int64_t *lda, scomplex *x, float *scale, float *cnorm, aocl_int64_t *info);
int clatrz_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work);
int clatzm_check(char *side, aocl_int64_t *m, aocl_int64_t *n, scomplex *v, aocl_int64_t *incv, scomplex *tau,
                 scomplex *c1, scomplex *c2, aocl_int64_t *ldc, scomplex *work);
int clauu2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int clauum_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int cpbcon_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab, float *anorm,
                 float *rcond, scomplex *work, float *rwork, aocl_int64_t *info);
int cpbequ_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab, float *s,
                 float *scond, float *amax, aocl_int64_t *info);
int cpbrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, scomplex *ab, aocl_int64_t *ldab,
                 scomplex *afb, aocl_int64_t *ldafb, scomplex *b, aocl_int64_t *ldb, scomplex *x,
                 aocl_int64_t *ldx, float *ferr, float *berr, scomplex *work, float *rwork,
                 aocl_int64_t *info);
int cpbstf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int cpbsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, scomplex *ab, aocl_int64_t *ldab,
                scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int cpbsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, scomplex *ab,
                 aocl_int64_t *ldab, scomplex *afb, aocl_int64_t *ldafb, char *equed, float *s, scomplex *b,
                 aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond, float *ferr, float *berr,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int cpbtf2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int cpbtrf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int cpbtrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, scomplex *ab, aocl_int64_t *ldab,
                 scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int cpftrf_check(char *transr, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *info);
int cpftri_check(char *transr, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *info);
int cpftrs_check(char *transr, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int cpocon_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *anorm, float *rcond,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int cpoequ_check(aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *s, float *scond, float *amax,
                 aocl_int64_t *info);
int cpoequb_check(aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *s, float *scond, float *amax,
                  aocl_int64_t *info);
int cporfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *af,
                 aocl_int64_t *ldaf, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr,
                 float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int cporfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, float *s, scomplex *b, aocl_int64_t *ldb, scomplex *x,
                  aocl_int64_t *ldx, float *rcond, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  scomplex *work, float *rwork, aocl_int64_t *info);
int cposv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int cposvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                 scomplex *af, aocl_int64_t *ldaf, char *equed, float *s, scomplex *b, aocl_int64_t *ldb,
                 scomplex *x, aocl_int64_t *ldx, float *rcond, float *ferr, float *berr, scomplex *work,
                 float *rwork, aocl_int64_t *info);
int cposvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, char *equed, float *s, scomplex *b, aocl_int64_t *ldb,
                  scomplex *x, aocl_int64_t *ldx, float *rcond, float *rpvgrw, float *berr,
                  aocl_int64_t *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  aocl_int64_t *nparams, float *params, scomplex *work, float *rwork, aocl_int64_t *info);
int cpotf2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int cpotrf_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int cpotri_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int cpotrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int cppcon_check(char *uplo, aocl_int64_t *n, scomplex *ap, float *anorm, float *rcond, scomplex *work,
                 float *rwork, aocl_int64_t *info);
int cppequ_check(char *uplo, aocl_int64_t *n, scomplex *ap, float *s, float *scond, float *amax,
                 aocl_int64_t *info);
int cpprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *afp, scomplex *b,
                 aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr, float *berr, scomplex *work,
                 float *rwork, aocl_int64_t *info);
int cppsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int cppsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *afp,
                 char *equed, float *s, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx,
                 float *rcond, float *ferr, float *berr, scomplex *work, float *rwork,
                 aocl_int64_t *info);
int cpptrf_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int64_t *info);
int cpptri_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int64_t *info);
int cpptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *b, aocl_int64_t *ldb,
                 aocl_int64_t *info);
int cpstf2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 float *tol, float *work, aocl_int64_t *info);
int cpstrf_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 float *tol, float *work, aocl_int64_t *info);
int cptcon_check(aocl_int64_t *n, float *d__, scomplex *e, float *anorm, float *rcond, float *rwork,
                 aocl_int64_t *info);
int cpteqr_check(char *compz, aocl_int64_t *n, float *d__, float *e, scomplex *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *info);
int cptrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, scomplex *e, float *df,
                 scomplex *ef, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr,
                 float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int cptsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, scomplex *e, scomplex *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int cptsvx_check(char *fact, aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, scomplex *e, float *df,
                 scomplex *ef, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int cpttrf_check(aocl_int64_t *n, float *d__, scomplex *e, aocl_int64_t *info);
int cpttrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, scomplex *e, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int cptts2_check(aocl_int64_t *iuplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, scomplex *e, scomplex *b,
                 aocl_int64_t *ldb);
int crot_check(aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx, scomplex *cy, aocl_int64_t *incy, float *c__,
               scomplex *s);
int cspcon_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int_t *ipiv, float *anorm, float *rcond,
                 scomplex *work, aocl_int64_t *info);
int cspmv_check(char *uplo, aocl_int64_t *n, scomplex *alpha, scomplex *ap, scomplex *x, aocl_int64_t *incx,
                scomplex *beta, scomplex *y, aocl_int64_t *incy);
int cspr_check(char *uplo, aocl_int64_t *n, scomplex *alpha, scomplex *x, aocl_int64_t *incx, scomplex *ap);
int csprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *afp, aocl_int_t *ipiv,
                 scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int cspsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, aocl_int_t *ipiv, scomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int cspsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, scomplex *afp,
                 aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int csptrf_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int_t *ipiv, aocl_int64_t *info);
int csptri_check(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int_t *ipiv, scomplex *work,
                 aocl_int64_t *info);
int csptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap, aocl_int_t *ipiv, scomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int csrscl_check(aocl_int64_t *n, float *sa, scomplex *sx, aocl_int64_t *incx);
int cstedc_check(char *compz, aocl_int64_t *n, float *d__, float *e, scomplex *z__, aocl_int64_t *ldz,
                 scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int cstegr_check(char *jobz, char *range, aocl_int64_t *n, float *d__, float *e, float *vl, float *vu,
                 aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, scomplex *z__,
                 aocl_int64_t *ldz, aocl_int64_t *isuppz, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int cstein_check(aocl_int64_t *n, float *d__, float *e, aocl_int64_t *m, float *w, aocl_int64_t *iblock,
                 aocl_int64_t *isplit, scomplex *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int cstemr_check(char *jobz, char *range, aocl_int64_t *n, float *d__, float *e, float *vl, float *vu,
                 aocl_int64_t *il, aocl_int64_t *iu, aocl_int64_t *m, float *w, scomplex *z__, aocl_int64_t *ldz,
                 aocl_int64_t *nzc, aocl_int64_t *isuppz, logical *tryrac, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int csteqr_check(char *compz, aocl_int64_t *n, float *d__, float *e, scomplex *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *info);
int csycon_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *anorm,
                 float *rcond, scomplex *work, aocl_int64_t *info);
int csycon_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      float *anorm, float *rcond, scomplex *work, aocl_int64_t *info);
int csyconv_check(char *uplo, char *way, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  scomplex *work, aocl_int64_t *info);
int csyequb_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, float *s, float *scond,
                  float *amax, scomplex *work, aocl_int64_t *info);
int csymv_check(char *uplo, aocl_int64_t *n, scomplex *alpha, scomplex *a, aocl_int64_t *lda, scomplex *x,
                aocl_int64_t *incx, scomplex *beta, scomplex *y, aocl_int64_t *incy);
int csyr_check(char *uplo, aocl_int64_t *n, scomplex *alpha, scomplex *x, aocl_int64_t *incx, scomplex *a,
               aocl_int64_t *lda);
int csyrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, scomplex *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int csyrfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *s, scomplex *b, aocl_int64_t *ldb,
                  scomplex *x, aocl_int64_t *ldx, float *rcond, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  scomplex *work, float *rwork, aocl_int64_t *info);
int csysv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int csysv_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                     aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *lwork,
                     aocl_int64_t *info);
int csysvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                 scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, scomplex *x,
                 aocl_int64_t *ldx, float *rcond, float *ferr, float *berr, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *info);
int csysvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                  scomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, float *s, scomplex *b,
                  aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *rcond, float *rpvgrw, float *berr,
                  aocl_int64_t *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  aocl_int64_t *nparams, float *params, scomplex *work, float *rwork, aocl_int64_t *info);
int csyswapr_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *i1, aocl_int64_t *i2);
int csytf2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int csytf2_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      aocl_int64_t *info);
int csytrf_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int csytrf_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int csytri_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                 aocl_int64_t *info);
int csytri2_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int csytri2x_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, scomplex *work,
                   aocl_int64_t *nb, aocl_int64_t *info);
int csytri_rook_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      scomplex *work, aocl_int64_t *info);
int csytrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int csytrs2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *info);
int csytrs_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ctbcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *kd, scomplex *ab,
                 aocl_int64_t *ldab, float *rcond, scomplex *work, float *rwork, aocl_int64_t *info);
int ctbrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 scomplex *ab, aocl_int64_t *ldab, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int ctbtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 scomplex *ab, aocl_int64_t *ldab, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ctfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, aocl_int64_t *m,
                aocl_int64_t *n, scomplex *alpha, scomplex *a, scomplex *b, aocl_int64_t *ldb);
int ctftri_check(char *transr, char *uplo, char *diag, aocl_int64_t *n, scomplex *a, aocl_int64_t *info);
int ctfttp_check(char *transr, char *uplo, aocl_int64_t *n, scomplex *arf, scomplex *ap, aocl_int64_t *info);
int ctfttr_check(char *transr, char *uplo, aocl_int64_t *n, scomplex *arf, scomplex *a, aocl_int64_t *lda,
                 aocl_int64_t *info);
int ctgevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, scomplex *s, aocl_int64_t *lds,
                 scomplex *p, aocl_int64_t *ldp, scomplex *vl, aocl_int64_t *ldvl, scomplex *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m, scomplex *work, float *rwork,
                 aocl_int64_t *info);
int ctgex2_check(logical *wantq, logical *wantz, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, scomplex *q, aocl_int64_t *ldq, scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *j1,
                 aocl_int64_t *info);
int ctgexc_check(logical *wantq, logical *wantz, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
                 aocl_int64_t *ldb, scomplex *q, aocl_int64_t *ldq, scomplex *z__, aocl_int64_t *ldz,
                 aocl_int64_t *ifst, aocl_int64_t *ilst, aocl_int64_t *info);
int ctgsen_check(aocl_int64_t *ijob, logical *wantq, logical *wantz, logical *select, aocl_int64_t *n,
                 scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *alpha,
                 scomplex *beta, scomplex *q, aocl_int64_t *ldq, scomplex *z__, aocl_int64_t *ldz, aocl_int64_t *m,
                 float *pl, float *pr, float *dif, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int ctgsja_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, aocl_int64_t *k,
                 aocl_int64_t *l, scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, float *tola,
                 float *tolb, float *alpha, float *beta, scomplex *u, aocl_int64_t *ldu, scomplex *v,
                 aocl_int64_t *ldv, scomplex *q, aocl_int64_t *ldq, scomplex *work, aocl_int64_t *ncycle,
                 aocl_int64_t *info);
int ctgsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 scomplex *b, aocl_int64_t *ldb, scomplex *vl, aocl_int64_t *ldvl, scomplex *vr,
                 aocl_int64_t *ldvr, float *s, float *dif, aocl_int64_t *mm, aocl_int64_t *m, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int ctgsy2_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 scomplex *b, aocl_int64_t *ldb, scomplex *c__, aocl_int64_t *ldc, scomplex *d__,
                 aocl_int64_t *ldd, scomplex *e, aocl_int64_t *lde, scomplex *f, aocl_int64_t *ldf, float *scale,
                 float *rdsum, float *rdscal, aocl_int64_t *info);
int ctgsyl_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 scomplex *b, aocl_int64_t *ldb, scomplex *c__, aocl_int64_t *ldc, scomplex *d__,
                 aocl_int64_t *ldd, scomplex *e, aocl_int64_t *lde, scomplex *f, aocl_int64_t *ldf, float *scale,
                 float *dif, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int ctpcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, scomplex *ap, float *rcond,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int ctpmqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                  aocl_int64_t *nb, scomplex *v, aocl_int64_t *ldv, scomplex *t, aocl_int64_t *ldt, scomplex *a,
                  aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *work, aocl_int64_t *info);
int ctpqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, aocl_int64_t *nb, scomplex *a, aocl_int64_t *lda,
                 scomplex *b, aocl_int64_t *ldb, scomplex *t, aocl_int64_t *ldt, scomplex *work,
                 aocl_int64_t *info);
int ctpqrt2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, scomplex *a, aocl_int64_t *lda, scomplex *b,
                  aocl_int64_t *ldb, scomplex *t, aocl_int64_t *ldt, aocl_int64_t *info);
int ctprfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, scomplex *v, aocl_int64_t *ldv, scomplex *t, aocl_int64_t *ldt,
                 scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *work,
                 aocl_int64_t *ldwork);
int ctprfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap,
                 scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr, float *berr,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int ctptri_check(char *uplo, char *diag, aocl_int64_t *n, scomplex *ap, aocl_int64_t *info);
int ctptrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *ap,
                 scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ctpttf_check(char *transr, char *uplo, aocl_int64_t *n, scomplex *ap, scomplex *arf, aocl_int64_t *info);
int ctpttr_check(char *uplo, aocl_int64_t *n, scomplex *ap, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int ctrcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                 float *rcond, scomplex *work, float *rwork, aocl_int64_t *info);
int ctrevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, scomplex *t, aocl_int64_t *ldt,
                 scomplex *vl, aocl_int64_t *ldvl, scomplex *vr, aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m,
                 scomplex *work, float *rwork, aocl_int64_t *info);
int ctrexc_check(char *compq, aocl_int64_t *n, scomplex *t, aocl_int64_t *ldt, scomplex *q, aocl_int64_t *ldq,
                 aocl_int64_t *ifst, aocl_int64_t *ilst, aocl_int64_t *info);
int ctrrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a,
                 aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx, float *ferr,
                 float *berr, scomplex *work, float *rwork, aocl_int64_t *info);
int ctrsen_check(char *job, char *compq, logical *select, aocl_int64_t *n, scomplex *t, aocl_int64_t *ldt,
                 scomplex *q, aocl_int64_t *ldq, scomplex *w, aocl_int64_t *m, float *s, float *sep,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int ctrsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, scomplex *t, aocl_int64_t *ldt,
                 scomplex *vl, aocl_int64_t *ldvl, scomplex *vr, aocl_int64_t *ldvr, float *s, float *sep,
                 aocl_int64_t *mm, aocl_int64_t *m, scomplex *work, aocl_int64_t *ldwork, float *rwork,
                 aocl_int64_t *info);
int ctrsyl_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, scomplex *a,
                 aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *c__, aocl_int64_t *ldc, float *scale,
                 aocl_int64_t *info);
int ctrsyl3_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, scomplex *a,
                  aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, scomplex *c__, aocl_int64_t *ldc, real *scale,
                  real *swork, aocl_int64_t *ldswork, aocl_int64_t *info);
int ctrti2_check(char *uplo, char *diag, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int ctrtri_check(char *uplo, char *diag, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int ctrtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, scomplex *a,
                 aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ctrttf_check(char *transr, char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *arf,
                 aocl_int64_t *info);
int ctrttp_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *ap, aocl_int64_t *info);
int ctzrqf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, aocl_int64_t *info);
int ctzrzf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cunbdb_check(char *trans, char *signs, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, scomplex *x11,
                 aocl_int64_t *ldx11, scomplex *x12, aocl_int64_t *ldx12, scomplex *x21, aocl_int64_t *ldx21,
                 scomplex *x22, aocl_int64_t *ldx22, float *theta, float *phi, scomplex *taup1,
                 scomplex *taup2, scomplex *tauq1, scomplex *tauq2, scomplex *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int cunbdb1_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, scomplex *x11, aocl_int64_t *ldx11, scomplex *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunbdb2_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, scomplex *x11, aocl_int64_t *ldx11, scomplex *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunbdb3_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, scomplex *x11, aocl_int64_t *ldx11, scomplex *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunbdb4_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, scomplex *x11, aocl_int64_t *ldx11, scomplex *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, scomplex *taup1, scomplex *taup2,
                  scomplex *tauq1, scomplex *phantom, scomplex *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);
int cunbdb5_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, scomplex *x1, aocl_int64_t *incx1, scomplex *x2,
                  aocl_int64_t *incx2, scomplex *q1, aocl_int64_t *ldq1, scomplex *q2, aocl_int64_t *ldq2,
                  scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunbdb6_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, scomplex *x1, aocl_int64_t *incx1, scomplex *x2,
                  aocl_int64_t *incx2, scomplex *q1, aocl_int64_t *ldq1, scomplex *q2, aocl_int64_t *ldq2,
                  scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cuncsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, scomplex *x11, aocl_int64_t *ldx11, scomplex *x12,
                 aocl_int64_t *ldx12, scomplex *x21, aocl_int64_t *ldx21, scomplex *x22, aocl_int64_t *ldx22,
                 float *theta, scomplex *u1, aocl_int64_t *ldu1, scomplex *u2, aocl_int64_t *ldu2,
                 scomplex *v1t, aocl_int64_t *ldv1t, scomplex *v2t, aocl_int64_t *ldv2t, scomplex *work,
                 aocl_int64_t *lwork, float *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *info);
int cuncsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q,
                     scomplex *x11, aocl_int64_t *ldx11, scomplex *x21, aocl_int64_t *ldx21, float *theta,
                     scomplex *u1, aocl_int64_t *ldu1, scomplex *u2, aocl_int64_t *ldu2, scomplex *v1t,
                     aocl_int64_t *ldv1t, scomplex *work, aocl_int64_t *lwork, float *rwork, aocl_int64_t *lrwork,
                     aocl_int64_t *iwork, aocl_int64_t *info);
int cung2l_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *info);
int cung2r_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *info);
int cungbr_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda,
                 scomplex *tau, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunghr_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cungl2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *info);
int cunglq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cungql_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cungqr_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cungr2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *info);
int cungrq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cungtr_check(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cunm2l_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *info);
int cunm2r_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *info);
int cunmbr_check(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                 scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunmhr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunml2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *info);
int cunmlq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cunmql_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cunmqr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cunmr2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *info);
int cunmr3_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                 scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                 scomplex *work, aocl_int64_t *info);
int cunmrq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cunmrz_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                 scomplex *a, aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int cunmtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, scomplex *a,
                 aocl_int64_t *lda, scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int cupgtr_check(char *uplo, aocl_int64_t *n, scomplex *ap, scomplex *tau, scomplex *q, aocl_int64_t *ldq,
                 scomplex *work, aocl_int64_t *info);
int cupmtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, scomplex *ap,
                 scomplex *tau, scomplex *c__, aocl_int64_t *ldc, scomplex *work, aocl_int64_t *info);
int dbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, aocl_int64_t *m,
                 aocl_int64_t *p, aocl_int64_t *q, double *theta, double *phi, double *u1, aocl_int64_t *ldu1,
                 double *u2, aocl_int64_t *ldu2, double *v1t, aocl_int64_t *ldv1t, double *v2t,
                 aocl_int64_t *ldv2t, double *b11d, double *b11e, double *b12d, double *b12e,
                 double *b21d, double *b21e, double *b22d, double *b22e, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dbdsdc_check(char *uplo, char *compq, aocl_int64_t *n, double *d__, double *e, double *u,
                 aocl_int64_t *ldu, double *vt, aocl_int64_t *ldvt, double *q, aocl_int64_t *iq, double *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dbdsqr_check(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru, aocl_int64_t *ncc, double *d__,
                 double *e, double *vt, aocl_int64_t *ldvt, double *u, aocl_int64_t *ldu, double *c__,
                 aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int ddisna_check(char *job, aocl_int64_t *m, aocl_int64_t *n, double *d__, double *sep, aocl_int64_t *info);
int dgbbrd_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ncc, aocl_int64_t *kl, aocl_int64_t *ku,
                 double *ab, aocl_int64_t *ldab, double *d__, double *e, double *q, aocl_int64_t *ldq,
                 double *pt, aocl_int64_t *ldpt, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *info);
int dgbcon_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, double *anorm, double *rcond, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dgbequ_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab, aocl_int64_t *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 aocl_int64_t *info);
int dgbequb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab, aocl_int64_t *ldab,
                  double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                  aocl_int64_t *info);
int dgbrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, double *ab,
                 aocl_int64_t *ldab, double *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, double *b, aocl_int64_t *ldb,
                 double *x, aocl_int64_t *ldx, double *ferr, double *berr, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dgbrfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  double *ab, aocl_int64_t *ldab, double *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                  double *r__, double *c__, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx,
                  double *rcond, double *berr, aocl_int64_t *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, aocl_int64_t *nparams, double *params, double *work,
                  aocl_int64_t *iwork, aocl_int64_t *info);
int dgbsv_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, double *ab, aocl_int64_t *ldab,
                aocl_int_t *ipiv, double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dgbsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                 double *ab, aocl_int64_t *ldab, double *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, char *equed,
                 double *r__, double *c__, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx,
                 double *rcond, double *ferr, double *berr, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dgbsvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  double *ab, aocl_int64_t *ldab, double *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                  char *equed, double *r__, double *c__, double *b, aocl_int64_t *ldb, double *x,
                  aocl_int64_t *ldx, double *rcond, double *rpvgrw, double *berr, aocl_int64_t *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, aocl_int64_t *nparams,
                  double *params, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dgbtf2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int dgbtrf_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int dgbtrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, double *ab,
                 aocl_int64_t *ldab, aocl_int_t *ipiv, double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dgebak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *scale,
                 aocl_int64_t *m, double *v, aocl_int64_t *ldv, aocl_int64_t *info);
int dgebal_check(char *job, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 double *scale, aocl_int64_t *info);
int dgebd2_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *d__, double *e,
                 double *tauq, double *taup, double *work, aocl_int64_t *info);
int dgebrd_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *d__, double *e,
                 double *tauq, double *taup, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dgecon_check(char *norm, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *anorm, double *rcond,
                 double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dgeequ_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, aocl_int64_t *info);
int dgeequb_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *r__, double *c__,
                  double *rowcnd, double *colcnd, double *amax, aocl_int64_t *info);
int dgees_check(char *jobvs, char *sort, L_fp select, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                aocl_int64_t *sdim, double *wr, double *wi, double *vs, aocl_int64_t *ldvs, double *work,
                aocl_int64_t *lwork, logical *bwork, aocl_int64_t *info);
int dgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, aocl_int64_t *sdim, double *wr, double *wi, double *vs, aocl_int64_t *ldvs,
                 double *rconde, double *rcondv, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, logical *bwork, aocl_int64_t *info);
int dgeev_check(char *jobvl, char *jobvr, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *wr,
                double *wi, double *vl, aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr, double *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
int dgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *wr, double *wi, double *vl, aocl_int64_t *ldvl, double *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *ilo, aocl_int64_t *ihi, double *scale, double *abnrm,
                 double *rconde, double *rcondv, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dgegs_check(char *jobvsl, char *jobvsr, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                aocl_int64_t *ldb, double *alphar, double *alphai, double *beta, double *vsl,
                aocl_int64_t *ldvsl, double *vsr, aocl_int64_t *ldvsr, double *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
int dgegv_check(char *jobvl, char *jobvr, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                aocl_int64_t *ldb, double *alphar, double *alphai, double *beta, double *vl,
                aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr, double *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
int dgehd2_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *info);
int dgehrd_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *a, aocl_int64_t *lda, double *tau,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dgejsv_check(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, aocl_int64_t *m,
                 aocl_int64_t *n, double *a, aocl_int64_t *lda, double *sva, double *u, aocl_int64_t *ldu,
                 double *v, aocl_int64_t *ldv, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dgelq2_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *info);
int dgelqf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgels_check(char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                double *b, aocl_int64_t *ldb, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dgelsd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, double *s, double *rcond, aocl_int64_t *rank, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int dgelss_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, double *s, double *rcond, aocl_int64_t *rank, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgelsx_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, double *rcond, aocl_int64_t *rank, double *work,
                 aocl_int64_t *info);
int dgelsy_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, double *rcond, aocl_int64_t *rank, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgemqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb,
                  double *v, aocl_int64_t *ldv, double *t, aocl_int64_t *ldt, double *c__, aocl_int64_t *ldc,
                  double *work, aocl_int64_t *info);
int dgeql2_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *info);
int dgeqlf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgeqp3_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *jpvt, double *tau,
                 double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dgeqpf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *jpvt, double *tau,
                 double *work, aocl_int64_t *info);
int dgeqr2_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *info);
int dgeqr2p_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                  aocl_int64_t *info);
int dgeqrf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgeqrfp_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int dgeqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, double *a, aocl_int64_t *lda, double *t,
                 aocl_int64_t *ldt, double *work, aocl_int64_t *info);
int dgeqrt2_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int dgeqrt3_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int dgerfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dgerfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                  double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *r__, double *c__, double *b,
                  aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond, double *berr,
                  aocl_int64_t *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  aocl_int64_t *nparams, double *params, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dgerq2_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *info);
int dgerqf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgesc2_check(aocl_int64_t *n, double *a, aocl_int64_t *lda, double *rhs, aocl_int_t *ipiv, aocl_int64_t *jpiv,
                 double *scale);
int dgesdd_check(char *jobz, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *s, double *u,
                 aocl_int64_t *ldu, double *vt, aocl_int64_t *ldvt, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dgesdd_fla_check(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                     double *s, double *u, aocl_int64_t *ldu, double *vt, aocl_int64_t *ldvt, double *work,
                     aocl_int64_t *lwork, aocl_int64_t *info);
int dgesv_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int dgesvd_check(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda,
                 double *s, double *u, aocl_int64_t *ldu, double *vt, aocl_int64_t *ldvt, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgesvj_check(char *joba, char *jobu, char *jobv, aocl_int64_t *m, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *sva, aocl_int64_t *mv, double *v, aocl_int64_t *ldv, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dgesvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                 double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, double *r__, double *c__,
                 double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond, double *ferr,
                 double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dgesvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda,
                  double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, double *r__, double *c__,
                  double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond, double *rpvgrw,
                  double *berr, aocl_int64_t *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, aocl_int64_t *nparams, double *params, double *work,
                  aocl_int64_t *iwork, aocl_int64_t *info);
int dgetc2_check(aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *jpiv, aocl_int64_t *info);
int dgetf2_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int dgetrf_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int dgetrfnp_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dgetrfnpi_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nfact, double *a, aocl_int64_t *lda, aocl_int64_t *info);
int dgetri_check(aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dgetrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dggbak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *lscale,
                 double *rscale, aocl_int64_t *m, double *v, aocl_int64_t *ldv, aocl_int64_t *info);
int dggbal_check(char *job, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb,
                 aocl_int64_t *ilo, aocl_int64_t *ihi, double *lscale, double *rscale, double *work,
                 aocl_int64_t *info);
int dgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, aocl_int64_t *n, double *a,
                aocl_int64_t *lda, double *b, aocl_int64_t *ldb, aocl_int64_t *sdim, double *alphar,
                double *alphai, double *beta, double *vsl, aocl_int64_t *ldvsl, double *vsr,
                aocl_int64_t *ldvsr, double *work, aocl_int64_t *lwork, logical *bwork, aocl_int64_t *info);
int dggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, aocl_int64_t *n,
                 double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb, aocl_int64_t *sdim, double *alphar,
                 double *alphai, double *beta, double *vsl, aocl_int64_t *ldvsl, double *vsr,
                 aocl_int64_t *ldvsr, double *rconde, double *rcondv, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, logical *bwork, aocl_int64_t *info);
int dggev_check(char *jobvl, char *jobvr, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                aocl_int64_t *ldb, double *alphar, double *alphai, double *beta, double *vl,
                aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr, double *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
int dggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *alphar, double *alphai,
                 double *beta, double *vl, aocl_int64_t *ldvl, double *vr, aocl_int64_t *ldvr, aocl_int64_t *ilo,
                 aocl_int64_t *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm,
                 double *rconde, double *rcondv, double *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 logical *bwork, aocl_int64_t *info);
int dggglm_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, double *d__, double *x, double *y, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dgghrd_check(char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *a,
                 aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *q, aocl_int64_t *ldq, double *z__,
                 aocl_int64_t *ldz, aocl_int64_t *info);
int dgglse_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb, double *c__, double *d__, double *x, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dggqrf_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, double *a, aocl_int64_t *lda, double *taua,
                 double *b, aocl_int64_t *ldb, double *taub, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dggrqf_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *taua,
                 double *b, aocl_int64_t *ldb, double *taub, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dggsvd_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, aocl_int64_t *k,
                 aocl_int64_t *l, double *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *alpha,
                 double *beta, double *u, aocl_int64_t *ldu, double *v, aocl_int64_t *ldv, double *q,
                 aocl_int64_t *ldq, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dggsvp_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *b, aocl_int64_t *ldb, double *tola, double *tolb, aocl_int64_t *k,
                 aocl_int64_t *l, double *u, aocl_int64_t *ldu, double *v, aocl_int64_t *ldv, double *q,
                 aocl_int64_t *ldq, aocl_int64_t *iwork, double *tau, double *work, aocl_int64_t *info);
int dgsvj0_check(char *jobv, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *d__,
                 double *sva, aocl_int64_t *mv, double *v, aocl_int64_t *ldv, double *eps, double *sfmin,
                 double *tol, aocl_int64_t *nsweep, double *work, aocl_int64_t *lwork, aocl_int64_t *info);
int dgsvj1_check(char *jobv, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *n1, double *a, aocl_int64_t *lda,
                 double *d__, double *sva, aocl_int64_t *mv, double *v, aocl_int64_t *ldv, double *eps,
                 double *sfmin, double *tol, aocl_int64_t *nsweep, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int dgtcon_check(char *norm, aocl_int64_t *n, double *dl, double *d__, double *du, double *du2,
                 aocl_int_t *ipiv, double *anorm, double *rcond, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int dgtrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *dl, double *d__, double *du,
                 double *dlf, double *df, double *duf, double *du2, aocl_int_t *ipiv, double *b,
                 aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *ferr, double *berr, double *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int dgtsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *dl, double *d__, double *du, double *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int dgtsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *dl, double *d__,
                 double *du, double *dlf, double *df, double *duf, double *du2, aocl_int_t *ipiv,
                 double *b, aocl_int64_t *ldb, double *x, aocl_int64_t *ldx, double *rcond, double *ferr,
                 double *berr, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int dgttrf_check(aocl_int64_t *n, double *dl, double *d__, double *du, double *du2, aocl_int_t *ipiv,
                 aocl_int64_t *info);
int dgttrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *dl, double *d__, double *du,
                 double *du2, aocl_int_t *ipiv, double *b, aocl_int64_t *ldb, aocl_int64_t *info);
int dgtts2_check(aocl_int64_t *itrans, aocl_int64_t *n, aocl_int64_t *nrhs, double *dl, double *d__, double *du,
                 double *du2, aocl_int_t *ipiv, double *b, aocl_int64_t *ldb);
int dhgeqz_check(char *job, char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 double *h__, aocl_int64_t *ldh, double *t, aocl_int64_t *ldt, double *alphar, double *alphai,
                 double *beta, double *q, aocl_int64_t *ldq, double *z__, aocl_int64_t *ldz, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int dhsein_check(char *side, char *eigsrc, char *initv, logical *select, aocl_int64_t *n, double *h__,
                 aocl_int64_t *ldh, double *wr, double *wi, double *vl, aocl_int64_t *ldvl, double *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m, double *work, aocl_int64_t *ifaill,
                 aocl_int64_t *ifailr, aocl_int64_t *info);
int dhseqr_check(char *job, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *h__,
                 aocl_int64_t *ldh, double *wr, double *wi, double *z__, aocl_int64_t *ldz, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
logical disnan_check(double *din);
int dla_gbamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *alpha,
                    double *ab, aocl_int64_t *ldab, double *x, aocl_int64_t *incx, double *beta, double *y,
                    aocl_int64_t *incy);
double dla_gbrcond_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *ab,
                         aocl_int64_t *ldab, double *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, aocl_int64_t *cmode,
                         double *c__, aocl_int64_t *info, double *work, aocl_int64_t *iwork);
int dla_gbrfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n, aocl_int64_t *kl,
                              aocl_int64_t *ku, aocl_int64_t *nrhs, double *ab, aocl_int64_t *ldab, double *afb,
                              aocl_int64_t *ldafb, aocl_int_t *ipiv, logical *colequ, double *c__,
                              double *b, aocl_int64_t *ldb, double *y, aocl_int64_t *ldy, double *berr_out__,
                              aocl_int64_t *n_norms__, double *err_bnds_norm__, double *err_bnds_comp__,
                              double *res, double *ayb, double *dy, double *y_tail__, double *rcond,
                              aocl_int64_t *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
double dla_gbrpvgrw_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *ncols, double *ab,
                          aocl_int64_t *ldab, double *afb, aocl_int64_t *ldafb);
int dla_geamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, double *alpha, double *a, aocl_int64_t *lda,
                    double *x, aocl_int64_t *incx, double *beta, double *y, aocl_int64_t *incy);
double dla_gercond_check(char *trans, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *af,
                         aocl_int64_t *ldaf, aocl_int_t *ipiv, aocl_int64_t *cmode, double *c__, aocl_int64_t *info,
                         double *work, aocl_int64_t *iwork);
int dla_gerfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n,
                              aocl_int64_t *nrhs, double *a, aocl_int64_t *lda, double *af, aocl_int64_t *ldaf,
                              aocl_int_t *ipiv, logical *colequ, double *c__, double *b, aocl_int64_t *ldb,
                              double *y, aocl_int64_t *ldy, double *berr_out__, aocl_int64_t *n_norms__,
                              double *errs_n__, double *errs_c__, double *res, double *ayb,
                              double *dy, double *y_tail__, double *rcond, aocl_int64_t *ithresh,
                              double *rthresh, double *dz_ub__, logical *ignore_cwise__,
                              aocl_int64_t *info);
double dla_gerpvgrw_check(aocl_int64_t *n, aocl_int64_t *ncols, double *a, aocl_int64_t *lda, double *af,
                          aocl_int64_t *ldaf);
int dla_lin_berr_check(aocl_int64_t *n, aocl_int64_t *nz, aocl_int64_t *nrhs, double *res, double *ayb,
                       double *berr);
double dla_porcond_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *af, aocl_int64_t *ldaf,
                         aocl_int64_t *cmode, double *c__, aocl_int64_t *info, double *work, aocl_int64_t *iwork);
int dla_porfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              double *a, aocl_int64_t *lda, double *af, aocl_int64_t *ldaf, logical *colequ,
                              double *c__, double *b, aocl_int64_t *ldb, double *y, aocl_int64_t *ldy,
                              double *berr_out__, aocl_int64_t *n_norms__, double *err_bnds_norm__,
                              double *err_bnds_comp__, double *res, double *ayb, double *dy,
                              double *y_tail__, double *rcond, aocl_int64_t *ithresh, double *rthresh,
                              double *dz_ub__, logical *ignore_cwise__, aocl_int64_t *info);
double dla_porpvgrw_check(char *uplo, aocl_int64_t *ncols, double *a, aocl_int64_t *lda, double *af,
                          aocl_int64_t *ldaf, double *work);
int dla_syamv_check(aocl_int64_t *uplo, aocl_int64_t *n, double *alpha, double *a, aocl_int64_t *lda, double *x,
                    aocl_int64_t *incx, double *beta, double *y, aocl_int64_t *incy);
double dla_syrcond_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *af, aocl_int64_t *ldaf,
                         aocl_int_t *ipiv, aocl_int64_t *cmode, double *c__, aocl_int64_t *info, double *work,
                         aocl_int64_t *iwork);
int dla_syrfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              double *a, aocl_int64_t *lda, double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv,
                              logical *colequ, double *c__, double *b, aocl_int64_t *ldb, double *y,
                              aocl_int64_t *ldy, double *berr_out__, aocl_int64_t *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, double *res,
                              double *ayb, double *dy, double *y_tail__, double *rcond,
                              aocl_int64_t *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
double dla_syrpvgrw_check(char *uplo, aocl_int64_t *n, aocl_int64_t *info, double *a, aocl_int64_t *lda,
                          double *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *work);
int dla_wwaddw_check(aocl_int64_t *n, double *x, double *y, double *w);
int dlabad_check(double *small_, double *large);
int dlabrd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, double *a, aocl_int64_t *lda, double *d__,
                 double *e, double *tauq, double *taup, double *x, aocl_int64_t *ldx, double *y,
                 aocl_int64_t *ldy);
int dlacn2_check(aocl_int64_t *n, double *v, double *x, aocl_int64_t *isgn, double *est, aocl_int64_t *kase,
                 aocl_int64_t *isave);
int dlacon_check(aocl_int64_t *n, double *v, double *x, aocl_int64_t *isgn, double *est, aocl_int64_t *kase);
int dlacpy_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *b,
                 aocl_int64_t *ldb);
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
int ilaver_check(aocl_int64_t *vers_major__, aocl_int64_t *vers_minor__, aocl_int64_t *vers_patch__);
int ilazlc_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda);
int ilazlr_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda);
int izmax1_check(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx);
int sbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, aocl_int64_t *m,
                 aocl_int64_t *p, aocl_int64_t *q, float *theta, float *phi, float *u1, aocl_int64_t *ldu1,
                 float *u2, aocl_int64_t *ldu2, float *v1t, aocl_int64_t *ldv1t, float *v2t, aocl_int64_t *ldv2t,
                 float *b11d, float *b11e, float *b12d, float *b12e, float *b21d, float *b21e,
                 float *b22d, float *b22e, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sbdsdc_check(char *uplo, char *compq, aocl_int64_t *n, float *d__, float *e, float *u, aocl_int64_t *ldu,
                 float *vt, aocl_int64_t *ldvt, float *q, aocl_int64_t *iq, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sbdsqr_check(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru, aocl_int64_t *ncc, float *d__,
                 float *e, float *vt, aocl_int64_t *ldvt, float *u, aocl_int64_t *ldu, float *c__,
                 aocl_int64_t *ldc, float *work, aocl_int64_t *info);
float scsum1_check(aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx);
int sdisna_check(char *job, aocl_int64_t *m, aocl_int64_t *n, float *d__, float *sep, aocl_int64_t *info);
int sgbbrd_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ncc, aocl_int64_t *kl, aocl_int64_t *ku,
                 float *ab, aocl_int64_t *ldab, float *d__, float *e, float *q, aocl_int64_t *ldq, float *pt,
                 aocl_int64_t *ldpt, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int sgbcon_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, float *anorm, float *rcond, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sgbequ_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int sgbequb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                  float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int sgbrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, float *ab,
                 aocl_int64_t *ldab, float *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, float *b, aocl_int64_t *ldb,
                 float *x, aocl_int64_t *ldx, float *ferr, float *berr, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sgbrfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  float *ab, aocl_int64_t *ldab, float *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, float *r__,
                  float *c__, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond,
                  float *berr, aocl_int64_t *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, aocl_int64_t *nparams, float *params, float *work,
                  aocl_int64_t *iwork, aocl_int64_t *info);
int sgbsv_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, float *ab, aocl_int64_t *ldab,
                aocl_int_t *ipiv, float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int sgbsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                 float *ab, aocl_int64_t *ldab, float *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, char *equed,
                 float *r__, float *c__, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx,
                 float *rcond, float *ferr, float *berr, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sgbsvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  float *ab, aocl_int64_t *ldab, float *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, char *equed,
                  float *r__, float *c__, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx,
                  float *rcond, float *rpvgrw, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sgbtf2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int sgbtrf_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int sgbtrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, float *ab,
                 aocl_int64_t *ldab, aocl_int_t *ipiv, float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int sgebak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *scale,
                 aocl_int64_t *m, float *v, aocl_int64_t *ldv, aocl_int64_t *info);
int sgebal_check(char *job, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 float *scale, aocl_int64_t *info);
int sgebd2_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *d__, float *e, float *tauq,
                 float *taup, float *work, aocl_int64_t *info);
int sgebrd_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *d__, float *e, float *tauq,
                 float *taup, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgecon_check(char *norm, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *anorm, float *rcond,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sgeequ_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *r__, float *c__,
                 float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int sgeequb_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *r__, float *c__,
                  float *rowcnd, float *colcnd, float *amax, aocl_int64_t *info);
int sgees_check(char *jobvs, char *sort, L_fp select, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                aocl_int64_t *sdim, float *wr, float *wi, float *vs, aocl_int64_t *ldvs, float *work,
                aocl_int64_t *lwork, logical *bwork, aocl_int64_t *info);
int sgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, aocl_int64_t *sdim, float *wr, float *wi, float *vs, aocl_int64_t *ldvs,
                 float *rconde, float *rcondv, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, logical *bwork, aocl_int64_t *info);
int sgeev_check(char *jobvl, char *jobvr, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *wr, float *wi,
                float *vl, aocl_int64_t *ldvl, float *vr, aocl_int64_t *ldvr, float *work, aocl_int64_t *lwork,
                aocl_int64_t *info);
int sgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, float *wr, float *wi, float *vl, aocl_int64_t *ldvl, float *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *ilo, aocl_int64_t *ihi, float *scale, float *abnrm,
                 float *rconde, float *rcondv, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sgegs_check(char *jobvsl, char *jobvsr, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                aocl_int64_t *ldb, float *alphar, float *alphai, float *beta, float *vsl, aocl_int64_t *ldvsl,
                float *vsr, aocl_int64_t *ldvsr, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgegv_check(char *jobvl, char *jobvr, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                aocl_int64_t *ldb, float *alphar, float *alphai, float *beta, float *vl, aocl_int64_t *ldvl,
                float *vr, aocl_int64_t *ldvr, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgehd2_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *info);
int sgehrd_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgejsv_check(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, aocl_int64_t *m,
                 aocl_int64_t *n, float *a, aocl_int64_t *lda, float *sva, float *u, aocl_int64_t *ldu, float *v,
                 aocl_int64_t *ldv, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int sgelq2_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *info);
int sgelqf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int sgels_check(char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                float *b, aocl_int64_t *ldb, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgelsd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, float *s, float *rcond, aocl_int64_t *rank, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int sgelss_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, float *s, float *rcond, aocl_int64_t *rank, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sgelsx_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, float *rcond, aocl_int64_t *rank, float *work,
                 aocl_int64_t *info);
int sgelsy_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, float *rcond, aocl_int64_t *rank, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int sgemqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb,
                  float *v, aocl_int64_t *ldv, float *t, aocl_int64_t *ldt, float *c__, aocl_int64_t *ldc,
                  float *work, aocl_int64_t *info);
int sgeql2_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *info);
int sgeqlf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int sgeqp3_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *jpvt, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgeqpf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *jpvt, float *tau,
                 float *work, aocl_int64_t *info);
int sgeqr2_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *info);
int sgeqr2p_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                  aocl_int64_t *info);
int sgeqrf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int sgeqrfp_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int sgeqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, float *a, aocl_int64_t *lda, float *t,
                 aocl_int64_t *ldt, float *work, aocl_int64_t *info);
int sgeqrt2_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int sgeqrt3_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int sgerfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sgerfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                  float *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *r__, float *c__, float *b,
                  aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond, float *berr,
                  aocl_int64_t *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  aocl_int64_t *nparams, float *params, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sgerq2_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *info);
int sgerqf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int sgesc2_check(aocl_int64_t *n, float *a, aocl_int64_t *lda, float *rhs, aocl_int_t *ipiv, aocl_int64_t *jpiv,
                 float *scale);
int sgesdd_check(char *jobz, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *s, float *u,
                 aocl_int64_t *ldu, float *vt, aocl_int64_t *ldvt, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int sgesv_check(aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int sgesvd_check(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *s,
                 float *u, aocl_int64_t *ldu, float *vt, aocl_int64_t *ldvt, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sgesvj_check(char *joba, char *jobu, char *jobv, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                 float *sva, aocl_int64_t *mv, float *v, aocl_int64_t *ldv, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sgesvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                 float *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, float *r__, float *c__,
                 float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond, float *ferr,
                 float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sgesvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                  float *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, float *r__, float *c__,
                  float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond, float *rpvgrw,
                  float *berr, aocl_int64_t *n_err_bnds__, float *err_bnds_norm__,
                  float *err_bnds_comp__, aocl_int64_t *nparams, float *params, float *work,
                  aocl_int64_t *iwork, aocl_int64_t *info);
int sgetc2_check(aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *jpiv, aocl_int64_t *info);
int sgetf2_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int sgetrf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int sgetrfnp_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int sgetrfnpi_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nfact, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int sgetri_check(aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sgetrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int sggbak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *lscale,
                 float *rscale, aocl_int64_t *m, float *v, aocl_int64_t *ldv, aocl_int64_t *info);
int sggbal_check(char *job, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb,
                 aocl_int64_t *ilo, aocl_int64_t *ihi, float *lscale, float *rscale, float *work,
                 aocl_int64_t *info);
int sgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, aocl_int64_t *n, float *a,
                aocl_int64_t *lda, float *b, aocl_int64_t *ldb, aocl_int64_t *sdim, float *alphar, float *alphai,
                float *beta, float *vsl, aocl_int64_t *ldvsl, float *vsr, aocl_int64_t *ldvsr, float *work,
                aocl_int64_t *lwork, logical *bwork, aocl_int64_t *info);
int sggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, aocl_int64_t *n,
                 float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb, aocl_int64_t *sdim, float *alphar,
                 float *alphai, float *beta, float *vsl, aocl_int64_t *ldvsl, float *vsr, aocl_int64_t *ldvsr,
                 float *rconde, float *rcondv, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, logical *bwork, aocl_int64_t *info);
int sggev_check(char *jobvl, char *jobvr, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                aocl_int64_t *ldb, float *alphar, float *alphai, float *beta, float *vl, aocl_int64_t *ldvl,
                float *vr, aocl_int64_t *ldvr, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *alphar, float *alphai, float *beta,
                 float *vl, aocl_int64_t *ldvl, float *vr, aocl_int64_t *ldvr, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 float *lscale, float *rscale, float *abnrm, float *bbnrm, float *rconde,
                 float *rcondv, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, logical *bwork,
                 aocl_int64_t *info);
int sggglm_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb,
                 float *d__, float *x, float *y, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgghrd_check(char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *q, aocl_int64_t *ldq, float *z__,
                 aocl_int64_t *ldz, aocl_int64_t *info);
int sgglse_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb,
                 float *c__, float *d__, float *x, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sggqrf_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, float *a, aocl_int64_t *lda, float *taua, float *b,
                 aocl_int64_t *ldb, float *taub, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sggrqf_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *taua, float *b,
                 aocl_int64_t *ldb, float *taub, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sggsvd_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, aocl_int64_t *k,
                 aocl_int64_t *l, float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *alpha,
                 float *beta, float *u, aocl_int64_t *ldu, float *v, aocl_int64_t *ldv, float *q,
                 aocl_int64_t *ldq, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sggsvp_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *tola, float *tolb, aocl_int64_t *k,
                 aocl_int64_t *l, float *u, aocl_int64_t *ldu, float *v, aocl_int64_t *ldv, float *q, aocl_int64_t *ldq,
                 aocl_int64_t *iwork, float *tau, float *work, aocl_int64_t *info);
int sgsvj0_check(char *jobv, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *d__, float *sva,
                 aocl_int64_t *mv, float *v, aocl_int64_t *ldv, float *eps, float *sfmin, float *tol,
                 aocl_int64_t *nsweep, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sgsvj1_check(char *jobv, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *n1, float *a, aocl_int64_t *lda,
                 float *d__, float *sva, aocl_int64_t *mv, float *v, aocl_int64_t *ldv, float *eps,
                 float *sfmin, float *tol, aocl_int64_t *nsweep, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sgtcon_check(char *norm, aocl_int64_t *n, float *dl, float *d__, float *du, float *du2,
                 aocl_int_t *ipiv, float *anorm, float *rcond, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sgtrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *dl, float *d__, float *du,
                 float *dlf, float *df, float *duf, float *du2, aocl_int_t *ipiv, float *b,
                 aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *ferr, float *berr, float *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int sgtsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, float *dl, float *d__, float *du, float *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int sgtsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *dl, float *d__,
                 float *du, float *dlf, float *df, float *duf, float *du2, aocl_int_t *ipiv, float *b,
                 aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond, float *ferr, float *berr,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sgttrf_check(aocl_int64_t *n, float *dl, float *d__, float *du, float *du2, aocl_int_t *ipiv,
                 aocl_int64_t *info);
int sgttrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, float *dl, float *d__, float *du,
                 float *du2, aocl_int_t *ipiv, float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int sgtts2_check(aocl_int64_t *itrans, aocl_int64_t *n, aocl_int64_t *nrhs, float *dl, float *d__, float *du,
                 float *du2, aocl_int_t *ipiv, float *b, aocl_int64_t *ldb);
int shgeqz_check(char *job, char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 float *h__, aocl_int64_t *ldh, float *t, aocl_int64_t *ldt, float *alphar, float *alphai,
                 float *beta, float *q, aocl_int64_t *ldq, float *z__, aocl_int64_t *ldz, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int shsein_check(char *side, char *eigsrc, char *initv, logical *select, aocl_int64_t *n, float *h__,
                 aocl_int64_t *ldh, float *wr, float *wi, float *vl, aocl_int64_t *ldvl, float *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m, float *work, aocl_int64_t *ifaill,
                 aocl_int64_t *ifailr, aocl_int64_t *info);
int shseqr_check(char *job, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *h__,
                 aocl_int64_t *ldh, float *wr, float *wi, float *z__, aocl_int64_t *ldz, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
logical sisnan_check(float *sin__);
int sla_gbamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *alpha,
                    float *ab, aocl_int64_t *ldab, float *x, aocl_int64_t *incx, float *beta, float *y,
                    aocl_int64_t *incy);
float sla_gbrcond_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, float *ab, aocl_int64_t *ldab,
                        float *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, aocl_int64_t *cmode, float *c__,
                        aocl_int64_t *info, float *work, aocl_int64_t *iwork);
int sla_gbrfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n, aocl_int64_t *kl,
                              aocl_int64_t *ku, aocl_int64_t *nrhs, float *ab, aocl_int64_t *ldab, float *afb,
                              aocl_int64_t *ldafb, aocl_int_t *ipiv, logical *colequ, float *c__, float *b,
                              aocl_int64_t *ldb, float *y, aocl_int64_t *ldy, float *berr_out__,
                              aocl_int64_t *n_norms__, float *err_bnds_norm__, float *err_bnds_comp__,
                              float *res, float *ayb, float *dy, float *y_tail__, float *rcond,
                              aocl_int64_t *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
float sla_gbrpvgrw_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *ncols, float *ab,
                         aocl_int64_t *ldab, float *afb, aocl_int64_t *ldafb);
int sla_geamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, float *alpha, float *a, aocl_int64_t *lda,
                    float *x, aocl_int64_t *incx, float *beta, float *y, aocl_int64_t *incy);
float sla_gercond_check(char *trans, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *af, aocl_int64_t *ldaf,
                        aocl_int_t *ipiv, aocl_int64_t *cmode, float *c__, aocl_int64_t *info, float *work,
                        aocl_int64_t *iwork);
int sla_gerfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n,
                              aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *af, aocl_int64_t *ldaf,
                              aocl_int_t *ipiv, logical *colequ, float *c__, float *b, aocl_int64_t *ldb,
                              float *y, aocl_int64_t *ldy, float *berr_out__, aocl_int64_t *n_norms__,
                              float *errs_n__, float *errs_c__, float *res, float *ayb, float *dy,
                              float *y_tail__, float *rcond, aocl_int64_t *ithresh, float *rthresh,
                              float *dz_ub__, logical *ignore_cwise__, aocl_int64_t *info);
float sla_gerpvgrw_check(aocl_int64_t *n, aocl_int64_t *ncols, float *a, aocl_int64_t *lda, float *af,
                         aocl_int64_t *ldaf);
int sla_lin_berr_check(aocl_int64_t *n, aocl_int64_t *nz, aocl_int64_t *nrhs, float *res, float *ayb, float *berr);
float sla_porcond_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *af, aocl_int64_t *ldaf,
                        aocl_int64_t *cmode, float *c__, aocl_int64_t *info, float *work, aocl_int64_t *iwork);
int sla_porfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a,
                              aocl_int64_t *lda, float *af, aocl_int64_t *ldaf, logical *colequ, float *c__,
                              float *b, aocl_int64_t *ldb, float *y, aocl_int64_t *ldy, float *berr_out__,
                              aocl_int64_t *n_norms__, float *err_bnds_norm__, float *err_bnds_comp__,
                              float *res, float *ayb, float *dy, float *y_tail__, float *rcond,
                              aocl_int64_t *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
float sla_porpvgrw_check(char *uplo, aocl_int64_t *ncols, float *a, aocl_int64_t *lda, float *af,
                         aocl_int64_t *ldaf, float *work);
int sla_syamv_check(aocl_int64_t *uplo, aocl_int64_t *n, float *alpha, float *a, aocl_int64_t *lda, float *x,
                    aocl_int64_t *incx, float *beta, float *y, aocl_int64_t *incy);
float sla_syrcond_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *af, aocl_int64_t *ldaf,
                        aocl_int_t *ipiv, aocl_int64_t *cmode, float *c__, aocl_int64_t *info, float *work,
                        aocl_int64_t *iwork);
int sla_syrfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a,
                              aocl_int64_t *lda, float *af, aocl_int64_t *ldaf, aocl_int_t *ipiv,
                              logical *colequ, float *c__, float *b, aocl_int64_t *ldb, float *y,
                              aocl_int64_t *ldy, float *berr_out__, aocl_int64_t *n_norms__,
                              float *err_bnds_norm__, float *err_bnds_comp__, float *res,
                              float *ayb, float *dy, float *y_tail__, float *rcond,
                              aocl_int64_t *ithresh, float *rthresh, float *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
float sla_syrpvgrw_check(char *uplo, aocl_int64_t *n, aocl_int64_t *info, float *a, aocl_int64_t *lda, float *af,
                         aocl_int64_t *ldaf, aocl_int_t *ipiv, float *work);
int sla_wwaddw_check(aocl_int64_t *n, float *x, float *y, float *w);
int slabad_check(float *small_, float *large);
int slabrd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, float *a, aocl_int64_t *lda, float *d__, float *e,
                 float *tauq, float *taup, float *x, aocl_int64_t *ldx, float *y, aocl_int64_t *ldy);
int slacn2_check(aocl_int64_t *n, float *v, float *x, aocl_int64_t *isgn, float *est, aocl_int64_t *kase,
                 aocl_int64_t *isave);
int slacon_check(aocl_int64_t *n, float *v, float *x, aocl_int64_t *isgn, float *est, aocl_int64_t *kase);
int slacpy_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb);
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
int slartgs_check(float *x, float *y, float *sigma, float *cs, float *sn);
int slartv_check(aocl_int64_t *n, float *x, aocl_int64_t *incx, float *y, aocl_int64_t *incy, float *c__, float *s,
                 aocl_int64_t *incc);
int slaruv_check(aocl_int64_t *iseed, aocl_int64_t *n, float *x);
int slarz_check(char *side, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, float *v, aocl_int64_t *incv, float *tau,
                float *c__, aocl_int64_t *ldc, float *work);
int slarzb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, float *v, aocl_int64_t *ldv, float *t, aocl_int64_t *ldt, float *c__,
                 aocl_int64_t *ldc, float *work, aocl_int64_t *ldwork);
int slarzt_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, float *v, aocl_int64_t *ldv,
                 float *tau, float *t, aocl_int64_t *ldt);
int slas2_check(float *f, float *g, float *h__, float *ssmin, float *ssmax);
int slascl_check(char *type__, aocl_int64_t *kl, aocl_int64_t *ku, float *cfrom, float *cto, aocl_int64_t *m,
                 aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int slascl2_check(aocl_int64_t *m, aocl_int64_t *n, float *d__, float *x, aocl_int64_t *ldx);
int slasd0_check(aocl_int64_t *n, aocl_int64_t *sqre, float *d__, float *e, float *u, aocl_int64_t *ldu, float *vt,
                 aocl_int64_t *ldvt, aocl_int64_t *smlsiz, aocl_int64_t *iwork, float *work, aocl_int64_t *info);
int slasd1_check(aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, float *d__, float *alpha, float *beta,
                 float *u, aocl_int64_t *ldu, float *vt, aocl_int64_t *ldvt, aocl_int64_t *idxq, aocl_int64_t *iwork,
                 float *work, aocl_int64_t *info);
int slasd2_check(aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *k, float *d__, float *z__,
                 float *alpha, float *beta, float *u, aocl_int64_t *ldu, float *vt, aocl_int64_t *ldvt,
                 float *dsigma, float *u2, aocl_int64_t *ldu2, float *vt2, aocl_int64_t *ldvt2, aocl_int64_t *idxp,
                 aocl_int64_t *idx, aocl_int64_t *idxc, aocl_int64_t *idxq, aocl_int64_t *coltyp, aocl_int64_t *info);
int slasd3_check(aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *k, float *d__, float *q,
                 aocl_int64_t *ldq, float *dsigma, float *u, aocl_int64_t *ldu, float *u2, aocl_int64_t *ldu2,
                 float *vt, aocl_int64_t *ldvt, float *vt2, aocl_int64_t *ldvt2, aocl_int64_t *idxc, aocl_int64_t *ctot,
                 float *z__, aocl_int64_t *info);
int slasd4_check(aocl_int64_t *n, aocl_int64_t *i__, float *d__, float *z__, float *delta, float *rho,
                 float *sigma, float *work, aocl_int64_t *info);
int slasd5_check(aocl_int64_t *i__, float *d__, float *z__, float *delta, float *rho, float *dsigma,
                 float *work);
int slasd6_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, float *d__, float *vf,
                 float *vl, float *alpha, float *beta, aocl_int64_t *idxq, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol, float *givnum, aocl_int64_t *ldgnum,
                 float *poles, float *difl, float *difr, float *z__, aocl_int64_t *k, float *c__,
                 float *s, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int slasd7_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *k, float *d__,
                 float *z__, float *zw, float *vf, float *vfw, float *vl, float *vlw, float *alpha,
                 float *beta, float *dsigma, aocl_int64_t *idx, aocl_int64_t *idxp, aocl_int64_t *idxq,
                 aocl_int64_t *perm, aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol, float *givnum,
                 aocl_int64_t *ldgnum, float *c__, float *s, aocl_int64_t *info);
int slasd8_check(aocl_int64_t *icompq, aocl_int64_t *k, float *d__, float *z__, float *vf, float *vl,
                 float *difl, float *difr, aocl_int64_t *lddifr, float *dsigma, float *work,
                 aocl_int64_t *info);
int slasda_check(aocl_int64_t *icompq, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *sqre, float *d__, float *e,
                 float *u, aocl_int64_t *ldu, float *vt, aocl_int64_t *k, float *difl, float *difr,
                 float *z__, float *poles, aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol,
                 aocl_int64_t *perm, float *givnum, float *c__, float *s, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int slasdq_check(char *uplo, aocl_int64_t *sqre, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru, aocl_int64_t *ncc,
                 float *d__, float *e, float *vt, aocl_int64_t *ldvt, float *u, aocl_int64_t *ldu, float *c__,
                 aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int slasdt_check(aocl_int64_t *n, aocl_int64_t *lvl, aocl_int64_t *nd, aocl_int64_t *inode, aocl_int64_t *ndiml,
                 aocl_int64_t *ndimr, aocl_int64_t *msub);
int slaset_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, float *alpha, float *beta, float *a,
                 aocl_int64_t *lda);
int slasq1_check(aocl_int64_t *n, float *d__, float *e, float *work, aocl_int64_t *info);
int slasq2_check(aocl_int64_t *n, float *z__, aocl_int64_t *info);
int slasq3_check(aocl_int64_t *i0, aocl_int64_t *n0, float *z__, aocl_int64_t *pp, float *dmin__, float *sigma,
                 float *desig, float *qmax, aocl_int64_t *nfail, aocl_int64_t *iter, aocl_int64_t *ndiv,
                 logical *ieee, aocl_int64_t *ttype, float *dmin1, float *dmin2, float *dn, float *dn1,
                 float *dn2, float *g, float *tau);
int slasq4_check(aocl_int64_t *i0, aocl_int64_t *n0, float *z__, aocl_int64_t *pp, aocl_int64_t *n0in, float *dmin__,
                 float *dmin1, float *dmin2, float *dn, float *dn1, float *dn2, float *tau,
                 aocl_int64_t *ttype, float *g);
int slasq5_check(aocl_int64_t *i0, aocl_int64_t *n0, float *z__, aocl_int64_t *pp, float *tau, float *sigma,
                 float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float *dnm2,
                 logical *ieee, float *eps);
int slasq6_check(aocl_int64_t *i0, aocl_int64_t *n0, float *z__, aocl_int64_t *pp, float *dmin__, float *dmin1,
                 float *dmin2, float *dn, float *dnm1, float *dnm2);
int slasr_check(char *side, char *pivot, char *direct, aocl_int64_t *m, aocl_int64_t *n, float *c__, float *s,
                float *a, aocl_int64_t *lda);
int slasrt_check(char *id, aocl_int64_t *n, float *d__, aocl_int64_t *info);
int slassq_check(aocl_int64_t *n, float *x, aocl_int64_t *incx, float *scale, float *sumsq);
int slasv2_check(float *f, float *g, float *h__, float *ssmin, float *ssmax, float *snr, float *csr,
                 float *snl, float *csl);
int slaswp_check(aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *k1, aocl_int64_t *k2, aocl_int_t *ipiv,
                 aocl_int64_t *incx);
int slasy2_check(logical *ltranl, logical *ltranr, aocl_int64_t *isgn, aocl_int64_t *n1, aocl_int64_t *n2,
                 float *tl, aocl_int64_t *ldtl, float *tr, aocl_int64_t *ldtr, float *b, aocl_int64_t *ldb,
                 float *scale, float *x, aocl_int64_t *ldx, float *xnorm, aocl_int64_t *info);
int slasyf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, float *a, aocl_int64_t *lda,
                 aocl_int_t *ipiv, float *w, aocl_int64_t *ldw, aocl_int64_t *info);
int slasyf_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, float *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, float *w, aocl_int64_t *ldw, aocl_int64_t *info);
int slatbs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, aocl_int64_t *kd,
                 float *ab, aocl_int64_t *ldab, float *x, float *scale, float *cnorm, aocl_int64_t *info);
int slatdf_check(aocl_int64_t *ijob, aocl_int64_t *n, float *z__, aocl_int64_t *ldz, float *rhs, float *rdsum,
                 float *rdscal, aocl_int_t *ipiv, aocl_int64_t *jpiv);
int slatps_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, float *ap, float *x,
                 float *scale, float *cnorm, aocl_int64_t *info);
int slatrd_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, float *a, aocl_int64_t *lda, float *e, float *tau,
                 float *w, aocl_int64_t *ldw);
int slatrs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, float *x, float *scale, float *cnorm, aocl_int64_t *info);
int slatrz_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, float *a, aocl_int64_t *lda, float *tau,
                 float *work);
int slatzm_check(char *side, aocl_int64_t *m, aocl_int64_t *n, float *v, aocl_int64_t *incv, float *tau, float *c1,
                 float *c2, aocl_int64_t *ldc, float *work);
int slauu2_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int slauum_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int sopgtr_check(char *uplo, aocl_int64_t *n, float *ap, float *tau, float *q, aocl_int64_t *ldq, float *work,
                 aocl_int64_t *info);
int sopmtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, float *ap, float *tau,
                 float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int sorbdb_check(char *trans, char *signs, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, float *x11,
                 aocl_int64_t *ldx11, float *x12, aocl_int64_t *ldx12, float *x21, aocl_int64_t *ldx21, float *x22,
                 aocl_int64_t *ldx22, float *theta, float *phi, float *taup1, float *taup2, float *tauq1,
                 float *tauq2, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorbdb1_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, float *x11, aocl_int64_t *ldx11, float *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorbdb2_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, float *x11, aocl_int64_t *ldx11, float *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorbdb3_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, float *x11, aocl_int64_t *ldx11, float *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorbdb4_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, float *x11, aocl_int64_t *ldx11, float *x21,
                  aocl_int64_t *ldx21, float *theta, float *phi, float *taup1, float *taup2,
                  float *tauq1, float *phantom, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorbdb5_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, float *x1, aocl_int64_t *incx1, float *x2,
                  aocl_int64_t *incx2, float *q1, aocl_int64_t *ldq1, float *q2, aocl_int64_t *ldq2, float *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int sorbdb6_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, float *x1, aocl_int64_t *incx1, float *x2,
                  aocl_int64_t *incx2, float *q1, aocl_int64_t *ldq1, float *q2, aocl_int64_t *ldq2, float *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int sorcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, float *x11, aocl_int64_t *ldx11, float *x12,
                 aocl_int64_t *ldx12, float *x21, aocl_int64_t *ldx21, float *x22, aocl_int64_t *ldx22,
                 float *theta, float *u1, aocl_int64_t *ldu1, float *u2, aocl_int64_t *ldu2, float *v1t,
                 aocl_int64_t *ldv1t, float *v2t, aocl_int64_t *ldv2t, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int sorcsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q,
                     float *x11, aocl_int64_t *ldx11, float *x21, aocl_int64_t *ldx21, float *theta,
                     float *u1, aocl_int64_t *ldu1, float *u2, aocl_int64_t *ldu2, float *v1t, aocl_int64_t *ldv1t,
                     float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int sorg2l_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *info);
int sorg2r_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *info);
int sorgbr_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorghr_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorgl2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *info);
int sorglq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorgql_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorgqr_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorgr2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *info);
int sorgrq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a, aocl_int64_t *lda, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int sorgtr_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int sorm2l_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int sorm2r_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int sormbr_check(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sormhr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 float *a, aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int sorml2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int sormlq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sormql_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sormqr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sormr2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int sormr3_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *info);
int sormrq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sormrz_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int sormtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, float *tau, float *c__, aocl_int64_t *ldc, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int spbcon_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab, float *anorm,
                 float *rcond, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int spbequ_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab, float *s,
                 float *scond, float *amax, aocl_int64_t *info);
int spbrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, float *ab, aocl_int64_t *ldab,
                 float *afb, aocl_int64_t *ldafb, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int spbstf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int spbsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, float *ab, aocl_int64_t *ldab,
                float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int spbsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, float *ab,
                 aocl_int64_t *ldab, float *afb, aocl_int64_t *ldafb, char *equed, float *s, float *b,
                 aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond, float *ferr, float *berr,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int spbtf2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int spbtrf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int spbtrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, float *ab, aocl_int64_t *ldab,
                 float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int spftrf_check(char *transr, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *info);
int spftri_check(char *transr, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *info);
int spftrs_check(char *transr, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, float *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int spocon_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *anorm, float *rcond,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int spoequ_check(aocl_int64_t *n, float *a, aocl_int64_t *lda, float *s, float *scond, float *amax,
                 aocl_int64_t *info);
int spoequb_check(aocl_int64_t *n, float *a, aocl_int64_t *lda, float *s, float *scond, float *amax,
                  aocl_int64_t *info);
int sporfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *af,
                 aocl_int64_t *ldaf, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *ferr,
                 float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sporfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                  float *af, aocl_int64_t *ldaf, float *s, float *b, aocl_int64_t *ldb, float *x,
                  aocl_int64_t *ldx, float *rcond, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sposv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int sposvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                 float *af, aocl_int64_t *ldaf, char *equed, float *s, float *b, aocl_int64_t *ldb, float *x,
                 aocl_int64_t *ldx, float *rcond, float *ferr, float *berr, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sposvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                  float *af, aocl_int64_t *ldaf, char *equed, float *s, float *b, aocl_int64_t *ldb, float *x,
                  aocl_int64_t *ldx, float *rcond, float *rpvgrw, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int spotf2_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int spotrf_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int spotri_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int spotrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int sppcon_check(char *uplo, aocl_int64_t *n, float *ap, float *anorm, float *rcond, float *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int sppequ_check(char *uplo, aocl_int64_t *n, float *ap, float *s, float *scond, float *amax,
                 aocl_int64_t *info);
int spprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, float *afp, float *b,
                 aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *ferr, float *berr, float *work,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int sppsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, float *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int sppsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, float *afp,
                 char *equed, float *s, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx,
                 float *rcond, float *ferr, float *berr, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int spptrf_check(char *uplo, aocl_int64_t *n, float *ap, aocl_int64_t *info);
int spptri_check(char *uplo, aocl_int64_t *n, float *ap, aocl_int64_t *info);
int spptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, float *b, aocl_int64_t *ldb,
                 aocl_int64_t *info);
int spstf2_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 float *tol, float *work, aocl_int64_t *info);
int spstrf_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 float *tol, float *work, aocl_int64_t *info);
int sptcon_check(aocl_int64_t *n, float *d__, float *e, float *anorm, float *rcond, float *work,
                 aocl_int64_t *info);
int spteqr_check(char *compz, aocl_int64_t *n, float *d__, float *e, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *info);
int sptrfs_check(aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, float *e, float *df, float *ef, float *b,
                 aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *ferr, float *berr, float *work,
                 aocl_int64_t *info);
int sptsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, float *e, float *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int sptsvx_check(char *fact, aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, float *e, float *df, float *ef,
                 float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond, float *ferr,
                 float *berr, float *work, aocl_int64_t *info);
int spttrf_check(aocl_int64_t *n, float *d__, float *e, aocl_int64_t *info);
int spttrs_check(aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, float *e, float *b, aocl_int64_t *ldb,
                 aocl_int64_t *info);
int sptts2_check(aocl_int64_t *n, aocl_int64_t *nrhs, float *d__, float *e, float *b, aocl_int64_t *ldb);
int srscl_check(aocl_int64_t *n, float *sa, float *sx, aocl_int64_t *incx);
int ssbev_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab, float *w,
                float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *info);
int ssbevd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab,
                 float *w, float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int ssbevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab,
                 aocl_int64_t *ldab, float *q, aocl_int64_t *ldq, float *vl, float *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int ssbgst_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, float *ab,
                 aocl_int64_t *ldab, float *bb, aocl_int64_t *ldbb, float *x, aocl_int64_t *ldx, float *work,
                 aocl_int64_t *info);
int ssbgv_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, float *ab,
                aocl_int64_t *ldab, float *bb, aocl_int64_t *ldbb, float *w, float *z__, aocl_int64_t *ldz,
                float *work, aocl_int64_t *info);
int ssbgvd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, float *ab,
                 aocl_int64_t *ldab, float *bb, aocl_int64_t *ldbb, float *w, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int ssbgvx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb,
                 float *ab, aocl_int64_t *ldab, float *bb, aocl_int64_t *ldbb, float *q, aocl_int64_t *ldq,
                 float *vl, float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m,
                 float *w, float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *iwork, aocl_int64_t *ifail,
                 aocl_int64_t *info);
int ssbtrd_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, float *ab, aocl_int64_t *ldab,
                 float *d__, float *e, float *q, aocl_int64_t *ldq, float *work, aocl_int64_t *info);
int ssfrk_check(char *transr, char *uplo, char *trans, aocl_int64_t *n, aocl_int64_t *k, float *alpha,
                float *a, aocl_int64_t *lda, float *beta, float *c__);
int sspcon_check(char *uplo, aocl_int64_t *n, float *ap, aocl_int_t *ipiv, float *anorm, float *rcond,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sspev_check(char *jobz, char *uplo, aocl_int64_t *n, float *ap, float *w, float *z__, aocl_int64_t *ldz,
                float *work, aocl_int64_t *info);
int sspevd_check(char *jobz, char *uplo, aocl_int64_t *n, float *ap, float *w, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int sspevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, float *ap, float *vl, float *vu,
                 aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, float *z__,
                 aocl_int64_t *ldz, float *work, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int sspgst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, float *ap, float *bp, aocl_int64_t *info);
int sspgv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, float *ap, float *bp, float *w,
                float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *info);
int sspgvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, float *ap, float *bp, float *w,
                 float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int sspgvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, float *ap,
                 float *bp, float *vl, float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol,
                 aocl_int64_t *m, float *w, float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int ssprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, float *afp, aocl_int_t *ipiv,
                 float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *ferr, float *berr,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int sspsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, aocl_int_t *ipiv, float *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int sspsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, float *afp,
                 aocl_int_t *ipiv, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond,
                 float *ferr, float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int ssptrd_check(char *uplo, aocl_int64_t *n, float *ap, float *d__, float *e, float *tau,
                 aocl_int64_t *info);
int ssptrf_check(char *uplo, aocl_int64_t *n, float *ap, aocl_int_t *ipiv, aocl_int64_t *info);
int ssptri_check(char *uplo, aocl_int64_t *n, float *ap, aocl_int_t *ipiv, float *work, aocl_int64_t *info);
int ssptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap, aocl_int_t *ipiv, float *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int sstebz_check(char *range, char *order, aocl_int64_t *n, float *vl, float *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, float *abstol, float *d__, float *e, aocl_int64_t *m, aocl_int64_t *nsplit,
                 float *w, aocl_int64_t *iblock, aocl_int64_t *isplit, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int sstedc_check(char *compz, aocl_int64_t *n, float *d__, float *e, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int sstegr_check(char *jobz, char *range, aocl_int64_t *n, float *d__, float *e, float *vl, float *vu,
                 aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, float *z__,
                 aocl_int64_t *ldz, aocl_int64_t *isuppz, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int sstein_check(aocl_int64_t *n, float *d__, float *e, aocl_int64_t *m, float *w, aocl_int64_t *iblock,
                 aocl_int64_t *isplit, float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int sstemr_check(char *jobz, char *range, aocl_int64_t *n, float *d__, float *e, float *vl, float *vu,
                 aocl_int64_t *il, aocl_int64_t *iu, aocl_int64_t *m, float *w, float *z__, aocl_int64_t *ldz,
                 aocl_int64_t *nzc, aocl_int64_t *isuppz, logical *tryrac, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int ssteqr_check(char *compz, aocl_int64_t *n, float *d__, float *e, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *info);
int ssterf_check(aocl_int64_t *n, float *d__, float *e, aocl_int64_t *info);
int sstev_check(char *jobz, aocl_int64_t *n, float *d__, float *e, float *z__, aocl_int64_t *ldz, float *work,
                aocl_int64_t *info);
int sstevd_check(char *jobz, aocl_int64_t *n, float *d__, float *e, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int sstevr_check(char *jobz, char *range, aocl_int64_t *n, float *d__, float *e, float *vl, float *vu,
                 aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, float *z__,
                 aocl_int64_t *ldz, aocl_int64_t *isuppz, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int sstevx_check(char *jobz, char *range, aocl_int64_t *n, float *d__, float *e, float *vl, float *vu,
                 aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, float *z__,
                 aocl_int64_t *ldz, float *work, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int ssycon_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *anorm,
                 float *rcond, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int ssycon_rook_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *anorm,
                      float *rcond, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int ssyconv_check(char *uplo, char *way, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  float *work, aocl_int64_t *info);
int ssyequb_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *s, float *scond,
                  float *amax, float *work, aocl_int64_t *info);
int ssyev_check(char *jobz, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *w, float *work,
                aocl_int64_t *lwork, aocl_int64_t *info);
int ssyevd_check(char *jobz, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *w, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int ssyevr_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *vl,
                 float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w,
                 float *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int ssyevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *vl,
                 float *vu, aocl_int64_t *il, aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w,
                 float *z__, aocl_int64_t *ldz, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int ssygs2_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int ssygst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int ssygv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                float *b, aocl_int64_t *ldb, float *w, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int ssygvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                 float *b, aocl_int64_t *ldb, float *w, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int ssygvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *vl, float *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, float *abstol, aocl_int64_t *m, float *w, float *z__, aocl_int64_t *ldz,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int ssyrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, float *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int ssyrfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                  float *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *s, float *b, aocl_int64_t *ldb,
                  float *x, aocl_int64_t *ldx, float *rcond, float *berr, aocl_int64_t *n_err_bnds__,
                  float *err_bnds_norm__, float *err_bnds_comp__, aocl_int64_t *nparams, float *params,
                  float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int ssysv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                float *b, aocl_int64_t *ldb, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int ssysv_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                     float *b, aocl_int64_t *ldb, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int ssysvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                 float *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, float *b, aocl_int64_t *ldb, float *x,
                 aocl_int64_t *ldx, float *rcond, float *ferr, float *berr, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int ssysvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda,
                  float *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, float *s, float *b,
                  aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *rcond, float *rpvgrw, float *berr,
                  aocl_int64_t *n_err_bnds__, float *err_bnds_norm__, float *err_bnds_comp__,
                  aocl_int64_t *nparams, float *params, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int ssyswapr_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *i1, aocl_int64_t *i2);
int ssytd2_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *d__, float *e, float *tau,
                 aocl_int64_t *info);
int ssytf2_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int ssytf2_rook_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int ssytrd_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *d__, float *e, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int ssytrf_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int ssytrf_rook_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *work,
                      aocl_int64_t *lwork, aocl_int64_t *info);
int ssytri_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *work,
                 aocl_int64_t *info);
int ssytri2_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int ssytri2x_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *work,
                   aocl_int64_t *nb, aocl_int64_t *info);
int ssytri_rook_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, float *work,
                      aocl_int64_t *info);
int ssytrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ssytrs2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  float *b, aocl_int64_t *ldb, float *work, aocl_int64_t *info);
int ssytrs_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, float *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int stbcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *kd, float *ab,
                 aocl_int64_t *ldab, float *rcond, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int stbrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 float *ab, aocl_int64_t *ldab, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx,
                 float *ferr, float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int stbtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 float *ab, aocl_int64_t *ldab, float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int stfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, aocl_int64_t *m,
                aocl_int64_t *n, float *alpha, float *a, float *b, aocl_int64_t *ldb);
int stftri_check(char *transr, char *uplo, char *diag, aocl_int64_t *n, float *a, aocl_int64_t *info);
int stfttp_check(char *transr, char *uplo, aocl_int64_t *n, float *arf, float *ap, aocl_int64_t *info);
int stfttr_check(char *transr, char *uplo, aocl_int64_t *n, float *arf, float *a, aocl_int64_t *lda,
                 aocl_int64_t *info);
int stgevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, float *s, aocl_int64_t *lds,
                 float *p, aocl_int64_t *ldp, float *vl, aocl_int64_t *ldvl, float *vr, aocl_int64_t *ldvr,
                 aocl_int64_t *mm, aocl_int64_t *m, float *work, aocl_int64_t *info);
int stgex2_check(logical *wantq, logical *wantz, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, float *q, aocl_int64_t *ldq, float *z__, aocl_int64_t *ldz, aocl_int64_t *j1,
                 aocl_int64_t *n1, aocl_int64_t *n2, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int stgexc_check(logical *wantq, logical *wantz, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, float *q, aocl_int64_t *ldq, float *z__, aocl_int64_t *ldz, aocl_int64_t *ifst,
                 aocl_int64_t *ilst, float *work, aocl_int64_t *lwork, aocl_int64_t *info);
int stgsen_check(aocl_int64_t *ijob, logical *wantq, logical *wantz, logical *select, aocl_int64_t *n,
                 float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *alphar, float *alphai,
                 float *beta, float *q, aocl_int64_t *ldq, float *z__, aocl_int64_t *ldz, aocl_int64_t *m,
                 float *pl, float *pr, float *dif, float *work, aocl_int64_t *lwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int stgsja_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, aocl_int64_t *k,
                 aocl_int64_t *l, float *a, aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *tola,
                 float *tolb, float *alpha, float *beta, float *u, aocl_int64_t *ldu, float *v,
                 aocl_int64_t *ldv, float *q, aocl_int64_t *ldq, float *work, aocl_int64_t *ncycle, aocl_int64_t *info);
int stgsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                 float *b, aocl_int64_t *ldb, float *vl, aocl_int64_t *ldvl, float *vr, aocl_int64_t *ldvr,
                 float *s, float *dif, aocl_int64_t *mm, aocl_int64_t *m, float *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int stgsy2_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                 float *b, aocl_int64_t *ldb, float *c__, aocl_int64_t *ldc, float *d__, aocl_int64_t *ldd,
                 float *e, aocl_int64_t *lde, float *f, aocl_int64_t *ldf, float *scale, float *rdsum,
                 float *rdscal, aocl_int64_t *iwork, aocl_int64_t *pq, aocl_int64_t *info);
int stgsyl_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                 float *b, aocl_int64_t *ldb, float *c__, aocl_int64_t *ldc, float *d__, aocl_int64_t *ldd,
                 float *e, aocl_int64_t *lde, float *f, aocl_int64_t *ldf, float *scale, float *dif,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int stpcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, float *ap, float *rcond,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int stpmqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                  aocl_int64_t *nb, float *v, aocl_int64_t *ldv, float *t, aocl_int64_t *ldt, float *a,
                  aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *work, aocl_int64_t *info);
int stpqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, aocl_int64_t *nb, float *a, aocl_int64_t *lda, float *b,
                 aocl_int64_t *ldb, float *t, aocl_int64_t *ldt, float *work, aocl_int64_t *info);
int stpqrt2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, float *a, aocl_int64_t *lda, float *b,
                  aocl_int64_t *ldb, float *t, aocl_int64_t *ldt, aocl_int64_t *info);
int stprfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, float *v, aocl_int64_t *ldv, float *t, aocl_int64_t *ldt, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *work, aocl_int64_t *ldwork);
int stprfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap,
                 float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *ferr, float *berr,
                 float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int stptri_check(char *uplo, char *diag, aocl_int64_t *n, float *ap, aocl_int64_t *info);
int stptrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, float *ap,
                 float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int stpttf_check(char *transr, char *uplo, aocl_int64_t *n, float *ap, float *arf, aocl_int64_t *info);
int stpttr_check(char *uplo, aocl_int64_t *n, float *ap, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int strcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, float *a, aocl_int64_t *lda,
                 float *rcond, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int strevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, float *t, aocl_int64_t *ldt,
                 float *vl, aocl_int64_t *ldvl, float *vr, aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m,
                 float *work, aocl_int64_t *info);
int strexc_check(char *compq, aocl_int64_t *n, float *t, aocl_int64_t *ldt, float *q, aocl_int64_t *ldq,
                 aocl_int64_t *ifst, aocl_int64_t *ilst, float *work, aocl_int64_t *info);
int strrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *x, aocl_int64_t *ldx, float *ferr,
                 float *berr, float *work, aocl_int64_t *iwork, aocl_int64_t *info);
int strsen_check(char *job, char *compq, logical *select, aocl_int64_t *n, float *t, aocl_int64_t *ldt,
                 float *q, aocl_int64_t *ldq, float *wr, float *wi, aocl_int64_t *m, float *s, float *sep,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int strsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, float *t, aocl_int64_t *ldt,
                 float *vl, aocl_int64_t *ldvl, float *vr, aocl_int64_t *ldvr, float *s, float *sep,
                 aocl_int64_t *mm, aocl_int64_t *m, float *work, aocl_int64_t *ldwork, aocl_int64_t *iwork,
                 aocl_int64_t *info);
int strsyl_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, float *c__, aocl_int64_t *ldc, float *scale,
                 aocl_int64_t *info);
int strsyl3_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, real *a,
                  aocl_int64_t *lda, real *b, aocl_int64_t *ldb, real *c__, aocl_int64_t *ldc, real *scale,
                  aocl_int64_t *iwork, aocl_int64_t *liwork, real *swork, aocl_int64_t *ldswork, aocl_int64_t *info);
int strti2_check(char *uplo, char *diag, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int strtri_check(char *uplo, char *diag, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info);
int strtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, float *a,
                 aocl_int64_t *lda, float *b, aocl_int64_t *ldb, aocl_int64_t *info);
int strttf_check(char *transr, char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *arf,
                 aocl_int64_t *info);
int strttp_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *ap, aocl_int64_t *info);
int stzrqf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, aocl_int64_t *info);
int stzrzf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *tau, float *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zbbcsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, aocl_int64_t *m,
                 aocl_int64_t *p, aocl_int64_t *q, double *theta, double *phi, dcomplex *u1, aocl_int64_t *ldu1,
                 dcomplex *u2, aocl_int64_t *ldu2, dcomplex *v1t, aocl_int64_t *ldv1t, dcomplex *v2t,
                 aocl_int64_t *ldv2t, double *b11d, double *b11e, double *b12d, double *b12e,
                 double *b21d, double *b21e, double *b22d, double *b22e, double *rwork,
                 aocl_int64_t *lrwork, aocl_int64_t *info);
int zbdsqr_check(char *uplo, aocl_int64_t *n, aocl_int64_t *ncvt, aocl_int64_t *nru, aocl_int64_t *ncc, double *d__,
                 double *e, dcomplex *vt, aocl_int64_t *ldvt, dcomplex *u, aocl_int64_t *ldu, dcomplex *c__,
                 aocl_int64_t *ldc, double *rwork, aocl_int64_t *info);
int zcgesv_check(aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, dcomplex *work, scomplex *swork,
                 double *rwork, aocl_int64_t *iter, aocl_int64_t *info);
int zcposv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, dcomplex *work, scomplex *swork,
                 double *rwork, aocl_int64_t *iter, aocl_int64_t *info);
int zdrscl_check(aocl_int64_t *n, double *sa, dcomplex *sx, aocl_int64_t *incx);
int zgbbrd_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ncc, aocl_int64_t *kl, aocl_int64_t *ku,
                 dcomplex *ab, aocl_int64_t *ldab, double *d__, double *e, dcomplex *q, aocl_int64_t *ldq,
                 dcomplex *pt, aocl_int64_t *ldpt, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 double *rwork, aocl_int64_t *info);
int zgbcon_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, double *anorm, double *rcond, dcomplex *work, double *rwork,
                 aocl_int64_t *info);
int zgbequ_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab, aocl_int64_t *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 aocl_int64_t *info);
int zgbequb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab, aocl_int64_t *ldab,
                  double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                  aocl_int64_t *info);
int zgbrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, dcomplex *ab,
                 aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zgbrfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  dcomplex *ab, aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                  double *r__, double *c__, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx,
                  double *rcond, double *berr, aocl_int64_t *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, aocl_int64_t *nparams, double *params, dcomplex *work,
                  double *rwork, aocl_int64_t *info);
int zgbsv_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, dcomplex *ab, aocl_int64_t *ldab,
                aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zgbsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                 dcomplex *ab, aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                 char *equed, double *r__, double *c__, dcomplex *b, aocl_int64_t *ldb, dcomplex *x,
                 aocl_int64_t *ldx, double *rcond, double *ferr, double *berr, dcomplex *work,
                 double *rwork, aocl_int64_t *info);
int zgbsvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs,
                  dcomplex *ab, aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv,
                  char *equed, double *r__, double *c__, dcomplex *b, aocl_int64_t *ldb, dcomplex *x,
                  aocl_int64_t *ldx, double *rcond, double *rpvgrw, double *berr, aocl_int64_t *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, aocl_int64_t *nparams,
                  double *params, dcomplex *work, double *rwork, aocl_int64_t *info);
int zgbtf2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int zgbtrf_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab, aocl_int64_t *ldab,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int zgbtrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *nrhs, dcomplex *ab,
                 aocl_int64_t *ldab, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zgebak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *scale,
                 aocl_int64_t *m, dcomplex *v, aocl_int64_t *ldv, aocl_int64_t *info);
int zgebal_check(char *job, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 double *scale, aocl_int64_t *info);
int zgebd2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *d__, double *e,
                 dcomplex *tauq, dcomplex *taup, dcomplex *work, aocl_int64_t *info);
int zgebrd_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *d__, double *e,
                 dcomplex *tauq, dcomplex *taup, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zgecon_check(char *norm, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *anorm, double *rcond,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zgeequ_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, aocl_int64_t *info);
int zgeequb_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *r__, double *c__,
                  double *rowcnd, double *colcnd, double *amax, aocl_int64_t *info);
int zgees_check(char *jobvs, char *sort, L_fp select, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                aocl_int64_t *sdim, dcomplex *w, dcomplex *vs, aocl_int64_t *ldvs, dcomplex *work,
                aocl_int64_t *lwork, double *rwork, logical *bwork, aocl_int64_t *info);
int zgeesx_check(char *jobvs, char *sort, L_fp select, char *sense, aocl_int64_t *n, dcomplex *a,
                 aocl_int64_t *lda, aocl_int64_t *sdim, dcomplex *w, dcomplex *vs, aocl_int64_t *ldvs,
                 double *rconde, double *rcondv, dcomplex *work, aocl_int64_t *lwork, double *rwork,
                 logical *bwork, aocl_int64_t *info);
int zgeev_check(char *jobvl, char *jobvr, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *w,
                dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr, aocl_int64_t *ldvr, dcomplex *work,
                aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zgeevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *w, dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *ilo, aocl_int64_t *ihi, double *scale, double *abnrm,
                 double *rconde, double *rcondv, dcomplex *work, aocl_int64_t *lwork, double *rwork,
                 aocl_int64_t *info);
int zgegs_check(char *jobvsl, char *jobvsr, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                aocl_int64_t *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vsl, aocl_int64_t *ldvsl,
                dcomplex *vsr, aocl_int64_t *ldvsr, dcomplex *work, aocl_int64_t *lwork, double *rwork,
                aocl_int64_t *info);
int zgegv_check(char *jobvl, char *jobvr, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                aocl_int64_t *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vl, aocl_int64_t *ldvl,
                dcomplex *vr, aocl_int64_t *ldvr, dcomplex *work, aocl_int64_t *lwork, double *rwork,
                aocl_int64_t *info);
int zgehd2_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *info);
int zgehrd_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zgelq2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *info);
int zgelqf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zgels_check(char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zgelsd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, double *s, double *rcond, aocl_int64_t *rank, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int zgelss_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, double *s, double *rcond, aocl_int64_t *rank, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zgelsx_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, double *rcond, aocl_int64_t *rank, dcomplex *work,
                 double *rwork, aocl_int64_t *info);
int zgelsy_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *jpvt, double *rcond, aocl_int64_t *rank, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zgemqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb,
                  dcomplex *v, aocl_int64_t *ldv, dcomplex *t, aocl_int64_t *ldt, dcomplex *c__, aocl_int64_t *ldc,
                  dcomplex *work, aocl_int64_t *info);
int zgeql2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *info);
int zgeqlf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zgeqp3_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zgeqpf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt, dcomplex *tau,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zgeqr2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *info);
int zgeqr2p_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                  aocl_int64_t *info);
int zgeqrf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zgeqrfp_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int zgeqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, dcomplex *a, aocl_int64_t *lda, dcomplex *t,
                 aocl_int64_t *ldt, dcomplex *work, aocl_int64_t *info);
int zgeqrt2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int zgeqrt3_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *t, aocl_int64_t *ldt,
                  aocl_int64_t *info);
int zgerfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zgerfsx_check(char *trans, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *r__, double *c__, dcomplex *b,
                  aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond, double *berr,
                  aocl_int64_t *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  aocl_int64_t *nparams, double *params, dcomplex *work, double *rwork, aocl_int64_t *info);
int zgerq2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *info);
int zgerqf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zgesc2_check(aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *rhs, aocl_int_t *ipiv, aocl_int64_t *jpiv,
                 double *scale);
int zgesdd_check(char *jobz, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *s,
                 dcomplex *u, aocl_int64_t *ldu, dcomplex *vt, aocl_int64_t *ldvt, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int zgesv_check(aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int zgesvd_check(char *jobu, char *jobvt, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 double *s, dcomplex *u, aocl_int64_t *ldu, dcomplex *vt, aocl_int64_t *ldvt, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zgesvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, double *r__, double *c__,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond, double *ferr,
                 double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zgesvxx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, double *r__, double *c__,
                  dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond,
                  double *rpvgrw, double *berr, aocl_int64_t *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, aocl_int64_t *nparams, double *params, dcomplex *work,
                  double *rwork, aocl_int64_t *info);
int zgetc2_check(aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *jpiv,
                 aocl_int64_t *info);
int zgetf2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int zgetrf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int zgetrfnp_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int zgetrfnpi_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nfact, dcomplex *a, aocl_int64_t *lda,
                    aocl_int64_t *info);
int zgetri_check(aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zgetrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zggbak_check(char *job, char *side, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, double *lscale,
                 double *rscale, aocl_int64_t *m, dcomplex *v, aocl_int64_t *ldv, aocl_int64_t *info);
int zggbal_check(char *job, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb,
                 aocl_int64_t *ilo, aocl_int64_t *ihi, double *lscale, double *rscale, double *work,
                 aocl_int64_t *info);
int zgges_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, aocl_int64_t *n, dcomplex *a,
                aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *sdim, dcomplex *alpha,
                dcomplex *beta, dcomplex *vsl, aocl_int64_t *ldvsl, dcomplex *vsr, aocl_int64_t *ldvsr,
                dcomplex *work, aocl_int64_t *lwork, double *rwork, logical *bwork, aocl_int64_t *info);
int zggesx_check(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, aocl_int64_t *n,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *sdim,
                 dcomplex *alpha, dcomplex *beta, dcomplex *vsl, aocl_int64_t *ldvsl, dcomplex *vsr,
                 aocl_int64_t *ldvsr, double *rconde, double *rcondv, dcomplex *work, aocl_int64_t *lwork,
                 double *rwork, aocl_int64_t *iwork, aocl_int64_t *liwork, logical *bwork, aocl_int64_t *info);
int zggev_check(char *jobvl, char *jobvr, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                aocl_int64_t *ldb, dcomplex *alpha, dcomplex *beta, dcomplex *vl, aocl_int64_t *ldvl,
                dcomplex *vr, aocl_int64_t *ldvr, dcomplex *work, aocl_int64_t *lwork, double *rwork,
                aocl_int64_t *info);
int zggevx_check(char *balanc, char *jobvl, char *jobvr, char *sense, aocl_int64_t *n, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *alpha, dcomplex *beta,
                 dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr, aocl_int64_t *ldvr, aocl_int64_t *ilo,
                 aocl_int64_t *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm,
                 double *rconde, double *rcondv, dcomplex *work, aocl_int64_t *lwork, double *rwork,
                 aocl_int64_t *iwork, logical *bwork, aocl_int64_t *info);
int zggglm_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *d__, dcomplex *x, dcomplex *y, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zgghrd_check(char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *q, aocl_int64_t *ldq, dcomplex *z__,
                 aocl_int64_t *ldz, aocl_int64_t *info);
int zgglse_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *c__, dcomplex *d__, dcomplex *x, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zggqrf_check(aocl_int64_t *n, aocl_int64_t *m, aocl_int64_t *p, dcomplex *a, aocl_int64_t *lda, dcomplex *taua,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *taub, dcomplex *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int zggrqf_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *taua,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *taub, dcomplex *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int zggsvd_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *p, aocl_int64_t *k,
                 aocl_int64_t *l, dcomplex *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, double *alpha,
                 double *beta, dcomplex *u, aocl_int64_t *ldu, dcomplex *v, aocl_int64_t *ldv, dcomplex *q,
                 aocl_int64_t *ldq, dcomplex *work, double *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int zggsvp_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, double *tola, double *tolb,
                 aocl_int64_t *k, aocl_int64_t *l, dcomplex *u, aocl_int64_t *ldu, dcomplex *v, aocl_int64_t *ldv,
                 dcomplex *q, aocl_int64_t *ldq, aocl_int64_t *iwork, double *rwork, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *info);
int zgtcon_check(char *norm, aocl_int64_t *n, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2,
                 aocl_int_t *ipiv, double *anorm, double *rcond, dcomplex *work, aocl_int64_t *info);
int zgtrfs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du,
                 dcomplex *dlf, dcomplex *df, dcomplex *duf, dcomplex *du2, aocl_int_t *ipiv,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zgtsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int zgtsvx_check(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *dl, dcomplex *d__,
                 dcomplex *du, dcomplex *dlf, dcomplex *df, dcomplex *duf, dcomplex *du2,
                 aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zgttrf_check(aocl_int64_t *n, dcomplex *dl, dcomplex *d__, dcomplex *du, dcomplex *du2,
                 aocl_int_t *ipiv, aocl_int64_t *info);
int zgttrs_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *dl, dcomplex *d__, dcomplex *du,
                 dcomplex *du2, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zgtts2_check(aocl_int64_t *itrans, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *dl, dcomplex *d__,
                 dcomplex *du, dcomplex *du2, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb);
int zhbev_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab,
                double *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, double *rwork,
                aocl_int64_t *info);
int zhbevd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab,
                 double *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, aocl_int64_t *lwork,
                 double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int zhbevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab,
                 aocl_int64_t *ldab, dcomplex *q, aocl_int64_t *ldq, double *vl, double *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w, dcomplex *z__, aocl_int64_t *ldz,
                 dcomplex *work, double *rwork, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int zhbgst_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, dcomplex *ab,
                 aocl_int64_t *ldab, dcomplex *bb, aocl_int64_t *ldbb, dcomplex *x, aocl_int64_t *ldx,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zhbgv_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, dcomplex *ab,
                aocl_int64_t *ldab, dcomplex *bb, aocl_int64_t *ldbb, double *w, dcomplex *z__, aocl_int64_t *ldz,
                dcomplex *work, double *rwork, aocl_int64_t *info);
int zhbgvd_check(char *jobz, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb, dcomplex *ab,
                 aocl_int64_t *ldab, dcomplex *bb, aocl_int64_t *ldbb, double *w, dcomplex *z__, aocl_int64_t *ldz,
                 dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int zhbgvx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, aocl_int64_t *ka, aocl_int64_t *kb,
                 dcomplex *ab, aocl_int64_t *ldab, dcomplex *bb, aocl_int64_t *ldbb, dcomplex *q,
                 aocl_int64_t *ldq, double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol,
                 aocl_int64_t *m, double *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, double *rwork,
                 aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int zhbtrd_check(char *vect, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab,
                 double *d__, double *e, dcomplex *q, aocl_int64_t *ldq, dcomplex *work, aocl_int64_t *info);
int zhecon_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *anorm,
                 double *rcond, dcomplex *work, aocl_int64_t *info);
int zhecon_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      double *anorm, double *rcond, dcomplex *work, aocl_int64_t *info);
int zheequb_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *s, double *scond,
                  double *amax, dcomplex *work, aocl_int64_t *info);
int zheev_check(char *jobz, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *w,
                dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zheevd_check(char *jobz, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *w,
                 dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int zheevr_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m,
                 double *w, dcomplex *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork,
                 aocl_int64_t *info);
int zheevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m,
                 double *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, aocl_int64_t *lwork,
                 double *rwork, aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int zhegs2_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int zhegst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int zhegv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                dcomplex *b, aocl_int64_t *ldb, double *w, dcomplex *work, aocl_int64_t *lwork, double *rwork,
                aocl_int64_t *info);
int zhegvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *b, aocl_int64_t *ldb, double *w, dcomplex *work, aocl_int64_t *lwork,
                 double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int zhegvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, double *vl, double *vu, aocl_int64_t *il,
                 aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w, dcomplex *z__, aocl_int64_t *ldz,
                 dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *iwork, aocl_int64_t *ifail,
                 aocl_int64_t *info);
int zherfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zherfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *s, dcomplex *b, aocl_int64_t *ldb,
                  dcomplex *x, aocl_int64_t *ldx, double *rcond, double *berr, aocl_int64_t *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, aocl_int64_t *nparams,
                  double *params, dcomplex *work, double *rwork, aocl_int64_t *info);
int zhesv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zhesv_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                     aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *lwork,
                     aocl_int64_t *info);
int zhesvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x,
                 aocl_int64_t *ldx, double *rcond, double *ferr, double *berr, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zhesvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, double *s, dcomplex *b,
                  aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond, double *rpvgrw,
                  double *berr, aocl_int64_t *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, aocl_int64_t *nparams, double *params, dcomplex *work,
                  double *rwork, aocl_int64_t *info);
int zheswapr_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *i1, aocl_int64_t *i2);
int zhetd2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *d__, double *e,
                 dcomplex *tau, aocl_int64_t *info);
int zhetf2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int zhetf2_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      aocl_int64_t *info);
int zhetrd_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *d__, double *e,
                 dcomplex *tau, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zhetrf_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zhetrf_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zhetri_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                 aocl_int64_t *info);
int zhetri2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int zhetri2x_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                   aocl_int64_t *nb, aocl_int64_t *info);
int zhetri_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      dcomplex *work, aocl_int64_t *info);
int zhetrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zhetrs2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *info);
int zhetrs_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zhfrk_check(char *transr, char *uplo, char *trans, aocl_int64_t *n, aocl_int64_t *k, double *alpha,
                dcomplex *a, aocl_int64_t *lda, double *beta, dcomplex *c__);
int zhgeqz_check(char *job, char *compq, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 dcomplex *h__, aocl_int64_t *ldh, dcomplex *t, aocl_int64_t *ldt, dcomplex *alpha,
                 dcomplex *beta, dcomplex *q, aocl_int64_t *ldq, dcomplex *z__, aocl_int64_t *ldz,
                 dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zhpcon_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int_t *ipiv, double *anorm, double *rcond,
                 dcomplex *work, aocl_int64_t *info);
int zhpev_check(char *jobz, char *uplo, aocl_int64_t *n, dcomplex *ap, double *w, dcomplex *z__,
                aocl_int64_t *ldz, dcomplex *work, double *rwork, aocl_int64_t *info);
int zhpevd_check(char *jobz, char *uplo, aocl_int64_t *n, dcomplex *ap, double *w, dcomplex *z__,
                 aocl_int64_t *ldz, dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *lrwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int zhpevx_check(char *jobz, char *range, char *uplo, aocl_int64_t *n, dcomplex *ap, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w,
                 dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, double *rwork, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int zhpgst_check(aocl_int64_t *itype, char *uplo, aocl_int64_t *n, dcomplex *ap, dcomplex *bp, aocl_int64_t *info);
int zhpgv_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, dcomplex *ap, dcomplex *bp,
                double *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, double *rwork,
                aocl_int64_t *info);
int zhpgvd_check(aocl_int64_t *itype, char *jobz, char *uplo, aocl_int64_t *n, dcomplex *ap, dcomplex *bp,
                 double *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, aocl_int64_t *lwork,
                 double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int zhpgvx_check(aocl_int64_t *itype, char *jobz, char *range, char *uplo, aocl_int64_t *n, dcomplex *ap,
                 dcomplex *bp, double *vl, double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol,
                 aocl_int64_t *m, double *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, double *rwork,
                 aocl_int64_t *iwork, aocl_int64_t *ifail, aocl_int64_t *info);
int zhprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *afp, aocl_int_t *ipiv,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zhpsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, aocl_int_t *ipiv, dcomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int zhpsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *afp,
                 aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zhptrd_check(char *uplo, aocl_int64_t *n, dcomplex *ap, double *d__, double *e, dcomplex *tau,
                 aocl_int64_t *info);
int zhptrf_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int_t *ipiv, aocl_int64_t *info);
int zhptri_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int_t *ipiv, dcomplex *work,
                 aocl_int64_t *info);
int zhptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, aocl_int_t *ipiv, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int zhsein_check(char *side, char *eigsrc, char *initv, logical *select, aocl_int64_t *n, dcomplex *h__,
                 aocl_int64_t *ldh, dcomplex *w, dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m, dcomplex *work, double *rwork,
                 aocl_int64_t *ifaill, aocl_int64_t *ifailr, aocl_int64_t *info);
int zhseqr_check(char *job, char *compz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, dcomplex *h__,
                 aocl_int64_t *ldh, dcomplex *w, dcomplex *z__, aocl_int64_t *ldz, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zla_gbamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, double *alpha,
                    dcomplex *ab, aocl_int64_t *ldab, dcomplex *x, aocl_int64_t *incx, double *beta,
                    double *y, aocl_int64_t *incy);
double zla_gbrcond_c_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab,
                           aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, double *c__,
                           logical *capply, aocl_int64_t *info, dcomplex *work, double *rwork);
double zla_gbrcond_x_check(char *trans, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab,
                           aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, dcomplex *x,
                           aocl_int64_t *info, dcomplex *work, double *rwork);
int zla_gbrfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n, aocl_int64_t *kl,
                              aocl_int64_t *ku, aocl_int64_t *nrhs, dcomplex *ab, aocl_int64_t *ldab,
                              dcomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, logical *colequ,
                              double *c__, dcomplex *b, aocl_int64_t *ldb, dcomplex *y, aocl_int64_t *ldy,
                              double *berr_out__, aocl_int64_t *n_norms__, double *err_bnds_norm__,
                              double *err_bnds_comp__, dcomplex *res, double *ayb, dcomplex *dy,
                              dcomplex *y_tail__, double *rcond, aocl_int64_t *ithresh, double *rthresh,
                              double *dz_ub__, logical *ignore_cwise__, aocl_int64_t *info);
double zla_gbrpvgrw_check(aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, aocl_int64_t *ncols, dcomplex *ab,
                          aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb);
int zla_geamv_check(aocl_int64_t *trans, aocl_int64_t *m, aocl_int64_t *n, double *alpha, dcomplex *a,
                    aocl_int64_t *lda, dcomplex *x, aocl_int64_t *incx, double *beta, double *y,
                    aocl_int64_t *incy);
double zla_gercond_c_check(char *trans, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, aocl_int_t *ipiv, double *c__, logical *capply,
                           aocl_int64_t *info, dcomplex *work, double *rwork);
double zla_gercond_x_check(char *trans, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *x, aocl_int64_t *info, dcomplex *work,
                           double *rwork);
int zla_gerfsx_extended_check(aocl_int64_t *prec_type__, aocl_int64_t *trans_type__, aocl_int64_t *n,
                              aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *af, aocl_int64_t *ldaf,
                              aocl_int_t *ipiv, logical *colequ, double *c__, dcomplex *b,
                              aocl_int64_t *ldb, dcomplex *y, aocl_int64_t *ldy, double *berr_out__,
                              aocl_int64_t *n_norms__, double *errs_n__, double *errs_c__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              aocl_int64_t *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
double zla_gerpvgrw_check(aocl_int64_t *n, aocl_int64_t *ncols, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                          aocl_int64_t *ldaf);
int zla_heamv_check(aocl_int64_t *uplo, aocl_int64_t *n, double *alpha, dcomplex *a, aocl_int64_t *lda,
                    dcomplex *x, aocl_int64_t *incx, double *beta, double *y, aocl_int64_t *incy);
double zla_hercond_c_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, aocl_int_t *ipiv, double *c__, logical *capply,
                           aocl_int64_t *info, dcomplex *work, double *rwork);
double zla_hercond_x_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *x, aocl_int64_t *info, dcomplex *work,
                           double *rwork);
int zla_herfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              dcomplex *a, aocl_int64_t *lda, dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv,
                              logical *colequ, double *c__, dcomplex *b, aocl_int64_t *ldb, dcomplex *y,
                              aocl_int64_t *ldy, double *berr_out__, aocl_int64_t *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              aocl_int64_t *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
double zla_herpvgrw_check(char *uplo, aocl_int64_t *n, aocl_int64_t *info, dcomplex *a, aocl_int64_t *lda,
                          dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *work);
int zla_lin_berr_check(aocl_int64_t *n, aocl_int64_t *nz, aocl_int64_t *nrhs, dcomplex *res, double *ayb,
                       double *berr);
double zla_porcond_c_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, double *c__, logical *capply, aocl_int64_t *info,
                           dcomplex *work, double *rwork);
double zla_porcond_x_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, dcomplex *x, aocl_int64_t *info, dcomplex *work,
                           double *rwork);
int zla_porfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              dcomplex *a, aocl_int64_t *lda, dcomplex *af, aocl_int64_t *ldaf,
                              logical *colequ, double *c__, dcomplex *b, aocl_int64_t *ldb, dcomplex *y,
                              aocl_int64_t *ldy, double *berr_out__, aocl_int64_t *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              aocl_int64_t *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
double zla_porpvgrw_check(char *uplo, aocl_int64_t *ncols, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                          aocl_int64_t *ldaf, double *work);
int zla_syamv_check(aocl_int64_t *uplo, aocl_int64_t *n, double *alpha, dcomplex *a, aocl_int64_t *lda,
                    dcomplex *x, aocl_int64_t *incx, double *beta, double *y, aocl_int64_t *incy);
double zla_syrcond_c_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, aocl_int_t *ipiv, double *c__, logical *capply,
                           aocl_int64_t *info, dcomplex *work, double *rwork);
double zla_syrcond_x_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                           aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *x, aocl_int64_t *info, dcomplex *work,
                           double *rwork);
int zla_syrfsx_extended_check(aocl_int64_t *prec_type__, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs,
                              dcomplex *a, aocl_int64_t *lda, dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv,
                              logical *colequ, double *c__, dcomplex *b, aocl_int64_t *ldb, dcomplex *y,
                              aocl_int64_t *ldy, double *berr_out__, aocl_int64_t *n_norms__,
                              double *err_bnds_norm__, double *err_bnds_comp__, dcomplex *res,
                              double *ayb, dcomplex *dy, dcomplex *y_tail__, double *rcond,
                              aocl_int64_t *ithresh, double *rthresh, double *dz_ub__,
                              logical *ignore_cwise__, aocl_int64_t *info);
double zla_syrpvgrw_check(char *uplo, aocl_int64_t *n, aocl_int64_t *info, dcomplex *a, aocl_int64_t *lda,
                          dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *work);
int zla_wwaddw_check(aocl_int64_t *n, dcomplex *x, dcomplex *y, dcomplex *w);
int zlabrd_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nb, dcomplex *a, aocl_int64_t *lda, double *d__,
                 double *e, dcomplex *tauq, dcomplex *taup, dcomplex *x, aocl_int64_t *ldx, dcomplex *y,
                 aocl_int64_t *ldy);
int zlacgv_check(aocl_int64_t *n, dcomplex *x, aocl_int64_t *incx);
int zlacn2_check(aocl_int64_t *n, dcomplex *v, dcomplex *x, double *est, aocl_int64_t *kase, aocl_int64_t *isave);
int zlacon_check(aocl_int64_t *n, dcomplex *v, dcomplex *x, double *est, aocl_int64_t *kase);
int zlacp2_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb);
int zlacpy_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb);
int zlacrm_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *b, aocl_int64_t *ldb,
                 dcomplex *c__, aocl_int64_t *ldc, double *rwork);
int zlacrt_check(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy, aocl_int64_t *incy,
                 dcomplex *c__, dcomplex *s);
void zladiv_check(dcomplex *ret_val, dcomplex *x, dcomplex *y);
int zlaed0_check(aocl_int64_t *qsiz, aocl_int64_t *n, double *d__, double *e, dcomplex *q, aocl_int64_t *ldq,
                 dcomplex *qstore, aocl_int64_t *ldqs, double *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int zlaed7_check(aocl_int64_t *n, aocl_int64_t *cutpnt, aocl_int64_t *qsiz, aocl_int64_t *tlvls, aocl_int64_t *curlvl,
                 aocl_int64_t *curpbm, double *d__, dcomplex *q, aocl_int64_t *ldq, double *rho,
                 aocl_int64_t *indxq, double *qstore, aocl_int64_t *qptr, aocl_int64_t *prmptr, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, double *givnum, dcomplex *work, double *rwork,
                 aocl_int64_t *iwork, aocl_int64_t *info);
int zlaed8_check(aocl_int64_t *k, aocl_int64_t *n, aocl_int64_t *qsiz, dcomplex *q, aocl_int64_t *ldq, double *d__,
                 double *rho, aocl_int64_t *cutpnt, double *z__, double *dlamda, dcomplex *q2,
                 aocl_int64_t *ldq2, double *w, aocl_int64_t *indxp, aocl_int64_t *indx, aocl_int64_t *indxq,
                 aocl_int64_t *perm, aocl_int64_t *givptr, aocl_int64_t *givcol, double *givnum, aocl_int64_t *info);
int zlaein_check(logical *rightv, logical *noinit, aocl_int64_t *n, dcomplex *h__, aocl_int64_t *ldh,
                 dcomplex *w, dcomplex *v, dcomplex *b, aocl_int64_t *ldb, double *rwork, double *eps3,
                 double *smlnum, aocl_int64_t *info);
int zlaesy_check(dcomplex *a, dcomplex *b, dcomplex *c__, dcomplex *rt1, dcomplex *rt2,
                 dcomplex *evscal, dcomplex *cs1, dcomplex *sn1);
int zlaev2_check(dcomplex *a, dcomplex *b, dcomplex *c__, double *rt1, double *rt2, double *cs1,
                 dcomplex *sn1);
int zlag2c_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, scomplex *sa, aocl_int64_t *ldsa,
                 aocl_int64_t *info);
int zlags2_check(logical *upper, double *a1, dcomplex *a2, double *a3, double *b1, dcomplex *b2,
                 double *b3, double *csu, dcomplex *snu, double *csv, dcomplex *snv, double *csq,
                 dcomplex *snq);
int zlagtm_check(char *trans, aocl_int64_t *n, aocl_int64_t *nrhs, double *alpha, dcomplex *dl, dcomplex *d__,
                 dcomplex *du, dcomplex *x, aocl_int64_t *ldx, double *beta, dcomplex *b, aocl_int64_t *ldb);
int zlahef_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, dcomplex *a, aocl_int64_t *lda,
                 aocl_int_t *ipiv, dcomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int zlahef_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, dcomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, dcomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int zlahqr_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 dcomplex *h__, aocl_int64_t *ldh, dcomplex *w, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 dcomplex *z__, aocl_int64_t *ldz, aocl_int64_t *info);
int zlahr2_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *t, aocl_int64_t *ldt, dcomplex *y, aocl_int64_t *ldy);
int zlahrd_check(aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *nb, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *t, aocl_int64_t *ldt, dcomplex *y, aocl_int64_t *ldy);
int zlaic1_check(aocl_int64_t *job, aocl_int64_t *j, dcomplex *x, double *sest, dcomplex *w, dcomplex *gamma,
                 double *sestpr, dcomplex *s, dcomplex *c__);
int zlals0_check(aocl_int64_t *icompq, aocl_int64_t *nl, aocl_int64_t *nr, aocl_int64_t *sqre, aocl_int64_t *nrhs,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *bx, aocl_int64_t *ldbx, aocl_int64_t *perm,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol, double *givnum, aocl_int64_t *ldgnum,
                 double *poles, double *difl, double *difr, double *z__, aocl_int64_t *k, double *c__,
                 double *s, double *rwork, aocl_int64_t *info);
int zlalsa_check(aocl_int64_t *icompq, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *bx, aocl_int64_t *ldbx, double *u, aocl_int64_t *ldu, double *vt,
                 aocl_int64_t *k, double *difl, double *difr, double *z__, double *poles,
                 aocl_int64_t *givptr, aocl_int64_t *givcol, aocl_int64_t *ldgcol, aocl_int64_t *perm, double *givnum,
                 double *c__, double *s, double *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
int zlalsd_check(char *uplo, aocl_int64_t *smlsiz, aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, double *e,
                 dcomplex *b, aocl_int64_t *ldb, double *rcond, aocl_int64_t *rank, dcomplex *work,
                 double *rwork, aocl_int64_t *iwork, aocl_int64_t *info);
double zlangb_check(char *norm, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab, aocl_int64_t *ldab,
                    double *work);
double zlange_check(char *norm, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *work);
double zlangt_check(char *norm, aocl_int64_t *n, dcomplex *dl, dcomplex *d__, dcomplex *du);
double zlanhb_check(char *norm, char *uplo, aocl_int64_t *n, aocl_int64_t *k, dcomplex *ab, aocl_int64_t *ldab,
                    double *work);
double zlanhe_check(char *norm, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *work);
double zlanhf_check(char *norm, char *transr, char *uplo, aocl_int64_t *n, dcomplex *a, double *work);
double zlanhp_check(char *norm, char *uplo, aocl_int64_t *n, dcomplex *ap, double *work);
double zlanhs_check(char *norm, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *work);
double zlanht_check(char *norm, aocl_int64_t *n, double *d__, dcomplex *e);
double zlansb_check(char *norm, char *uplo, aocl_int64_t *n, aocl_int64_t *k, dcomplex *ab, aocl_int64_t *ldab,
                    double *work);
double zlansp_check(char *norm, char *uplo, aocl_int64_t *n, dcomplex *ap, double *work);
double zlansy_check(char *norm, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *work);
double zlantb_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *k, dcomplex *ab,
                    aocl_int64_t *ldab, double *work);
double zlantp_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, dcomplex *ap, double *work);
double zlantr_check(char *norm, char *uplo, char *diag, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a,
                    aocl_int64_t *lda, double *work);
int zlapll_check(aocl_int64_t *n, dcomplex *x, aocl_int64_t *incx, dcomplex *y, aocl_int64_t *incy, double *ssmin);
int zlapmr_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, dcomplex *x, aocl_int64_t *ldx, aocl_int64_t *k);
int zlapmt_check(logical *forwrd, aocl_int64_t *m, aocl_int64_t *n, dcomplex *x, aocl_int64_t *ldx, aocl_int64_t *k);
int zlaqgb_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *kl, aocl_int64_t *ku, dcomplex *ab, aocl_int64_t *ldab,
                 double *r__, double *c__, double *rowcnd, double *colcnd, double *amax,
                 char *equed);
int zlaqge_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *r__, double *c__,
                 double *rowcnd, double *colcnd, double *amax, char *equed);
int zlaqhb_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab, double *s,
                 double *scond, double *amax, char *equed);
int zlaqhe_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *s, double *scond,
                 double *amax, char *equed);
int zlaqhp_check(char *uplo, aocl_int64_t *n, dcomplex *ap, double *s, double *scond, double *amax,
                 char *equed);
int zlaqp2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt,
                 dcomplex *tau, double *vn1, double *vn2, dcomplex *work);
int zlaqps_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *offset, aocl_int64_t *nb, aocl_int64_t *kb, dcomplex *a,
                 aocl_int64_t *lda, aocl_int64_t *jpvt, dcomplex *tau, double *vn1, double *vn2,
                 dcomplex *auxv, dcomplex *f, aocl_int64_t *ldf);
int zlaqr0_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 dcomplex *h__, aocl_int64_t *ldh, dcomplex *w, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zlaqr1_check(aocl_int64_t *n, dcomplex *h__, aocl_int64_t *ldh, dcomplex *s1, dcomplex *s2, dcomplex *v);
int zlaqr2_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, dcomplex *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 dcomplex *z__, aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, dcomplex *sh, dcomplex *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, dcomplex *t, aocl_int64_t *ldt, aocl_int64_t *nv, dcomplex *wv,
                 aocl_int64_t *ldwv, dcomplex *work, aocl_int64_t *lwork);
int zlaqr3_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ktop, aocl_int64_t *kbot,
                 aocl_int64_t *nw, dcomplex *h__, aocl_int64_t *ldh, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 dcomplex *z__, aocl_int64_t *ldz, aocl_int64_t *ns, aocl_int64_t *nd, dcomplex *sh, dcomplex *v,
                 aocl_int64_t *ldv, aocl_int64_t *nh, dcomplex *t, aocl_int64_t *ldt, aocl_int64_t *nv, dcomplex *wv,
                 aocl_int64_t *ldwv, dcomplex *work, aocl_int64_t *lwork);
int zlaqr4_check(logical *wantt, logical *wantz, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 dcomplex *h__, aocl_int64_t *ldh, dcomplex *w, aocl_int64_t *iloz, aocl_int64_t *ihiz,
                 dcomplex *z__, aocl_int64_t *ldz, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zlaqr5_check(logical *wantt, logical *wantz, aocl_int64_t *kacc22, aocl_int64_t *n, aocl_int64_t *ktop,
                 aocl_int64_t *kbot, aocl_int64_t *nshfts, dcomplex *s, dcomplex *h__, aocl_int64_t *ldh,
                 aocl_int64_t *iloz, aocl_int64_t *ihiz, dcomplex *z__, aocl_int64_t *ldz, dcomplex *v,
                 aocl_int64_t *ldv, dcomplex *u, aocl_int64_t *ldu, aocl_int64_t *nv, dcomplex *wv, aocl_int64_t *ldwv,
                 aocl_int64_t *nh, dcomplex *wh, aocl_int64_t *ldwh);
int zlaqsb_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab, double *s,
                 double *scond, double *amax, char *equed);
int zlaqsp_check(char *uplo, aocl_int64_t *n, dcomplex *ap, double *s, double *scond, double *amax,
                 char *equed);
int zlaqsy_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *s, double *scond,
                 double *amax, char *equed);
int zlar1v_check(aocl_int64_t *n, aocl_int64_t *b1, aocl_int64_t *bn, double *lambda, double *d__, double *l,
                 double *ld, double *lld, double *pivmin, double *gaptol, dcomplex *z__,
                 logical *wantnc, aocl_int64_t *negcnt, double *ztz, double *mingma, aocl_int64_t *r__,
                 aocl_int64_t *isuppz, double *nrminv, double *resid, double *rqcorr, double *work);
int zlar2v_check(aocl_int64_t *n, dcomplex *x, dcomplex *y, dcomplex *z__, aocl_int64_t *incx, double *c__,
                 dcomplex *s, aocl_int64_t *incc);
int zlarcm_check(aocl_int64_t *m, aocl_int64_t *n, double *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb,
                 dcomplex *c__, aocl_int64_t *ldc, double *rwork);
int zlarf_check(char *side, aocl_int64_t *m, aocl_int64_t *n, dcomplex *v, aocl_int64_t *incv, dcomplex *tau,
                dcomplex *c__, aocl_int64_t *ldc, dcomplex *work);
int zlarfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, dcomplex *v, aocl_int64_t *ldv, dcomplex *t, aocl_int64_t *ldt, dcomplex *c__,
                 aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *ldwork);
int zlarfg_check(aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx, dcomplex *tau);
int zlarfgp_check(aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx, dcomplex *tau);
int zlarft_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, dcomplex *v, aocl_int64_t *ldv,
                 dcomplex *tau, dcomplex *t, aocl_int64_t *ldt);
int zlarfx_check(char *side, aocl_int64_t *m, aocl_int64_t *n, dcomplex *v, dcomplex *tau, dcomplex *c__,
                 aocl_int64_t *ldc, dcomplex *work);
int zlargv_check(aocl_int64_t *n, dcomplex *x, aocl_int64_t *incx, dcomplex *y, aocl_int64_t *incy, double *c__,
                 aocl_int64_t *incc);
int zlarnv_check(aocl_int64_t *idist, aocl_int64_t *iseed, aocl_int64_t *n, dcomplex *x);
int zlarrv_check(aocl_int64_t *n, double *vl, double *vu, double *d__, double *l, double *pivmin,
                 aocl_int64_t *isplit, aocl_int64_t *m, aocl_int64_t *dol, aocl_int64_t *dou, double *minrgp,
                 double *rtol1, double *rtol2, double *w, double *werr, double *wgap,
                 aocl_int64_t *iblock, aocl_int64_t *indexw, double *gers, dcomplex *z__, aocl_int64_t *ldz,
                 aocl_int64_t *isuppz, double *work, aocl_int64_t *iwork, aocl_int64_t *info);
int zlarscl2_check(aocl_int64_t *m, aocl_int64_t *n, double *d__, dcomplex *x, aocl_int64_t *ldx);
int zlartg_check(dcomplex *f, dcomplex *g, double *cs, dcomplex *sn, dcomplex *r__);
int zlartv_check(aocl_int64_t *n, dcomplex *x, aocl_int64_t *incx, dcomplex *y, aocl_int64_t *incy, double *c__,
                 dcomplex *s, aocl_int64_t *incc);
int zlarz_check(char *side, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, dcomplex *v, aocl_int64_t *incv,
                dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work);
int zlarzb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, dcomplex *v, aocl_int64_t *ldv, dcomplex *t, aocl_int64_t *ldt,
                 dcomplex *c__, aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *ldwork);
int zlarzt_check(char *direct, char *storev, aocl_int64_t *n, aocl_int64_t *k, dcomplex *v, aocl_int64_t *ldv,
                 dcomplex *tau, dcomplex *t, aocl_int64_t *ldt);
int zlascl_check(char *type__, aocl_int64_t *kl, aocl_int64_t *ku, double *cfrom, double *cto, aocl_int64_t *m,
                 aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int zlascl2_check(aocl_int64_t *m, aocl_int64_t *n, double *d__, dcomplex *x, aocl_int64_t *ldx);
int zlaset_check(char *uplo, aocl_int64_t *m, aocl_int64_t *n, dcomplex *alpha, dcomplex *beta, dcomplex *a,
                 aocl_int64_t *lda);
int zlasr_check(char *side, char *pivot, char *direct, aocl_int64_t *m, aocl_int64_t *n, double *c__,
                double *s, dcomplex *a, aocl_int64_t *lda);
int zlassq_check(aocl_int64_t *n, dcomplex *x, aocl_int64_t *incx, double *scale, double *sumsq);
int zlaswp_check(aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *k1, aocl_int64_t *k2, aocl_int_t *ipiv,
                 aocl_int64_t *incx);
int zlasyf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, dcomplex *a, aocl_int64_t *lda,
                 aocl_int_t *ipiv, dcomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int zlasyf_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, aocl_int64_t *kb, dcomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, dcomplex *w, aocl_int64_t *ldw, aocl_int64_t *info);
int zlat2c_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, scomplex *sa, aocl_int64_t *ldsa,
                 aocl_int64_t *info);
int zlatbs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, aocl_int64_t *kd,
                 dcomplex *ab, aocl_int64_t *ldab, dcomplex *x, double *scale, double *cnorm,
                 aocl_int64_t *info);
int zlatdf_check(aocl_int64_t *ijob, aocl_int64_t *n, dcomplex *z__, aocl_int64_t *ldz, dcomplex *rhs,
                 double *rdsum, double *rdscal, aocl_int_t *ipiv, aocl_int64_t *jpiv);
int zlatps_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, dcomplex *ap,
                 dcomplex *x, double *scale, double *cnorm, aocl_int64_t *info);
int zlatrd_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nb, dcomplex *a, aocl_int64_t *lda, double *e,
                 dcomplex *tau, dcomplex *w, aocl_int64_t *ldw);
int zlatrs_check(char *uplo, char *trans, char *diag, char *normin, aocl_int64_t *n, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *x, double *scale, double *cnorm, aocl_int64_t *info);
int zlatrz_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work);
int zlatzm_check(char *side, aocl_int64_t *m, aocl_int64_t *n, dcomplex *v, aocl_int64_t *incv, dcomplex *tau,
                 dcomplex *c1, dcomplex *c2, aocl_int64_t *ldc, dcomplex *work);
int zlauu2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int zlauum_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int zpbcon_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab, double *anorm,
                 double *rcond, dcomplex *work, double *rwork, aocl_int64_t *info);
int zpbequ_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab, double *s,
                 double *scond, double *amax, aocl_int64_t *info);
int zpbrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, dcomplex *ab, aocl_int64_t *ldab,
                 dcomplex *afb, aocl_int64_t *ldafb, dcomplex *b, aocl_int64_t *ldb, dcomplex *x,
                 aocl_int64_t *ldx, double *ferr, double *berr, dcomplex *work, double *rwork,
                 aocl_int64_t *info);
int zpbstf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int zpbsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, dcomplex *ab, aocl_int64_t *ldab,
                dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zpbsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, dcomplex *ab,
                 aocl_int64_t *ldab, dcomplex *afb, aocl_int64_t *ldafb, char *equed, double *s, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zpbtf2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int zpbtrf_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab, aocl_int64_t *ldab, aocl_int64_t *info);
int zpbtrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs, dcomplex *ab, aocl_int64_t *ldab,
                 dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zpftrf_check(char *transr, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *info);
int zpftri_check(char *transr, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *info);
int zpftrs_check(char *transr, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int zpocon_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *anorm, double *rcond,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zpoequ_check(aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *s, double *scond, double *amax,
                 aocl_int64_t *info);
int zpoequb_check(aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *s, double *scond, double *amax,
                  aocl_int64_t *info);
int zporfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                 aocl_int64_t *ldaf, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr,
                 double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zporfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, double *s, dcomplex *b, aocl_int64_t *ldb, dcomplex *x,
                  aocl_int64_t *ldx, double *rcond, double *berr, aocl_int64_t *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, aocl_int64_t *nparams,
                  double *params, dcomplex *work, double *rwork, aocl_int64_t *info);
int zposv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int zposvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *af, aocl_int64_t *ldaf, char *equed, double *s, dcomplex *b, aocl_int64_t *ldb,
                 dcomplex *x, aocl_int64_t *ldx, double *rcond, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zposvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, char *equed, double *s, dcomplex *b, aocl_int64_t *ldb,
                  dcomplex *x, aocl_int64_t *ldx, double *rcond, double *rpvgrw, double *berr,
                  aocl_int64_t *n_err_bnds__, double *err_bnds_norm__, double *err_bnds_comp__,
                  aocl_int64_t *nparams, double *params, dcomplex *work, double *rwork, aocl_int64_t *info);
int zpotf2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int zpotrf_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int zpotri_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int zpotrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int zppcon_check(char *uplo, aocl_int64_t *n, dcomplex *ap, double *anorm, double *rcond, dcomplex *work,
                 double *rwork, aocl_int64_t *info);
int zppequ_check(char *uplo, aocl_int64_t *n, dcomplex *ap, double *s, double *scond, double *amax,
                 aocl_int64_t *info);
int zpprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *afp, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zppsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int zppsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *afp,
                 char *equed, double *s, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx,
                 double *rcond, double *ferr, double *berr, dcomplex *work, double *rwork,
                 aocl_int64_t *info);
int zpptrf_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int64_t *info);
int zpptri_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int64_t *info);
int zpptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *b, aocl_int64_t *ldb,
                 aocl_int64_t *info);
int zpstf2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 double *tol, double *work, aocl_int64_t *info);
int zpstrf_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *piv, aocl_int64_t *rank,
                 double *tol, double *work, aocl_int64_t *info);
int zptcon_check(aocl_int64_t *n, double *d__, dcomplex *e, double *anorm, double *rcond, double *rwork,
                 aocl_int64_t *info);
int zpteqr_check(char *compz, aocl_int64_t *n, double *d__, double *e, dcomplex *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *info);
int zptrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, dcomplex *e, double *df,
                 dcomplex *ef, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr,
                 double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zptsv_check(aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, dcomplex *e, dcomplex *b, aocl_int64_t *ldb,
                aocl_int64_t *info);
int zptsvx_check(char *fact, aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, dcomplex *e, double *df,
                 dcomplex *ef, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zpttrf_check(aocl_int64_t *n, double *d__, dcomplex *e, aocl_int64_t *info);
int zpttrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, dcomplex *e, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int zptts2_check(aocl_int64_t *iuplo, aocl_int64_t *n, aocl_int64_t *nrhs, double *d__, dcomplex *e, dcomplex *b,
                 aocl_int64_t *ldb);
int zrot_check(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy, aocl_int64_t *incy, double *c__,
               dcomplex *s);
int zspcon_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int_t *ipiv, double *anorm, double *rcond,
                 dcomplex *work, aocl_int64_t *info);
int zspmv_check(char *uplo, aocl_int64_t *n, dcomplex *alpha, dcomplex *ap, dcomplex *x, aocl_int64_t *incx,
                dcomplex *beta, dcomplex *y, aocl_int64_t *incy);
int zspr_check(char *uplo, aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx, dcomplex *ap);
int zsprfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *afp, aocl_int_t *ipiv,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int zspsv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, aocl_int_t *ipiv, dcomplex *b,
                aocl_int64_t *ldb, aocl_int64_t *info);
int zspsvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, dcomplex *afp,
                 aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zsptrf_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int_t *ipiv, aocl_int64_t *info);
int zsptri_check(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int_t *ipiv, dcomplex *work,
                 aocl_int64_t *info);
int zsptrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap, aocl_int_t *ipiv, dcomplex *b,
                 aocl_int64_t *ldb, aocl_int64_t *info);
int zstedc_check(char *compz, aocl_int64_t *n, double *d__, double *e, dcomplex *z__, aocl_int64_t *ldz,
                 dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork,
                 aocl_int64_t *liwork, aocl_int64_t *info);
int zstegr_check(char *jobz, char *range, aocl_int64_t *n, double *d__, double *e, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, double *abstol, aocl_int64_t *m, double *w,
                 dcomplex *z__, aocl_int64_t *ldz, aocl_int64_t *isuppz, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int zstein_check(aocl_int64_t *n, double *d__, double *e, aocl_int64_t *m, double *w, aocl_int64_t *iblock,
                 aocl_int64_t *isplit, dcomplex *z__, aocl_int64_t *ldz, double *work, aocl_int64_t *iwork,
                 aocl_int64_t *ifail, aocl_int64_t *info);
int zstemr_check(char *jobz, char *range, aocl_int64_t *n, double *d__, double *e, double *vl,
                 double *vu, aocl_int64_t *il, aocl_int64_t *iu, aocl_int64_t *m, double *w, dcomplex *z__,
                 aocl_int64_t *ldz, aocl_int64_t *nzc, aocl_int64_t *isuppz, logical *tryrac, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int zsteqr_check(char *compz, aocl_int64_t *n, double *d__, double *e, dcomplex *z__, aocl_int64_t *ldz,
                 double *work, aocl_int64_t *info);
int zsycon_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, double *anorm,
                 double *rcond, dcomplex *work, aocl_int64_t *info);
int zsycon_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      double *anorm, double *rcond, dcomplex *work, aocl_int64_t *info);
int zsyconv_check(char *uplo, char *way, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  dcomplex *work, aocl_int64_t *info);
int zsyequb_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *s, double *scond,
                  double *amax, dcomplex *work, aocl_int64_t *info);
int zsymv_check(char *uplo, aocl_int64_t *n, dcomplex *alpha, dcomplex *a, aocl_int64_t *lda, dcomplex *x,
                aocl_int64_t *incx, dcomplex *beta, dcomplex *y, aocl_int64_t *incy);
int zsyr_check(char *uplo, aocl_int64_t *n, dcomplex *alpha, dcomplex *x, aocl_int64_t *incx, dcomplex *a,
               aocl_int64_t *lda);
int zsyrfs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, dcomplex *af,
                 aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int zsyrfsx_check(char *uplo, char *equed, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, double *s, dcomplex *b, aocl_int64_t *ldb,
                  dcomplex *x, aocl_int64_t *ldx, double *rcond, double *berr, aocl_int64_t *n_err_bnds__,
                  double *err_bnds_norm__, double *err_bnds_comp__, aocl_int64_t *nparams,
                  double *params, dcomplex *work, double *rwork, aocl_int64_t *info);
int zsysv_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zsysv_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                     aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *lwork,
                     aocl_int64_t *info);
int zsysvx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, dcomplex *x,
                 aocl_int64_t *ldx, double *rcond, double *ferr, double *berr, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *info);
int zsysvxx_check(char *fact, char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                  dcomplex *af, aocl_int64_t *ldaf, aocl_int_t *ipiv, char *equed, double *s, dcomplex *b,
                  aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *rcond, double *rpvgrw,
                  double *berr, aocl_int64_t *n_err_bnds__, double *err_bnds_norm__,
                  double *err_bnds_comp__, aocl_int64_t *nparams, double *params, dcomplex *work,
                  double *rwork, aocl_int64_t *info);
int zsyswapr_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *i1, aocl_int64_t *i2);
int zsytf2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info);
int zsytf2_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      aocl_int64_t *info);
int zsytrf_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zsytrf_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zsytri_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                 aocl_int64_t *info);
int zsytri2_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                  aocl_int64_t *lwork, aocl_int64_t *info);
int zsytri2x_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv, dcomplex *work,
                   aocl_int64_t *nb, aocl_int64_t *info);
int zsytri_rook_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                      dcomplex *work, aocl_int64_t *info);
int zsytrs_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                 dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int zsytrs2_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda, aocl_int_t *ipiv,
                  dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *info);
int zsytrs_rook_check(char *uplo, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a, aocl_int64_t *lda,
                      aocl_int_t *ipiv, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ztbcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, aocl_int64_t *kd, dcomplex *ab,
                 aocl_int64_t *ldab, double *rcond, dcomplex *work, double *rwork, aocl_int64_t *info);
int ztbrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 dcomplex *ab, aocl_int64_t *ldab, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx,
                 double *ferr, double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int ztbtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *kd, aocl_int64_t *nrhs,
                 dcomplex *ab, aocl_int64_t *ldab, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ztfsm_check(char *transr, char *side, char *uplo, char *trans, char *diag, aocl_int64_t *m,
                aocl_int64_t *n, dcomplex *alpha, dcomplex *a, dcomplex *b, aocl_int64_t *ldb);
int ztftri_check(char *transr, char *uplo, char *diag, aocl_int64_t *n, dcomplex *a, aocl_int64_t *info);
int ztfttp_check(char *transr, char *uplo, aocl_int64_t *n, dcomplex *arf, dcomplex *ap, aocl_int64_t *info);
int ztfttr_check(char *transr, char *uplo, aocl_int64_t *n, dcomplex *arf, dcomplex *a, aocl_int64_t *lda,
                 aocl_int64_t *info);
int ztgevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, dcomplex *s, aocl_int64_t *lds,
                 dcomplex *p, aocl_int64_t *ldp, dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr,
                 aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m, dcomplex *work, double *rwork,
                 aocl_int64_t *info);
int ztgex2_check(logical *wantq, logical *wantz, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *q, aocl_int64_t *ldq, dcomplex *z__, aocl_int64_t *ldz, aocl_int64_t *j1,
                 aocl_int64_t *info);
int ztgexc_check(logical *wantq, logical *wantz, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                 aocl_int64_t *ldb, dcomplex *q, aocl_int64_t *ldq, dcomplex *z__, aocl_int64_t *ldz,
                 aocl_int64_t *ifst, aocl_int64_t *ilst, aocl_int64_t *info);
int ztgsen_check(aocl_int64_t *ijob, logical *wantq, logical *wantz, logical *select, aocl_int64_t *n,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *alpha,
                 dcomplex *beta, dcomplex *q, aocl_int64_t *ldq, dcomplex *z__, aocl_int64_t *ldz, aocl_int64_t *m,
                 double *pl, double *pr, double *dif, dcomplex *work, aocl_int64_t *lwork,
                 aocl_int64_t *iwork, aocl_int64_t *liwork, aocl_int64_t *info);
int ztgsja_check(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *n, aocl_int64_t *k,
                 aocl_int64_t *l, dcomplex *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, double *tola,
                 double *tolb, double *alpha, double *beta, dcomplex *u, aocl_int64_t *ldu, dcomplex *v,
                 aocl_int64_t *ldv, dcomplex *q, aocl_int64_t *ldq, dcomplex *work, aocl_int64_t *ncycle,
                 aocl_int64_t *info);
int ztgsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr,
                 aocl_int64_t *ldvr, double *s, double *dif, aocl_int64_t *mm, aocl_int64_t *m, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int ztgsy2_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *c__, aocl_int64_t *ldc, dcomplex *d__,
                 aocl_int64_t *ldd, dcomplex *e, aocl_int64_t *lde, dcomplex *f, aocl_int64_t *ldf, double *scale,
                 double *rdsum, double *rdscal, aocl_int64_t *info);
int ztgsyl_check(char *trans, aocl_int64_t *ijob, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *c__, aocl_int64_t *ldc, dcomplex *d__,
                 aocl_int64_t *ldd, dcomplex *e, aocl_int64_t *lde, dcomplex *f, aocl_int64_t *ldf, double *scale,
                 double *dif, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *iwork, aocl_int64_t *info);
int ztpcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, dcomplex *ap, double *rcond,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int ztpmqrt_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                  aocl_int64_t *nb, dcomplex *v, aocl_int64_t *ldv, dcomplex *t, aocl_int64_t *ldt, dcomplex *a,
                  aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *work, aocl_int64_t *info);
int ztpqrt_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, aocl_int64_t *nb, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *t, aocl_int64_t *ldt, dcomplex *work,
                 aocl_int64_t *info);
int ztpqrt2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *l, dcomplex *a, aocl_int64_t *lda, dcomplex *b,
                  aocl_int64_t *ldb, dcomplex *t, aocl_int64_t *ldt, aocl_int64_t *info);
int ztprfb_check(char *side, char *trans, char *direct, char *storev, aocl_int64_t *m, aocl_int64_t *n,
                 aocl_int64_t *k, aocl_int64_t *l, dcomplex *v, aocl_int64_t *ldv, dcomplex *t, aocl_int64_t *ldt,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *work,
                 aocl_int64_t *ldwork);
int ztprfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap,
                 dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr, double *berr,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int ztptri_check(char *uplo, char *diag, aocl_int64_t *n, dcomplex *ap, aocl_int64_t *info);
int ztptrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *ap,
                 dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ztpttf_check(char *transr, char *uplo, aocl_int64_t *n, dcomplex *ap, dcomplex *arf, aocl_int64_t *info);
int ztpttr_check(char *uplo, aocl_int64_t *n, dcomplex *ap, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int ztrcon_check(char *norm, char *uplo, char *diag, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                 double *rcond, dcomplex *work, double *rwork, aocl_int64_t *info);
int ztrevc_check(char *side, char *howmny, logical *select, aocl_int64_t *n, dcomplex *t, aocl_int64_t *ldt,
                 dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr, aocl_int64_t *ldvr, aocl_int64_t *mm, aocl_int64_t *m,
                 dcomplex *work, double *rwork, aocl_int64_t *info);
int ztrexc_check(char *compq, aocl_int64_t *n, dcomplex *t, aocl_int64_t *ldt, dcomplex *q, aocl_int64_t *ldq,
                 aocl_int64_t *ifst, aocl_int64_t *ilst, aocl_int64_t *info);
int ztrrfs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *x, aocl_int64_t *ldx, double *ferr,
                 double *berr, dcomplex *work, double *rwork, aocl_int64_t *info);
int ztrsen_check(char *job, char *compq, logical *select, aocl_int64_t *n, dcomplex *t, aocl_int64_t *ldt,
                 dcomplex *q, aocl_int64_t *ldq, dcomplex *w, aocl_int64_t *m, double *s, double *sep,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int ztrsna_check(char *job, char *howmny, logical *select, aocl_int64_t *n, dcomplex *t, aocl_int64_t *ldt,
                 dcomplex *vl, aocl_int64_t *ldvl, dcomplex *vr, aocl_int64_t *ldvr, double *s, double *sep,
                 aocl_int64_t *mm, aocl_int64_t *m, dcomplex *work, aocl_int64_t *ldwork, double *rwork,
                 aocl_int64_t *info);
int ztrsyl_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *c__, aocl_int64_t *ldc,
                 double *scale, aocl_int64_t *info);
int ztrsyl3_check(char *trana, char *tranb, aocl_int64_t *isgn, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a,
                  aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, dcomplex *c__, aocl_int64_t *ldc,
                  doublereal *scale, doublereal *swork, aocl_int64_t *ldswork, aocl_int64_t *info);
int ztrti2_check(char *uplo, char *diag, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int ztrtri_check(char *uplo, char *diag, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *info);
int ztrtrs_check(char *uplo, char *trans, char *diag, aocl_int64_t *n, aocl_int64_t *nrhs, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *b, aocl_int64_t *ldb, aocl_int64_t *info);
int ztrttf_check(char *transr, char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *arf,
                 aocl_int64_t *info);
int ztrttp_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *ap, aocl_int64_t *info);
int ztzrqf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, aocl_int64_t *info);
int ztzrzf_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zunbdb_check(char *trans, char *signs, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, dcomplex *x11,
                 aocl_int64_t *ldx11, dcomplex *x12, aocl_int64_t *ldx12, dcomplex *x21, aocl_int64_t *ldx21,
                 dcomplex *x22, aocl_int64_t *ldx22, double *theta, double *phi, dcomplex *taup1,
                 dcomplex *taup2, dcomplex *tauq1, dcomplex *tauq2, dcomplex *work, aocl_int64_t *lwork,
                 aocl_int64_t *info);
int zunbdb1_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, dcomplex *x11, aocl_int64_t *ldx11, dcomplex *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunbdb2_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, dcomplex *x11, aocl_int64_t *ldx11, dcomplex *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunbdb3_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, dcomplex *x11, aocl_int64_t *ldx11, dcomplex *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunbdb4_check(aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, dcomplex *x11, aocl_int64_t *ldx11, dcomplex *x21,
                  aocl_int64_t *ldx21, double *theta, double *phi, dcomplex *taup1, dcomplex *taup2,
                  dcomplex *tauq1, dcomplex *phantom, dcomplex *work, aocl_int64_t *lwork,
                  aocl_int64_t *info);
int zunbdb5_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, dcomplex *x1, aocl_int64_t *incx1, dcomplex *x2,
                  aocl_int64_t *incx2, dcomplex *q1, aocl_int64_t *ldq1, dcomplex *q2, aocl_int64_t *ldq2,
                  dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunbdb6_check(aocl_int64_t *m1, aocl_int64_t *m2, aocl_int64_t *n, dcomplex *x1, aocl_int64_t *incx1, dcomplex *x2,
                  aocl_int64_t *incx2, dcomplex *q1, aocl_int64_t *ldq1, dcomplex *q2, aocl_int64_t *ldq2,
                  dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zuncsd_check(char *jobu1, char *jobu2, char *jobv1t, char *jobv2t, char *trans, char *signs,
                 aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q, dcomplex *x11, aocl_int64_t *ldx11, dcomplex *x12,
                 aocl_int64_t *ldx12, dcomplex *x21, aocl_int64_t *ldx21, dcomplex *x22, aocl_int64_t *ldx22,
                 double *theta, dcomplex *u1, aocl_int64_t *ldu1, dcomplex *u2, aocl_int64_t *ldu2,
                 dcomplex *v1t, aocl_int64_t *ldv1t, dcomplex *v2t, aocl_int64_t *ldv2t, dcomplex *work,
                 aocl_int64_t *lwork, double *rwork, aocl_int64_t *lrwork, aocl_int64_t *iwork, aocl_int64_t *info);
int zuncsd2by1_check(char *jobu1, char *jobu2, char *jobv1t, aocl_int64_t *m, aocl_int64_t *p, aocl_int64_t *q,
                     dcomplex *x11, aocl_int64_t *ldx11, dcomplex *x21, aocl_int64_t *ldx21, double *theta,
                     dcomplex *u1, aocl_int64_t *ldu1, dcomplex *u2, aocl_int64_t *ldu2, dcomplex *v1t,
                     aocl_int64_t *ldv1t, dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *lrwork,
                     aocl_int64_t *iwork, aocl_int64_t *info);
int zung2l_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *info);
int zung2r_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *info);
int zungbr_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *tau, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunghr_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zungl2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *info);
int zunglq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zungql_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zungqr_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zungr2_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *info);
int zungrq_check(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zungtr_check(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zunm2l_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *info);
int zunm2r_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *info);
int zunmbr_check(char *vect, char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunmhr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunml2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *info);
int zunmlq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zunmql_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zunmqr_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zunmr2_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *info);
int zunmr3_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc,
                 dcomplex *work, aocl_int64_t *info);
int zunmrq_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zunmrz_check(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, aocl_int64_t *l,
                 dcomplex *a, aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc,
                 dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info);
int zunmtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, dcomplex *a,
                 aocl_int64_t *lda, dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work,
                 aocl_int64_t *lwork, aocl_int64_t *info);
int zupgtr_check(char *uplo, aocl_int64_t *n, dcomplex *ap, dcomplex *tau, dcomplex *q, aocl_int64_t *ldq,
                 dcomplex *work, aocl_int64_t *info);
int zupmtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, dcomplex *ap,
                 dcomplex *tau, dcomplex *c__, aocl_int64_t *ldc, dcomplex *work, aocl_int64_t *info);
#endif