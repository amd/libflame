/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_common.h
 *  @brief Defines function declarations of validate functions of APIs to be
 *         used in test suite.
 *  */

#ifndef VALIDATE_COMMON_H
#define VALIDATE_COMMON_H

void validate_geqrf(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    void *T_test, integer datatype, double err_thresh);

void validate_gerq2(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    void *T_test, integer datatype, double err_thresh);

void validate_gerqf(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    void *T_test, integer datatype, double err_thresh);

void validate_gelqf(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    void *T_test, integer datatype, double err_thresh);

void validate_gesdd(char *tst_api, char *jobz, integer m, integer n, void *A, void *A_test,
                    integer lda, void *s, void *s_in, void *U, integer ldu, void *V, integer ldvt,
                    integer datatype, double err_thresh, char imatrix, void *scal);

void validate_gesvd(char *tst_api, char *jobu, char *jobvt, integer m, integer n, void *A,
                    void *A_test, integer lda, void *s, void *s_test, void *U, integer ldu, void *V,
                    integer ldvt, integer datatype, double err_thresh, FILE *g_ext_fptr,
                    char imatrix, void *scal);

void validate_getrf(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    integer *IPIV, integer datatype, double err_thresh, char imatrix);
void validate_getrf_internal(integer m_A, integer n_A, void *A, void *A_test, integer lda,
                             integer *IPIV, integer datatype, char imatrix, double *resid1,
                             double *resid2);

void validate_getri(char *tst_api, integer m_A, integer n_A, void *A, void *A_inv, integer lda,
                    integer *IPIV, integer datatype, double err_thresh, char imatrix);

void validate_getrs(char *tst_api, char *trans, integer n, integer nrhs, void *A, integer lda,
                    void *B, integer ldb, void *X, integer datatype, double err_thresh,
                    char imatrix, void *scal);

void validate_orgqr(char *tst_api, integer m, integer n, void *A, integer lda, void *Q, void *R,
                    integer datatype, double err_thresh, char imatrix);

void validate_potrf(char *tst_api, char *uplo, integer m, void *A, void *A_test, integer lda,
                    integer datatype, double err_thresh);

void validate_potrs(char *tst_api, integer n, integer nrhs, void *A, integer lda, void *X, void *B,
                    integer ldb, integer datatype, double err_thresh, char imatrix);

void validate_syevd(char *tst_api, char *jobz, integer n, void *A, void *A_test, integer lda,
                    void *w, integer datatype, double err_thresh);

void validate_geevx(char *tst_api, char *jobvl, char *jobvr, char *sense, char *balanc, integer m,
                    void *A, void *A_test, integer lda, void *VL, integer ldvl, void *VR,
                    integer ldvr, void *w, void *wr, void *wi, void *scale, void *abnrm,
                    void *rconde, void *rcondv, integer datatype, char imatrix, void *scal,
                    double err_thresh, void *wr_in, void *wi_in);

void validate_geev(char *tst_api, char *jobvl, char *jobvr, integer m, void *A, void *A_test,
                   integer lda, void *VL, integer ldvt, void *VR, integer ldvr, void *w, void *wr,
                   void *wi, integer datatype, char imatrix, void *scal, double err_thresh,
                   void *wr_in, void *wi_in);

void validate_geqp3(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    integer *jpvt, void *T_test, integer datatype, double err_thresh, char imatrix);

void validate_gesv(char *tst_api, integer n, integer nrhs, void *A, integer lda, void *B,
                   integer ldb, void *X, integer datatype, double err_thresh, char imatrix,
                   void *scal);

void validate_ggev(char *tst_api, char *jobvl, char *jobvr, integer n, void *A, integer lda,
                   void *B, integer ldb, void *alpha, void *alphar, void *alphai, void *beta,
                   void *VL, integer ldvl, void *VR, integer ldvr, integer datatype,
                   double err_thresh);

void validate_ggevx(char *tst_api, char *balanc, char *jobvl, char *jobvr, char *sense, integer n,
                    void *A, integer lda, void *B, integer ldb, void *alpha, void *alphar,
                    void *alphai, void *beta, void *VL, integer ldvl, void *VR, integer ldvr,
                    integer datatype, double err_thresh);

/* This function will validate STEDC() output eigenvectors and orthogonal
   matrices only if compz != N, as output will not be generated
   if compz = N.*/
void validate_stedc(char *tst_api, char compz, integer n, void *D_test, void *Z_input, void *Z,
                    integer ldz, integer datatype, double err_thresh);
void validate_eigen_values(char *tst_api, integer datatype, integer n, void *alpha, void *alphar,
                           void *alphai, void *beta, void *alphae, void *alphaer, void *alphaei,
                           void *betae, double *resid1, double *resid2, double *resid3);
void validate_hgeqz_comp_n(char *tst_api, integer datatype, integer n, void *H_test, void *H_ntest,
                           integer ldh, void *T_test, void *T_ntest, integer ldt, void *alpha,
                           void *alphar, void *alphai, void *beta, void *alphan, void *alphanr,
                           void *alphani, void *betan, double err_thresh);
void validate_hgeqz_eigen_values(char *tst_api, integer datatype, integer n, void *alpha,
                                 void *alphar, void *alphai, void *beta, void *alphae,
                                 void *alphaer, void *alphaei, void *betae, double err_thresh);
void validate_hgeqz(char *tst_api, char *job, char *compq, char *compz, integer n, void *H,
                    void *H_test, void *A, integer ldh, void *T, void *T_test, void *B, integer ldt,
                    void *Q, void *Q_test, void *Q_A, integer ldq, void *Z, void *Z_test, void *Z_A,
                    integer ldz, integer datatype, double err_thresh, char imatrix);

void validate_hseqr(char *tst_api, char *job, char *compz, integer n, void *H, void *H_test,
                    integer ldh, void *Z, void *Z_test, integer ldz, void *wr, void *wr_in,
                    void *wi, void *wi_in, void *w, integer datatype, double err_thresh,
                    integer *ilo, integer *ihi, char imatrix, void *scal_H);

void validate_spffrt2(char *tst_api, integer n, integer ncolm, void *A, void *AP, integer datatype,
                      double err_thresh);

void validate_spffrtx(char *tst_api, integer n, integer ncolm, void *A, void *AP, integer datatype,
                      double err_thresh);

void validate_gghrd_int(char *tst_api, char *compq, char *compz, integer n, void *A, void *A_test,
                        integer lda, void *B, void *B_test, integer ldb, void *Q, void *Q_test,
                        integer ldq, void *Z, void *Z_test, integer ldz, integer datatype,
                        double *resid);
void validate_gghrd(char *tst_api, char *compq, char *compz, integer n, void *A, void *A_test,
                    integer lda, void *B, void *B_test, integer ldb, void *Q, void *Q_test,
                    integer ldq, void *Z, void *Z_test, integer ldz, integer datatype,
                    double err_thresh);

void validate_gehrd(char *tst_api, integer n, integer ilo, integer ihi, void *A, void *A_test,
                    integer lda, void *tau, integer datatype, double err_thresh);

void validate_rot(char *tst_api, integer datatype, integer n, void *cx, void *cx_test, integer incx,
                  void *cy, void *cy_test, integer incy, void *c, void *s, double err_thresh);

void validate_lartg(char *tst_api, integer datatype, void *f, void *g, void *r, void *c, void *s,
                    double err_thresh);

void validate_gels(char *tst_api, char *trans, integer m, integer n, integer nrhs, void *A,
                   integer lda, void *B, integer ldb, void *x, integer datatype, double err_thresh,
                   char imatrix);

void validate_larfg(char *tst_api, integer datatype, integer n, integer incx, integer x_length,
                    void *x, void *v, void *tau, double err_thresh);
void validate_larf(char *tst_api, integer datatype, char side, integer m, integer n, void *v,
                   integer incv, void *c__, integer ldc__, void *c__out, integer ldc__out,
                   void *tau, double err_thresh);

void validate_gtsv(char *tst_api, integer datatype, integer n, integer nrhs, void *B, integer ldb,
                   void *X, void *Xact, integer ldx, void *dl, void *d, void *du, void *dl_save,
                   void *d_save, void *du_save, void *scal, char imatrix, double err_thresh);

void validate_syev(char *tst_api, char *jobz, char *range, integer n, void *A, void *A_test,
                   integer lda, integer il, integer iu, void *L, void *lambda, void *ifail,
                   integer datatype, double err_thresh, char imatrix, void *scal);

void validate_gesvdx(char *tst_api, char *jobu, char *jobvt, char range, integer m, integer n,
                     void *A, void *A_test, integer lda, void *vl, void *vu, integer il, integer iu,
                     integer ns, void *s, void *s_test, void *U, integer ldu, void *V, integer ldvt,
                     integer datatype, double err_thresh, FILE *g_ext_fptr, void *scal,
                     char imatrix);
/* This function validates LU factorization output by reconstructing into input band storage
   matrix.*/
void validate_gbtrf(char *tst_api, integer m_A, integer n_A, integer kl, integer ku, void *AB,
                    void *AB_test, integer ldab, integer *IPIV, integer datatype,
                    double err_thresh);

void validate_gelsd(char *tst_api, integer m, integer n, integer NRHS, void *A, integer lda,
                    void *B, integer ldb, void *S, void *X, void *rcond, integer *rank,
                    integer datatype, double err_thresh, char imatrix);

void validate_sytrf(char *tst_api, char *uplo, integer n, integer lda, void *A_res,
                    integer datatype, integer *ipiv, double err_thresh, void *A);

void validate_hetrf(char *tst_api, char *uplo, integer n, integer lda, void *A_res,
                    integer datatype, integer *ipiv, double err_thresh, void *A);
void validate_hetrf_rook(char *tst_api, char *uplo, integer n, integer lda, void *A_res,
                         integer datatype, integer *ipiv, double err_thresh, void *A);

/* GGEV API validate case for JOBVL= JOBVR = N */
void validate_ggev_EVs(char *tst_api, integer m, void *alpha, void *alphar, void *alphai,
                       void *beta, void *alpha_copy, void *alphar_copy, void *alphai_copy,
                       void *beta_copy, integer datatype, double err_thresh);

void validate_sygvd(char *tst_api, integer itype, char *jobz, char *range, char *uplo, integer n,
                    void *A, void *A_test, integer lda, void *B, void *B_test, integer ldb,
                    integer il, integer iu, void *lambda_orig, void *lambda_out, void *ifail,
                    integer datatype, double err_thresh, char imatrix, void *scal);

/* Validate function for lange */
void validate_lange(char *tst_api, integer datatype, char norm_type, integer m, integer n,
                    integer lda, void *A, void *result, double err_thresh);

void validate_gecon(char *tst_api, integer datatype, char norm, integer n, void *A, void *A_save,
                    integer lda, double err_thresh, char imatrix_char);

void validate_getrfnpi(char *tst_api, integer m_A, integer n_A, integer nfact, void *A,
                       void *A_test, integer lda, integer *IPIV, integer datatype,
                       double err_thresh, char imatrix);
void validate_hetri_rook(char *tst_api, char uplo, integer n, void *A, void *A_inv, integer lda,
                         integer *ipiv, integer datatype, double err_thresh, char imatrix);
void validate_ormqr(char *tst_api, char side, char trans, integer m, integer n, integer k, void *A,
                    integer lda, void *C, void *Tau, integer ldc, void *C_test,
                    integer datatype, double err_thresh, char imatrix);

#endif // VALIDATE_COMMON_H
