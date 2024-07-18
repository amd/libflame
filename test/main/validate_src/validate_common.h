/******************************************************************************
 * Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_common.h
 *  @brief Defines function declarations of validate functions of APIs to be
 *         used in test suite.
 *  */

#ifndef VALIDATE_COMMON_H
#define VALIDATE_COMMON_H

void validate_geqrf(integer m_A, integer n_A, void *A, void *A_test, integer lda, void *T_test,
                    integer datatype, double *residual, integer *info);

void validate_gerq2(integer m_A, integer n_A, void *A, void *A_test, integer lda, void *T_test,
                    integer datatype, double *residual, integer *info);

void validate_gerqf(integer m_A, integer n_A, void *A, void *A_test, integer lda, void *T_test,
                    integer datatype, double *residual, integer *info);

void validate_gelqf(integer m_A, integer n_A, void *A, void *A_test, integer lda, void *T_test,
                    integer datatype, double *residual, integer *info);

void validate_gesdd(char *jobz, integer m, integer n, void *A, void *A_test, integer lda, void *s,
                    void *U, integer ldu, void *V, integer ldvt, integer datatype, double *residual,
                    integer *info);

void validate_gesvd(char *jobu, char *jobvt, integer m, integer n, void *A, void *A_test,
                    integer lda, void *s, void *s_test, void *U, integer ldu, void *V, integer ldvt,
                    integer datatype, double *residual, integer *info, FILE *g_ext_fptr,
                    char imatrix, void *scal);

void validate_getrf(integer m_A, integer n_A, void *A, void *A_test, integer lda, integer *IPIV,
                    integer datatype, double *residual, integer *info, char imatrix);

void validate_getri(integer m_A, integer n_A, void *A, void *A_inv, integer lda, integer *IPIV,
                    integer datatype, double *residual, integer *info);

void validate_getrs(char *trans, integer n, integer nrhs, void *A, integer lda, void *B,
                    integer ldb, void *X, integer datatype, double *residual, integer *info);

void validate_orgqr(integer m, integer n, void *A, integer lda, void *Q, void *R, void *work,
                    integer datatype, double *residual, integer *info);

void validate_potrf(char *uplo, integer m, void *A, void *A_test, integer lda, integer datatype,
                    double *residual, integer *info);

void validate_potrs(integer n, integer nrhs, void *A, integer lda, void *X, void *B, integer ldb,
                    integer datatype, double *residual, integer *info);

void validate_syevd(char *jobz, integer n, void *A, void *A_test, integer lda, void *w,
                    integer datatype, double *residual, integer *info);

void validate_geevx(char *jobvl, char *jobvr, char *sense, char *balanc, integer m, void *A,
                    void *A_test, integer lda, void *VL, integer ldvl, void *VR, integer ldvr,
                    void *w, void *wr, void *wi, void *scale, void *abnrm, void *rconde,
                    void *rcondv, integer datatype, double *residual, integer *info, void *wr_in,
                    void *wi_in);

void validate_geev(char *jobvl, char *jobvr, integer m, void *A, void *A_test, integer lda,
                   void *VL, integer ldvt, void *VR, integer ldvr, void *w, void *wr, void *wi,
                   integer datatype, char imatrix, void *scal, double *residual, integer *info,
                   void *wr_in, void *wi_in);

void validate_geqp3(integer m_A, integer n_A, void *A, void *A_test, integer lda, integer *jpvt,
                    void *T_test, integer datatype, double *residual, integer *info);

void validate_gesv(integer n, integer nrhs, void *A, integer lda, void *B, integer ldb, void *X,
                   integer datatype, double *residual, integer *info);

void validate_ggev(char *jobvl, char *jobvr, integer n, void *A, integer lda, void *B, integer ldb,
                   void *alpha, void *alphar, void *alphai, void *beta, void *VL, integer ldvl,
                   void *VR, integer ldvr, integer datatype, double *residual, integer *info);

void validate_ggevx(char *balanc, char *jobvl, char *jobvr, char *sense, integer n, void *A,
                    integer lda, void *B, integer ldb, void *alpha, void *alphar, void *alphai,
                    void *beta, void *VL, integer ldvl, void *VR, integer ldvr, integer datatype,
                    double *residual, integer *info);

/* This function will validate STEDC() output eigenvectors and orthogonal
   matrices only if compz != N, as output will not be generated
   if compz = N.*/
void validate_stedc(char compz, integer n, void *D_test, void *Z_input, void *Z, integer ldz,
                    integer datatype, double *residual, integer *info);

void validate_hgeqz(char *job, char *compq, char *compz, integer n, void *H, void *H_test, void *A,
                    integer ldh, void *T, void *T_test, void *B, integer ldt, void *Q, void *Q_test,
                    void *Q_A, integer ldq, void *Z, void *Z_test, void *Z_A, integer ldz,
                    integer datatype, double *residual, integer *info);

void validate_hseqr(char *job, char *compz, integer n, void *H, void *H_test, integer ldh, void *Z,
                    void *Z_test, integer ldz, void *wr, void *wr_in, void *wi, void *wi_in,
                    void *w, integer datatype, double *residual, integer *info, integer *ilo,
                    integer *ihi);

void validate_spffrt2(integer n, integer ncolm, void *A, void *AP, integer datatype,
                      double *residual);

void validate_spffrtx(integer n, integer ncolm, void *A, void *AP, integer datatype,
                      double *residual);

void validate_gghrd(char *compq, char *compz, integer n, void *A, void *A_test, integer lda,
                    void *B, void *B_test, integer ldb, void *Q, void *Q_test, integer ldq, void *Z,
                    void *Z_test, integer ldz, integer datatype, double *residual, integer *info);

void validate_gehrd(integer n, integer ilo, integer ihi, void *A, void *A_test, integer lda,
                    void *tau, integer datatype, double *residual, integer *info);

void validate_rot(integer datatype, integer n, void *cx, void *cx_test, integer incx, void *cy,
                  void *cy_test, integer incy, void *c, void *s, double *residual);

void validate_lartg(integer datatype, void *f, void *g, void *r, void *c, void *s,
                    double *residual);

void validate_gels(char *trans, integer m, integer n, integer nrhs, void *A, integer lda, void *B,
                   integer ldb, void *x, integer datatype, double *residual, integer *info,
                   char imatrix);

void validate_larfg(integer datatype, integer n, integer incx, integer x_length, void *x, void *v,
                    void *tau, double *residual);
void validate_larf(integer datatype, char side, integer m, integer n, void *v, integer incv,
                   void *c__, integer ldc__, void *c__out, integer ldc__out, void *tau,
                   double *residual);

void validate_gtsv(integer datatype, integer n, integer nrhs, void *B, integer ldb, void *X,
                   void *Xact, integer ldx, void *dl, void *d, void *du, void *dl_save,
                   void *d_save, void *du_save, integer info, void *scal, char imatrix,
                   double *residual);

void validate_syev(char *jobz, char *range, integer n, void *A, void *A_test, integer lda,
                   integer il, integer iu, void *L, void *lambda, void *ifail, integer datatype,
                   double *residual, char imatrix, void *scal);

void validate_gesvdx(char *jobu, char *jobvt, char range, integer m, integer n, void *A,
                     void *A_test, integer lda, void *vl, void *vu, integer il, integer iu,
                     integer ns, void *s, void *s_test, void *U, integer ldu, void *V, integer ldvt,
                     integer datatype, double *residual, integer *info, FILE *g_ext_fptr,
                     void *scal, char imatrix);
/* This function validates LU factorization output by reconstructing into input band storage
   matrix.*/
void validate_gbtrf(integer m_A, integer n_A, integer kl, integer ku, void *AB, void *AB_test,
                    integer ldab, integer *IPIV, integer datatype, double *residual, integer *info);

void validate_gelsd(integer m, integer n, integer NRHS, void *A, integer lda, void *B, integer ldb,
                    void *S, void *X, void *rcond, integer *rank, integer datatype,
                    double *residual, char imatrix);

void validate_sytrf(char *uplo, integer n, integer lda, void *A_res, integer datatype,
                    integer *ipiv, double *residual, integer *info, void *A);

/* GGEV API validate case for JOBVL= JOBVR = N */
void validate_ggev_EVs(integer m, void *alpha, void *alphar, void *alphai, void *beta,
                       void *alpha_copy, void *alphar_copy, void *alphai_copy, void *beta_copy,
                       integer datatype, double *residual);
#endif // VALIDATE_COMMON_H
