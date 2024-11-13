/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#ifdef __cplusplus
extern "C" {
#endif
void invoke_cpp_gbtrf(integer datatype, integer *m, integer *n, integer *kl, integer *ku, void *ab,
                      integer *ldab, integer *ipiv, integer *info);
void invoke_cpp_gbtrs(integer datatype, char *trans, integer *n, integer *kl, integer *ku,
                      integer *nrhs, void *ab, integer *ldab, integer *ipiv, void *b, integer *ldb,
                      integer *info);
void invoke_cpp_geev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                    void *wr, void *wi, void *w, void *vl, integer *ldvl, void *vr, integer *ldvr,
                    void *work, integer *lwork, void *rwork, integer *info);
void invoke_cpp_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                  void *a, integer *lda, void *wr, void *wi, void *w, void *vl, integer *ldvl,
                  void *vr, integer *ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                  void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                  integer *iwork, integer *info);
void invoke_cpp_gehrd(integer datatype, integer *n, integer *ilo, integer *ihi, void *A, integer *lda,
                  void *tau, void *work, integer *lwork, integer *info);
void invoke_cpp_gelqf(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                  void *work, integer *lwork, integer *info);
void invoke_cpp_gels(integer datatype, char *trans, integer *m, integer *n, integer *nrhs, void *A,
                     integer *lda, void *B, integer *ldb, void *work, integer *lwork, integer *info);
void invoke_cpp_gelsd(integer datatype, integer *m, integer *n, integer *nrhs, void *a, integer *lda,
                      void *b, integer *ldb, void *s, void *rcond, integer *rank, void *work,
                      integer *lwork, void *rwork, integer *iwork, integer *info);
void invoke_cpp_gelss(integer datatype, integer *m, integer *n, integer *nrhs, void *A, integer *lda,
                      void *B, integer *ldb, void *s, void *rcond, integer *rank, void *work,
                      integer *lwork, void *rwork, integer *info);
void invoke_cpp_geqp3(integer datatype, integer *m, integer *n, void *a, integer *lda, integer* jpvt, void *tau,
                      void *work, integer *lwork, void *rwork, integer *info);
void invoke_cpp_geqrf(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                    void *work, integer *lwork, integer *info);
void invoke_cpp_gesdd(integer datatype, char *jobz, integer *m, integer *n, void *a, integer *lda,
                    void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                    integer *lwork, void *rwork, integer *iwork, integer *info);
void invoke_cpp_gesv(integer datatype, integer *n, integer *nrhs, void *a, integer *lda, integer *ipiv,
                     void *b, integer *ldb, integer *info);
void invoke_cpp_gesvd(integer datatype, char *jobu, char *jobvt, integer *m, integer *n, void *a,
                      integer *lda, void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                      integer *lwork, void *rwork, integer *info);
void invoke_cpp_gesvdx(integer datatype, char *jobu, char *jobvt, char *range, integer *m, integer *n,
                    void *a, integer *lda, void *vl, void *vu, integer *il, integer *iu, integer *ns,
                    void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                    integer *lwork, integer *iwork, void *rwork, integer *info);
void invoke_cpp_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv,
                    integer *info);
void invoke_cpp_getri(integer datatype, integer *n, void *a, integer *lda, integer *ipiv, void *work,
                      integer *lwork, integer *info);
void invoke_cpp_getrs(integer datatype, char *trans, integer *n, integer *nrhs, void *a, integer *lda,
                      integer *ipiv, void *b, integer *ldb, integer *info);
void invoke_cpp_ggev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                 void *b, integer *ldb, void *alpha, void *alphar, void *alphai, void *beta,
                 void *vl, integer *ldvl, void *vr, integer *ldvr, void *work, integer *lwork,
                 void *rwork, integer *info);
void invoke_cpp_ggevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                      void *a, integer *lda, void *b, integer *ldb, void *alpha, void *alphar,
                      void *alphai, void *beta, void *vl, integer *ldvl, void *vr, integer *ldvr,
                      integer *ilo, integer *ihi, void *lscale, void *rscale, void *abnrm, void *bbnrm,
                      void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                      integer *iwork, integer *bwork, integer *info);
void invoke_cpp_gghrd(integer datatype, char *compq, char *compz, integer *n, integer *ilo,
                      integer *ihi, void *a, integer *lda, void *b, integer *ldb, void *q, integer *ldq,
                      void *z, integer *ldz, integer *info);
void invoke_cpp_gtsv(integer datatype, integer *n, integer *nrhs, void *dl, void *d, void *du, void *b,
                     integer *ldb, integer *info);
void invoke_cpp_hetrf_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                           integer *ipiv, void *work, integer *lwork, integer *info);
void invoke_cpp_hgeqz(integer datatype, char *job, char *compq, char *compz, integer *n, integer *ilo,
                      integer *ihi, void *h, integer *ldh, void *t, integer *ldt, void *alpha,
                      void *alphar, void *alphai, void *beta, void *q, integer *ldq, void *z,
                      integer *ldz, void *work, integer *lwork, void *rwork, integer *info);
void invoke_cpp_hseqr(integer datatype, char *job, char *compz, integer *n, integer *ilo, integer *ihi,
                      void *h, integer *ldh, void *w, void *wr, void *wi, void *z, integer *ldz,
                      void *work, integer *lwork, integer *info);
void invoke_cpp_lartg(integer datatype, void *f, void *g, void *c, void *s, void *r);
#ifdef __cplusplus
}
#endif