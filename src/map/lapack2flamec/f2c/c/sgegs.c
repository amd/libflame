/* ../netlib/sgegs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static real c_b36 = 0.f;
static real c_b37 = 1.f;
/* > \brief <b> SGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for GE mat rices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGEGS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgegs.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgegs.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgegs.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR, */
/* ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, */
/* LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVSL, JOBVSR */
/* INTEGER INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/* $ B( LDB, * ), BETA( * ), VSL( LDVSL, * ), */
/* $ VSR( LDVSR, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine SGGES. */
/* > */
/* > SGEGS computes the eigenvalues, real Schur form, and, optionally, */
/* > left and or/right Schur vectors of a real matrix pair (A,B). */
/* > Given two square matrices A and B, the generalized real Schur */
/* > factorization has the form */
/* > */
/* > A = Q*S*Z**T, B = Q*T*Z**T */
/* > */
/* > where Q and Z are orthogonal matrices, T is upper triangular, and S */
/* > is an upper quasi-triangular matrix with 1-by-1 and 2-by-2 diagonal */
/* > blocks, the 2-by-2 blocks corresponding to scomplex conjugate pairs */
/* > of eigenvalues of (A,B). The columns of Q are the left Schur vectors */
/* > and the columns of Z are the right Schur vectors. */
/* > */
/* > If only the eigenvalues of (A,B) are needed, the driver routine */
/* > SGEGV should be used instead. See SGEGV for a description of the */
/* > eigenvalues of the generalized nonsymmetric eigenvalue problem */
/* > (GNEP). */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVSL */
/* > \verbatim */
/* > JOBVSL is CHARACTER*1 */
/* > = 'N': do not compute the left Schur vectors;
 */
/* > = 'V': compute the left Schur vectors (returned in VSL). */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVSR */
/* > \verbatim */
/* > JOBVSR is CHARACTER*1 */
/* > = 'N': do not compute the right Schur vectors;
 */
/* > = 'V': compute the right Schur vectors (returned in VSR). */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A, B, VSL, and VSR. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, N) */
/* > On entry, the matrix A. */
/* > On exit, the upper quasi-triangular matrix S from the */
/* > generalized real Schur factorization. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB, N) */
/* > On entry, the matrix B. */
/* > On exit, the upper triangular matrix T from the generalized */
/* > real Schur factorization. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* > ALPHAR is REAL array, dimension (N) */
/* > The real parts of each scalar alpha defining an eigenvalue */
/* > of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is REAL array, dimension (N) */
/* > The imaginary parts of each scalar alpha defining an */
/* > eigenvalue of GNEP. If ALPHAI(j) is zero, then the j-th */
/* > eigenvalue is real;
if positive, then the j-th and (j+1)-st */
/* > eigenvalues are a scomplex conjugate pair, with */
/* > ALPHAI(j+1) = -ALPHAI(j). */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is REAL array, dimension (N) */
/* > The scalars beta that define the eigenvalues of GNEP. */
/* > Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and */
/* > beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* > pair (A,B), in one of the forms lambda = alpha/beta or */
/* > mu = beta/alpha. Since either lambda or mu may overflow, */
/* > they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* > VSL is REAL array, dimension (LDVSL,N) */
/* > If JOBVSL = 'V', the matrix of left Schur vectors Q. */
/* > Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* > LDVSL is INTEGER */
/* > The leading dimension of the matrix VSL. LDVSL >=1, and */
/* > if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* > VSR is REAL array, dimension (LDVSR,N) */
/* > If JOBVSR = 'V', the matrix of right Schur vectors Z. */
/* > Not referenced if JOBVSR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSR */
/* > \verbatim */
/* > LDVSR is INTEGER */
/* > The leading dimension of the matrix VSR. LDVSR >= 1, and */
/* > if JOBVSR = 'V', LDVSR >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,4*N). */
/* > For good performance, LWORK must generally be larger. */
/* > To compute the optimal value of LWORK, call ILAENV to get */
/* > blocksizes (for SGEQRF, SORMQR, and SORGQR.) Then compute: */
/* > NB -- MAX of the blocksizes for SGEQRF, SORMQR, and SORGQR */
/* > The optimal LWORK is 2*N + N*(NB+1). */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > = 1,...,N: */
/* > The QZ iteration failed. (A,B) are not in Schur */
/* > form, but ALPHAR(j), ALPHAI(j), and BETA(j) should */
/* > be correct for j=INFO+1,...,N. */
/* > > N: errors that usually indicate LAPACK problems: */
/* > =N+1: error return from SGGBAL */
/* > =N+2: error return from SGEQRF */
/* > =N+3: error return from SORMQR */
/* > =N+4: error return from SORGQR */
/* > =N+5: error return from SGGHRD */
/* > =N+6: error return from SHGEQZ (other than failed */
/* > iteration) */
/* > =N+7: error return from SGGBAK (computing VSL) */
/* > =N+8: error return from SGGBAK (computing VSR) */
/* > =N+9: error return from SLASCL (various places) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realGEeigen */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sgegs_(char *jobvsl, char *jobvsr, aocl_int_t *n, real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *alphar, real *alphai, real *beta, real *vsl, aocl_int_t *ldvsl, real *vsr, aocl_int_t *ldvsr, real *work, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgegs(jobvsl, jobvsr, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t ldvsl_64 = *ldvsl;
    aocl_int64_t ldvsr_64 = *ldvsr;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgegs(jobvsl, jobvsr, &n_64, a, &lda_64, b, &ldb_64, alphar, alphai, beta, vsl, &ldvsl_64, vsr, &ldvsr_64, work, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sgegs(char *jobvsl, char *jobvsr, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *b,
            aocl_int64_t *ldb, real *alphar, real *alphai, real *beta, real *vsl,
            aocl_int64_t *ldvsl, real *vsr, aocl_int64_t *ldvsr, real *work, aocl_int64_t *lwork,
            aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgegs inputs: jobvsl %c ,jobvsr %c ,n %" FLA_IS ",lda %" FLA_IS
                      ",ldb %" FLA_IS ",ldvsl %" FLA_IS ",ldvsr %" FLA_IS ",lwork %" FLA_IS "",
                      *jobvsl, *jobvsr, *n, *lda, *ldb, *ldvsl, *ldvsr, *lwork);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, vsr_dim1, vsr_offset,
        i__1, i__2;
    /* Local variables */
    aocl_int64_t nb, nb1, nb2, nb3, ihi, ilo;
    real eps, anrm, bnrm;
    aocl_int64_t itau, lopt;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t ileft, iinfo, icols;
    logical ilvsl;
    aocl_int64_t iwork;
    logical ilvsr;
    aocl_int64_t irows;
    logical ilascl, ilbscl;
    real safmin;
    real bignum;
    aocl_int64_t ijobvl, iright;
    aocl_int64_t ijobvr;
    real anrmto;
    aocl_int64_t lwkmin;
    real bnrmto;
    real smlnum;
    aocl_int64_t lwkopt;
    logical lquery;
    /* -- LAPACK driver routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    vsl_dim1 = *ldvsl;
    vsl_offset = 1 + vsl_dim1;
    vsl -= vsl_offset;
    vsr_dim1 = *ldvsr;
    vsr_offset = 1 + vsr_dim1;
    vsr -= vsr_offset;
    --work;
    /* Function Body */
    if(lsame_(jobvsl, "N", 1, 1))
    {
        ijobvl = 1;
        ilvsl = FALSE_;
    }
    else if(lsame_(jobvsl, "V", 1, 1))
    {
        ijobvl = 2;
        ilvsl = TRUE_;
    }
    else
    {
        ijobvl = -1;
        ilvsl = FALSE_;
    }
    if(lsame_(jobvsr, "N", 1, 1))
    {
        ijobvr = 1;
        ilvsr = FALSE_;
    }
    else if(lsame_(jobvsr, "V", 1, 1))
    {
        ijobvr = 2;
        ilvsr = TRUE_;
    }
    else
    {
        ijobvr = -1;
        ilvsr = FALSE_;
    }
    /* Test the input arguments */
    /* Computing MAX */
    i__1 = *n << 2;
    lwkmin = fla_max(i__1, 1);
    lwkopt = lwkmin;
    work[1] = (real)lwkopt;
    lquery = *lwork == -1;
    *info = 0;
    if(ijobvl <= 0)
    {
        *info = -1;
    }
    else if(ijobvr <= 0)
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldvsl < 1 || ilvsl && *ldvsl < *n)
    {
        *info = -12;
    }
    else if(*ldvsr < 1 || ilvsr && *ldvsr < *n)
    {
        *info = -14;
    }
    else if(*lwork < lwkmin && !lquery)
    {
        *info = -16;
    }
    if(*info == 0)
    {
        nb1 = aocl_lapack_ilaenv(&c__1, "SGEQRF", " ", n, n, &c_n1, &c_n1);
        nb2 = aocl_lapack_ilaenv(&c__1, "SORMQR", " ", n, n, n, &c_n1);
        nb3 = aocl_lapack_ilaenv(&c__1, "SORGQR", " ", n, n, n, &c_n1);
        /* Computing MAX */
        i__1 = fla_max(nb1, nb2);
        nb = fla_max(i__1, nb3);
        lopt = (*n << 1) + *n * (nb + 1);
        work[1] = (real)lopt;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SGEGS ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants */
    eps = slamch_("E") * slamch_("B");
    safmin = slamch_("S");
    smlnum = *n * safmin / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = aocl_lapack_slange("M", n, n, &a[a_offset], lda, &work[1]);
    ilascl = FALSE_;
    if(anrm > 0.f && anrm < smlnum)
    {
        anrmto = smlnum;
        ilascl = TRUE_;
    }
    else if(anrm > bignum)
    {
        anrmto = bignum;
        ilascl = TRUE_;
    }
    if(ilascl)
    {
        aocl_lapack_slascl("G", &c_n1, &c_n1, &anrm, &anrmto, n, n, &a[a_offset], lda, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 9;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    /* Scale B if max element outside range [SMLNUM,BIGNUM] */
    bnrm = aocl_lapack_slange("M", n, n, &b[b_offset], ldb, &work[1]);
    ilbscl = FALSE_;
    if(bnrm > 0.f && bnrm < smlnum)
    {
        bnrmto = smlnum;
        ilbscl = TRUE_;
    }
    else if(bnrm > bignum)
    {
        bnrmto = bignum;
        ilbscl = TRUE_;
    }
    if(ilbscl)
    {
        aocl_lapack_slascl("G", &c_n1, &c_n1, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 9;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    /* Permute the matrix to make it more nearly triangular */
    /* Workspace layout: (2*N words -- "work..." not actually used) */
    /* left_permutation, right_permutation, work... */
    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;
    aocl_lapack_sggbal("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[ileft],
                       &work[iright], &work[iwork], &iinfo);
    if(iinfo != 0)
    {
        *info = *n + 1;
        goto L10;
    }
    /* Reduce B to triangular form, and initialize VSL and/or VSR */
    /* Workspace layout: ("work..." must have at least N words) */
    /* left_permutation, right_permutation, tau, work... */
    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    aocl_lapack_sgeqrf(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[iwork],
                       &i__1, &iinfo);
    if(iinfo >= 0)
    {
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer)work[iwork] + iwork - 1; // , expr subst
        lwkopt = fla_max(i__1, i__2);
    }
    if(iinfo != 0)
    {
        *info = *n + 2;
        goto L10;
    }
    i__1 = *lwork + 1 - iwork;
    aocl_lapack_sormqr("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &work[itau],
                       &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &iinfo);
    if(iinfo >= 0)
    {
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer)work[iwork] + iwork - 1; // , expr subst
        lwkopt = fla_max(i__1, i__2);
    }
    if(iinfo != 0)
    {
        *info = *n + 3;
        goto L10;
    }
    if(ilvsl)
    {
        aocl_lapack_slaset("Full", n, n, &c_b36, &c_b37, &vsl[vsl_offset], ldvsl);
        i__1 = irows - 1;
        i__2 = irows - 1;
        aocl_lapack_slacpy("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb,
                           &vsl[ilo + 1 + ilo * vsl_dim1], ldvsl);
        i__1 = *lwork + 1 - iwork;
        aocl_lapack_sorgqr(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &work[itau],
                           &work[iwork], &i__1, &iinfo);
        if(iinfo >= 0)
        {
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = (integer)work[iwork] + iwork - 1; // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        if(iinfo != 0)
        {
            *info = *n + 4;
            goto L10;
        }
    }
    if(ilvsr)
    {
        aocl_lapack_slaset("Full", n, n, &c_b36, &c_b37, &vsr[vsr_offset], ldvsr);
    }
    /* Reduce to generalized Hessenberg form */
    aocl_lapack_sgghrd(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb,
                       &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &iinfo);
    if(iinfo != 0)
    {
        *info = *n + 5;
        goto L10;
    }
    /* Perform QZ algorithm, computing Schur vectors if desired */
    /* Workspace layout: ("work..." must have at least 1 word) */
    /* left_permutation, right_permutation, work... */
    iwork = itau;
    i__1 = *lwork + 1 - iwork;
    aocl_lapack_shgeqz("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb,
                       &alphar[1], &alphai[1], &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset],
                       ldvsr, &work[iwork], &i__1, &iinfo);
    if(iinfo >= 0)
    {
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer)work[iwork] + iwork - 1; // , expr subst
        lwkopt = fla_max(i__1, i__2);
    }
    if(iinfo != 0)
    {
        if(iinfo > 0 && iinfo <= *n)
        {
            *info = iinfo;
        }
        else if(iinfo > *n && iinfo <= *n << 1)
        {
            *info = iinfo - *n;
        }
        else
        {
            *info = *n + 6;
        }
        goto L10;
    }
    /* Apply permutation to VSL and VSR */
    if(ilvsl)
    {
        aocl_lapack_sggbak("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n,
                           &vsl[vsl_offset], ldvsl, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 7;
            goto L10;
        }
    }
    if(ilvsr)
    {
        aocl_lapack_sggbak("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n,
                           &vsr[vsr_offset], ldvsr, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 8;
            goto L10;
        }
    }
    /* Undo scaling */
    if(ilascl)
    {
        aocl_lapack_slascl("H", &c_n1, &c_n1, &anrmto, &anrm, n, n, &a[a_offset], lda, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 9;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        aocl_lapack_slascl("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphar[1], n, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 9;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        aocl_lapack_slascl("G", &c_n1, &c_n1, &anrmto, &anrm, n, &c__1, &alphai[1], n, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 9;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    if(ilbscl)
    {
        aocl_lapack_slascl("U", &c_n1, &c_n1, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 9;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        aocl_lapack_slascl("G", &c_n1, &c_n1, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 9;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
L10:
    work[1] = (real)lwkopt;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGEGS */
}
/* sgegs_ */
