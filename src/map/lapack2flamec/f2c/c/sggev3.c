/* ./sggev3.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static real c_b36 = 0.f;
static real c_b37 = 1.f;
/* > \brief <b> SGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for GE mat rices (blocked algorithm)</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGGEV3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggev3.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggev3.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggev3.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, */
/* $ ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, */
/* $ INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVL, JOBVR */
/* INTEGER INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/* $ B( LDB, * ), BETA( * ), VL( LDVL, * ), */
/* $ VR( LDVR, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGEV3 computes for a pair of N-by-N real nonsymmetric matrices (A,B) */
/* > the generalized eigenvalues, and optionally, the left and/or right */
/* > generalized eigenvectors. */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar */
/* > lambda or a ratio alpha/beta = lambda, such that A - lambda*B is */
/* > singular. It is usually represented as the pair (alpha,beta), as */
/* > there is a reasonable interpretation for beta=0, and even for both */
/* > being zero. */
/* > */
/* > The right eigenvector v(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* > */
/* > A * v(j) = lambda(j) * B * v(j). */
/* > */
/* > The left eigenvector u(j) corresponding to the eigenvalue lambda(j) */
/* > of (A,B) satisfies */
/* > */
/* > u(j)**H * A = lambda(j) * u(j)**H * B . */
/* > */
/* > where u(j)**H is the conjugate-transpose of u(j). */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVL */
/* > \verbatim */
/* > JOBVL is CHARACTER*1 */
/* > = 'N': do not compute the left generalized eigenvectors;
 */
/* > = 'V': compute the left generalized eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* > JOBVR is CHARACTER*1 */
/* > = 'N': do not compute the right generalized eigenvectors;
 */
/* > = 'V': compute the right generalized eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A, B, VL, and VR. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA, N) */
/* > On entry, the matrix A in the pair (A,B). */
/* > On exit, A has been overwritten. */
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
/* > On entry, the matrix B in the pair (A,B). */
/* > On exit, B has been overwritten. */
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
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* > ALPHAI is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is REAL array, dimension (N) */
/* > On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will */
/* > be the generalized eigenvalues. If ALPHAI(j) is zero, then */
/* > the j-th eigenvalue is real;
if positive, then the j-th and */
/* > (j+1)-st eigenvalues are a complex conjugate pair, with */
/* > ALPHAI(j+1) negative. */
/* > */
/* > Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j) */
/* > may easily over- or underflow, and BETA(j) may even be zero. */
/* > Thus, the user should avoid naively computing the ratio */
/* > alpha/beta. However, ALPHAR and ALPHAI will be always less */
/* > than and usually comparable with norm(A) in magnitude, and */
/* > BETA always less than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* > VL is REAL array, dimension (LDVL,N) */
/* > If JOBVL = 'V', the left eigenvectors u(j) are stored one */
/* > after another in the columns of VL, in the same order as */
/* > their eigenvalues. If the j-th eigenvalue is real, then */
/* > u(j) = VL(:,j), the j-th column of VL. If the j-th and */
/* > (j+1)-th eigenvalues form a complex conjugate pair, then */
/* > u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1). */
/* > Each eigenvector is scaled so the largest component has */
/* > f2c_abs(real part)+f2c_abs(imag. part)=1. */
/* > Not referenced if JOBVL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVL */
/* > \verbatim */
/* > LDVL is INTEGER */
/* > The leading dimension of the matrix VL. LDVL >= 1, and */
/* > if JOBVL = 'V', LDVL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VR */
/* > \verbatim */
/* > VR is REAL array, dimension (LDVR,N) */
/* > If JOBVR = 'V', the right eigenvectors v(j) are stored one */
/* > after another in the columns of VR, in the same order as */
/* > their eigenvalues. If the j-th eigenvalue is real, then */
/* > v(j) = VR(:,j), the j-th column of VR. If the j-th and */
/* > (j+1)-th eigenvalues form a complex conjugate pair, then */
/* > v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1). */
/* > Each eigenvector is scaled so the largest component has */
/* > f2c_abs(real part)+f2c_abs(imag. part)=1. */
/* > Not referenced if JOBVR = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVR */
/* > \verbatim */
/* > LDVR is INTEGER */
/* > The leading dimension of the matrix VR. LDVR >= 1, and */
/* > if JOBVR = 'V', LDVR >= N. */
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
/* > The QZ iteration failed. No eigenvectors have been */
/* > calculated, but ALPHAR(j), ALPHAI(j), and BETA(j) */
/* > should be correct for j=INFO+1,...,N. */
/* > > N: =N+1: other than QZ iteration failed in SLAQZ0. */
/* > =N+2: error return from STGEVC. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup ggev3 */
/* ===================================================================== */
/* Subroutine */
void sggev3_(char *jobvl, char *jobvr, integer *n, real *a, integer *lda, real *b, integer *ldb,
             real *alphar, real *alphai, real *beta, real *vl, integer *ldvl, real *vr,
             integer *ldvr, real *work, integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sggev3 inputs: jobvl %c, jobvr %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS ", ldvl %" FLA_IS ", ldvr %" FLA_IS "",
                      *jobvl, *jobvr, *n, *lda, *ldb, *ldvl, *ldvr);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer jc, in, jr, ihi, ilo;
    real eps;
    logical ilv;
    real anrm, bnrm;
    integer ierr, itau;
    real temp;
    logical ilvl, ilvr;
    integer iwrk;
    extern logical lsame_(char *, char *, integer, integer);
    integer ileft, icols, irows;
    extern /* Subroutine */
        void
        sgghd3_(char *, char *, integer *, integer *, integer *, real *, integer *, real *,
                integer *, real *, integer *, real *, integer *, real *, integer *, integer *),
        sggbak_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *,
                integer *, integer *),
        sggbal_(char *, integer *, real *, integer *, real *, integer *, integer *, integer *,
                real *, real *, real *, integer *);
    logical ilascl, ilbscl;
    extern real slamch_(char *), slange_(char *, integer *, integer *, real *, integer *, real *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical ldumma[1];
    char chtemp[1];
    real bignum;
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *);
    integer ijobvl, iright;
    extern /* Subroutine */
        void
        sgeqrf_(integer *, integer *, real *, integer *, real *, real *, integer *, integer *);
    integer ijobvr;
    extern /* Subroutine */
        void
        slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *),
        slaset_(char *, integer *, integer *, real *, real *, real *, integer *),
        stgevc_(char *, char *, logical *, integer *, real *, integer *, real *, integer *, real *,
                integer *, real *, integer *, integer *, integer *, real *, integer *);
    real anrmto, bnrmto, smlnum;
    extern /* Subroutine */
        void
        shgeqz_(char *, char *, char *, integer *, integer *, integer *, real *, integer *, real *,
                integer *, real *, real *, real *, real *, integer *, real *, integer *, real *,
                integer *, integer *);
    extern /* Subroutine */
        void
        sorgqr_(integer *, integer *, integer *, real *, integer *, real *, real *, integer *,
                integer *);
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
        void
        sormqr_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *,
                integer *, real *, integer *, integer *);
    extern real sroundup_lwork(integer *);
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Local Arrays .. */
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
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    /* Function Body */
    if(lsame_(jobvl, "N", 1, 1))
    {
        ijobvl = 1;
        ilvl = FALSE_;
    }
    else if(lsame_(jobvl, "V", 1, 1))
    {
        ijobvl = 2;
        ilvl = TRUE_;
    }
    else
    {
        ijobvl = -1;
        ilvl = FALSE_;
    }
    if(lsame_(jobvr, "N", 1, 1))
    {
        ijobvr = 1;
        ilvr = FALSE_;
    }
    else if(lsame_(jobvr, "V", 1, 1))
    {
        ijobvr = 2;
        ilvr = TRUE_;
    }
    else
    {
        ijobvr = -1;
        ilvr = FALSE_;
    }
    ilv = ilvl || ilvr;
    /* Test the input arguments */
    *info = 0;
    lquery = *lwork == -1;
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
    else if(*ldvl < 1 || ilvl && *ldvl < *n)
    {
        *info = -12;
    }
    else if(*ldvr < 1 || ilvr && *ldvr < *n)
    {
        *info = -14;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 3; // , expr subst
        if(*lwork < fla_max(i__1, i__2) && !lquery)
        {
            *info = -16;
        }
    }
    /* Compute workspace */
    if(*info == 0)
    {
        sgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
        /* Computing MAX */
        i__1 = 1, i__2 = *n << 3;
        i__1 = fla_max(i__1, i__2);
        i__2 = *n * 3 + (integer)work[1]; // ; expr subst
        lwkopt = fla_max(i__1, i__2);
        sormqr_("L", "T", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], lda, &work[1], &c_n1,
                &ierr);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = *n * 3 + (integer)work[1]; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        sgghd3_(jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset],
                ldvl, &vr[vr_offset], ldvr, &work[1], &c_n1, &ierr);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = *n * 3 + (integer)work[1]; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        if(ilvl)
        {
            sorgqr_(n, n, n, &vl[vl_offset], ldvl, &work[1], &work[1], &c_n1, &ierr);
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = *n * 3 + (integer)work[1]; // , expr subst
            lwkopt = fla_max(i__1, i__2);
            shgeqz_("S", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], ldb,
                    &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], ldvr,
                    &work[1], &c_n1, &ierr);
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = (*n << 1) + (integer)work[1]; // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        else
        {
            shgeqz_("E", jobvl, jobvr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], ldb,
                    &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], ldvr,
                    &work[1], &c_n1, &ierr);
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = (*n << 1) + (integer)work[1]; // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        work[1] = sroundup_lwork(&lwkopt);
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGGEV3", &i__1, (ftnlen)6);
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
    eps = slamch_("P");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = slange_("M", n, n, &a[a_offset], lda, &work[1]);
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
        slascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &ierr);
    }
    /* Scale B if max element outside range [SMLNUM,BIGNUM] */
    bnrm = slange_("M", n, n, &b[b_offset], ldb, &work[1]);
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
        slascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &ierr);
    }
    /* Permute the matrices A, B to isolate eigenvalues if possible */
    ileft = 1;
    iright = *n + 1;
    iwrk = iright + *n;
    sggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[ileft], &work[iright],
            &work[iwrk], &ierr);
    /* Reduce B to triangular form (QR decomposition of B) */
    irows = ihi + 1 - ilo;
    if(ilv)
    {
        icols = *n + 1 - ilo;
    }
    else
    {
        icols = irows;
    }
    itau = iwrk;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    sgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[iwrk], &i__1, &ierr);
    /* Apply the orthogonal transformation to matrix A */
    i__1 = *lwork + 1 - iwrk;
    sormqr_("L", "T", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &work[itau],
            &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &ierr);
    /* Initialize VL */
    if(ilvl)
    {
        slaset_("Full", n, n, &c_b36, &c_b37, &vl[vl_offset], ldvl);
        if(irows > 1)
        {
            i__1 = irows - 1;
            i__2 = irows - 1;
            slacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb,
                    &vl[ilo + 1 + ilo * vl_dim1], ldvl);
        }
        i__1 = *lwork + 1 - iwrk;
        sorgqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[itau], &work[iwrk],
                &i__1, &ierr);
    }
    /* Initialize VR */
    if(ilvr)
    {
        slaset_("Full", n, n, &c_b36, &c_b37, &vr[vr_offset], ldvr);
    }
    /* Reduce to generalized Hessenberg form */
    if(ilv)
    {
        /* Eigenvectors requested -- work on whole matrix. */
        i__1 = *lwork + 1 - iwrk;
        sgghd3_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset],
                ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1, &ierr);
    }
    else
    {
        i__1 = *lwork + 1 - iwrk;
        sgghd3_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda,
                &b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr,
                &work[iwrk], &i__1, &ierr);
    }
    /* Perform QZ algorithm (Compute eigenvalues, and optionally, the */
    /* Schur forms and Schur vectors) */
    iwrk = itau;
    if(ilv)
    {
        *(unsigned char *)chtemp = 'S';
    }
    else
    {
        *(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwrk;
    shgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &alphar[1],
            &alphai[1], &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1,
            &ierr);
    if(ierr != 0)
    {
        if(ierr > 0 && ierr <= *n)
        {
            *info = ierr;
        }
        else if(ierr > *n && ierr <= *n << 1)
        {
            *info = ierr - *n;
        }
        else
        {
            *info = *n + 1;
        }
        goto L110;
    }
    /* Compute Eigenvectors */
    if(ilv)
    {
        if(ilvl)
        {
            if(ilvr)
            {
                *(unsigned char *)chtemp = 'B';
            }
            else
            {
                *(unsigned char *)chtemp = 'L';
            }
        }
        else
        {
            *(unsigned char *)chtemp = 'R';
        }
        stgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl,
                &vr[vr_offset], ldvr, n, &in, &work[iwrk], &ierr);
        if(ierr != 0)
        {
            *info = *n + 2;
            goto L110;
        }
        /* Undo balancing on VL and VR and normalization */
        if(ilvl)
        {
            sggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vl[vl_offset], ldvl,
                    &ierr);
            i__1 = *n;
            for(jc = 1; jc <= i__1; ++jc)
            {
                if(alphai[jc] < 0.f)
                {
                    goto L50;
                }
                temp = 0.f;
                if(alphai[jc] == 0.f)
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        /* Computing MAX */
                        r__2 = temp;
                        r__3 = (r__1 = vl[jr + jc * vl_dim1], f2c_abs(r__1)); // , expr subst
                        temp = fla_max(r__2, r__3);
                        /* L10: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        /* Computing MAX */
                        r__3 = temp;
                        r__4
                            = (r__1 = vl[jr + jc * vl_dim1], f2c_abs(r__1))
                              + (r__2 = vl[jr + (jc + 1) * vl_dim1], f2c_abs(r__2)); // , expr subst
                        temp = fla_max(r__3, r__4);
                        /* L20: */
                    }
                }
                if(temp < smlnum)
                {
                    goto L50;
                }
                temp = 1.f / temp;
                if(alphai[jc] == 0.f)
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        vl[jr + jc * vl_dim1] *= temp;
                        /* L30: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        vl[jr + jc * vl_dim1] *= temp;
                        vl[jr + (jc + 1) * vl_dim1] *= temp;
                        /* L40: */
                    }
                }
            L50:;
            }
        }
        if(ilvr)
        {
            sggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &vr[vr_offset], ldvr,
                    &ierr);
            i__1 = *n;
            for(jc = 1; jc <= i__1; ++jc)
            {
                if(alphai[jc] < 0.f)
                {
                    goto L100;
                }
                temp = 0.f;
                if(alphai[jc] == 0.f)
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        /* Computing MAX */
                        r__2 = temp;
                        r__3 = (r__1 = vr[jr + jc * vr_dim1], f2c_abs(r__1)); // , expr subst
                        temp = fla_max(r__2, r__3);
                        /* L60: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        /* Computing MAX */
                        r__3 = temp;
                        r__4
                            = (r__1 = vr[jr + jc * vr_dim1], f2c_abs(r__1))
                              + (r__2 = vr[jr + (jc + 1) * vr_dim1], f2c_abs(r__2)); // , expr subst
                        temp = fla_max(r__3, r__4);
                        /* L70: */
                    }
                }
                if(temp < smlnum)
                {
                    goto L100;
                }
                temp = 1.f / temp;
                if(alphai[jc] == 0.f)
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        vr[jr + jc * vr_dim1] *= temp;
                        /* L80: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(jr = 1; jr <= i__2; ++jr)
                    {
                        vr[jr + jc * vr_dim1] *= temp;
                        vr[jr + (jc + 1) * vr_dim1] *= temp;
                        /* L90: */
                    }
                }
            L100:;
            }
        }
        /* End of eigenvector calculation */
    }
/* Undo scaling if necessary */
L110:
    if(ilascl)
    {
        slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphar[1], n, &ierr);
        slascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alphai[1], n, &ierr);
    }
    if(ilbscl)
    {
        slascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &ierr);
    }
    work[1] = sroundup_lwork(&lwkopt);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGGEV3 */
}
/* sggev3_ */
