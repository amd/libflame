/* ./cggev.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {0.f, 0.f};
static complex c_b2 = {1.f, 0.f};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
/* > \brief <b> CGGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for GE matr ices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGGEV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggev.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggev.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggev.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, */
/* VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVL, JOBVR */
/* INTEGER INFO, LDA, LDB, LDVL, LDVR, LWORK, N */
/* .. */
/* .. Array Arguments .. */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/* $ BETA( * ), VL( LDVL, * ), VR( LDVR, * ), */
/* $ WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGEV computes for a pair of N-by-N complex nonsymmetric matrices */
/* > (A,B), the generalized eigenvalues, and optionally, the left and/or */
/* > right generalized eigenvectors. */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar */
/* > lambda or a ratio alpha/beta = lambda, such that A - lambda*B is */
/* > singular. It is usually represented as the pair (alpha,beta), as */
/* > there is a reasonable interpretation for beta=0, and even for both */
/* > being zero. */
/* > */
/* > The right generalized eigenvector v(j) corresponding to the */
/* > generalized eigenvalue lambda(j) of (A,B) satisfies */
/* > */
/* > A * v(j) = lambda(j) * B * v(j). */
/* > */
/* > The left generalized eigenvector u(j) corresponding to the */
/* > generalized eigenvalues lambda(j) of (A,B) satisfies */
/* > */
/* > u(j)**H * A = lambda(j) * u(j)**H * B */
/* > */
/* > where u(j)**H is the conjugate-transpose of u(j). */
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
/* > A is COMPLEX array, dimension (LDA, N) */
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
/* > B is COMPLEX array, dimension (LDB, N) */
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
/* > \param[out] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is COMPLEX array, dimension (N) */
/* > On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the */
/* > generalized eigenvalues. */
/* > */
/* > Note: the quotients ALPHA(j)/BETA(j) may easily over- or */
/* > underflow, and BETA(j) may even be zero. Thus, the user */
/* > should avoid naively computing the ratio alpha/beta. */
/* > However, ALPHA will be always less than and usually */
/* > comparable with norm(A) in magnitude, and BETA always less */
/* > than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* > VL is COMPLEX array, dimension (LDVL,N) */
/* > If JOBVL = 'V', the left generalized eigenvectors u(j) are */
/* > stored one after another in the columns of VL, in the same */
/* > order as their eigenvalues. */
/* > Each eigenvector is scaled so the largest component has */
/* > f2c_abs(real part) + f2c_abs(imag. part) = 1. */
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
/* > VR is COMPLEX array, dimension (LDVR,N) */
/* > If JOBVR = 'V', the right generalized eigenvectors v(j) are */
/* > stored one after another in the columns of VR, in the same */
/* > order as their eigenvalues. */
/* > Each eigenvector is scaled so the largest component has */
/* > f2c_abs(real part) + f2c_abs(imag. part) = 1. */
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
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,2*N). */
/* > For good performance, LWORK must generally be larger. */
/* > */
/* > If LWORK = -1, then a workspace query is assumed;
the routine */
/* > only calculates the optimal size of the WORK array, returns */
/* > this value as the first entry of the WORK array, and no error */
/* > message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (8*N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > =1,...,N: */
/* > The QZ iteration failed. No eigenvectors have been */
/* > calculated, but ALPHA(j) and BETA(j) should be */
/* > correct for j=INFO+1,...,N. */
/* > > N: =N+1: other then QZ iteration failed in CHGEQZ, */
/* > =N+2: error return from CTGEVC. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup ggev */
/* ===================================================================== */
/* Subroutine */
void cggev_(char *jobvl, char *jobvr, integer *n, complex *a, integer *lda, complex *b,
            integer *ldb, complex *alpha, complex *beta, complex *vl, integer *ldvl, complex *vr,
            integer *ldvr, complex *work, integer *lwork, real *rwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cggev inputs: jobvl %c, jobvr %c, n %lld, lda %lld, ldb %lld, ldvl %lld, ldvr %lld, "
             "lwork %lld",
             *jobvl, *jobvr, *n, *lda, *ldb, *ldvl, *ldvr, *lwork);
#else
    snprintf(buffer, 256,
             "cggev inputs: jobvl %c, jobvr %c, n %d, lda %d, ldb %d, ldvl %d, ldvr %d, lwork %d",
             *jobvl, *jobvr, *n, *lda, *ldb, *ldvl, *ldvr, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2,
        i__3, i__4;
    real r__1, r__2, r__3, r__4;
    complex q__1;
    /* Builtin functions */
    double sqrt(doublereal), r_imag(complex *);
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
    integer ileft, icols, irwrk, irows;
    extern /* Subroutine */
        void
        cggbak_(char *, char *, integer *, integer *, integer *, real *, real *, integer *,
                complex *, integer *, integer *),
        cggbal_(char *, integer *, complex *, integer *, complex *, integer *, integer *, integer *,
                real *, real *, real *, integer *);
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
        void
        cgghrd_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                integer *, complex *, integer *, complex *, integer *, integer *),
        clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *,
                integer *, integer *);
    logical ilascl, ilbscl;
    extern /* Subroutine */
        void
        cgeqrf_(integer *, integer *, complex *, integer *, complex *, complex *, integer *,
                integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *),
        claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *),
        ctgevc_(char *, char *, logical *, integer *, complex *, integer *, complex *, integer *,
                complex *, integer *, complex *, integer *, integer *, integer *, complex *, real *,
                integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical ldumma[1];
    char chtemp[1];
    real bignum;
    extern /* Subroutine */
        void
        chgeqz_(char *, char *, char *, integer *, integer *, integer *, complex *, integer *,
                complex *, integer *, complex *, complex *, complex *, integer *, complex *,
                integer *, complex *, integer *, real *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    integer ijobvl, iright, ijobvr;
    extern /* Subroutine */
        void
        cungqr_(integer *, integer *, integer *, complex *, integer *, complex *, complex *,
                integer *, integer *);
    real anrmto;
    integer lwkmin;
    real bnrmto;
    extern /* Subroutine */
        void
        cunmqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                complex *, integer *, complex *, integer *, integer *);
    real smlnum;
    integer lwkopt;
    logical lquery;
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
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
    --alpha;
    --beta;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    --work;
    --rwork;
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
        *info = -11;
    }
    else if(*ldvr < 1 || ilvr && *ldvr < *n)
    {
        *info = -13;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV. The workspace is */
    /* computed assuming ILO = 1 and IHI = N, the worst case.) */
    if(*info == 0)
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        lwkmin = fla_max(i__1, i__2);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n + *n * ilaenv_(&c__1, "CGEQRF", " ", n, &c__1, n, &c__0); // , expr subst
        lwkopt = fla_max(i__1, i__2);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = *n + *n * ilaenv_(&c__1, "CUNMQR", " ", n, &c__1, n, &c__0); // , expr subst
        lwkopt = fla_max(i__1, i__2);
        if(ilvl)
        {
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = *n + *n * ilaenv_(&c__1, "CUNGQR", " ", n, &c__1, n, &c_n1); // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        r__1 = sroundup_lwork(&lwkopt);
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
        if(*lwork < lwkmin && !lquery)
        {
            *info = -15;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGGEV ", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Get machine constants */
    eps = slamch_("E") * slamch_("B");
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = clange_("M", n, n, &a[a_offset], lda, &rwork[1]);
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
        clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &ierr);
    }
    /* Scale B if max element outside range [SMLNUM,BIGNUM] */
    bnrm = clange_("M", n, n, &b[b_offset], ldb, &rwork[1]);
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
        clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &ierr);
    }
    /* Permute the matrices A, B to isolate eigenvalues if possible */
    /* (Real Workspace: need 6*N) */
    ileft = 1;
    iright = *n + 1;
    irwrk = iright + *n;
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright],
            &rwork[irwrk], &ierr);
    /* Reduce B to triangular form (QR decomposition of B) */
    /* (Complex Workspace: need N, prefer N*NB) */
    irows = ihi + 1 - ilo;
    if(ilv)
    {
        icols = *n + 1 - ilo;
    }
    else
    {
        icols = irows;
    }
    itau = 1;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[iwrk], &i__1, &ierr);
    /* Apply the orthogonal transformation to matrix A */
    /* (Complex Workspace: need N, prefer N*NB) */
    i__1 = *lwork + 1 - iwrk;
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &work[itau],
            &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &ierr);
    /* Initialize VL */
    /* (Complex Workspace: need N, prefer N*NB) */
    if(ilvl)
    {
        claset_("Full", n, n, &c_b1, &c_b2, &vl[vl_offset], ldvl);
        if(irows > 1)
        {
            i__1 = irows - 1;
            i__2 = irows - 1;
            clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb,
                    &vl[ilo + 1 + ilo * vl_dim1], ldvl);
        }
        i__1 = *lwork + 1 - iwrk;
        cungqr_(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[itau], &work[iwrk],
                &i__1, &ierr);
    }
    /* Initialize VR */
    if(ilvr)
    {
        claset_("Full", n, n, &c_b1, &c_b2, &vr[vr_offset], ldvr);
    }
    /* Reduce to generalized Hessenberg form */
    if(ilv)
    {
        /* Eigenvectors requested -- work on whole matrix. */
        cgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset],
                ldvl, &vr[vr_offset], ldvr, &ierr);
    }
    else
    {
        cgghrd_("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda,
                &b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &ierr);
    }
    /* Perform QZ algorithm (Compute eigenvalues, and optionally, the */
    /* Schur form and Schur vectors) */
    /* (Complex Workspace: need N) */
    /* (Real Workspace: need N) */
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
    chgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &alpha[1],
            &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &work[iwrk], &i__1, &rwork[irwrk],
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
        goto L70;
    }
    /* Compute Eigenvectors */
    /* (Real Workspace: need 2*N) */
    /* (Complex Workspace: need 2*N) */
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
        ctgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl,
                &vr[vr_offset], ldvr, n, &in, &work[iwrk], &rwork[irwrk], &ierr);
        if(ierr != 0)
        {
            *info = *n + 2;
            goto L70;
        }
        /* Undo balancing on VL and VR and normalization */
        /* (Workspace: none needed) */
        if(ilvl)
        {
            cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &vl[vl_offset], ldvl,
                    &ierr);
            i__1 = *n;
            for(jc = 1; jc <= i__1; ++jc)
            {
                temp = 0.f;
                i__2 = *n;
                for(jr = 1; jr <= i__2; ++jr)
                {
                    /* Computing MAX */
                    i__3 = jr + jc * vl_dim1;
                    r__3 = temp;
                    r__4 = (r__1 = vl[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&vl[jr + jc * vl_dim1]), f2c_abs(r__2)); // , expr subst
                    temp = fla_max(r__3, r__4);
                    /* L10: */
                }
                if(temp < smlnum)
                {
                    goto L30;
                }
                temp = 1.f / temp;
                i__2 = *n;
                for(jr = 1; jr <= i__2; ++jr)
                {
                    i__3 = jr + jc * vl_dim1;
                    i__4 = jr + jc * vl_dim1;
                    q__1.r = temp * vl[i__4].r;
                    q__1.i = temp * vl[i__4].i; // , expr subst
                    vl[i__3].r = q__1.r;
                    vl[i__3].i = q__1.i; // , expr subst
                    /* L20: */
                }
            L30:;
            }
        }
        if(ilvr)
        {
            cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &vr[vr_offset], ldvr,
                    &ierr);
            i__1 = *n;
            for(jc = 1; jc <= i__1; ++jc)
            {
                temp = 0.f;
                i__2 = *n;
                for(jr = 1; jr <= i__2; ++jr)
                {
                    /* Computing MAX */
                    i__3 = jr + jc * vr_dim1;
                    r__3 = temp;
                    r__4 = (r__1 = vr[i__3].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&vr[jr + jc * vr_dim1]), f2c_abs(r__2)); // , expr subst
                    temp = fla_max(r__3, r__4);
                    /* L40: */
                }
                if(temp < smlnum)
                {
                    goto L60;
                }
                temp = 1.f / temp;
                i__2 = *n;
                for(jr = 1; jr <= i__2; ++jr)
                {
                    i__3 = jr + jc * vr_dim1;
                    i__4 = jr + jc * vr_dim1;
                    q__1.r = temp * vr[i__4].r;
                    q__1.i = temp * vr[i__4].i; // , expr subst
                    vr[i__3].r = q__1.r;
                    vr[i__3].i = q__1.i; // , expr subst
                    /* L50: */
                }
            L60:;
            }
        }
    }
/* Undo scaling if necessary */
L70:
    if(ilascl)
    {
        clascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &ierr);
    }
    if(ilbscl)
    {
        clascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &ierr);
    }
    r__1 = sroundup_lwork(&lwkopt);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGGEV */
}
/* cggev_ */
