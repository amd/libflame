/* ./cgges3.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {0.f, 0.f};
static complex c_b2 = {1.f, 0.f};
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
/* > \brief <b> CGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of
 * Schur vectors for GE matrices (blocked algorithm)</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGGES3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgges3.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgges3.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgges3.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, */
/* $ LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, */
/* $ WORK, LWORK, RWORK, BWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVSL, JOBVSR, SORT */
/* INTEGER INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL BWORK( * ) */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/* $ BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), */
/* $ WORK( * ) */
/* .. */
/* .. Function Arguments .. */
/* LOGICAL SELCTG */
/* EXTERNAL SELCTG */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGES3 computes for a pair of N-by-N complex nonsymmetric matrices */
/* > (A,B), the generalized eigenvalues, the generalized complex Schur */
/* > form (S, T), and optionally left and/or right Schur vectors (VSL */
/* > and VSR). This gives the generalized Schur factorization */
/* > */
/* > (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H ) */
/* > */
/* > where (VSR)**H is the conjugate-transpose of VSR. */
/* > */
/* > Optionally, it also orders the eigenvalues so that a selected cluster */
/* > of eigenvalues appears in the leading diagonal blocks of the upper */
/* > triangular matrix S and the upper triangular matrix T. The leading */
/* > columns of VSL and VSR then form an unitary basis for the */
/* > corresponding left and right eigenspaces (deflating subspaces). */
/* > */
/* > (If only the generalized eigenvalues are needed, use the driver */
/* > CGGEV instead, which is faster.) */
/* > */
/* > A generalized eigenvalue for a pair of matrices (A,B) is a scalar w */
/* > or a ratio alpha/beta = w, such that A - w*B is singular. It is */
/* > usually represented as the pair (alpha,beta), as there is a */
/* > reasonable interpretation for beta=0, and even for both being zero. */
/* > */
/* > A pair of matrices (S,T) is in generalized complex Schur form if S */
/* > and T are upper triangular and, in addition, the diagonal elements */
/* > of T are non-negative real numbers. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVSL */
/* > \verbatim */
/* > JOBVSL is CHARACTER*1 */
/* > = 'N': do not compute the left Schur vectors;
 */
/* > = 'V': compute the left Schur vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVSR */
/* > \verbatim */
/* > JOBVSR is CHARACTER*1 */
/* > = 'N': do not compute the right Schur vectors;
 */
/* > = 'V': compute the right Schur vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] SORT */
/* > \verbatim */
/* > SORT is CHARACTER*1 */
/* > Specifies whether or not to order the eigenvalues on the */
/* > diagonal of the generalized Schur form. */
/* > = 'N': Eigenvalues are not ordered;
 */
/* > = 'S': Eigenvalues are ordered (see SELCTG). */
/* > \endverbatim */
/* > */
/* > \param[in] SELCTG */
/* > \verbatim */
/* > SELCTG is a LOGICAL FUNCTION of two COMPLEX arguments */
/* > SELCTG must be declared EXTERNAL in the calling subroutine. */
/* > If SORT = 'N', SELCTG is not referenced. */
/* > If SORT = 'S', SELCTG is used to select eigenvalues to sort */
/* > to the top left of the Schur form. */
/* > An eigenvalue ALPHA(j)/BETA(j) is selected if */
/* > SELCTG(ALPHA(j),BETA(j)) is true. */
/* > */
/* > Note that a selected complex eigenvalue may no longer satisfy */
/* > SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since */
/* > ordering may change the value of complex eigenvalues */
/* > (especially if the eigenvalue is ill-conditioned), in this */
/* > case INFO is set to N+2 (See INFO below). */
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
/* > A is COMPLEX array, dimension (LDA, N) */
/* > On entry, the first of the pair of matrices. */
/* > On exit, A has been overwritten by its generalized Schur */
/* > form S. */
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
/* > On entry, the second of the pair of matrices. */
/* > On exit, B has been overwritten by its generalized Schur */
/* > form T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] SDIM */
/* > \verbatim */
/* > SDIM is INTEGER */
/* > If SORT = 'N', SDIM = 0. */
/* > If SORT = 'S', SDIM = number of eigenvalues (after sorting) */
/* > for which SELCTG is true. */
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
/* > generalized eigenvalues. ALPHA(j), j=1,...,N and BETA(j), */
/* > j=1,...,N are the diagonals of the complex Schur form (A,B) */
/* > output by CGGES3. The BETA(j) will be non-negative real. */
/* > */
/* > Note: the quotients ALPHA(j)/BETA(j) may easily over- or */
/* > underflow, and BETA(j) may even be zero. Thus, the user */
/* > should avoid naively computing the ratio alpha/beta. */
/* > However, ALPHA will be always less than and usually */
/* > comparable with norm(A) in magnitude, and BETA always less */
/* > than and usually comparable with norm(B). */
/* > \endverbatim */
/* > */
/* > \param[out] VSL */
/* > \verbatim */
/* > VSL is COMPLEX array, dimension (LDVSL,N) */
/* > If JOBVSL = 'V', VSL will contain the left Schur vectors. */
/* > Not referenced if JOBVSL = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVSL */
/* > \verbatim */
/* > LDVSL is INTEGER */
/* > The leading dimension of the matrix VSL. LDVSL >= 1, and */
/* > if JOBVSL = 'V', LDVSL >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VSR */
/* > \verbatim */
/* > VSR is COMPLEX array, dimension (LDVSR,N) */
/* > If JOBVSR = 'V', VSR will contain the right Schur vectors. */
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
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
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
/* > \param[out] BWORK */
/* > \verbatim */
/* > BWORK is LOGICAL array, dimension (N) */
/* > Not referenced if SORT = 'N'. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > =1,...,N: */
/* > The QZ iteration failed. (A,B) are not in Schur */
/* > form, but ALPHA(j) and BETA(j) should be correct for */
/* > j=INFO+1,...,N. */
/* > > N: =N+1: other than QZ iteration failed in CLAQZ0 */
/* > =N+2: after reordering, roundoff changed values of */
/* > some complex eigenvalues so that leading */
/* > eigenvalues in the Generalized Schur form no */
/* > longer satisfy SELCTG=.TRUE. This could also */
/* > be caused due to scaling. */
/* > =N+3: reordering failed in CTGSEN. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gges3 */
/* ===================================================================== */
/* Subroutine */
void cgges3_(char *jobvsl, char *jobvsr, char *sort, L_fp2 selctg, integer *n, complex *a,
             integer *lda, complex *b, integer *ldb, integer *sdim, complex *alpha, complex *beta,
             complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, complex *work,
             integer *lwork, real *rwork, logical *bwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cgges3 inputs: jobvsl %c, jobvsr %c, sort %c, n %lld, lda %lld, ldb %lld, ldvsl "
             "%lld, ldvsr %lld, lwork %lld",
             *jobvsl, *jobvsr, *sort, *n, *lda, *ldb, *ldvsl, *ldvsr, *lwork);
#else
    snprintf(buffer, 256,
             "cgges3 inputs: jobvsl %c, jobvsr %c, sort %c, n %d, lda %d, ldb %d, ldvsl %d, ldvsr "
             "%d, lwork %d",
             *jobvsl, *jobvsr, *sort, *n, *lda, *ldb, *ldvsl, *ldvsr, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, vsr_dim1, vsr_offset, i__1,
        i__2;
    complex q__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    real dif[2];
    integer ihi, ilo;
    real eps, anrm, bnrm;
    integer idum[1], ierr, itau, iwrk;
    real pvsl, pvsr;
    extern logical lsame_(char *, char *, integer, integer);
    integer ileft, icols;
    logical cursl, ilvsl, ilvsr;
    integer irwrk;
    extern /* Subroutine */
        void
        cgghd3_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                integer *, complex *, integer *, complex *, integer *, complex *, integer *,
                integer *);
    integer irows;
    extern /* Subroutine */
        void
        cggbak_(char *, char *, integer *, integer *, integer *, real *, real *, integer *,
                complex *, integer *, integer *),
        cggbal_(char *, integer *, complex *, integer *, complex *, integer *, integer *, integer *,
                real *, real *, real *, integer *);
    extern real clange_(char *, integer *, integer *, complex *, integer *, real *);
    extern /* Subroutine */
        void
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
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    extern /* Subroutine */
        void
        chgeqz_(char *, char *, char *, integer *, integer *, integer *, complex *, integer *,
                complex *, integer *, complex *, complex *, complex *, integer *, complex *,
                integer *, complex *, integer *, real *, integer *),
        ctgsen_(integer *, logical *, logical *, logical *, integer *, complex *, integer *,
                complex *, integer *, complex *, complex *, complex *, integer *, complex *,
                integer *, integer *, real *, real *, real *, complex *, integer *, integer *,
                integer *, integer *);
    integer ijobvl, iright, ijobvr;
    real anrmto, bnrmto;
    logical lastsl;
    extern /* Subroutine */
        void
        cungqr_(integer *, integer *, integer *, complex *, integer *, complex *, complex *,
                integer *, integer *),
        cunmqr_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                complex *, integer *, complex *, integer *, integer *);
    real smlnum;
    logical wantst, lquery;
    integer lwkopt;
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Function Arguments .. */
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
    --alpha;
    --beta;
    vsl_dim1 = *ldvsl;
    vsl_offset = 1 + vsl_dim1;
    vsl -= vsl_offset;
    vsr_dim1 = *ldvsr;
    vsr_offset = 1 + vsr_dim1;
    vsr -= vsr_offset;
    --work;
    --rwork;
    --bwork;
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
    wantst = lsame_(sort, "S", 1, 1);
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
    else if(!wantst && !lsame_(sort, "N", 1, 1))
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -9;
    }
    else if(*ldvsl < 1 || ilvsl && *ldvsl < *n)
    {
        *info = -14;
    }
    else if(*ldvsr < 1 || ilvsr && *ldvsr < *n)
    {
        *info = -16;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        if(*lwork < fla_max(i__1, i__2) && !lquery)
        {
            *info = -18;
        }
    }
    /* Compute workspace */
    if(*info == 0)
    {
        cgeqrf_(n, n, &b[b_offset], ldb, &work[1], &work[1], &c_n1, &ierr);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n + (integer)work[1].r; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        cunmqr_("L", "C", n, n, n, &b[b_offset], ldb, &work[1], &a[a_offset], lda, &work[1], &c_n1,
                &ierr);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = *n + (integer)work[1].r; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        if(ilvsl)
        {
            cungqr_(n, n, n, &vsl[vsl_offset], ldvsl, &work[1], &work[1], &c_n1, &ierr);
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = *n + (integer)work[1].r; // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        cgghd3_(jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset],
                ldvsl, &vsr[vsr_offset], ldvsr, &work[1], &c_n1, &ierr);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = *n + (integer)work[1].r; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        chgeqz_("S", jobvsl, jobvsr, n, &c__1, n, &a[a_offset], lda, &b[b_offset], ldb, &alpha[1],
                &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[1], &c_n1,
                &rwork[1], &ierr);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer)work[1].r; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        if(wantst)
        {
            ctgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[b_offset], ldb,
                    &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim,
                    &pvsl, &pvsr, dif, &work[1], &c_n1, idum, &c__1, &ierr);
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = (integer)work[1].r; // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        q__1.r = (real)lwkopt;
        q__1.i = 0.f; // , expr subst
        work[1].r = q__1.r;
        work[1].i = q__1.i; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGGES3", &i__1, (ftnlen)6);
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
        *sdim = 0;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Get machine constants */
    eps = slamch_("P");
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
    /* Permute the matrix to make it more nearly triangular */
    ileft = 1;
    iright = *n + 1;
    irwrk = iright + *n;
    cggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright],
            &rwork[irwrk], &ierr);
    /* Reduce B to triangular form (QR decomposition of B) */
    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    cgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[iwrk], &i__1, &ierr);
    /* Apply the orthogonal transformation to matrix A */
    i__1 = *lwork + 1 - iwrk;
    cunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &work[itau],
            &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &ierr);
    /* Initialize VSL */
    if(ilvsl)
    {
        claset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl);
        if(irows > 1)
        {
            i__1 = irows - 1;
            i__2 = irows - 1;
            clacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb,
                    &vsl[ilo + 1 + ilo * vsl_dim1], ldvsl);
        }
        i__1 = *lwork + 1 - iwrk;
        cungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &work[itau], &work[iwrk],
                &i__1, &ierr);
    }
    /* Initialize VSR */
    if(ilvsr)
    {
        claset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr);
    }
    /* Reduce to generalized Hessenberg form */
    i__1 = *lwork + 1 - iwrk;
    cgghd3_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset],
            ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1, &ierr);
    *sdim = 0;
    /* Perform QZ algorithm, computing Schur vectors if desired */
    iwrk = itau;
    i__1 = *lwork + 1 - iwrk;
    chgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &alpha[1],
            &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, &work[iwrk], &i__1,
            &rwork[irwrk], &ierr);
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
        goto L30;
    }
    /* Sort eigenvalues ALPHA/BETA if desired */
    if(wantst)
    {
        /* Undo scaling on eigenvalues before selecting */
        if(ilascl)
        {
            clascl_("G", &c__0, &c__0, &anrm, &anrmto, n, &c__1, &alpha[1], n, &ierr);
        }
        if(ilbscl)
        {
            clascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, &c__1, &beta[1], n, &ierr);
        }
        /* Select eigenvalues */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
            /* L10: */
        }
        i__1 = *lwork - iwrk + 1;
        ctgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[b_offset], ldb,
                &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl,
                &pvsr, dif, &work[iwrk], &i__1, idum, &c__1, &ierr);
        if(ierr == 1)
        {
            *info = *n + 3;
        }
    }
    /* Apply back-permutation to VSL and VSR */
    if(ilvsl)
    {
        cggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &vsl[vsl_offset], ldvsl,
                &ierr);
    }
    if(ilvsr)
    {
        cggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &vsr[vsr_offset], ldvsr,
                &ierr);
    }
    /* Undo scaling */
    if(ilascl)
    {
        clascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &ierr);
        clascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &ierr);
    }
    if(ilbscl)
    {
        clascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &ierr);
        clascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &ierr);
    }
    if(wantst)
    {
        /* Check if reordering is correct */
        lastsl = TRUE_;
        *sdim = 0;
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            cursl = (*selctg)(&alpha[i__], &beta[i__]);
            if(cursl)
            {
                ++(*sdim);
            }
            if(cursl && !lastsl)
            {
                *info = *n + 2;
            }
            lastsl = cursl;
            /* L20: */
        }
    }
L30:
    q__1.r = (real)lwkopt;
    q__1.i = 0.f; // , expr subst
    work[1].r = q__1.r;
    work[1].i = q__1.i; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGGES3 */
}
/* cgges3_ */
