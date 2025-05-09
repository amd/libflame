/* ./zgges.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 = {0., 0.};
static doublecomplex c_b2 = {1., 0.};
static integer c__1 = 1;
static integer c__0 = 0;
static integer c_n1 = -1;
/* > \brief <b> ZGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
 * vectors f or GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGGES + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgges.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgges.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgges.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, */
/* SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, */
/* LWORK, RWORK, BWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBVSL, JOBVSR, SORT */
/* INTEGER INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM */
/* .. */
/* .. Array Arguments .. */
/* LOGICAL BWORK( * ) */
/* DOUBLE PRECISION RWORK( * ) */
/* COMPLEX*16 A( LDA, * ), ALPHA( * ), B( LDB, * ), */
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
/* > ZGGES computes for a pair of N-by-N complex nonsymmetric matrices */
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
/* > ZGGEV instead, which is faster.) */
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
/* > SELCTG is a LOGICAL FUNCTION of two COMPLEX*16 arguments */
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
/* > A is COMPLEX*16 array, dimension (LDA, N) */
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
/* > B is COMPLEX*16 array, dimension (LDB, N) */
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
/* > ALPHA is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is COMPLEX*16 array, dimension (N) */
/* > On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the */
/* > generalized eigenvalues. ALPHA(j), j=1,...,N and BETA(j), */
/* > j=1,...,N are the diagonals of the complex Schur form (A,B) */
/* > output by ZGGES. The BETA(j) will be non-negative real. */
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
/* > VSL is COMPLEX*16 array, dimension (LDVSL,N) */
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
/* > VSR is COMPLEX*16 array, dimension (LDVSR,N) */
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
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* > RWORK is DOUBLE PRECISION array, dimension (8*N) */
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
/* > > N: =N+1: other than QZ iteration failed in ZHGEQZ */
/* > =N+2: after reordering, roundoff changed values of */
/* > some complex eigenvalues so that leading */
/* > eigenvalues in the Generalized Schur form no */
/* > longer satisfy SELCTG=.TRUE. This could also */
/* > be caused due to scaling. */
/* > =N+3: reordering failed in ZTGSEN. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gges */
/* ===================================================================== */
/* Subroutine */
void zgges_(char *jobvsl, char *jobvsr, char *sort, L_fpz2 selctg, integer *n, doublecomplex *a,
            integer *lda, doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha,
            doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr,
            integer *ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork,
            integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgges inputs: jobvsl %c, jobvsr %c, sort %c, n %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS ", sdim %" FLA_IS ", ldvsl %" FLA_IS ", ldvsr %" FLA_IS "",
                      *jobvsl, *jobvsr, *sort, *n, *lda, *ldb, *sdim, *ldvsl, *ldvsr);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vsl_dim1, vsl_offset, vsr_dim1, vsr_offset, i__1,
        i__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    doublereal dif[2];
    integer ihi, ilo;
    doublereal eps, anrm, bnrm;
    integer idum[1], ierr, itau, iwrk;
    doublereal pvsl, pvsr;
    extern logical lsame_(char *, char *, integer, integer);
    integer ileft, icols;
    logical cursl, ilvsl, ilvsr;
    integer irwrk, irows;
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
        void
        zggbak_(char *, char *, integer *, integer *, integer *, doublereal *, doublereal *,
                integer *, doublecomplex *, integer *, integer *),
        zggbal_(char *, integer *, doublecomplex *, integer *, doublecomplex *, integer *,
                integer *, integer *, doublereal *, doublereal *, doublereal *, integer *);
    logical ilascl, ilbscl;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, integer *,
                              doublereal *);
    doublereal bignum;
    integer ijobvl, iright;
    extern /* Subroutine */
        void
        zgghrd_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *, integer *,
                integer *),
        zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublecomplex *, integer *, integer *);
    integer ijobvr;
    extern /* Subroutine */
        void
        zgeqrf_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *,
                integer *, integer *);
    doublereal anrmto;
    integer lwkmin;
    logical lastsl;
    doublereal bnrmto;
    extern /* Subroutine */
        void
        zlacpy_(char *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                integer *),
        zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
                integer *),
        zhgeqz_(char *, char *, char *, integer *, integer *, integer *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
                integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *,
                integer *),
        ztgsen_(integer *, logical *, logical *, logical *, integer *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
                integer *, doublecomplex *, integer *, integer *, doublereal *, doublereal *,
                doublereal *, doublecomplex *, integer *, integer *, integer *, integer *);
    doublereal smlnum;
    logical wantst, lquery;
    integer lwkopt;
    extern /* Subroutine */
        void
        zungqr_(integer *, integer *, integer *, doublecomplex *, integer *, doublecomplex *,
                doublecomplex *, integer *, integer *),
        zunmqr_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *,
                doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, integer *);
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
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV.) */
    if(*info == 0)
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        lwkmin = fla_max(i__1, i__2);
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n + *n * ilaenv_(&c__1, "ZGEQRF", " ", n, &c__1, n, &c__0); // , expr subst
        lwkopt = fla_max(i__1, i__2);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = *n + *n * ilaenv_(&c__1, "ZUNMQR", " ", n, &c__1, n, &c_n1); // , expr subst
        lwkopt = fla_max(i__1, i__2);
        if(ilvsl)
        {
            /* Computing MAX */
            i__1 = lwkopt;
            i__2 = *n + *n * ilaenv_(&c__1, "ZUNGQR", " ", n, &c__1, n, &c_n1); // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        work[1].r = (doublereal)lwkopt;
        work[1].i = 0.; // , expr subst
        if(*lwork < lwkmin && !lquery)
        {
            *info = -18;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGGES ", &i__1, (ftnlen)6);
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
        *sdim = 0;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants */
    eps = dlamch_("P");
    smlnum = dlamch_("S");
    bignum = 1. / smlnum;
    smlnum = sqrt(smlnum) / eps;
    bignum = 1. / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = zlange_("M", n, n, &a[a_offset], lda, &rwork[1]);
    ilascl = FALSE_;
    if(anrm > 0. && anrm < smlnum)
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
        zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, n, &a[a_offset], lda, &ierr);
    }
    /* Scale B if max element outside range [SMLNUM,BIGNUM] */
    bnrm = zlange_("M", n, n, &b[b_offset], ldb, &rwork[1]);
    ilbscl = FALSE_;
    if(bnrm > 0. && bnrm < smlnum)
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
        zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, n, &b[b_offset], ldb, &ierr);
    }
    /* Permute the matrix to make it more nearly triangular */
    /* (Real Workspace: need 6*N) */
    ileft = 1;
    iright = *n + 1;
    irwrk = iright + *n;
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright],
            &rwork[irwrk], &ierr);
    /* Reduce B to triangular form (QR decomposition of B) */
    /* (Complex Workspace: need N, prefer N*NB) */
    irows = ihi + 1 - ilo;
    icols = *n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
    i__1 = *lwork + 1 - iwrk;
    zgeqrf_(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[iwrk], &i__1, &ierr);
    /* Apply the orthogonal transformation to matrix A */
    /* (Complex Workspace: need N, prefer N*NB) */
    i__1 = *lwork + 1 - iwrk;
    zunmqr_("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &work[itau],
            &a[ilo + ilo * a_dim1], lda, &work[iwrk], &i__1, &ierr);
    /* Initialize VSL */
    /* (Complex Workspace: need N, prefer N*NB) */
    if(ilvsl)
    {
        zlaset_("Full", n, n, &c_b1, &c_b2, &vsl[vsl_offset], ldvsl);
        if(irows > 1)
        {
            i__1 = irows - 1;
            i__2 = irows - 1;
            zlacpy_("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb,
                    &vsl[ilo + 1 + ilo * vsl_dim1], ldvsl);
        }
        i__1 = *lwork + 1 - iwrk;
        zungqr_(&irows, &irows, &irows, &vsl[ilo + ilo * vsl_dim1], ldvsl, &work[itau], &work[iwrk],
                &i__1, &ierr);
    }
    /* Initialize VSR */
    if(ilvsr)
    {
        zlaset_("Full", n, n, &c_b1, &c_b2, &vsr[vsr_offset], ldvsr);
    }
    /* Reduce to generalized Hessenberg form */
    /* (Workspace: none needed) */
    zgghrd_(jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset],
            ldvsl, &vsr[vsr_offset], ldvsr, &ierr);
    *sdim = 0;
    /* Perform QZ algorithm, computing Schur vectors if desired */
    /* (Complex Workspace: need N) */
    /* (Real Workspace: need N) */
    iwrk = itau;
    i__1 = *lwork + 1 - iwrk;
    zhgeqz_("S", jobvsl, jobvsr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb, &alpha[1],
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
    /* (Workspace: none needed) */
    if(wantst)
    {
        /* Undo scaling on eigenvalues before selecting */
        if(ilascl)
        {
            zlascl_("G", &c__0, &c__0, &anrm, &anrmto, n, &c__1, &alpha[1], n, &ierr);
        }
        if(ilbscl)
        {
            zlascl_("G", &c__0, &c__0, &bnrm, &bnrmto, n, &c__1, &beta[1], n, &ierr);
        }
        /* Select eigenvalues */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            bwork[i__] = (*selctg)(&alpha[i__], &beta[i__]);
            /* L10: */
        }
        i__1 = *lwork - iwrk + 1;
        ztgsen_(&c__0, &ilvsl, &ilvsr, &bwork[1], n, &a[a_offset], lda, &b[b_offset], ldb,
                &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr, sdim, &pvsl,
                &pvsr, dif, &work[iwrk], &i__1, idum, &c__1, &ierr);
        if(ierr == 1)
        {
            *info = *n + 3;
        }
    }
    /* Apply back-permutation to VSL and VSR */
    /* (Workspace: none needed) */
    if(ilvsl)
    {
        zggbak_("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &vsl[vsl_offset], ldvsl,
                &ierr);
    }
    if(ilvsr)
    {
        zggbak_("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n, &vsr[vsr_offset], ldvsr,
                &ierr);
    }
    /* Undo scaling */
    if(ilascl)
    {
        zlascl_("U", &c__0, &c__0, &anrmto, &anrm, n, n, &a[a_offset], lda, &ierr);
        zlascl_("G", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &alpha[1], n, &ierr);
    }
    if(ilbscl)
    {
        zlascl_("U", &c__0, &c__0, &bnrmto, &bnrm, n, n, &b[b_offset], ldb, &ierr);
        zlascl_("G", &c__0, &c__0, &bnrmto, &bnrm, n, &c__1, &beta[1], n, &ierr);
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
    work[1].r = (doublereal)lwkopt;
    work[1].i = 0.; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGGES */
}
/* zgges_ */
