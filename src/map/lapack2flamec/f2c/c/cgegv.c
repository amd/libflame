/* ../netlib/cgegv.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {0.f, 0.f};
static scomplex c_b2 = {1.f, 0.f};
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static real c_b29 = 1.f;
/* > \brief <b> CGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors
 * for GE mat rices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEGV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgegv.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgegv.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgegv.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, */
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
/* > This routine is deprecated and has been replaced by routine CGGEV. */
/* > */
/* > CGEGV computes the eigenvalues and, optionally, the left and/or right */
/* > eigenvectors of a scomplex matrix pair (A,B). */
/* > Given two square matrices A and B, */
/* > the generalized nonsymmetric eigenvalue problem (GNEP) is to find the */
/* > eigenvalues lambda and corresponding (non-zero) eigenvectors x such */
/* > that */
/* > A*x = lambda*B*x. */
/* > */
/* > An alternate form is to find the eigenvalues mu and corresponding */
/* > eigenvectors y such that */
/* > mu*A*y = B*y. */
/* > */
/* > These two forms are equivalent with mu = 1/lambda and x = y if */
/* > neither lambda nor mu is zero. In order to deal with the case that */
/* > lambda or mu is zero or small, two values alpha and beta are returned */
/* > for each eigenvalue, such that lambda = alpha/beta and */
/* > mu = beta/alpha. */
/* > */
/* > The vectors x and y in the above equations are right eigenvectors of */
/* > the matrix pair (A,B). Vectors u and v satisfying */
/* > u**H*A = lambda*u**H*B or mu*v**H*A = v**H*B */
/* > are left eigenvectors of (A,B). */
/* > */
/* > Note: this routine performs "full balancing" on A and B */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBVL */
/* > \verbatim */
/* > JOBVL is CHARACTER*1 */
/* > = 'N': do not compute the left generalized eigenvectors;
 */
/* > = 'V': compute the left generalized eigenvectors (returned */
/* > in VL). */
/* > \endverbatim */
/* > */
/* > \param[in] JOBVR */
/* > \verbatim */
/* > JOBVR is CHARACTER*1 */
/* > = 'N': do not compute the right generalized eigenvectors;
 */
/* > = 'V': compute the right generalized eigenvectors (returned */
/* > in VR). */
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
/* > On entry, the matrix A. */
/* > If JOBVL = 'V' or JOBVR = 'V', then on exit A */
/* > contains the Schur form of A from the generalized Schur */
/* > factorization of the pair (A,B) after balancing. If no */
/* > eigenvectors were computed, then only the diagonal elements */
/* > of the Schur form will be correct. See CGGHRD and CHGEQZ */
/* > for details. */
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
/* > On entry, the matrix B. */
/* > If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the */
/* > upper triangular matrix obtained from B in the generalized */
/* > Schur factorization of the pair (A,B) after balancing. */
/* > If no eigenvectors were computed, then only the diagonal */
/* > elements of B will be correct. See CGGHRD and CHGEQZ for */
/* > details. */
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
/* > The scomplex scalars alpha that define the eigenvalues of */
/* > GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is COMPLEX array, dimension (N) */
/* > The scomplex scalars beta that define the eigenvalues of GNEP. */
/* > */
/* > Together, the quantities alpha = ALPHA(j) and beta = BETA(j) */
/* > represent the j-th eigenvalue of the matrix pair (A,B), in */
/* > one of the forms lambda = alpha/beta or mu = beta/alpha. */
/* > Since either lambda or mu may overflow, they should not, */
/* > in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[out] VL */
/* > \verbatim */
/* > VL is COMPLEX array, dimension (LDVL,N) */
/* > If JOBVL = 'V', the left eigenvectors u(j) are stored */
/* > in the columns of VL, in the same order as their eigenvalues. */
/* > Each eigenvector is scaled so that its largest component has */
/* > f2c_abs(real part) + f2c_abs(imag. part) = 1, except for eigenvectors */
/* > corresponding to an eigenvalue with alpha = beta = 0, which */
/* > are set to zero. */
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
/* > If JOBVR = 'V', the right eigenvectors x(j) are stored */
/* > in the columns of VR, in the same order as their eigenvalues. */
/* > Each eigenvector is scaled so that its largest component has */
/* > f2c_abs(real part) + f2c_abs(imag. part) = 1, except for eigenvectors */
/* > corresponding to an eigenvalue with alpha = beta = 0, which */
/* > are set to zero. */
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
/* > To compute the optimal value of LWORK, call ILAENV to get */
/* > blocksizes (for CGEQRF, CUNMQR, and CUNGQR.) Then compute: */
/* > NB -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR;
 */
/* > The optimal LWORK is MAX( 2*N, N*(NB+1) ). */
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
/* > > N: errors that usually indicate LAPACK problems: */
/* > =N+1: error return from CGGBAL */
/* > =N+2: error return from CGEQRF */
/* > =N+3: error return from CUNMQR */
/* > =N+4: error return from CUNGQR */
/* > =N+5: error return from CGGHRD */
/* > =N+6: error return from CHGEQZ (other than failed */
/* > iteration) */
/* > =N+7: error return from CTGEVC */
/* > =N+8: error return from CGGBAK (computing VL) */
/* > =N+9: error return from CGGBAK (computing VR) */
/* > =N+10: error return from CLASCL (various calls) */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexGEeigen */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Balancing */
/* > --------- */
/* > */
/* > This driver calls CGGBAL to both permute and scale rows and columns */
/* > of A and B. The permutations PL and PR are chosen so that PL*A*PR */
/* > and PL*B*R will be upper triangular except for the diagonal blocks */
/* > A(i:j,i:j) and B(i:j,i:j), with i and j as close together as */
/* > possible. The diagonal scaling matrices DL and DR are chosen so */
/* > that the pair DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to */
/* > one (except for the elements that start out zero.) */
/* > */
/* > After the eigenvalues and eigenvectors of the balanced matrices */
/* > have been computed, CGGBAK transforms the eigenvectors back to what */
/* > they would have been (in perfect arithmetic) if they had not been */
/* > balanced. */
/* > */
/* > Contents of A and B on Exit */
/* > -------- -- - --- - -- ---- */
/* > */
/* > If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or */
/* > both), then on exit the arrays A and B will contain the scomplex Schur */
/* > form[*] of the "balanced" versions of A and B. If no eigenvectors */
/* > are computed, then only the diagonal blocks will be correct. */
/* > */
/* > [*] In other words, upper triangular form. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgegv_(char *jobvl, char *jobvr, aocl_int_t *n, scomplex *a, aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, aocl_int_t *ldvl, scomplex *vr, aocl_int_t *ldvr, scomplex *work, aocl_int_t *lwork, real *rwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgegv(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t ldvl_64 = *ldvl;
    aocl_int64_t ldvr_64 = *ldvr;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgegv(jobvl, jobvr, &n_64, a, &lda_64, b, &ldb_64, alpha, beta, vl, &ldvl_64, vr, &ldvr_64, work, &lwork_64, rwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgegv(char *jobvl, char *jobvr, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, scomplex *b,
            aocl_int64_t *ldb, scomplex *alpha, scomplex *beta, scomplex *vl, aocl_int64_t *ldvl,
            scomplex *vr, aocl_int64_t *ldvr, scomplex *work, aocl_int64_t *lwork, real *rwork,
            aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cgegv inputs: jobvl %c, jobvr %c, n %lld, lda %lld, ldb %lld, ldvl %lld, ldvr %lld, "
             "lwork %lld",
             *jobvl, *jobvr, *n, *lda, *ldb, *ldvl, *ldvr, *lwork);
#else
    snprintf(buffer, 256,
             "cgegv inputs: jobvl %c, jobvr %c, n %d, lda %d, ldb %d, ldvl %d, ldvr %d, lwork %d",
             *jobvl, *jobvr, *n, *lda, *ldb, *ldvl, *ldvr, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1,
        i__2, i__3, i__4;
    real r__1, r__2, r__3, r__4;
    scomplex q__1, q__2;
    /* Builtin functions */
    double r_imag(scomplex *);
    /* Local variables */
    aocl_int64_t jc, nb, in, jr, nb1, nb2, nb3, ihi, ilo;
    real eps;
    logical ilv;
    real absb, anrm, bnrm;
    aocl_int64_t itau;
    real temp;
    logical ilvl, ilvr;
    aocl_int64_t lopt;
    real anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t ileft, iinfo, icols, iwork, irows;
    real salfai;
    real salfar;
    extern real slamch_(char *);
    real safmin;
    real safmax;
    char chtemp[1];
    logical ldumma[1];
    aocl_int64_t ijobvl, iright;
    logical ilimit;
    aocl_int64_t ijobvr;
    aocl_int64_t lwkmin;
    aocl_int64_t irwork, lwkopt;
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
    /* Computing MAX */
    i__1 = *n << 1;
    lwkmin = fla_max(i__1, 1);
    lwkopt = lwkmin;
    work[1].real = (real)lwkopt;
    work[1].imag = 0.f; // , expr subst
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
    else if(*ldvl < 1 || ilvl && *ldvl < *n)
    {
        *info = -11;
    }
    else if(*ldvr < 1 || ilvr && *ldvr < *n)
    {
        *info = -13;
    }
    else if(*lwork < lwkmin && !lquery)
    {
        *info = -15;
    }
    if(*info == 0)
    {
        nb1 = aocl_lapack_ilaenv(&c__1, "CGEQRF", " ", n, n, &c_n1, &c_n1);
        nb2 = aocl_lapack_ilaenv(&c__1, "CUNMQR", " ", n, n, n, &c_n1);
        nb3 = aocl_lapack_ilaenv(&c__1, "CUNGQR", " ", n, n, n, &c_n1);
        /* Computing MAX */
        i__1 = fla_max(nb1, nb2);
        nb = fla_max(i__1, nb3);
        /* Computing MAX */
        i__1 = *n << 1;
        i__2 = *n * (nb + 1); // , expr subst
        lopt = fla_max(i__1, i__2);
        work[1].real = (real)lopt;
        work[1].imag = 0.f; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGEGV ", &i__1, (ftnlen)6);
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
    safmin = slamch_("S");
    safmin += safmin;
    safmax = 1.f / safmin;
    /* Scale A */
    anrm = aocl_lapack_clange("M", n, n, &a[a_offset], lda, &rwork[1]);
    anrm1 = anrm;
    anrm2 = 1.f;
    if(anrm < 1.f)
    {
        if(safmax * anrm < 1.f)
        {
            anrm1 = safmin;
            anrm2 = safmax * anrm;
        }
    }
    if(anrm > 0.f)
    {
        aocl_lapack_clascl("G", &c_n1, &c_n1, &anrm, &c_b29, n, n, &a[a_offset], lda, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 10;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    /* Scale B */
    bnrm = aocl_lapack_clange("M", n, n, &b[b_offset], ldb, &rwork[1]);
    bnrm1 = bnrm;
    bnrm2 = 1.f;
    if(bnrm < 1.f)
    {
        if(safmax * bnrm < 1.f)
        {
            bnrm1 = safmin;
            bnrm2 = safmax * bnrm;
        }
    }
    if(bnrm > 0.f)
    {
        aocl_lapack_clascl("G", &c_n1, &c_n1, &bnrm, &c_b29, n, n, &b[b_offset], ldb, &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 10;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    /* Permute the matrix to make it more nearly triangular */
    /* Also "balance" the matrix. */
    ileft = 1;
    iright = *n + 1;
    irwork = iright + *n;
    aocl_lapack_cggbal("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[ileft],
                       &rwork[iright], &rwork[irwork], &iinfo);
    if(iinfo != 0)
    {
        *info = *n + 1;
        goto L80;
    }
    /* Reduce B to triangular form, and initialize VL and/or VR */
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
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    aocl_lapack_cgeqrf(&irows, &icols, &b[ilo + ilo * b_dim1], ldb, &work[itau], &work[iwork],
                       &i__1, &iinfo);
    if(iinfo >= 0)
    {
        /* Computing MAX */
        i__3 = iwork;
        i__1 = lwkopt;
        i__2 = (integer)work[i__3].real + iwork - 1; // , expr subst
        lwkopt = fla_max(i__1, i__2);
    }
    if(iinfo != 0)
    {
        *info = *n + 2;
        goto L80;
    }
    i__1 = *lwork + 1 - iwork;
    aocl_lapack_cunmqr("L", "C", &irows, &icols, &irows, &b[ilo + ilo * b_dim1], ldb, &work[itau],
                       &a[ilo + ilo * a_dim1], lda, &work[iwork], &i__1, &iinfo);
    if(iinfo >= 0)
    {
        /* Computing MAX */
        i__3 = iwork;
        i__1 = lwkopt;
        i__2 = (integer)work[i__3].real + iwork - 1; // , expr subst
        lwkopt = fla_max(i__1, i__2);
    }
    if(iinfo != 0)
    {
        *info = *n + 3;
        goto L80;
    }
    if(ilvl)
    {
        aocl_lapack_claset("Full", n, n, &c_b1, &c_b2, &vl[vl_offset], ldvl);
        i__1 = irows - 1;
        i__2 = irows - 1;
        aocl_lapack_clacpy("L", &i__1, &i__2, &b[ilo + 1 + ilo * b_dim1], ldb,
                           &vl[ilo + 1 + ilo * vl_dim1], ldvl);
        i__1 = *lwork + 1 - iwork;
        aocl_lapack_cungqr(&irows, &irows, &irows, &vl[ilo + ilo * vl_dim1], ldvl, &work[itau],
                           &work[iwork], &i__1, &iinfo);
        if(iinfo >= 0)
        {
            /* Computing MAX */
            i__3 = iwork;
            i__1 = lwkopt;
            i__2 = (integer)work[i__3].real + iwork - 1; // , expr subst
            lwkopt = fla_max(i__1, i__2);
        }
        if(iinfo != 0)
        {
            *info = *n + 4;
            goto L80;
        }
    }
    if(ilvr)
    {
        aocl_lapack_claset("Full", n, n, &c_b1, &c_b2, &vr[vr_offset], ldvr);
    }
    /* Reduce to generalized Hessenberg form */
    if(ilv)
    {
        /* Eigenvectors requested -- work on whole matrix. */
        aocl_lapack_cgghrd(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb,
                           &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &iinfo);
    }
    else
    {
        aocl_lapack_cgghrd("N", "N", &irows, &c__1, &irows, &a[ilo + ilo * a_dim1], lda,
                           &b[ilo + ilo * b_dim1], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr,
                           &iinfo);
    }
    if(iinfo != 0)
    {
        *info = *n + 5;
        goto L80;
    }
    /* Perform QZ algorithm */
    iwork = itau;
    if(ilv)
    {
        *(unsigned char *)chtemp = 'S';
    }
    else
    {
        *(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwork;
    aocl_lapack_chgeqz(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], ldb,
                       &alpha[1], &beta[1], &vl[vl_offset], ldvl, &vr[vr_offset], ldvr,
                       &work[iwork], &i__1, &rwork[irwork], &iinfo);
    if(iinfo >= 0)
    {
        /* Computing MAX */
        i__3 = iwork;
        i__1 = lwkopt;
        i__2 = (integer)work[i__3].real + iwork - 1; // , expr subst
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
        goto L80;
    }
    if(ilv)
    {
        /* Compute Eigenvectors */
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
        aocl_lapack_ctgevc(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb,
                           &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[iwork],
                           &rwork[irwork], &iinfo);
        if(iinfo != 0)
        {
            *info = *n + 7;
            goto L80;
        }
        /* Undo balancing on VL and VR, rescale */
        if(ilvl)
        {
            aocl_lapack_cggbak("P", "L", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
                               &vl[vl_offset], ldvl, &iinfo);
            if(iinfo != 0)
            {
                *info = *n + 8;
                goto L80;
            }
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
                    r__4 = (r__1 = vl[i__3].real, f2c_abs(r__1))
                           + (r__2 = r_imag(&vl[jr + jc * vl_dim1]), f2c_abs(r__2)); // , expr subst
                    temp = fla_max(r__3, r__4);
                    /* L10: */
                }
                if(temp < safmin)
                {
                    goto L30;
                }
                temp = 1.f / temp;
                i__2 = *n;
                for(jr = 1; jr <= i__2; ++jr)
                {
                    i__3 = jr + jc * vl_dim1;
                    i__4 = jr + jc * vl_dim1;
                    q__1.real = temp * vl[i__4].real;
                    q__1.imag = temp * vl[i__4].imag; // , expr subst
                    vl[i__3].real = q__1.real;
                    vl[i__3].imag = q__1.imag; // , expr subst
                    /* L20: */
                }
            L30:;
            }
        }
        if(ilvr)
        {
            aocl_lapack_cggbak("P", "R", n, &ilo, &ihi, &rwork[ileft], &rwork[iright], n,
                               &vr[vr_offset], ldvr, &iinfo);
            if(iinfo != 0)
            {
                *info = *n + 9;
                goto L80;
            }
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
                    r__4 = (r__1 = vr[i__3].real, f2c_abs(r__1))
                           + (r__2 = r_imag(&vr[jr + jc * vr_dim1]), f2c_abs(r__2)); // , expr subst
                    temp = fla_max(r__3, r__4);
                    /* L40: */
                }
                if(temp < safmin)
                {
                    goto L60;
                }
                temp = 1.f / temp;
                i__2 = *n;
                for(jr = 1; jr <= i__2; ++jr)
                {
                    i__3 = jr + jc * vr_dim1;
                    i__4 = jr + jc * vr_dim1;
                    q__1.real = temp * vr[i__4].real;
                    q__1.imag = temp * vr[i__4].imag; // , expr subst
                    vr[i__3].real = q__1.real;
                    vr[i__3].imag = q__1.imag; // , expr subst
                    /* L50: */
                }
            L60:;
            }
        }
        /* End of eigenvector calculation */
    }
    /* Undo scaling in alpha, beta */
    /* Note: this does not give the alpha and beta for the unscaled */
    /* problem. */
    /* Un-scaling is limited to avoid underflow in alpha and beta */
    /* if they are significant. */
    i__1 = *n;
    for(jc = 1; jc <= i__1; ++jc)
    {
        i__2 = jc;
        absar = (r__1 = alpha[i__2].real, f2c_abs(r__1));
        absai = (r__1 = r_imag(&alpha[jc]), f2c_abs(r__1));
        i__2 = jc;
        absb = (r__1 = beta[i__2].real, f2c_abs(r__1));
        i__2 = jc;
        salfar = anrm * alpha[i__2].real;
        salfai = anrm * r_imag(&alpha[jc]);
        i__2 = jc;
        sbeta = bnrm * beta[i__2].real;
        ilimit = FALSE_;
        scale = 1.f;
        /* Check for significant underflow in imaginary part of ALPHA */
        /* Computing MAX */
        r__1 = safmin, r__2 = eps * absar;
        r__1 = fla_max(r__1, r__2);
        r__2 = eps * absb; // ; expr subst
        if(f2c_abs(salfai) < safmin && absai >= fla_max(r__1, r__2))
        {
            ilimit = TRUE_;
            /* Computing MAX */
            r__1 = safmin;
            r__2 = anrm2 * absai; // , expr subst
            scale = safmin / anrm1 / fla_max(r__1, r__2);
        }
        /* Check for significant underflow in real part of ALPHA */
        /* Computing MAX */
        r__1 = safmin, r__2 = eps * absai;
        r__1 = fla_max(r__1, r__2);
        r__2 = eps * absb; // ; expr subst
        if(f2c_abs(salfar) < safmin && absar >= fla_max(r__1, r__2))
        {
            ilimit = TRUE_;
            /* Computing MAX */
            /* Computing MAX */
            r__3 = safmin;
            r__4 = anrm2 * absar; // , expr subst
            r__1 = scale;
            r__2 = safmin / anrm1 / fla_max(r__3, r__4); // , expr subst
            scale = fla_max(r__1, r__2);
        }
        /* Check for significant underflow in BETA */
        /* Computing MAX */
        r__1 = safmin, r__2 = eps * absar;
        r__1 = fla_max(r__1, r__2);
        r__2 = eps * absai; // ; expr subst
        if(f2c_abs(sbeta) < safmin && absb >= fla_max(r__1, r__2))
        {
            ilimit = TRUE_;
            /* Computing MAX */
            /* Computing MAX */
            r__3 = safmin;
            r__4 = bnrm2 * absb; // , expr subst
            r__1 = scale;
            r__2 = safmin / bnrm1 / fla_max(r__3, r__4); // , expr subst
            scale = fla_max(r__1, r__2);
        }
        /* Check for possible overflow when limiting scaling */
        if(ilimit)
        {
            /* Computing MAX */
            r__1 = f2c_abs(salfar), r__2 = f2c_abs(salfai);
            r__1 = fla_max(r__1, r__2);
            r__2 = f2c_abs(sbeta); // ; expr subst
            temp = scale * safmin * fla_max(r__1, r__2);
            if(temp > 1.f)
            {
                scale /= temp;
            }
            if(scale < 1.f)
            {
                ilimit = FALSE_;
            }
        }
        /* Recompute un-scaled ALPHA, BETA if necessary. */
        if(ilimit)
        {
            i__2 = jc;
            salfar = scale * alpha[i__2].real * anrm;
            salfai = scale * r_imag(&alpha[jc]) * anrm;
            i__2 = jc;
            q__2.real = scale * beta[i__2].real;
            q__2.imag = scale * beta[i__2].imag; // , expr subst
            q__1.real = bnrm * q__2.real;
            q__1.imag = bnrm * q__2.imag; // , expr subst
            sbeta = q__1.real;
        }
        i__2 = jc;
        q__1.real = salfar;
        q__1.imag = salfai; // , expr subst
        alpha[i__2].real = q__1.real;
        alpha[i__2].imag = q__1.imag; // , expr subst
        i__2 = jc;
        beta[i__2].real = sbeta;
        beta[i__2].imag = 0.f; // , expr subst
        /* L70: */
    }
L80:
    work[1].real = (real)lwkopt;
    work[1].imag = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGEGV */
}
/* cgegv_ */
