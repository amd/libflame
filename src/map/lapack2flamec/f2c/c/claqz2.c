/* ./claqz2.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {0.f, 0.f};
static complex c_b2 = {1.f, 0.f};
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;
/* > \brief \b CLAQZ2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLAQZ2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/CLAQZ2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/CLAQZ2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/CLAQZ2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLAQZ2( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, */
/* $ LDB, Q, LDQ, Z, LDZ, NS, ND, ALPHA, BETA, QC, LDQC, ZC, LDZC, */
/* $ WORK, LWORK, RWORK, REC, INFO ) */
/* IMPLICIT NONE */
/* Arguments */
/* LOGICAL ILSCHUR, ILQ, ILZ */
/* INTEGER N, ILO, IHI, NW, LDA, LDB, LDQ, LDZ, */
/* $ LDQC, LDZC, LWORK, REC */
/* COMPLEX A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ Z( LDZ, * ), ALPHA( * ), BETA( * ) */
/* INTEGER, INTENT( OUT ) NS, ND, INFO */
/* COMPLEX QC( LDQC, * ), ZC( LDZC, * ), WORK( * ) */
/* REAL RWORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLAQZ2 performs AED */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] ILSCHUR */
/* > \verbatim */
/* > ILSCHUR is LOGICAL */
/* > Determines whether or not to update the full Schur form */
/* > \endverbatim */
/* > */
/* > \param[in] ILQ */
/* > \verbatim */
/* > ILQ is LOGICAL */
/* > Determines whether or not to update the matrix Q */
/* > \endverbatim */
/* > */
/* > \param[in] ILZ */
/* > \verbatim */
/* > ILZ is LOGICAL */
/* > Determines whether or not to update the matrix Z */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrices A, B, Q, and Z. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ILO */
/* > \verbatim */
/* > ILO is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] IHI */
/* > \verbatim */
/* > IHI is INTEGER */
/* > ILO and IHI mark the rows and columns of (A,B) which */
/* > are to be normalized */
/* > \endverbatim */
/* > */
/* > \param[in] NW */
/* > \verbatim */
/* > NW is INTEGER */
/* > The desired size of the deflation window. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is COMPLEX array, dimension (LDZ, N) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] NS */
/* > \verbatim */
/* > NS is INTEGER */
/* > The number of unconverged eigenvalues available to */
/* > use as shifts. */
/* > \endverbatim */
/* > */
/* > \param[out] ND */
/* > \verbatim */
/* > ND is INTEGER */
/* > The number of converged eigenvalues found. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX array, dimension (N) */
/* > Each scalar alpha defining an eigenvalue */
/* > of GNEP. */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is COMPLEX array, dimension (N) */
/* > The scalars beta that define the eigenvalues of GNEP. */
/* > Together, the quantities alpha = ALPHA(j) and */
/* > beta = BETA(j) represent the j-th eigenvalue of the matrix */
/* > pair (A,B), in one of the forms lambda = alpha/beta or */
/* > mu = beta/alpha. Since either lambda or mu may overflow, */
/* > they should not, in general, be computed. */
/* > \endverbatim */
/* > */
/* > \param[in,out] QC */
/* > \verbatim */
/* > QC is COMPLEX array, dimension (LDQC, NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDQC */
/* > \verbatim */
/* > LDQC is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in,out] ZC */
/* > \verbatim */
/* > ZC is COMPLEX array, dimension (LDZC, NW) */
/* > \endverbatim */
/* > */
/* > \param[in] LDZC */
/* > \verbatim */
/* > LDZ is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO >= 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
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
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[in] REC */
/* > \verbatim */
/* > REC is INTEGER */
/* > REC indicates the current recursion level. Should be set */
/* > to 0 on first call. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Thijs Steel, KU Leuven, KU Leuven */
/* > \date May 2020 */
/* > \ingroup laqz2 */
/* > */
/* ===================================================================== */
/* Subroutine */
void claqz2_(logical *ilschur, logical *ilq, logical *ilz, integer *n, integer *ilo, integer *ihi,
             integer *nw, complex *a, integer *lda, complex *b, integer *ldb, complex *q,
             integer *ldq, complex *z__, integer *ldz, integer *ns, integer *nd, complex *alpha,
             complex *beta, complex *qc, integer *ldqc, complex *zc, integer *ldzc, complex *work,
             integer *lwork, real *rwork, integer *rec, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF(
        "claqz2 inputs: n %" FLA_IS ",ilo %" FLA_IS ",ihi %" FLA_IS ",nw %" FLA_IS ",lda %" FLA_IS
        ",ldb %" FLA_IS ",ldq %" FLA_IS ",ldz %" FLA_IS
        ",ldqc %" FLA_IS ",ldzc %" FLA_IS ",lwork %" FLA_IS ",rec %" FLA_IS "",
        *n, *ilo, *ihi, *nw, *lda, *ldb, *ldq, *ldz, *ldqc, *ldzc, *lwork, *rec);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, z_offset, qc_dim1,
        qc_offset, zc_dim1, zc_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double c_abs(complex *);
    void r_cnjg(complex *, complex *);
    /* Local variables */
    integer lworkreq, i__, j, k;
    complex s;
    real c1;
    integer k2;
    complex s1;
    integer jw;
    real ulp;
    integer ctgexc_info__, ifst;
    complex temp;
    extern /* Subroutine */
        void
        crot_(integer *, complex *, integer *, complex *, integer *, real *, complex *);
    integer ilst;
    extern /* Subroutine */
        void
        cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, complex *, integer *);
    complex atemp;           
    integer kwbot;
    real tempr;
    integer kwtop;
    extern /* Subroutine */
        void
        claqz0_(char *, char *, char *, integer *, integer *, integer *, complex *, integer *,
                complex *, integer *, complex *, complex *, complex *, integer *, complex *,
                integer *, complex *, integer *, real *, integer *, integer *);
    integer qz_small_info__;
    extern /* Subroutine */
        void
        claqz1_(logical *, logical *, integer *, integer *, integer *, integer *, complex *,
                integer *, complex *, integer *, integer *, integer *, complex *, integer *,
                integer *, integer *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *),
        claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *);
    real safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        ctgexc_(logical *, logical *, integer *, complex *, integer *, complex *, integer *,
                complex *, integer *, complex *, integer *, integer *, integer *, integer *),
        clartg_(complex *, complex *, real *, complex *, complex *);
    integer istopm;
    real smlnum;
    integer istartm;
    /* Arguments */
    /* Parameters */
    /* Local Scalars */
    /* External Functions */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --alpha;
    --beta;
    qc_dim1 = *ldqc;
    qc_offset = 1 + qc_dim1;
    qc -= qc_offset;
    zc_dim1 = *ldzc;
    zc_offset = 1 + zc_dim1;
    zc -= zc_offset;
    --work;
    --rwork;
    /* Function Body */
    *info = 0;
    /* Set up deflation window */
    /* Computing MIN */
    i__1 = *nw;
    i__2 = *ihi - *ilo + 1; // , expr subst
    jw = fla_min(i__1, i__2);
    kwtop = *ihi - jw + 1;
    if(kwtop == *ilo)
    {
        s.r = 0.f;
        s.i = 0.f; // , expr subst
    }
    else
    {
        i__1 = kwtop + (kwtop - 1) * a_dim1;
        s.r = a[i__1].r;
        s.i = a[i__1].i; // , expr subst
    }
    /* Determine required workspace */
    ifst = 1;
    ilst = jw;
    i__1 = *rec + 1;
    claqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda,
            &b[kwtop + kwtop * b_dim1], ldb, &alpha[1], &beta[1], &qc[qc_offset], ldqc,
            &zc[zc_offset], ldzc, &work[1], &c_n1, &rwork[1], &i__1, &qz_small_info__);
    /* Computing 2nd power */
    i__1 = jw;
    lworkreq = (integer)work[1].r + (i__1 * i__1 << 1);
    /* Computing MAX */
    /* Computing 2nd power */
    i__3 = *nw;
    i__1 = lworkreq, i__2 = *n * *nw;
    i__1 = fla_max(i__1, i__2);
    i__2 = (i__3 * i__3 << 1) + *n; // ; expr subst
    lworkreq = fla_max(i__1, i__2);
    if(*lwork == -1)
    {
        /* workspace query, quick return */
        work[1].r = (real)lworkreq;
        work[1].i = 0.f; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(*lwork < lworkreq)
    {
        *info = -26;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CLAQZ2", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants */
    safmin = slamch_("SAFE MINIMUM");
    ulp = slamch_("PRECISION");
    smlnum = safmin * ((real)(*n) / ulp);
    if(*ihi == kwtop)
    {
        /* 1 by 1 deflation window, just try a regular deflation */
        i__1 = kwtop;
        i__2 = kwtop + kwtop * a_dim1;
        alpha[i__1].r = a[i__2].r;
        alpha[i__1].i = a[i__2].i; // , expr subst
        i__1 = kwtop;
        i__2 = kwtop + kwtop * b_dim1;
        beta[i__1].r = b[i__2].r;
        beta[i__1].i = b[i__2].i; // , expr subst
        *ns = 1;
        *nd = 0;
        /* Computing MAX */
        r__1 = smlnum;
        r__2 = ulp * c_abs(&a[kwtop + kwtop * a_dim1]); // , expr subst
        if(c_abs(&s) <= fla_max(r__1, r__2))
        {
            *ns = 0;
            *nd = 1;
            if(kwtop > *ilo)
            {
                i__1 = kwtop + (kwtop - 1) * a_dim1;
                a[i__1].r = 0.f;
                a[i__1].i = 0.f; // , expr subst
            }
        }
    }
    /* Store window in case of convergence failure */
    clacpy_("ALL", &jw, &jw, &a[kwtop + kwtop * a_dim1], lda, &work[1], &jw);
    /* Computing 2nd power */
    i__1 = jw;
    clacpy_("ALL", &jw, &jw, &b[kwtop + kwtop * b_dim1], ldb, &work[i__1 * i__1 + 1], &jw);
    /* Transform window to real schur form */
    claset_("FULL", &jw, &jw, &c_b1, &c_b2, &qc[qc_offset], ldqc);
    claset_("FULL", &jw, &jw, &c_b1, &c_b2, &zc[zc_offset], ldzc);
    /* Computing 2nd power */
    i__1 = jw;
    /* Computing 2nd power */
    i__3 = jw;
    i__2 = *lwork - (i__3 * i__3 << 1);
    i__4 = *rec + 1;
    claqz0_("S", "V", "V", &jw, &c__1, &jw, &a[kwtop + kwtop * a_dim1], lda,
            &b[kwtop + kwtop * b_dim1], ldb, &alpha[1], &beta[1], &qc[qc_offset], ldqc,
            &zc[zc_offset], ldzc, &work[(i__1 * i__1 << 1) + 1], &i__2, &rwork[1], &i__4,
            &qz_small_info__);
    if(qz_small_info__ != 0)
    {
        /* Convergence failure, restore the window and exit */
        *nd = 0;
        *ns = jw - qz_small_info__;
        clacpy_("ALL", &jw, &jw, &work[1], &jw, &a[kwtop + kwtop * a_dim1], lda);
        /* Computing 2nd power */
        i__1 = jw;
        clacpy_("ALL", &jw, &jw, &work[i__1 * i__1 + 1], &jw, &b[kwtop + kwtop * b_dim1], ldb);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Deflation detection loop */
    if(kwtop == *ilo || s.r == 0.f && s.i == 0.f)
    {
        kwbot = kwtop - 1;
    }
    else
    {
        kwbot = *ihi;
        k = 1;
        k2 = 1;
        while(k <= jw)
        {
            /* Try to deflate eigenvalue */
            tempr = c_abs(&a[kwbot + kwbot * a_dim1]);
            if(tempr == 0.f)
            {
                tempr = c_abs(&s);
            }
            i__1 = (kwbot - kwtop + 1) * qc_dim1 + 1;
            q__1.r = s.r * qc[i__1].r - s.i * qc[i__1].i;
            q__1.i = s.r * qc[i__1].i + s.i * qc[i__1].r; // , expr subst
            /* Computing MAX */
            r__1 = ulp * tempr;
            if(c_abs(&q__1) <= fla_max(r__1, smlnum))
            {
                /* Deflatable */
                --kwbot;
            }
            else
            {
                /* Not deflatable, move out of the way */
                ifst = kwbot - kwtop + 1;
                ilst = k2;
                ctgexc_(&c_true, &c_true, &jw, &a[kwtop + kwtop * a_dim1], lda,
                        &b[kwtop + kwtop * b_dim1], ldb, &qc[qc_offset], ldqc, &zc[zc_offset], ldzc,
                        &ifst, &ilst, &ctgexc_info__);
                ++k2;
            }
            ++k;
        }
    }
    /* Store eigenvalues */
    *nd = *ihi - kwbot;
    *ns = jw - *nd;
    k = kwtop;
    while(k <= *ihi)
    {
        i__1 = k;
        i__2 = k + k * a_dim1;
        alpha[i__1].r = a[i__2].r;
        alpha[i__1].i = a[i__2].i; // , expr subst
        i__1 = k;
        i__2 = k + k * b_dim1;
        beta[i__1].r = b[i__2].r;
        beta[i__1].i = b[i__2].i; // , expr subst
        ++k;
    }
    if(kwtop != *ilo && (s.r != 0.f || s.i != 0.f))
    {
        /* Reflect spike back, this will create optimally packed bulges */
        /* A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 ) *CONJG( QC( 1, */
        /* $ 1:JW-ND ) ) */
        i__1 = kwtop + (kwtop - 1) * a_dim1;
        atemp.r = a[i__1].r;
        atemp.i = a[i__1].i; // , expr subst
        j = 1;
        i__1 = kwbot;
        for(i__ = kwtop; i__ <= i__1; ++i__)
        {
            i__2 = i__ + (kwtop - 1) * a_dim1;
            r_cnjg(&q__2, &qc[j * qc_dim1 + 1]);
            q__1.r = atemp.r * q__2.r - atemp.i * q__2.i;
            q__1.i = atemp.r * q__2.i + atemp.i * q__2.r; // , expr subst
            a[i__2].r = q__1.r;
            a[i__2].i = q__1.i; // , expr subst
            ++j;
        }
        i__1 = kwtop;
        for(k = kwbot - 1; k >= i__1; --k)
        {
            clartg_(&a[k + (kwtop - 1) * a_dim1], &a[k + 1 + (kwtop - 1) * a_dim1], &c1, &s1,
                    &temp);
            i__2 = k + (kwtop - 1) * a_dim1;
            a[i__2].r = temp.r;
            a[i__2].i = temp.i; // , expr subst
            i__2 = k + 1 + (kwtop - 1) * a_dim1;
            a[i__2].r = 0.f;
            a[i__2].i = 0.f; // , expr subst
            /* Computing MAX */
            i__2 = kwtop;
            i__3 = k - 1; // , expr subst
            k2 = fla_max(i__2, i__3);
            i__2 = *ihi - k2 + 1;
            crot_(&i__2, &a[k + k2 * a_dim1], lda, &a[k + 1 + k2 * a_dim1], lda, &c1, &s1);
            i__2 = *ihi - (k - 1) + 1;
            crot_(&i__2, &b[k + (k - 1) * b_dim1], ldb, &b[k + 1 + (k - 1) * b_dim1], ldb, &c1,
                  &s1);
            r_cnjg(&q__1, &s1);
            crot_(&jw, &qc[(k - kwtop + 1) * qc_dim1 + 1], &c__1,
                  &qc[(k + 1 - kwtop + 1) * qc_dim1 + 1], &c__1, &c1, &q__1);
        }
        /* Chase bulges down */
        istartm = kwtop;
        istopm = *ihi;
        k = kwbot - 1;
        while(k >= kwtop)
        {
            /* Move bulge down and remove it */
            i__1 = kwbot - 1;
            for(k2 = k; k2 <= i__1; ++k2)
            {
                i__2 = kwtop + jw - 1;
                claqz1_(&c_true, &c_true, &k2, &kwtop, &i__2, &kwbot, &a[a_offset], lda,
                        &b[b_offset], ldb, &jw, &kwtop, &qc[qc_offset], ldqc, &jw, &kwtop,
                        &zc[zc_offset], ldzc);
            }
            --k;
        }
    }
    /* Apply Qc and Zc to rest of the matrix */
    if(*ilschur)
    {
        istartm = 1;
        istopm = *n;
    }
    else
    {
        istartm = *ilo;
        istopm = *ihi;
    }
    if(istopm - *ihi > 0)
    {
        i__1 = istopm - *ihi;
        cgemm_("C", "N", &jw, &i__1, &jw, &c_b2, &qc[qc_offset], ldqc,
               &a[kwtop + (*ihi + 1) * a_dim1], lda, &c_b1, &work[1], &jw);
        i__1 = istopm - *ihi;
        clacpy_("ALL", &jw, &i__1, &work[1], &jw, &a[kwtop + (*ihi + 1) * a_dim1], lda);
        i__1 = istopm - *ihi;
        cgemm_("C", "N", &jw, &i__1, &jw, &c_b2, &qc[qc_offset], ldqc,
               &b[kwtop + (*ihi + 1) * b_dim1], ldb, &c_b1, &work[1], &jw);
        i__1 = istopm - *ihi;
        clacpy_("ALL", &jw, &i__1, &work[1], &jw, &b[kwtop + (*ihi + 1) * b_dim1], ldb);
    }
    if(*ilq)
    {
        cgemm_("N", "N", n, &jw, &jw, &c_b2, &q[kwtop * q_dim1 + 1], ldq, &qc[qc_offset], ldqc,
               &c_b1, &work[1], n);
        clacpy_("ALL", n, &jw, &work[1], n, &q[kwtop * q_dim1 + 1], ldq);
    }
    if(kwtop - 1 - istartm + 1 > 0)
    {
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        cgemm_("N", "N", &i__1, &jw, &jw, &c_b2, &a[istartm + kwtop * a_dim1], lda, &zc[zc_offset],
               ldzc, &c_b1, &work[1], &i__2);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        clacpy_("ALL", &i__1, &jw, &work[1], &i__2, &a[istartm + kwtop * a_dim1], lda);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        cgemm_("N", "N", &i__1, &jw, &jw, &c_b2, &b[istartm + kwtop * b_dim1], ldb, &zc[zc_offset],
               ldzc, &c_b1, &work[1], &i__2);
        i__1 = kwtop - istartm;
        i__2 = kwtop - istartm;
        clacpy_("ALL", &i__1, &jw, &work[1], &i__2, &b[istartm + kwtop * b_dim1], ldb);
    }
    if(*ilz)
    {
        cgemm_("N", "N", n, &jw, &jw, &c_b2, &z__[kwtop * z_dim1 + 1], ldz, &zc[zc_offset], ldzc,
               &c_b1, &work[1], n);
        clacpy_("ALL", n, &jw, &work[1], n, &z__[kwtop * z_dim1 + 1], ldz);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* claqz2_ */
