/* ../netlib/v3.9.0/cggsvp3.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {0.f, 0.f};
static complex c_b2 = {1.f, 0.f};
static integer c_n1 = -1;
/* > \brief \b CGGSVP3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGGSVP3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggsvp3
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggsvp3
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggsvp3
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/* TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/* IWORK, RWORK, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBQ, JOBU, JOBV */
/* INTEGER INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK */
/* REAL TOLA, TOLB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL RWORK( * ) */
/* COMPLEX A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGSVP3 computes unitary matrices U, V and Q such that */
/* > */
/* > N-K-L K L */
/* > U**H*A*Q = K ( 0 A12 A13 ) if M-K-L >= 0;
 */
/* > L ( 0 0 A23 ) */
/* > M-K-L ( 0 0 0 ) */
/* > */
/* > N-K-L K L */
/* > = K ( 0 A12 A13 ) if M-K-L < 0;
 */
/* > M-K ( 0 0 A23 ) */
/* > */
/* > N-K-L K L */
/* > V**H*B*Q = L ( 0 0 B13 ) */
/* > P-L ( 0 0 0 ) */
/* > */
/* > where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular */
/* > upper triangular;
A23 is L-by-L upper triangular if M-K-L >= 0, */
/* > otherwise A23 is (M-K)-by-L upper trapezoidal. K+L = the effective */
/* > numerical rank of the (M+P)-by-N matrix (A**H,B**H)**H. */
/* > */
/* > This decomposition is the preprocessing step for computing the */
/* > Generalized Singular Value Decomposition (GSVD), see subroutine */
/* > CGGSVD3. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU */
/* > \verbatim */
/* > JOBU is CHARACTER*1 */
/* > = 'U': Unitary matrix U is computed;
 */
/* > = 'N': U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* > JOBV is CHARACTER*1 */
/* > = 'V': Unitary matrix V is computed;
 */
/* > = 'N': V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* > JOBQ is CHARACTER*1 */
/* > = 'Q': Unitary matrix Q is computed;
 */
/* > = 'N': Q is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows of the matrix B. P >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, A contains the triangular (or trapezoidal) matrix */
/* > described in the Purpose section. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,N) */
/* > On entry, the P-by-N matrix B. */
/* > On exit, B contains the triangular matrix described in */
/* > the Purpose section. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in] TOLA */
/* > \verbatim */
/* > TOLA is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] TOLB */
/* > \verbatim */
/* > TOLB is REAL */
/* > */
/* > TOLA and TOLB are the thresholds to determine the effective */
/* > numerical rank of matrix B and a subblock of A. Generally, */
/* > they are set to */
/* > TOLA = MAX(M,N)*norm(A)*MACHEPS, */
/* > TOLB = MAX(P,N)*norm(B)*MACHEPS. */
/* > The size of TOLA and TOLB may affect the size of backward */
/* > errors of the decomposition. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* > K is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] L */
/* > \verbatim */
/* > L is INTEGER */
/* > */
/* > On exit, K and L specify the dimension of the subblocks */
/* > described in Purpose section. */
/* > K + L = effective numerical rank of (A**H,B**H)**H. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* > U is COMPLEX array, dimension (LDU,M) */
/* > If JOBU = 'U', U contains the unitary matrix U. */
/* > If JOBU = 'N', U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* > LDU is INTEGER */
/* > The leading dimension of the array U. LDU >= fla_max(1,M) if */
/* > JOBU = 'U';
LDU >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* > V is COMPLEX array, dimension (LDV,P) */
/* > If JOBV = 'V', V contains the unitary matrix V. */
/* > If JOBV = 'N', V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV is INTEGER */
/* > The leading dimension of the array V. LDV >= fla_max(1,P) if */
/* > JOBV = 'V';
LDV >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* > Q is COMPLEX array, dimension (LDQ,N) */
/* > If JOBQ = 'Q', Q contains the unitary matrix Q. */
/* > If JOBQ = 'N', Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* > LDQ is INTEGER */
/* > The leading dimension of the array Q. LDQ >= fla_max(1,N) if */
/* > JOBQ = 'Q';
LDQ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (N) */
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
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date August 2015 */
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The subroutine uses LAPACK subroutine CGEQP3 for the QR factorization */
/* > with column pivoting to detect the effective numerical rank of the */
/* > a matrix. It may be replaced by a better rank determination strategy. */
/* > */
/* > CGGSVP3 replaces the deprecated subroutine CGGSVP. */
/* > */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void cggsvp3_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, complex *a,
              integer *lda, complex *b, integer *ldb, real *tola, real *tolb, integer *k,
              integer *l, complex *u, integer *ldu, complex *v, integer *ldv, complex *q,
              integer *ldq, integer *iwork, real *rwork, complex *tau, complex *work,
              integer *lwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cggsvp3 inputs: jobu %c, jobv %c, jobq %c, m %lld, p %lld, n %lld, lda %lld, ldb "
             "%lld, ldu %lld, ldv %lld, ldq %lld, lwork %lld",
             *jobu, *jobv, *jobq, *m, *p, *n, *lda, *ldb, *ldu, *ldv, *ldq, *lwork);
#else
    snprintf(buffer, 256,
             "cggsvp3 inputs: jobu %c, jobv %c, jobq %c, m %d, p %d, n %d, lda %d, ldb %d, ldu %d, "
             "ldv %d, ldq %d, lwork %d",
             *jobu, *jobv, *jobq, *m, *p, *n, *lda, *ldb, *ldu, *ldv, *ldq, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, u_offset, v_dim1,
        v_offset, i__1, i__2, i__3;
    complex q__1;
    /* Builtin functions */
    double c_abs(complex *);
    /* Local variables */
    integer i__, j;
    extern logical lsame_(char *, char *, integer, integer);
    logical wantq, wantu, wantv;
    extern /* Subroutine */
        void
        cgeqp3_(integer *, integer *, complex *, integer *, integer *, complex *, complex *,
                integer *, real *, integer *),
        cgeqr2_(integer *, integer *, complex *, integer *, complex *, complex *, integer *),
        cgerq2_(integer *, integer *, complex *, integer *, complex *, complex *, integer *),
        cung2r_(integer *, integer *, integer *, complex *, integer *, complex *, complex *,
                integer *),
        cunm2r_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                complex *, integer *, complex *, integer *),
        cunmr2_(char *, char *, integer *, integer *, integer *, complex *, integer *, complex *,
                complex *, integer *, complex *, integer *),
        clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *),
        claset_(char *, integer *, integer *, complex *, complex *, complex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        clapmt_(logical *, integer *, integer *, complex *, integer *, integer *);
    logical forwrd;
    integer lwkopt;
    logical lquery;
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* August 2015 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --iwork;
    --rwork;
    --tau;
    --work;
    /* Function Body */
    wantu = lsame_(jobu, "U", 1, 1);
    wantv = lsame_(jobv, "V", 1, 1);
    wantq = lsame_(jobq, "Q", 1, 1);
    forwrd = TRUE_;
    lquery = *lwork == -1;
    lwkopt = 1;
    /* Test the input arguments */
    *info = 0;
    if(!(wantu || lsame_(jobu, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(wantv || lsame_(jobv, "N", 1, 1)))
    {
        *info = -2;
    }
    else if(!(wantq || lsame_(jobq, "N", 1, 1)))
    {
        *info = -3;
    }
    else if(*m < 0)
    {
        *info = -4;
    }
    else if(*p < 0)
    {
        *info = -5;
    }
    else if(*n < 0)
    {
        *info = -6;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -8;
    }
    else if(*ldb < fla_max(1, *p))
    {
        *info = -10;
    }
    else if(*ldu < 1 || wantu && *ldu < *m)
    {
        *info = -16;
    }
    else if(*ldv < 1 || wantv && *ldv < *p)
    {
        *info = -18;
    }
    else if(*ldq < 1 || wantq && *ldq < *n)
    {
        *info = -20;
    }
    else if(*lwork < 1 && !lquery)
    {
        *info = -24;
    }
    /* Compute workspace */
    if(*info == 0)
    {
        cgeqp3_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], &c_n1, &rwork[1], info);
        lwkopt = (integer)work[1].r;
        if(wantv)
        {
            lwkopt = fla_max(lwkopt, *p);
        }
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = fla_min(*n, *p); // , expr subst
        lwkopt = fla_max(i__1, i__2);
        lwkopt = fla_max(lwkopt, *m);
        if(wantq)
        {
            lwkopt = fla_max(lwkopt, *n);
        }
        cgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], &c_n1, &rwork[1], info);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer)work[1].r; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        lwkopt = fla_max(1, lwkopt);
        q__1.r = (real)lwkopt;
        q__1.i = 0.f; // , expr subst
        work[1].r = q__1.r;
        work[1].i = q__1.i; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGGSVP3", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* QR with column pivoting of B: B*P = V*( S11 S12 ) */
    /* ( 0 0 ) */
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        iwork[i__] = 0;
        /* L10: */
    }
    cgeqp3_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], lwork, &rwork[1], info);
    /* Update A := A*P */
    clapmt_(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);
    /* Determine the effective rank of matrix B. */
    *l = 0;
    i__1 = fla_min(*p, *n);
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(c_abs(&b[i__ + i__ * b_dim1]) > *tolb)
        {
            ++(*l);
        }
        /* L20: */
    }
    if(wantv)
    {
        /* Copy the details of V, and form V. */
        claset_("Full", p, p, &c_b1, &c_b1, &v[v_offset], ldv);
        if(*p > 1)
        {
            i__1 = *p - 1;
            clacpy_("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], ldv);
        }
        i__1 = fla_min(*p, *n);
        cung2r_(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
    }
    /* Clean up B */
    i__1 = *l - 1;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *l;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * b_dim1;
            b[i__3].r = 0.f;
            b[i__3].i = 0.f; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    if(*p > *l)
    {
        i__1 = *p - *l;
        claset_("Full", &i__1, n, &c_b1, &c_b1, &b[*l + 1 + b_dim1], ldb);
    }
    if(wantq)
    {
        /* Set Q = I and Update Q := Q*P */
        claset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq);
        clapmt_(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
    }
    if(*p >= *l && *n != *l)
    {
        /* RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z */
        cgerq2_(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);
        /* Update A := A*Z**H */
        cunmr2_("Right", "Conjugate transpose", m, n, l, &b[b_offset], ldb, &tau[1], &a[a_offset],
                lda, &work[1], info);
        if(wantq)
        {
            /* Update Q := Q*Z**H */
            cunmr2_("Right", "Conjugate transpose", n, n, l, &b[b_offset], ldb, &tau[1],
                    &q[q_offset], ldq, &work[1], info);
        }
        /* Clean up B */
        i__1 = *n - *l;
        claset_("Full", l, &i__1, &c_b1, &c_b1, &b[b_offset], ldb);
        i__1 = *n;
        for(j = *n - *l + 1; j <= i__1; ++j)
        {
            i__2 = *l;
            for(i__ = j - *n + *l + 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * b_dim1;
                b[i__3].r = 0.f;
                b[i__3].i = 0.f; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
    }
    /* Let N-L L */
    /* A = ( A11 A12 ) M, */
    /* then the following does the complete QR decomposition of A11: */
    /* A11 = U*( 0 T12 )*P1**H */
    /* ( 0 0 ) */
    i__1 = *n - *l;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        iwork[i__] = 0;
        /* L70: */
    }
    i__1 = *n - *l;
    cgeqp3_(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], lwork, &rwork[1], info);
    /* Determine the effective rank of A11 */
    *k = 0;
    /* Computing MIN */
    i__2 = *m;
    i__3 = *n - *l; // , expr subst
    i__1 = fla_min(i__2, i__3);
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(c_abs(&a[i__ + i__ * a_dim1]) > *tola)
        {
            ++(*k);
        }
        /* L80: */
    }
    /* Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N ) */
    /* Computing MIN */
    i__2 = *m;
    i__3 = *n - *l; // , expr subst
    i__1 = fla_min(i__2, i__3);
    cunm2r_("Left", "Conjugate transpose", m, l, &i__1, &a[a_offset], lda, &tau[1],
            &a[(*n - *l + 1) * a_dim1 + 1], lda, &work[1], info);
    if(wantu)
    {
        /* Copy the details of U, and form U */
        claset_("Full", m, m, &c_b1, &c_b1, &u[u_offset], ldu);
        if(*m > 1)
        {
            i__1 = *m - 1;
            i__2 = *n - *l;
            clacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2], ldu);
        }
        /* Computing MIN */
        i__2 = *m;
        i__3 = *n - *l; // , expr subst
        i__1 = fla_min(i__2, i__3);
        cung2r_(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
    }
    if(wantq)
    {
        /* Update Q( 1:N, 1:N-L ) = Q( 1:N, 1:N-L )*P1 */
        i__1 = *n - *l;
        clapmt_(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
    }
    /* Clean up A: set the strictly lower triangular part of */
    /* A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */
    i__1 = *k - 1;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *k;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * a_dim1;
            a[i__3].r = 0.f;
            a[i__3].i = 0.f; // , expr subst
            /* L90: */
        }
        /* L100: */
    }
    if(*m > *k)
    {
        i__1 = *m - *k;
        i__2 = *n - *l;
        claset_("Full", &i__1, &i__2, &c_b1, &c_b1, &a[*k + 1 + a_dim1], lda);
    }
    if(*n - *l > *k)
    {
        /* RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */
        i__1 = *n - *l;
        cgerq2_(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);
        if(wantq)
        {
            /* Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H */
            i__1 = *n - *l;
            cunmr2_("Right", "Conjugate transpose", n, &i__1, k, &a[a_offset], lda, &tau[1],
                    &q[q_offset], ldq, &work[1], info);
        }
        /* Clean up A */
        i__1 = *n - *l - *k;
        claset_("Full", k, &i__1, &c_b1, &c_b1, &a[a_offset], lda);
        i__1 = *n - *l;
        for(j = *n - *l - *k + 1; j <= i__1; ++j)
        {
            i__2 = *k;
            for(i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * a_dim1;
                a[i__3].r = 0.f;
                a[i__3].i = 0.f; // , expr subst
                /* L110: */
            }
            /* L120: */
        }
    }
    if(*m > *k)
    {
        /* QR factorization of A( K+1:M,N-L+1:N ) */
        i__1 = *m - *k;
        cgeqr2_(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &work[1], info);
        if(wantu)
        {
            /* Update U(:,K+1:M) := U(:,K+1:M)*U1 */
            i__1 = *m - *k;
            /* Computing MIN */
            i__3 = *m - *k;
            i__2 = fla_min(i__3, *l);
            cunm2r_("Right", "No transpose", m, &i__1, &i__2, &a[*k + 1 + (*n - *l + 1) * a_dim1],
                    lda, &tau[1], &u[(*k + 1) * u_dim1 + 1], ldu, &work[1], info);
        }
        /* Clean up */
        i__1 = *n;
        for(j = *n - *l + 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * a_dim1;
                a[i__3].r = 0.f;
                a[i__3].i = 0.f; // , expr subst
                /* L130: */
            }
            /* L140: */
        }
    }
    q__1.r = (real)lwkopt;
    q__1.i = 0.f; // , expr subst
    work[1].r = q__1.r;
    work[1].i = q__1.i; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGGSVP3 */
}
/* cggsvp3_ */
