/* ./sggsvp3.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c_n1 = -1;
static real c_b14 = 0.f;
static real c_b24 = 1.f;
/* > \brief \b SGGSVP3 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGGSVP3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggsvp3
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggsvp3
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggsvp3
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/* TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/* IWORK, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBQ, JOBU, JOBV */
/* INTEGER INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK */
/* REAL TOLA, TOLB */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* REAL A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGSVP3 computes orthogonal matrices U, V and Q such that */
/* > */
/* > N-K-L K L */
/* > U**T*A*Q = K ( 0 A12 A13 ) if M-K-L >= 0;
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
/* > V**T*B*Q = L ( 0 0 B13 ) */
/* > P-L ( 0 0 0 ) */
/* > */
/* > where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular */
/* > upper triangular;
A23 is L-by-L upper triangular if M-K-L >= 0, */
/* > otherwise A23 is (M-K)-by-L upper trapezoidal. K+L = the effective */
/* > numerical rank of the (M+P)-by-N matrix (A**T,B**T)**T. */
/* > */
/* > This decomposition is the preprocessing step for computing the */
/* > Generalized Singular Value Decomposition (GSVD), see subroutine */
/* > SGGSVD3. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU */
/* > \verbatim */
/* > JOBU is CHARACTER*1 */
/* > = 'U': Orthogonal matrix U is computed;
 */
/* > = 'N': U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* > JOBV is CHARACTER*1 */
/* > = 'V': Orthogonal matrix V is computed;
 */
/* > = 'N': V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* > JOBQ is CHARACTER*1 */
/* > = 'Q': Orthogonal matrix Q is computed;
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
/* > A is REAL array, dimension (LDA,N) */
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
/* > B is REAL array, dimension (LDB,N) */
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
/* > K + L = effective numerical rank of (A**T,B**T)**T. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* > U is REAL array, dimension (LDU,M) */
/* > If JOBU = 'U', U contains the orthogonal matrix U. */
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
/* > V is REAL array, dimension (LDV,P) */
/* > If JOBV = 'V', V contains the orthogonal matrix V. */
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
/* > Q is REAL array, dimension (LDQ,N) */
/* > If JOBQ = 'Q', Q contains the orthogonal matrix Q. */
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
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is REAL array, dimension (N) */
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
/* > \ingroup ggsvp3 */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The subroutine uses LAPACK subroutine SGEQP3 for the QR factorization */
/* > with column pivoting to detect the effective numerical rank of the */
/* > a matrix. It may be replaced by a better rank determination strategy. */
/* > */
/* > SGGSVP3 replaces the deprecated subroutine SGGSVP. */
/* > */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sggsvp3_(char *jobu, char *jobv, char *jobq, aocl_int_t *m, aocl_int_t *p, aocl_int_t *n,
              real *a, aocl_int_t *lda, real *b, aocl_int_t *ldb, real *tola, real *tolb,
              aocl_int_t *k, aocl_int_t *l, real *u, aocl_int_t *ldu, real *v, aocl_int_t *ldv,
              real *q, aocl_int_t *ldq, aocl_int_t *iwork, real *tau, real *work, aocl_int_t *lwork,
              aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv,
                        q, ldq, iwork, tau, work, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t p_64 = *p;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t k_64 = *k;
    aocl_int64_t l_64 = *l;
    aocl_int64_t ldu_64 = *ldu;
    aocl_int64_t ldv_64 = *ldv;
    aocl_int64_t ldq_64 = *ldq;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sggsvp3(jobu, jobv, jobq, &m_64, &p_64, &n_64, a, &lda_64, b, &ldb_64, tola, tolb,
                        &k_64, &l_64, u, &ldu_64, v, &ldv_64, q, &ldq_64, iwork, tau, work,
                        &lwork_64, &info_64);

    *k = (aocl_int_t)k_64;
    *l = (aocl_int_t)l_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sggsvp3(char *jobu, char *jobv, char *jobq, aocl_int64_t *m, aocl_int64_t *p,
                         aocl_int64_t *n, real *a, aocl_int64_t *lda, real *b, aocl_int64_t *ldb,
                         real *tola, real *tolb, aocl_int64_t *k, aocl_int64_t *l, real *u,
                         aocl_int64_t *ldu, real *v, aocl_int64_t *ldv, real *q, aocl_int64_t *ldq,
                         aocl_int_t *iwork, real *tau, real *work, aocl_int64_t *lwork,
                         aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sggsvp3 inputs: jobu %c ,jobv %c ,jobq %c ,m %" FLA_IS ",p %" FLA_IS
                      ",n %" FLA_IS ",lda %" FLA_IS ",ldb %" FLA_IS ",l %" FLA_IS ",ldu %" FLA_IS
                      ",ldv %" FLA_IS ",ldq %" FLA_IS ",lwork %" FLA_IS "",
                      *jobu, *jobv, *jobq, *m, *p, *n, *lda, *ldb, *l, *ldu, *ldv, *ldq, *lwork);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, u_offset, v_dim1,
        v_offset, i__1, i__2, i__3;
    real r__1;
    /* Local variables */
    aocl_int64_t i__, j;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical wantq, wantu, wantv;
    logical forwrd;
    aocl_int64_t lwkopt;
    logical lquery;
    /* -- LAPACK computational routine -- */
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
        aocl_lapack_sgeqp3(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], &c_n1, info);
        lwkopt = (integer)work[1];
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
        aocl_lapack_sgeqp3(m, n, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], &c_n1, info);
        /* Computing MAX */
        i__1 = lwkopt;
        i__2 = (integer)work[1]; // , expr subst
        lwkopt = fla_max(i__1, i__2);
        lwkopt = fla_max(1, lwkopt);
        work[1] = aocl_lapack_sroundup_lwork(&lwkopt);
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SGGSVP3", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
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
    aocl_lapack_sgeqp3(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], lwork, info);
    /* Update A := A*P */
    aocl_lapack_slapmt(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);
    /* Determine the effective rank of matrix B. */
    *l = 0;
    i__1 = fla_min(*p, *n);
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if((r__1 = b[i__ + i__ * b_dim1], f2c_abs(r__1)) > *tolb)
        {
            ++(*l);
        }
        /* L20: */
    }
    if(wantv)
    {
        /* Copy the details of V, and form V. */
        aocl_lapack_slaset("Full", p, p, &c_b14, &c_b14, &v[v_offset], ldv);
        if(*p > 1)
        {
            i__1 = *p - 1;
            aocl_lapack_slacpy("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], ldv);
        }
        i__1 = fla_min(*p, *n);
        aocl_lapack_sorg2r(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
    }
    /* Clean up B */
    i__1 = *l - 1;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *l;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            b[i__ + j * b_dim1] = 0.f;
            /* L30: */
        }
        /* L40: */
    }
    if(*p > *l)
    {
        i__1 = *p - *l;
        aocl_lapack_slaset("Full", &i__1, n, &c_b14, &c_b14, &b[*l + 1 + b_dim1], ldb);
    }
    if(wantq)
    {
        /* Set Q = I and Update Q := Q*P */
        aocl_lapack_slaset("Full", n, n, &c_b14, &c_b24, &q[q_offset], ldq);
        aocl_lapack_slapmt(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
    }
    if(*p >= *l && *n != *l)
    {
        /* RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z */
        aocl_lapack_sgerq2(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);
        /* Update A := A*Z**T */
        aocl_lapack_sormr2("Right", "Transpose", m, n, l, &b[b_offset], ldb, &tau[1], &a[a_offset],
                           lda, &work[1], info);
        if(wantq)
        {
            /* Update Q := Q*Z**T */
            aocl_lapack_sormr2("Right", "Transpose", n, n, l, &b[b_offset], ldb, &tau[1],
                               &q[q_offset], ldq, &work[1], info);
        }
        /* Clean up B */
        i__1 = *n - *l;
        aocl_lapack_slaset("Full", l, &i__1, &c_b14, &c_b14, &b[b_offset], ldb);
        i__1 = *n;
        for(j = *n - *l + 1; j <= i__1; ++j)
        {
            i__2 = *l;
            for(i__ = j - *n + *l + 1; i__ <= i__2; ++i__)
            {
                b[i__ + j * b_dim1] = 0.f;
                /* L50: */
            }
            /* L60: */
        }
    }
    /* Let N-L L */
    /* A = ( A11 A12 ) M, */
    /* then the following does the complete QR decomposition of A11: */
    /* A11 = U*( 0 T12 )*P1**T */
    /* ( 0 0 ) */
    i__1 = *n - *l;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        iwork[i__] = 0;
        /* L70: */
    }
    i__1 = *n - *l;
    aocl_lapack_sgeqp3(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], lwork, info);
    /* Determine the effective rank of A11 */
    *k = 0;
    /* Computing MIN */
    i__2 = *m;
    i__3 = *n - *l; // , expr subst
    i__1 = fla_min(i__2, i__3);
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if((r__1 = a[i__ + i__ * a_dim1], f2c_abs(r__1)) > *tola)
        {
            ++(*k);
        }
        /* L80: */
    }
    /* Update A12 := U**T*A12, where A12 = A( 1:M, N-L+1:N ) */
    /* Computing MIN */
    i__2 = *m;
    i__3 = *n - *l; // , expr subst
    i__1 = fla_min(i__2, i__3);
    aocl_lapack_sorm2r("Left", "Transpose", m, l, &i__1, &a[a_offset], lda, &tau[1],
                       &a[(*n - *l + 1) * a_dim1 + 1], lda, &work[1], info);
    if(wantu)
    {
        /* Copy the details of U, and form U */
        aocl_lapack_slaset("Full", m, m, &c_b14, &c_b14, &u[u_offset], ldu);
        if(*m > 1)
        {
            i__1 = *m - 1;
            i__2 = *n - *l;
            aocl_lapack_slacpy("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2], ldu);
        }
        /* Computing MIN */
        i__2 = *m;
        i__3 = *n - *l; // , expr subst
        i__1 = fla_min(i__2, i__3);
        aocl_lapack_sorg2r(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
    }
    if(wantq)
    {
        /* Update Q( 1:N, 1:N-L ) = Q( 1:N, 1:N-L )*P1 */
        i__1 = *n - *l;
        aocl_lapack_slapmt(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
    }
    /* Clean up A: set the strictly lower triangular part of */
    /* A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */
    i__1 = *k - 1;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *k;
        for(i__ = j + 1; i__ <= i__2; ++i__)
        {
            a[i__ + j * a_dim1] = 0.f;
            /* L90: */
        }
        /* L100: */
    }
    if(*m > *k)
    {
        i__1 = *m - *k;
        i__2 = *n - *l;
        aocl_lapack_slaset("Full", &i__1, &i__2, &c_b14, &c_b14, &a[*k + 1 + a_dim1], lda);
    }
    if(*n - *l > *k)
    {
        /* RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */
        i__1 = *n - *l;
        aocl_lapack_sgerq2(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);
        if(wantq)
        {
            /* Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**T */
            i__1 = *n - *l;
            aocl_lapack_sormr2("Right", "Transpose", n, &i__1, k, &a[a_offset], lda, &tau[1],
                               &q[q_offset], ldq, &work[1], info);
        }
        /* Clean up A */
        i__1 = *n - *l - *k;
        aocl_lapack_slaset("Full", k, &i__1, &c_b14, &c_b14, &a[a_offset], lda);
        i__1 = *n - *l;
        for(j = *n - *l - *k + 1; j <= i__1; ++j)
        {
            i__2 = *k;
            for(i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__)
            {
                a[i__ + j * a_dim1] = 0.f;
                /* L110: */
            }
            /* L120: */
        }
    }
    if(*m > *k)
    {
        /* QR factorization of A( K+1:M,N-L+1:N ) */
        i__1 = *m - *k;
        aocl_lapack_sgeqr2(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &work[1],
                           info);
        if(wantu)
        {
            /* Update U(:,K+1:M) := U(:,K+1:M)*U1 */
            i__1 = *m - *k;
            /* Computing MIN */
            i__3 = *m - *k;
            i__2 = fla_min(i__3, *l);
            aocl_lapack_sorm2r("Right", "No transpose", m, &i__1, &i__2,
                               &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1],
                               &u[(*k + 1) * u_dim1 + 1], ldu, &work[1], info);
        }
        /* Clean up */
        i__1 = *n;
        for(j = *n - *l + 1; j <= i__1; ++j)
        {
            i__2 = *m;
            for(i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__)
            {
                a[i__ + j * a_dim1] = 0.f;
                /* L130: */
            }
            /* L140: */
        }
    }
    work[1] = aocl_lapack_sroundup_lwork(&lwkopt);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGGSVP3 */
}
/* sggsvp3_ */
