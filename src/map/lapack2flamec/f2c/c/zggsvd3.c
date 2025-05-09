/* ../netlib/v3.9.0/zggsvd3.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
static integer c__1 = 1;
/* > \brief <b> ZGGSVD3 computes the singular value decomposition (SVD) for OTHER matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGGSVD3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggsvd3
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggsvd3
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggsvd3
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, */
/* LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, */
/* LWORK, RWORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBQ, JOBU, JOBV */
/* INTEGER INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IWORK( * ) */
/* DOUBLE PRECISION ALPHA( * ), BETA( * ), RWORK( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ U( LDU, * ), V( LDV, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGSVD3 computes the generalized singular value decomposition (GSVD) */
/* > of an M-by-N complex matrix A and P-by-N complex matrix B: */
/* > */
/* > U**H*A*Q = D1*( 0 R ), V**H*B*Q = D2*( 0 R ) */
/* > */
/* > where U, V and Q are unitary matrices. */
/* > Let K+L = the effective numerical rank of the */
/* > matrix (A**H,B**H)**H, then R is a (K+L)-by-(K+L) nonsingular upper */
/* > triangular matrix, D1 and D2 are M-by-(K+L) and P-by-(K+L) "diagonal" */
/* > matrices and of the following structures, respectively: */
/* > */
/* > If M-K-L >= 0, */
/* > */
/* > K L */
/* > D1 = K ( I 0 ) */
/* > L ( 0 C ) */
/* > M-K-L ( 0 0 ) */
/* > */
/* > K L */
/* > D2 = L ( 0 S ) */
/* > P-L ( 0 0 ) */
/* > */
/* > N-K-L K L */
/* > ( 0 R ) = K ( 0 R11 R12 ) */
/* > L ( 0 0 R22 ) */
/* > where */
/* > */
/* > C = diag( ALPHA(K+1), ... , ALPHA(K+L) ), */
/* > S = diag( BETA(K+1), ... , BETA(K+L) ), */
/* > C**2 + S**2 = I. */
/* > */
/* > R is stored in A(1:K+L,N-K-L+1:N) on exit. */
/* > */
/* > If M-K-L < 0, */
/* > */
/* > K M-K K+L-M */
/* > D1 = K ( I 0 0 ) */
/* > M-K ( 0 C 0 ) */
/* > */
/* > K M-K K+L-M */
/* > D2 = M-K ( 0 S 0 ) */
/* > K+L-M ( 0 0 I ) */
/* > P-L ( 0 0 0 ) */
/* > */
/* > N-K-L K M-K K+L-M */
/* > ( 0 R ) = K ( 0 R11 R12 R13 ) */
/* > M-K ( 0 0 R22 R23 ) */
/* > K+L-M ( 0 0 0 R33 ) */
/* > */
/* > where */
/* > */
/* > C = diag( ALPHA(K+1), ... , ALPHA(M) ), */
/* > S = diag( BETA(K+1), ... , BETA(M) ), */
/* > C**2 + S**2 = I. */
/* > */
/* > (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored */
/* > ( 0 R22 R23 ) */
/* > in B(M-K+1:L,N+M-K-L+1:N) on exit. */
/* > */
/* > The routine computes C, S, R, and optionally the unitary */
/* > transformation matrices U, V and Q. */
/* > */
/* > In particular, if B is an N-by-N nonsingular matrix, then the GSVD of */
/* > A and B implicitly gives the SVD of A*inv(B): */
/* > A*inv(B) = U*(D1*inv(D2))*V**H. */
/* > If ( A**H,B**H)**H has orthonormal columns, then the GSVD of A and B is also */
/* > equal to the CS decomposition of A and B. Furthermore, the GSVD can */
/* > be used to derive the solution of the eigenvalue problem: */
/* > A**H*A x = lambda* B**H*B x. */
/* > In some literature, the GSVD of A and B is presented in the form */
/* > U**H*A*X = ( 0 D1 ), V**H*B*X = ( 0 D2 ) */
/* > where U and V are orthogonal and X is nonsingular, and D1 and D2 are */
/* > ``diagonal''. The former GSVD form can be converted to the latter */
/* > form by taking the nonsingular matrix X as */
/* > */
/* > X = Q*( I 0 ) */
/* > ( 0 inv(R) ) */
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
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* > P is INTEGER */
/* > The number of rows of the matrix B. P >= 0. */
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
/* > described in Purpose. */
/* > K + L = effective numerical rank of (A**H,B**H)**H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, A contains the triangular matrix R, or part of R. */
/* > See Purpose for details. */
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
/* > B is COMPLEX*16 array, dimension (LDB,N) */
/* > On entry, the P-by-N matrix B. */
/* > On exit, B contains part of the triangular matrix R if */
/* > M-K-L < 0. See Purpose for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,P). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* > ALPHA is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* > BETA is DOUBLE PRECISION array, dimension (N) */
/* > */
/* > On exit, ALPHA and BETA contain the generalized singular */
/* > value pairs of A and B;
 */
/* > ALPHA(1:K) = 1, */
/* > BETA(1:K) = 0, */
/* > and if M-K-L >= 0, */
/* > ALPHA(K+1:K+L) = C, */
/* > BETA(K+1:K+L) = S, */
/* > or if M-K-L < 0, */
/* > ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0 */
/* > BETA(K+1:M) =S, BETA(M+1:K+L) =1 */
/* > and */
/* > ALPHA(K+L+1:N) = 0 */
/* > BETA(K+L+1:N) = 0 */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* > U is COMPLEX*16 array, dimension (LDU,M) */
/* > If JOBU = 'U', U contains the M-by-M unitary matrix U. */
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
/* > V is COMPLEX*16 array, dimension (LDV,P) */
/* > If JOBV = 'V', V contains the P-by-P unitary matrix V. */
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
/* > Q is COMPLEX*16 array, dimension (LDQ,N) */
/* > If JOBQ = 'Q', Q contains the N-by-N unitary matrix Q. */
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
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
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
/* > RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
/* > On exit, IWORK stores the sorting information. More */
/* > precisely, the following loop will sort ALPHA */
/* > for I = K+1, fla_min(M,K+L) */
/* > swap ALPHA(I) and ALPHA(IWORK(I)) */
/* > endfor */
/* > such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit. */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > > 0: if INFO = 1, the Jacobi-type procedure failed to */
/* > converge. For further details, see subroutine ZTGSJA. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > TOLA DOUBLE PRECISION */
/* > TOLB DOUBLE PRECISION */
/* > TOLA and TOLB are the thresholds to determine the effective */
/* > rank of (A**H,B**H)**H. Generally, they are set to */
/* > TOLA = MAX(M,N)*norm(A)*MACHEPS, */
/* > TOLB = MAX(P,N)*norm(B)*MACHEPS. */
/* > The size of TOLA and TOLB may affect the size of backward */
/* > errors of the decomposition. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date August 2015 */
/* > \ingroup complex16GEsing */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Ming Gu and Huan Ren, Computer Science Division, University of */
/* > California at Berkeley, USA */
/* > */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > ZGGSVD3 replaces the deprecated subroutine ZGGSVD. */
/* > */
/* ===================================================================== */
/* Subroutine */
void zggsvd3_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k,
              integer *l, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
              doublereal *alpha, doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v,
              integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, integer *lwork,
              doublereal *rwork, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zggsvd3 inputs: jobu %c, jobv %c, jobq %c, m %" FLA_IS ", n %" FLA_IS
                      ", p %" FLA_IS ", k %" FLA_IS ", l %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS
                      ", ldu %" FLA_IS ", ldv %" FLA_IS ", ldq %" FLA_IS "",
                      *jobu, *jobv, *jobq, *m, *n, *p, *k, *l, *lda, *ldb, *ldu, *ldv, *ldq);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, u_offset, v_dim1,
        v_offset, i__1, i__2;
    doublecomplex z__1;
    /* Local variables */
    integer i__, j;
    doublereal ulp;
    integer ibnd;
    doublereal tola;
    integer isub;
    doublereal tolb, unfl, temp, smax;
    extern logical lsame_(char *, char *, integer, integer);
    doublereal anorm, bnorm;
    extern /* Subroutine */
        void
        dcopy_(integer *, doublereal *, integer *, doublereal *, integer *);
    logical wantq, wantu, wantv;
    extern doublereal dlamch_(char *);
    integer ncycle;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, integer *,
                              doublereal *);
    extern /* Subroutine */
        void
        ztgsja_(char *, char *, char *, integer *, integer *, integer *, integer *, integer *,
                doublecomplex *, integer *, doublecomplex *, integer *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublecomplex *, integer *, doublecomplex *, integer *,
                doublecomplex *, integer *, doublecomplex *, integer *, integer *);
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
        void
        zggsvp3_(char *, char *, char *, integer *, integer *, integer *, doublecomplex *,
                 integer *, doublecomplex *, integer *, doublereal *, doublereal *, integer *,
                 integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublecomplex *,
                 integer *, integer *, doublereal *, doublecomplex *, doublecomplex *, integer *,
                 integer *);
    /* -- LAPACK driver routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* August 2015 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Decode and test the input parameters */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --alpha;
    --beta;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --work;
    --rwork;
    --iwork;
    /* Function Body */
    wantu = lsame_(jobu, "U", 1, 1);
    wantv = lsame_(jobv, "V", 1, 1);
    wantq = lsame_(jobq, "Q", 1, 1);
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
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*p < 0)
    {
        *info = -6;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -10;
    }
    else if(*ldb < fla_max(1, *p))
    {
        *info = -12;
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
        zggsvp3_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], ldb, &tola, &tolb, k,
                 l, &u[u_offset], ldu, &v[v_offset], ldv, &q[q_offset], ldq, &iwork[1], &rwork[1],
                 &work[1], &work[1], &c_n1, info);
        lwkopt = *n + (integer)work[1].r;
        /* Computing MAX */
        i__1 = *n << 1;
        lwkopt = fla_max(i__1, lwkopt);
        lwkopt = fla_max(1, lwkopt);
        z__1.r = (doublereal)lwkopt;
        z__1.i = 0.; // , expr subst
        work[1].r = z__1.r;
        work[1].i = z__1.i; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGGSVD3", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Compute the Frobenius norm of matrices A and B */
    anorm = zlange_("1", m, n, &a[a_offset], lda, &rwork[1]);
    bnorm = zlange_("1", p, n, &b[b_offset], ldb, &rwork[1]);
    /* Get machine precision and set up threshold for determining */
    /* the effective numerical rank of the matrices A and B. */
    ulp = dlamch_("Precision");
    unfl = dlamch_("Safe Minimum");
    tola = fla_max(*m, *n) * fla_max(anorm, unfl) * ulp;
    tolb = fla_max(*p, *n) * fla_max(bnorm, unfl) * ulp;
    i__1 = *lwork - *n;
    zggsvp3_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], ldb, &tola, &tolb, k, l,
             &u[u_offset], ldu, &v[v_offset], ldv, &q[q_offset], ldq, &iwork[1], &rwork[1],
             &work[1], &work[*n + 1], &i__1, info);
    /* Compute the GSVD of two upper "triangular" matrices */
    ztgsja_(jobu, jobv, jobq, m, p, n, k, l, &a[a_offset], lda, &b[b_offset], ldb, &tola, &tolb,
            &alpha[1], &beta[1], &u[u_offset], ldu, &v[v_offset], ldv, &q[q_offset], ldq, &work[1],
            &ncycle, info);
    /* Sort the singular values and store the pivot indices in IWORK */
    /* Copy ALPHA to RWORK, then sort ALPHA in RWORK */
    dcopy_(n, &alpha[1], &c__1, &rwork[1], &c__1);
    /* Computing MIN */
    i__1 = *l;
    i__2 = *m - *k; // , expr subst
    ibnd = fla_min(i__1, i__2);
    i__1 = ibnd;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Scan for largest ALPHA(K+I) */
        isub = i__;
        smax = rwork[*k + i__];
        i__2 = ibnd;
        for(j = i__ + 1; j <= i__2; ++j)
        {
            temp = rwork[*k + j];
            if(temp > smax)
            {
                isub = j;
                smax = temp;
            }
            /* L10: */
        }
        if(isub != i__)
        {
            rwork[*k + isub] = rwork[*k + i__];
            rwork[*k + i__] = smax;
            iwork[*k + i__] = *k + isub;
        }
        else
        {
            iwork[*k + i__] = *k + i__;
        }
        /* L20: */
    }
    z__1.r = (doublereal)lwkopt;
    z__1.i = 0.; // , expr subst
    work[1].r = z__1.r;
    work[1].i = z__1.i; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGGSVD3 */
}
/* zggsvd3_ */
