/* ztgsja.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 = {0., 0.};
static doublecomplex c_b2 = {1., 0.};
static integer c__1 = 1;
static doublereal c_b39 = -1.;
static doublereal c_b42 = 1.;
/* > \brief \b ZTGSJA */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZTGSJA + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsja.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsja.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsja.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, */
/* LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, */
/* Q, LDQ, WORK, NCYCLE, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER JOBQ, JOBU, JOBV */
/* INTEGER INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, */
/* $ NCYCLE, P */
/* DOUBLE PRECISION TOLA, TOLB */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION ALPHA( * ), BETA( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/* $ U( LDU, * ), V( LDV, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTGSJA computes the generalized singular value decomposition (GSVD) */
/* > of two complex upper triangular (or trapezoidal) matrices A and B. */
/* > */
/* > On entry, it is assumed that matrices A and B have the following */
/* > forms, which may be obtained by the preprocessing subroutine ZGGSVP */
/* > from a general M-by-N matrix A and P-by-N matrix B: */
/* > */
/* > N-K-L K L */
/* > A = K ( 0 A12 A13 ) if M-K-L >= 0;
 */
/* > L ( 0 0 A23 ) */
/* > M-K-L ( 0 0 0 ) */
/* > */
/* > N-K-L K L */
/* > A = K ( 0 A12 A13 ) if M-K-L < 0;
 */
/* > M-K ( 0 0 A23 ) */
/* > */
/* > N-K-L K L */
/* > B = L ( 0 0 B13 ) */
/* > P-L ( 0 0 0 ) */
/* > */
/* > where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular */
/* > upper triangular;
A23 is L-by-L upper triangular if M-K-L >= 0, */
/* > otherwise A23 is (M-K)-by-L upper trapezoidal. */
/* > */
/* > On exit, */
/* > */
/* > U**H *A*Q = D1*( 0 R ), V**H *B*Q = D2*( 0 R ), */
/* > */
/* > where U, V and Q are unitary matrices. */
/* > R is a nonsingular upper triangular matrix, and D1 */
/* > and D2 are ``diagonal'' matrices, which are of the following */
/* > structures: */
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
/* > ( 0 R ) = K ( 0 R11 R12 ) K */
/* > L ( 0 0 R22 ) L */
/* > */
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
/* > C = diag( ALPHA(K+1), ... , ALPHA(M) ), */
/* > S = diag( BETA(K+1), ... , BETA(M) ), */
/* > C**2 + S**2 = I. */
/* > */
/* > R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored */
/* > ( 0 R22 R23 ) */
/* > in B(M-K+1:L,N+M-K-L+1:N) on exit. */
/* > */
/* > The computation of the unitary transformation matrices U, V or Q */
/* > is optional. These matrices may either be formed explicitly, or they */
/* > may be postmultiplied into input matrices U1, V1, or Q1. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] JOBU */
/* > \verbatim */
/* > JOBU is CHARACTER*1 */
/* > = 'U': U must contain a unitary matrix U1 on entry, and */
/* > the product U1*U is returned;
 */
/* > = 'I': U is initialized to the unit matrix, and the */
/* > unitary matrix U is returned;
 */
/* > = 'N': U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* > JOBV is CHARACTER*1 */
/* > = 'V': V must contain a unitary matrix V1 on entry, and */
/* > the product V1*V is returned;
 */
/* > = 'I': V is initialized to the unit matrix, and the */
/* > unitary matrix V is returned;
 */
/* > = 'N': V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* > JOBQ is CHARACTER*1 */
/* > = 'Q': Q must contain a unitary matrix Q1 on entry, and */
/* > the product Q1*Q is returned;
 */
/* > = 'I': Q is initialized to the unit matrix, and the */
/* > unitary matrix Q is returned;
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
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* > L is INTEGER */
/* > */
/* > K and L specify the subblocks in the input matrices A and B: */
/* > A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,,N-L+1:N) */
/* > of A and B, whose GSVD is going to be computed by ZTGSJA. */
/* > See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular */
/* > matrix R or part of R. See Purpose for details. */
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
/* > On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains */
/* > a part of R. See Purpose for details. */
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
/* > TOLA is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] TOLB */
/* > \verbatim */
/* > TOLB is DOUBLE PRECISION */
/* > */
/* > TOLA and TOLB are the convergence criteria for the Jacobi- */
/* > Kogbetliantz iteration procedure. Generally, they are the */
/* > same as used in the preprocessing step, say */
/* > TOLA = MAX(M,N)*norm(A)*MAZHEPS, */
/* > TOLB = MAX(P,N)*norm(B)*MAZHEPS. */
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
/* > ALPHA(K+1:K+L) = diag(C), */
/* > BETA(K+1:K+L) = diag(S), */
/* > or if M-K-L < 0, */
/* > ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0 */
/* > BETA(K+1:M) = S, BETA(M+1:K+L) = 1. */
/* > Furthermore, if K+L < N, */
/* > ALPHA(K+L+1:N) = 0 and */
/* > BETA(K+L+1:N) = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* > U is COMPLEX*16 array, dimension (LDU,M) */
/* > On entry, if JOBU = 'U', U must contain a matrix U1 (usually */
/* > the unitary matrix returned by ZGGSVP). */
/* > On exit, */
/* > if JOBU = 'I', U contains the unitary matrix U;
 */
/* > if JOBU = 'U', U contains the product U1*U. */
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
/* > \param[in,out] V */
/* > \verbatim */
/* > V is COMPLEX*16 array, dimension (LDV,P) */
/* > On entry, if JOBV = 'V', V must contain a matrix V1 (usually */
/* > the unitary matrix returned by ZGGSVP). */
/* > On exit, */
/* > if JOBV = 'I', V contains the unitary matrix V;
 */
/* > if JOBV = 'V', V contains the product V1*V. */
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
/* > \param[in,out] Q */
/* > \verbatim */
/* > Q is COMPLEX*16 array, dimension (LDQ,N) */
/* > On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually */
/* > the unitary matrix returned by ZGGSVP). */
/* > On exit, */
/* > if JOBQ = 'I', Q contains the unitary matrix Q;
 */
/* > if JOBQ = 'Q', Q contains the product Q1*Q. */
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
/* > WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] NCYCLE */
/* > \verbatim */
/* > NCYCLE is INTEGER */
/* > The number of cycles required for convergence. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value. */
/* > = 1: the procedure does not converge after MAXIT cycles. */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > MAXIT INTEGER */
/* > MAXIT specifies the total loops that the iterative procedure */
/* > may take. If after MAXIT cycles, the routine fails to */
/* > converge, we return INFO = 1. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complex16OTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > ZTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce */
/* > fla_min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L */
/* > matrix B13 to the form: */
/* > */
/* > U1**H *A13*Q1 = C1*R1;
V1**H *B13*Q1 = S1*R1, */
/* > */
/* > where U1, V1 and Q1 are unitary matrix. */
/* > C1 and S1 are diagonal matrices satisfying */
/* > */
/* > C1**2 + S1**2 = I, */
/* > */
/* > and R1 is an L-by-L nonsingular upper triangular matrix. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void ztgsja_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k,
             integer *l, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb,
             doublereal *tola, doublereal *tolb, doublereal *alpha, doublereal *beta,
             doublecomplex *u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q,
             integer *ldq, doublecomplex *work, integer *ncycle, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ztgsja inputs: jobu %c, jobv %c, jobq %c, m %" FLA_IS ", p %" FLA_IS
                      ", n %" FLA_IS ", k %" FLA_IS ", l %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS
                      ", ldu %" FLA_IS ", ldv %" FLA_IS ", ldq %" FLA_IS "",
                      *jobu, *jobv, *jobq, *m, *p, *n, *k, *l, *lda, *ldb, *ldu, *ldv, *ldq);
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, u_offset, v_dim1,
        v_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;
    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublereal a1, b1, a3, b3;
    doublecomplex a2, b2;
    doublereal csq, csu, csv;
    doublecomplex snq;
    doublereal rwk;
    doublecomplex snu, snv;
    extern /* Subroutine */
        void
        zrot_(integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *,
              doublecomplex *);
    doublereal gamma;
    extern logical lsame_(char *, char *, integer, integer);
    logical initq, initu, initv, wantq, upper;
    doublereal error, ssmin;
    logical wantu, wantv;
    extern /* Subroutine */
        void
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        zlags2_(logical *, doublereal *, doublecomplex *, doublereal *, doublereal *,
                doublecomplex *, doublereal *, doublereal *, doublecomplex *, doublereal *,
                doublecomplex *, doublereal *, doublecomplex *);
    integer kcycle;
    extern /* Subroutine */
        void
        dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len),
        zdscal_(integer *, doublereal *, doublecomplex *, integer *),
        zlapll_(integer *, doublecomplex *, integer *, doublecomplex *, integer *, doublereal *),
        zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
                integer *);
    doublereal hugenum;
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
    /* Function Body */
    hugenum = 1.7976931348623157e308;
    initu = lsame_(jobu, "I", 1, 1);
    wantu = initu || lsame_(jobu, "U", 1, 1);
    initv = lsame_(jobv, "I", 1, 1);
    wantv = initv || lsame_(jobv, "V", 1, 1);
    initq = lsame_(jobq, "I", 1, 1);
    wantq = initq || lsame_(jobq, "Q", 1, 1);
    *info = 0;
    if(!(initu || wantu || lsame_(jobu, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(initv || wantv || lsame_(jobv, "N", 1, 1)))
    {
        *info = -2;
    }
    else if(!(initq || wantq || lsame_(jobq, "N", 1, 1)))
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
        *info = -10;
    }
    else if(*ldb < fla_max(1, *p))
    {
        *info = -12;
    }
    else if(*ldu < 1 || wantu && *ldu < *m)
    {
        *info = -18;
    }
    else if(*ldv < 1 || wantv && *ldv < *p)
    {
        *info = -20;
    }
    else if(*ldq < 1 || wantq && *ldq < *n)
    {
        *info = -22;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTGSJA", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Initialize U, V and Q, if necessary */
    if(initu)
    {
        zlaset_("Full", m, m, &c_b1, &c_b2, &u[u_offset], ldu);
    }
    if(initv)
    {
        zlaset_("Full", p, p, &c_b1, &c_b2, &v[v_offset], ldv);
    }
    if(initq)
    {
        zlaset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq);
    }
    /* Loop until convergence */
    upper = FALSE_;
    for(kcycle = 1; kcycle <= 40; ++kcycle)
    {
        upper = !upper;
        i__1 = *l - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = *l;
            for(j = i__ + 1; j <= i__2; ++j)
            {
                a1 = 0.;
                a2.r = 0.;
                a2.i = 0.; // , expr subst
                a3 = 0.;
                if(*k + i__ <= *m)
                {
                    i__3 = *k + i__ + (*n - *l + i__) * a_dim1;
                    a1 = a[i__3].r;
                }
                if(*k + j <= *m)
                {
                    i__3 = *k + j + (*n - *l + j) * a_dim1;
                    a3 = a[i__3].r;
                }
                i__3 = i__ + (*n - *l + i__) * b_dim1;
                b1 = b[i__3].r;
                i__3 = j + (*n - *l + j) * b_dim1;
                b3 = b[i__3].r;
                if(upper)
                {
                    if(*k + i__ <= *m)
                    {
                        i__3 = *k + i__ + (*n - *l + j) * a_dim1;
                        a2.r = a[i__3].r;
                        a2.i = a[i__3].i; // , expr subst
                    }
                    i__3 = i__ + (*n - *l + j) * b_dim1;
                    b2.r = b[i__3].r;
                    b2.i = b[i__3].i; // , expr subst
                }
                else
                {
                    if(*k + j <= *m)
                    {
                        i__3 = *k + j + (*n - *l + i__) * a_dim1;
                        a2.r = a[i__3].r;
                        a2.i = a[i__3].i; // , expr subst
                    }
                    i__3 = j + (*n - *l + i__) * b_dim1;
                    b2.r = b[i__3].r;
                    b2.i = b[i__3].i; // , expr subst
                }
                zlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &csv, &snv, &csq, &snq);
                /* Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A */
                if(*k + j <= *m)
                {
                    d_cnjg(&z__1, &snu);
                    zrot_(l, &a[*k + j + (*n - *l + 1) * a_dim1], lda,
                          &a[*k + i__ + (*n - *l + 1) * a_dim1], lda, &csu, &z__1);
                }
                /* Update I-th and J-th rows of matrix B: V**H *B */
                d_cnjg(&z__1, &snv);
                zrot_(l, &b[j + (*n - *l + 1) * b_dim1], ldb, &b[i__ + (*n - *l + 1) * b_dim1], ldb,
                      &csv, &z__1);
                /* Update (N-L+I)-th and (N-L+J)-th columns of matrices */
                /* A and B: A*Q and B*Q */
                /* Computing MIN */
                i__4 = *k + *l;
                i__3 = fla_min(i__4, *m);
                zrot_(&i__3, &a[(*n - *l + j) * a_dim1 + 1], &c__1,
                      &a[(*n - *l + i__) * a_dim1 + 1], &c__1, &csq, &snq);
                zrot_(l, &b[(*n - *l + j) * b_dim1 + 1], &c__1, &b[(*n - *l + i__) * b_dim1 + 1],
                      &c__1, &csq, &snq);
                if(upper)
                {
                    if(*k + i__ <= *m)
                    {
                        i__3 = *k + i__ + (*n - *l + j) * a_dim1;
                        a[i__3].r = 0.;
                        a[i__3].i = 0.; // , expr subst
                    }
                    i__3 = i__ + (*n - *l + j) * b_dim1;
                    b[i__3].r = 0.;
                    b[i__3].i = 0.; // , expr subst
                }
                else
                {
                    if(*k + j <= *m)
                    {
                        i__3 = *k + j + (*n - *l + i__) * a_dim1;
                        a[i__3].r = 0.;
                        a[i__3].i = 0.; // , expr subst
                    }
                    i__3 = j + (*n - *l + i__) * b_dim1;
                    b[i__3].r = 0.;
                    b[i__3].i = 0.; // , expr subst
                }
                /* Ensure that the diagonal elements of A and B are real. */
                if(*k + i__ <= *m)
                {
                    i__3 = *k + i__ + (*n - *l + i__) * a_dim1;
                    i__4 = *k + i__ + (*n - *l + i__) * a_dim1;
                    d__1 = a[i__4].r;
                    a[i__3].r = d__1;
                    a[i__3].i = 0.; // , expr subst
                }
                if(*k + j <= *m)
                {
                    i__3 = *k + j + (*n - *l + j) * a_dim1;
                    i__4 = *k + j + (*n - *l + j) * a_dim1;
                    d__1 = a[i__4].r;
                    a[i__3].r = d__1;
                    a[i__3].i = 0.; // , expr subst
                }
                i__3 = i__ + (*n - *l + i__) * b_dim1;
                i__4 = i__ + (*n - *l + i__) * b_dim1;
                d__1 = b[i__4].r;
                b[i__3].r = d__1;
                b[i__3].i = 0.; // , expr subst
                i__3 = j + (*n - *l + j) * b_dim1;
                i__4 = j + (*n - *l + j) * b_dim1;
                d__1 = b[i__4].r;
                b[i__3].r = d__1;
                b[i__3].i = 0.; // , expr subst
                /* Update unitary matrices U, V, Q, if desired. */
                if(wantu && *k + j <= *m)
                {
                    zrot_(m, &u[(*k + j) * u_dim1 + 1], &c__1, &u[(*k + i__) * u_dim1 + 1], &c__1,
                          &csu, &snu);
                }
                if(wantv)
                {
                    zrot_(p, &v[j * v_dim1 + 1], &c__1, &v[i__ * v_dim1 + 1], &c__1, &csv, &snv);
                }
                if(wantq)
                {
                    zrot_(n, &q[(*n - *l + j) * q_dim1 + 1], &c__1,
                          &q[(*n - *l + i__) * q_dim1 + 1], &c__1, &csq, &snq);
                }
                /* L10: */
            }
            /* L20: */
        }
        if(!upper)
        {
            /* The matrices A13 and B13 were lower triangular at the start */
            /* of the cycle, and are now upper triangular. */
            /* Convergence test: test the parallelism of the corresponding */
            /* rows of A and B. */
            error = 0.;
            /* Computing MIN */
            i__2 = *l;
            i__3 = *m - *k; // , expr subst
            i__1 = fla_min(i__2, i__3);
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = *l - i__ + 1;
                zcopy_(&i__2, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda, &work[1], &c__1);
                i__2 = *l - i__ + 1;
                zcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &work[*l + 1], &c__1);
                i__2 = *l - i__ + 1;
                zlapll_(&i__2, &work[1], &c__1, &work[*l + 1], &c__1, &ssmin);
                error = fla_max(error, ssmin);
                /* L30: */
            }
            if(f2c_dabs(error) <= fla_min(*tola, *tolb))
            {
                goto L50;
            }
        }
        /* End of cycle loop */
        /* L40: */
    }
    /* The algorithm has not converged after MAXIT cycles. */
    *info = 1;
    goto L100;
L50: /* If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged. */
    /* Compute the generalized singular value pairs (ALPHA, BETA), and */
    /* set the triangular matrix R to array A. */
    i__1 = *k;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        alpha[i__] = 1.;
        beta[i__] = 0.;
        /* L60: */
    }
    /* Computing MIN */
    i__2 = *l;
    i__3 = *m - *k; // , expr subst
    i__1 = fla_min(i__2, i__3);
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = *k + i__ + (*n - *l + i__) * a_dim1;
        a1 = a[i__2].r;
        i__2 = i__ + (*n - *l + i__) * b_dim1;
        b1 = b[i__2].r;
        gamma = b1 / a1;
        if(gamma <= hugenum && gamma >= -hugenum)
        {
            if(gamma < 0.)
            {
                i__2 = *l - i__ + 1;
                zdscal_(&i__2, &c_b39, &b[i__ + (*n - *l + i__) * b_dim1], ldb);
                if(wantv)
                {
                    zdscal_(p, &c_b39, &v[i__ * v_dim1 + 1], &c__1);
                }
            }
            d__1 = f2c_dabs(gamma);
            dlartg_(&d__1, &c_b42, &beta[*k + i__], &alpha[*k + i__], &rwk);
            if(alpha[*k + i__] >= beta[*k + i__])
            {
                i__2 = *l - i__ + 1;
                d__1 = 1. / alpha[*k + i__];
                zdscal_(&i__2, &d__1, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda);
            }
            else
            {
                i__2 = *l - i__ + 1;
                d__1 = 1. / beta[*k + i__];
                zdscal_(&i__2, &d__1, &b[i__ + (*n - *l + i__) * b_dim1], ldb);
                i__2 = *l - i__ + 1;
                zcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb,
                       &a[*k + i__ + (*n - *l + i__) * a_dim1], lda);
            }
        }
        else
        {
            alpha[*k + i__] = 0.;
            beta[*k + i__] = 1.;
            i__2 = *l - i__ + 1;
            zcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb,
                   &a[*k + i__ + (*n - *l + i__) * a_dim1], lda);
        }
        /* L70: */
    }
    /* Post-assignment */
    i__1 = *k + *l;
    for(i__ = *m + 1; i__ <= i__1; ++i__)
    {
        alpha[i__] = 0.;
        beta[i__] = 1.;
        /* L80: */
    }
    if(*k + *l < *n)
    {
        i__1 = *n;
        for(i__ = *k + *l + 1; i__ <= i__1; ++i__)
        {
            alpha[i__] = 0.;
            beta[i__] = 0.;
            /* L90: */
        }
    }
L100:
    *ncycle = kcycle;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZTGSJA */
}
/* ztgsja_ */
