/* ./zgelsy.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublecomplex c_b1 = {0., 0.};
static doublecomplex c_b2 = {1., 0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static integer c__2 = 2;
/* > \brief <b> ZGELSY solves overdetermined or underdetermined systems for GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZGELSY + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelsy.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelsy.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelsy.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/* WORK, LWORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/* DOUBLE PRECISION RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* DOUBLE PRECISION RWORK( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGELSY computes the minimum-norm solution to a complex linear least */
/* > squares problem: */
/* > minimize || A * X - B || */
/* > using a complete orthogonal factorization of A. A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call;
they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/* > matrix X. */
/* > */
/* > The routine first computes a QR factorization with column pivoting: */
/* > A * P = Q * [ R11 R12 ] */
/* > [ 0 R22 ] */
/* > with R11 defined as the largest leading submatrix whose estimated */
/* > condition number is less than 1/RCOND. The order of R11, RANK, */
/* > is the effective rank of A. */
/* > */
/* > Then, R22 is considered to be negligible, and R12 is annihilated */
/* > by unitary transformations from the right, arriving at the */
/* > complete orthogonal factorization: */
/* > A * P = Q * [ T11 0 ] * Z */
/* > [ 0 0 ] */
/* > The minimum-norm solution is then */
/* > X = P * Z**H [ inv(T11)*Q1**H*B ] */
/* > [ 0 ] */
/* > where Q1 consists of the first RANK columns of Q. */
/* > */
/* > This routine is basically identical to the original xGELSX except */
/* > three differences: */
/* > o The permutation of matrix B (the right hand side) is faster and */
/* > more simple. */
/* > o The call to the subroutine xGEQPF has been substituted by the */
/* > the call to the subroutine xGEQP3. This subroutine is a Blas-3 */
/* > version of the QR factorization with column pivoting. */
/* > o Matrix B (the right hand side) is updated with Blas-3. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of */
/* > columns of matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, A has been overwritten by details of its */
/* > complete orthogonal factorization. */
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
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > On entry, the M-by-NRHS right hand side matrix B. */
/* > On exit, the N-by-NRHS solution matrix X. */
/* > If M = 0 or N = 0, B is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,M,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] JPVT */
/* > \verbatim */
/* > JPVT is INTEGER array, dimension (N) */
/* > On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted */
/* > to the front of AP, otherwise column i is a free column. */
/* > On exit, if JPVT(i) = k, then the i-th column of A*P */
/* > was the k-th column of A. */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* > RCOND is DOUBLE PRECISION */
/* > RCOND is used to determine the effective rank of A, which */
/* > is defined as the order of the largest leading triangular */
/* > submatrix R11 in the QR factorization with pivoting of A, */
/* > whose estimated condition number < 1/RCOND. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* > RANK is INTEGER */
/* > The effective rank of A, i.e., the order of the submatrix */
/* > R11. This is the same as the order of the submatrix T11 */
/* > in the complete orthogonal factorization of A. */
/* > If NRHS = 0, RANK = 0 on output. */
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
/* > The unblocked strategy requires that: */
/* > LWORK >= MN + MAX( 2*MN, N+1, MN+NRHS ) */
/* > where MN = fla_min(M,N). */
/* > The block algorithm requires that: */
/* > LWORK >= MN + MAX( 2*MN, NB*(N+1), MN+MN*NB, MN+NB*NRHS ) */
/* > where NB is an upper bound on the blocksize returned */
/* > by ILAENV for the routines ZGEQP3, ZTZRZF, CTZRQF, ZUNMQR, */
/* > and ZUNMRZ. */
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
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gelsy */
/* > \par Contributors: */
/* ================== */
/* > */
/* > A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA \n */
/* > E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* > G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain \n */
/* > */
/* ===================================================================== */
/* Subroutine */
void zgelsy_(integer *m, integer *n, integer *nrhs, doublecomplex *a, integer *lda,
             doublecomplex *b, integer *ldb, integer *jpvt, doublereal *rcond, integer *rank,
             doublecomplex *work, integer *lwork, doublereal *rwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgelsy inputs: m %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *m, *n, *nrhs, *lda, *ldb);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    /* Local variables */
    integer i__, j;
    doublecomplex c1, c2, s1, s2;
    integer nb, mn, nb1, nb2, nb3, nb4;
    doublereal anrm, bnrm, smin, smax;
    integer iascl, ibscl, ismin, ismax;
    doublereal wsize;
    extern /* Subroutine */
        void
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        ztrsm_(char *, char *, char *, char *, integer *, integer *, doublecomplex *,
               doublecomplex *, integer *, doublecomplex *, integer *),
        zlaic1_(integer *, integer *, doublecomplex *, doublereal *, doublecomplex *,
                doublecomplex *, doublereal *, doublecomplex *, doublecomplex *),
        zgeqp3_(integer *, integer *, doublecomplex *, integer *, integer *, doublecomplex *,
                doublecomplex *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, integer *, integer *);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, integer *,
                              doublereal *);
    doublereal bignum;
    extern /* Subroutine */
        void
        zlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublecomplex *, integer *, integer *),
        zlaset_(char *, integer *, integer *, doublecomplex *, doublecomplex *, doublecomplex *,
                integer *);
    doublereal sminpr, smaxpr, smlnum;
    integer lwkopt;
    logical lquery;
    extern /* Subroutine */
        void
        zunmqr_(char *, char *, integer *, integer *, integer *, doublecomplex *, integer *,
                doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *, integer *),
        zunmrz_(char *, char *, integer *, integer *, integer *, integer *, doublecomplex *,
                integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *,
                integer *),
        ztzrzf_(integer *, integer *, doublecomplex *, integer *, doublecomplex *, doublecomplex *,
                integer *, integer *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --jpvt;
    --work;
    --rwork;
    /* Function Body */
    mn = fla_min(*m, *n);
    ismin = mn + 1;
    ismax = (mn << 1) + 1;
    /* Test the input arguments. */
    *info = 0;
    nb1 = ilaenv_(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1);
    nb2 = ilaenv_(&c__1, "ZGERQF", " ", m, n, &c_n1, &c_n1);
    nb3 = ilaenv_(&c__1, "ZUNMQR", " ", m, n, nrhs, &c_n1);
    nb4 = ilaenv_(&c__1, "ZUNMRQ", " ", m, n, nrhs, &c_n1);
    /* Computing MAX */
    i__1 = fla_max(nb1, nb2);
    i__1 = fla_max(i__1, nb3); // , expr subst
    nb = fla_max(i__1, nb4);
    /* Computing MAX */
    i__1 = 1, i__2 = mn + (*n << 1) + nb * (*n + 1);
    i__1 = fla_max(i__1, i__2);
    i__2 = (mn << 1) + nb * *nrhs; // ; expr subst
    lwkopt = fla_max(i__1, i__2);
    z__1.r = (doublereal)lwkopt;
    z__1.i = 0.; // , expr subst
    work[1].r = z__1.r;
    work[1].i = z__1.i; // , expr subst
    lquery = *lwork == -1;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*nrhs < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = fla_max(1, *m);
        if(*ldb < fla_max(i__1, *n))
        {
            *info = -7;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            i__1 = mn << 1, i__2 = *n + 1;
            i__1 = fla_max(i__1, i__2);
            i__2 = mn + *nrhs; // ; expr subst
            if(*lwork < mn + fla_max(i__1, i__2) && !lquery)
            {
                *info = -12;
            }
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGELSY", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    /* Computing MIN */
    i__1 = fla_min(*m, *n);
    if(fla_min(i__1, *nrhs) == 0)
    {
        *rank = 0;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine parameters */
    smlnum = dlamch_("S") / dlamch_("P");
    bignum = 1. / smlnum;
    /* Scale A, B if max entries outside range [SMLNUM,BIGNUM] */
    anrm = zlange_("M", m, n, &a[a_offset], lda, &rwork[1]);
    iascl = 0;
    if(anrm > 0. && anrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        zlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info);
        iascl = 1;
    }
    else if(anrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        zlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info);
        iascl = 2;
    }
    else if(anrm == 0.)
    {
        /* Matrix all zero. Return zero solution. */
        i__1 = fla_max(*m, *n);
        zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        *rank = 0;
        goto L70;
    }
    bnrm = zlange_("M", m, nrhs, &b[b_offset], ldb, &rwork[1]);
    ibscl = 0;
    if(bnrm > 0. && bnrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        zlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb, info);
        ibscl = 1;
    }
    else if(bnrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        zlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb, info);
        ibscl = 2;
    }
    /* Compute QR factorization with column pivoting of A: */
    /* A * P = Q * R */
    i__1 = *lwork - mn;
    zgeqp3_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], &i__1, &rwork[1], info);
    i__1 = mn + 1;
    wsize = mn + work[i__1].r;
    /* complex workspace: MN+NB*(N+1). real workspace 2*N. */
    /* Details of Householder rotations stored in WORK(1:MN). */
    /* Determine RANK using incremental condition estimation */
    i__1 = ismin;
    work[i__1].r = 1.;
    work[i__1].i = 0.; // , expr subst
    i__1 = ismax;
    work[i__1].r = 1.;
    work[i__1].i = 0.; // , expr subst
    smax = z_abs(&a[a_dim1 + 1]);
    smin = smax;
    if(z_abs(&a[a_dim1 + 1]) == 0.)
    {
        *rank = 0;
        i__1 = fla_max(*m, *n);
        zlaset_("F", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        goto L70;
    }
    else
    {
        *rank = 1;
    }
L10:
    if(*rank < mn)
    {
        i__ = *rank + 1;
        zlaic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[i__ + i__ * a_dim1],
                &sminpr, &s1, &c1);
        zlaic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[i__ + i__ * a_dim1],
                &smaxpr, &s2, &c2);
        if(smaxpr * *rcond <= sminpr)
        {
            i__1 = *rank;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                i__2 = ismin + i__ - 1;
                i__3 = ismin + i__ - 1;
                z__1.r = s1.r * work[i__3].r - s1.i * work[i__3].i;
                z__1.i = s1.r * work[i__3].i + s1.i * work[i__3].r; // , expr subst
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
                i__2 = ismax + i__ - 1;
                i__3 = ismax + i__ - 1;
                z__1.r = s2.r * work[i__3].r - s2.i * work[i__3].i;
                z__1.i = s2.r * work[i__3].i + s2.i * work[i__3].r; // , expr subst
                work[i__2].r = z__1.r;
                work[i__2].i = z__1.i; // , expr subst
                /* L20: */
            }
            i__1 = ismin + *rank;
            work[i__1].r = c1.r;
            work[i__1].i = c1.i; // , expr subst
            i__1 = ismax + *rank;
            work[i__1].r = c2.r;
            work[i__1].i = c2.i; // , expr subst
            smin = sminpr;
            smax = smaxpr;
            ++(*rank);
            goto L10;
        }
    }
    /* complex workspace: 3*MN. */
    /* Logically partition R = [ R11 R12 ] */
    /* [ 0 R22 ] */
    /* where R11 = R(1:RANK,1:RANK) */
    /* [R11,R12] = [ T11, 0 ] * Y */
    if(*rank < *n)
    {
        i__1 = *lwork - (mn << 1);
        ztzrzf_(rank, n, &a[a_offset], lda, &work[mn + 1], &work[(mn << 1) + 1], &i__1, info);
    }
    /* complex workspace: 2*MN. */
    /* Details of Householder rotations stored in WORK(MN+1:2*MN) */
    /* B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS) */
    i__1 = *lwork - (mn << 1);
    zunmqr_("Left", "Conjugate transpose", m, nrhs, &mn, &a[a_offset], lda, &work[1], &b[b_offset],
            ldb, &work[(mn << 1) + 1], &i__1, info);
    /* Computing MAX */
    i__1 = (mn << 1) + 1;
    d__1 = wsize;
    d__2 = (mn << 1) + work[i__1].r; // , expr subst
    wsize = fla_max(d__1, d__2);
    /* complex workspace: 2*MN+NB*NRHS. */
    /* B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */
    ztrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b2, &a[a_offset], lda,
           &b[b_offset], ldb);
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *n;
        for(i__ = *rank + 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * b_dim1;
            b[i__3].r = 0.;
            b[i__3].i = 0.; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    /* B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS) */
    if(*rank < *n)
    {
        i__1 = *n - *rank;
        i__2 = *lwork - (mn << 1);
        zunmrz_("Left", "Conjugate transpose", n, nrhs, rank, &i__1, &a[a_offset], lda,
                &work[mn + 1], &b[b_offset], ldb, &work[(mn << 1) + 1], &i__2, info);
    }
    /* complex workspace: 2*MN+NRHS. */
    /* B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = jpvt[i__];
            i__4 = i__ + j * b_dim1;
            work[i__3].r = b[i__4].r;
            work[i__3].i = b[i__4].i; // , expr subst
            /* L50: */
        }
        zcopy_(n, &work[1], &c__1, &b[j * b_dim1 + 1], &c__1);
        /* L60: */
    }
    /* complex workspace: N. */
    /* Undo scaling */
    if(iascl == 1)
    {
        zlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb, info);
        zlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], lda, info);
    }
    else if(iascl == 2)
    {
        zlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb, info);
        zlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], lda, info);
    }
    if(ibscl == 1)
    {
        zlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb, info);
    }
    else if(ibscl == 2)
    {
        zlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb, info);
    }
L70:
    z__1.r = (doublereal)lwkopt;
    z__1.i = 0.; // , expr subst
    work[1].r = z__1.r;
    work[1].i = z__1.i; // , expr subst
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZGELSY */
}
/* zgelsy_ */
