/* ../netlib/dgelsx.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__0 = 0;
static doublereal c_b13 = 0.;
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b36 = 1.;
/* > \brief <b> DGELSX solves overdetermined or underdetermined systems for GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DGELSX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelsx.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelsx.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelsx.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, */
/* WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, M, N, NRHS, RANK */
/* DOUBLE PRECISION RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER JPVT( * ) */
/* DOUBLE PRECISION A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > This routine is deprecated and has been replaced by routine DGELSY. */
/* > */
/* > DGELSX computes the minimum-norm solution to a real linear least */
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
/* > by orthogonal transformations from the right, arriving at the */
/* > complete orthogonal factorization: */
/* > A * P = Q * [ T11 0 ] * Z */
/* > [ 0 0 ] */
/* > The minimum-norm solution is then */
/* > X = P * Z**T [ inv(T11)*Q1**T*B ] */
/* > [ 0 ] */
/* > where Q1 consists of the first RANK columns of Q. */
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
/* > A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* > On entry, the M-by-NRHS right hand side matrix B. */
/* > On exit, the N-by-NRHS solution matrix X. */
/* > If m >= n and RANK = n, the residual sum-of-squares for */
/* > the solution in the i-th column is given by the sum of */
/* > squares of elements N+1:M in that column. */
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
/* > On entry, if JPVT(i) .ne. 0, the i-th column of A is an */
/* > initial column, otherwise it is a free column. Before */
/* > the QR factorization of A, all initial columns are */
/* > permuted to the leading positions;
only the remaining */
/* > free columns are moved as a result of column pivoting */
/* > during the factorization. */
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
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension */
/* > (fla_max( fla_min(M,N)+3*N, 2*fla_min(M,N)+NRHS )), */
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
/* > \date November 2011 */
/* > \ingroup doubleGEsolve */
/* ===================================================================== */
/* Subroutine */
void dgelsx_(integer *m, integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *b,
             integer *ldb, integer *jpvt, doublereal *rcond, integer *rank, doublereal *work,
             integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgelsx inputs: m %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS "",
                      *m, *n, *nrhs, *lda, *ldb);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    doublereal d__1;
    /* Local variables */
    integer i__, j, k;
    doublereal c1, c2, s1, s2, t1, t2;
    integer mn;
    doublereal anrm, bnrm, smin, smax;
    integer iascl, ibscl, ismin, ismax;
    extern /* Subroutine */
        void
        dtrsm_(char *, char *, char *, char *, integer *, integer *, doublereal *, doublereal *,
               integer *, doublereal *, integer *),
        dlaic1_(integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *,
                doublereal *, doublereal *, doublereal *),
        dorm2r_(char *, char *, integer *, integer *, integer *, doublereal *, integer *,
                doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *),
        dlange_(char *, integer *, integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */
        void
        dlascl_(char *, integer *, integer *, doublereal *, doublereal *, integer *, integer *,
                doublereal *, integer *, integer *),
        dgeqpf_(integer *, integer *, doublereal *, integer *, integer *, doublereal *,
                doublereal *, integer *),
        dlaset_(char *, integer *, integer *, doublereal *, doublereal *, doublereal *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    doublereal bignum;
    extern /* Subroutine */
        void
        dlatzm_(char *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *,
                doublereal *, integer *, doublereal *);
    doublereal sminpr, smaxpr, smlnum;
    extern /* Subroutine */
        void
        dtzrqf_(integer *, integer *, doublereal *, integer *, doublereal *, integer *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
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
    /* Function Body */
    mn = fla_min(*m, *n);
    ismin = mn + 1;
    ismax = (mn << 1) + 1;
    /* Test the input arguments. */
    *info = 0;
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
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("DGELSX", &i__1, (ftnlen)6);
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
    /* Scale A, B if max elements outside range [SMLNUM,BIGNUM] */
    anrm = dlange_("M", m, n, &a[a_offset], lda, &work[1]);
    iascl = 0;
    if(anrm > 0. && anrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        dlascl_("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info);
        iascl = 1;
    }
    else if(anrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        dlascl_("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info);
        iascl = 2;
    }
    else if(anrm == 0.)
    {
        /* Matrix all zero. Return zero solution. */
        i__1 = fla_max(*m, *n);
        dlaset_("F", &i__1, nrhs, &c_b13, &c_b13, &b[b_offset], ldb);
        *rank = 0;
        goto L100;
    }
    bnrm = dlange_("M", m, nrhs, &b[b_offset], ldb, &work[1]);
    ibscl = 0;
    if(bnrm > 0. && bnrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        dlascl_("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb, info);
        ibscl = 1;
    }
    else if(bnrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        dlascl_("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb, info);
        ibscl = 2;
    }
    /* Compute QR factorization with column pivoting of A: */
    /* A * P = Q * R */
    dgeqpf_(m, n, &a[a_offset], lda, &jpvt[1], &work[1], &work[mn + 1], info);
    /* workspace 3*N. Details of Householder rotations stored */
    /* in WORK(1:MN). */
    /* Determine RANK using incremental condition estimation */
    work[ismin] = 1.;
    work[ismax] = 1.;
    smax = (d__1 = a[a_dim1 + 1], f2c_dabs(d__1));
    smin = smax;
    if((d__1 = a[a_dim1 + 1], f2c_dabs(d__1)) == 0.)
    {
        *rank = 0;
        i__1 = fla_max(*m, *n);
        dlaset_("F", &i__1, nrhs, &c_b13, &c_b13, &b[b_offset], ldb);
        goto L100;
    }
    else
    {
        *rank = 1;
    }
L10:
    if(*rank < mn)
    {
        i__ = *rank + 1;
        dlaic1_(&c__2, rank, &work[ismin], &smin, &a[i__ * a_dim1 + 1], &a[i__ + i__ * a_dim1],
                &sminpr, &s1, &c1);
        dlaic1_(&c__1, rank, &work[ismax], &smax, &a[i__ * a_dim1 + 1], &a[i__ + i__ * a_dim1],
                &smaxpr, &s2, &c2);
        if(smaxpr * *rcond <= sminpr)
        {
            i__1 = *rank;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                work[ismin + i__ - 1] = s1 * work[ismin + i__ - 1];
                work[ismax + i__ - 1] = s2 * work[ismax + i__ - 1];
                /* L20: */
            }
            work[ismin + *rank] = c1;
            work[ismax + *rank] = c2;
            smin = sminpr;
            smax = smaxpr;
            ++(*rank);
            goto L10;
        }
    }
    /* Logically partition R = [ R11 R12 ] */
    /* [ 0 R22 ] */
    /* where R11 = R(1:RANK,1:RANK) */
    /* [R11,R12] = [ T11, 0 ] * Y */
    if(*rank < *n)
    {
        dtzrqf_(rank, n, &a[a_offset], lda, &work[mn + 1], info);
    }
    /* Details of Householder rotations stored in WORK(MN+1:2*MN) */
    /* B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */
    dorm2r_("Left", "Transpose", m, nrhs, &mn, &a[a_offset], lda, &work[1], &b[b_offset], ldb,
            &work[(mn << 1) + 1], info);
    /* workspace NRHS */
    /* B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS) */
    dtrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b36, &a[a_offset], lda,
           &b[b_offset], ldb);
    i__1 = *n;
    for(i__ = *rank + 1; i__ <= i__1; ++i__)
    {
        i__2 = *nrhs;
        for(j = 1; j <= i__2; ++j)
        {
            b[i__ + j * b_dim1] = 0.;
            /* L30: */
        }
        /* L40: */
    }
    /* B(1:N,1:NRHS) := Y**T * B(1:N,1:NRHS) */
    if(*rank < *n)
    {
        i__1 = *rank;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = *n - *rank + 1;
            dlatzm_("Left", &i__2, nrhs, &a[i__ + (*rank + 1) * a_dim1], lda, &work[mn + i__],
                    &b[i__ + b_dim1], &b[*rank + 1 + b_dim1], ldb, &work[(mn << 1) + 1]);
            /* L50: */
        }
    }
    /* workspace NRHS */
    /* B(1:N,1:NRHS) := P * B(1:N,1:NRHS) */
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            work[(mn << 1) + i__] = 1.;
            /* L60: */
        }
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(work[(mn << 1) + i__] == 1.)
            {
                if(jpvt[i__] != i__)
                {
                    k = i__;
                    t1 = b[k + j * b_dim1];
                    t2 = b[jpvt[k] + j * b_dim1];
                L70:
                    b[jpvt[k] + j * b_dim1] = t1;
                    work[(mn << 1) + k] = 0.;
                    t1 = t2;
                    k = jpvt[k];
                    t2 = b[jpvt[k] + j * b_dim1];
                    if(jpvt[k] != i__)
                    {
                        goto L70;
                    }
                    b[i__ + j * b_dim1] = t1;
                    work[(mn << 1) + k] = 0.;
                }
            }
            /* L80: */
        }
        /* L90: */
    }
    /* Undo scaling */
    if(iascl == 1)
    {
        dlascl_("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb, info);
        dlascl_("U", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[a_offset], lda, info);
    }
    else if(iascl == 2)
    {
        dlascl_("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb, info);
        dlascl_("U", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[a_offset], lda, info);
    }
    if(ibscl == 1)
    {
        dlascl_("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb, info);
    }
    else if(ibscl == 2)
    {
        dlascl_("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb, info);
    }
L100:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DGELSX */
}
/* dgelsx_ */
