/* ./sgelst.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static real c_b12 = 0.f;
static aocl_int64_t c__2 = 2;
static aocl_int64_t c__0 = 0;
/* > \brief <b> SGELST solves overdetermined or underdetermined systems for GE matrices using QR or
 * LQ factori zation with compact WY representation of Q.</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGELST + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelst.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelst.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelst.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGELST( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, */
/* INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDA, LDB, LWORK, M, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGELST solves overdetermined or underdetermined real linear systems */
/* > involving an M-by-N matrix A, or its transpose, using a QR or LQ */
/* > factorization of A with compact WY representation of Q. */
/* > It is assumed that A has full rank. */
/* > */
/* > The following options are provided: */
/* > */
/* > 1. If TRANS = 'N' and m >= n: find the least squares solution of */
/* > an overdetermined system, i.e., solve the least squares problem */
/* > minimize || B - A*X ||. */
/* > */
/* > 2. If TRANS = 'N' and m < n: find the minimum norm solution of */
/* > an underdetermined system A * X = B. */
/* > */
/* > 3. If TRANS = 'T' and m >= n: find the minimum norm solution of */
/* > an underdetermined system A**T * X = B. */
/* > */
/* > 4. If TRANS = 'T' and m < n: find the least squares solution of */
/* > an overdetermined system, i.e., solve the least squares problem */
/* > minimize || B - A**T * X ||. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call;
they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution */
/* > matrix X. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': the linear system involves A;
 */
/* > = 'T': the linear system involves A**T. */
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
/* > The number of columns of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of */
/* > columns of the matrices B and X. NRHS >=0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, */
/* > if M >= N, A is overwritten by details of its QR */
/* > factorization as returned by SGEQRT;
 */
/* > if M < N, A is overwritten by details of its LQ */
/* > factorization as returned by SGELQT. */
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
/* > B is REAL array, dimension (LDB,NRHS) */
/* > On entry, the matrix B of right hand side vectors, stored */
/* > columnwise;
B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/* > if TRANS = 'T'. */
/* > On exit, if INFO = 0, B is overwritten by the solution */
/* > vectors, stored columnwise: */
/* > if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/* > squares solution vectors;
the residual sum of squares for the */
/* > solution in each column is given by the sum of squares of */
/* > elements N+1 to M in that column;
 */
/* > if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* > minimum norm solution vectors;
 */
/* > if TRANS = 'T' and m >= n, rows 1 to M of B contain the */
/* > minimum norm solution vectors;
 */
/* > if TRANS = 'T' and m < n, rows 1 to M of B contain the */
/* > least squares solution vectors;
the residual sum of squares */
/* > for the solution in each column is given by the sum of */
/* > squares of elements M+1 to N in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= MAX(1,M,N). */
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
/* > LWORK >= fla_max( 1, MN + fla_max( MN, NRHS ) ). */
/* > For optimal performance, */
/* > LWORK >= fla_max( 1, (MN + fla_max( MN, NRHS ))*NB ). */
/* > where MN = fla_min(M,N) and NB is the optimum block size. */
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
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the i-th diagonal element of the */
/* > triangular factor of A is zero, so that A does not have */
/* > full rank;
the least squares solution could not be */
/* > computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gelst */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > November 2022, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sgelst_(char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *nrhs, real *a, aocl_int_t *lda,
             real *b, aocl_int_t *ldb, real *work, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgelst(trans, &m_64, &n_64, &nrhs_64, a, &lda_64, b, &ldb_64, work, &lwork_64,
                       &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sgelst(char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, real *a,
                        aocl_int64_t *lda, real *b, aocl_int64_t *ldb, real *work,
                        aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgelst inputs: trans %c, m %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS
                      ", lda %" FLA_IS ", ldb %" FLA_IS "",
                      *trans, *m, *n, *nrhs, *lda, *ldb);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    aocl_int64_t i__, j, nb, mn;
    real anrm, bnrm;
    aocl_int64_t brow;
    logical tpsd;
    aocl_int64_t iascl, ibscl;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t nbmin;
    real rwork[1];
    aocl_int64_t lwopt;
    aocl_int64_t scllen;
    real bignum;
    aocl_int64_t mnnrhs;
    real smlnum;
    logical lquery;
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --work;
    /* Function Body */
    *info = 0;
    mn = fla_min(*m, *n);
    lquery = *lwork == -1;
    if(!(lsame_(trans, "N", 1, 1) || lsame_(trans, "T", 1, 1)))
    {
        *info = -1;
    }
    else if(*m < 0)
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*nrhs < 0)
    {
        *info = -4;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -6;
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = fla_max(1, *m);
        if(*ldb < fla_max(i__1, *n))
        {
            *info = -8;
        }
        else /* if(complicated condition) */
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = mn + fla_max(mn, *nrhs); // , expr subst
            if(*lwork < fla_max(i__1, i__2) && !lquery)
            {
                *info = -10;
            }
        }
    }
    /* Figure out optimal block size and optimal workspace size */
    if(*info == 0 || *info == -10)
    {
        tpsd = TRUE_;
        if(lsame_(trans, "N", 1, 1))
        {
            tpsd = FALSE_;
        }
        nb = aocl_lapack_ilaenv(&c__1, "SGELST", " ", m, n, &c_n1, &c_n1);
        mnnrhs = fla_max(mn, *nrhs);
        /* Computing MAX */
        i__1 = 1;
        i__2 = (mn + mnnrhs) * nb; // , expr subst
        lwopt = fla_max(i__1, i__2);
        work[1] = aocl_lapack_sroundup_lwork(&lwopt);
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SGELST ", &i__1, (ftnlen)6);
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
        i__1 = fla_max(*m, *n);
        aocl_lapack_slaset("Full", &i__1, nrhs, &c_b12, &c_b12, &b[b_offset], ldb);
        work[1] = aocl_lapack_sroundup_lwork(&lwopt);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* *GEQRT and *GELQT routines cannot accept NB larger than fla_min(M,N) */
    if(nb > mn)
    {
        nb = mn;
    }
    /* Determine the block size from the supplied LWORK */
    /* ( at this stage we know that LWORK >= (minimum required workspace, */
    /* but it may be less than optimal) */
    /* Computing MIN */
    i__1 = nb;
    i__2 = *lwork / (mn + mnnrhs); // , expr subst
    nb = fla_min(i__1, i__2);
    /* The minimum value of NB, when blocked code is used */
    /* Computing MAX */
    i__1 = 2;
    i__2 = aocl_lapack_ilaenv(&c__2, "SGELST", " ", m, n, &c_n1, &c_n1); // , expr subst
    nbmin = fla_max(i__1, i__2);
    if(nb < nbmin)
    {
        nb = 1;
    }
    /* Get machine parameters */
    smlnum = slamch_("S") / slamch_("P");
    bignum = 1.f / smlnum;
    /* Scale A, B if max element outside range [SMLNUM,BIGNUM] */
    anrm = aocl_lapack_slange("M", m, n, &a[a_offset], lda, rwork);
    iascl = 0;
    if(anrm > 0.f && anrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        aocl_lapack_slascl("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info);
        iascl = 1;
    }
    else if(anrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        aocl_lapack_slascl("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info);
        iascl = 2;
    }
    else if(anrm == 0.f)
    {
        /* Matrix all zero. Return zero solution. */
        i__1 = fla_max(*m, *n);
        aocl_lapack_slaset("Full", &i__1, nrhs, &c_b12, &c_b12, &b[b_offset], ldb);
        work[1] = aocl_lapack_sroundup_lwork(&lwopt);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    brow = *m;
    if(tpsd)
    {
        brow = *n;
    }
    bnrm = aocl_lapack_slange("M", &brow, nrhs, &b[b_offset], ldb, rwork);
    ibscl = 0;
    if(bnrm > 0.f && bnrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        aocl_lapack_slascl("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 1;
    }
    else if(bnrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        aocl_lapack_slascl("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 2;
    }
    if(*m >= *n)
    {
        /* M > N: */
        /* Compute the blocked QR factorization of A, */
        /* using the compact WY representation of Q, */
        /* workspace at least N, optimally N*NB. */
        aocl_lapack_sgeqrt(m, n, &nb, &a[a_offset], lda, &work[1], &nb, &work[mn * nb + 1], info);
        if(!tpsd)
        {
            /* M > N, A is not transposed: */
            /* Overdetermined system of equations, */
            /* least-squares problem, min || A * X - B ||. */
            /* Compute B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            aocl_lapack_sgemqrt("Left", "Transpose", m, nrhs, n, &nb, &a[a_offset], lda, &work[1],
                                &nb, &b[b_offset], ldb, &work[mn * nb + 1], info);
            /* Compute B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */
            aocl_lapack_strtrs("Upper", "No transpose", "Non-unit", n, nrhs, &a[a_offset], lda,
                               &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            scllen = *n;
        }
        else
        {
            /* M > N, A is transposed: */
            /* Underdetermined system of equations, */
            /* minimum norm solution of A**T * X = B. */
            /* Compute B := inv(R**T) * B in two row blocks of B. */
            /* Block 1: B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */
            aocl_lapack_strtrs("Upper", "Transpose", "Non-unit", n, nrhs, &a[a_offset], lda,
                               &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* Block 2: Zero out all rows below the N-th row in B: */
            /* B(N+1:M,1:NRHS) = ZERO */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *m;
                for(i__ = *n + 1; i__ <= i__2; ++i__)
                {
                    b[i__ + j * b_dim1] = 0.f;
                }
            }
            /* Compute B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            aocl_lapack_sgemqrt("Left", "No transpose", m, nrhs, n, &nb, &a[a_offset], lda,
                                &work[1], &nb, &b[b_offset], ldb, &work[mn * nb + 1], info);
            scllen = *m;
        }
    }
    else
    {
        /* M < N: */
        /* Compute the blocked LQ factorization of A, */
        /* using the compact WY representation of Q, */
        /* workspace at least M, optimally M*NB. */
        aocl_lapack_sgelqt(m, n, &nb, &a[a_offset], lda, &work[1], &nb, &work[mn * nb + 1], info);
        if(!tpsd)
        {
            /* M < N, A is not transposed: */
            /* Underdetermined system of equations, */
            /* minimum norm solution of A * X = B. */
            /* Compute B := inv(L) * B in two row blocks of B. */
            /* Block 1: B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */
            aocl_lapack_strtrs("Lower", "No transpose", "Non-unit", m, nrhs, &a[a_offset], lda,
                               &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* Block 2: Zero out all rows below the M-th row in B: */
            /* B(M+1:N,1:NRHS) = ZERO */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *n;
                for(i__ = *m + 1; i__ <= i__2; ++i__)
                {
                    b[i__ + j * b_dim1] = 0.f;
                }
            }
            /* Compute B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            aocl_lapack_sgemlqt("Left", "Transpose", n, nrhs, m, &nb, &a[a_offset], lda, &work[1],
                                &nb, &b[b_offset], ldb, &work[mn * nb + 1], info);
            scllen = *n;
        }
        else
        {
            /* M < N, A is transposed: */
            /* Overdetermined system of equations, */
            /* least-squares problem, min || A**T * X - B ||. */
            /* Compute B(1:N,1:NRHS) := Q * B(1:N,1:NRHS), */
            /* using the compact WY representation of Q, */
            /* workspace at least NRHS, optimally NRHS*NB. */
            aocl_lapack_sgemlqt("Left", "No transpose", n, nrhs, m, &nb, &a[a_offset], lda,
                                &work[1], &nb, &b[b_offset], ldb, &work[mn * nb + 1], info);
            /* Compute B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */
            aocl_lapack_strtrs("Lower", "Transpose", "Non-unit", m, nrhs, &a[a_offset], lda,
                               &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            scllen = *m;
        }
    }
    /* Undo scaling */
    if(iascl == 1)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
    else if(iascl == 2)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
    if(ibscl == 1)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
    else if(ibscl == 2)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
    work[1] = aocl_lapack_sroundup_lwork(&lwopt);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGELST */
}
/* sgelst_ */
