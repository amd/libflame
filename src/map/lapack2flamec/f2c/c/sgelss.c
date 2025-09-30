/* ./sgelss.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__6 = 6;
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c__1 = 1;
static aocl_int64_t c__0 = 0;
static real c_b50 = 0.f;
static real c_b83 = 1.f;
/* > \brief <b> SGELSS solves overdetermined or underdetermined systems for GE matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGELSS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelss.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelss.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelss.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, */
/* WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LDB, LWORK, M, N, NRHS, RANK */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), B( LDB, * ), S( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGELSS computes the minimum norm solution to a real linear least */
/* > squares problem: */
/* > */
/* > Minimize 2-norm(| b - A*x |). */
/* > */
/* > using the singular value decomposition (SVD) of A. A is an M-by-N */
/* > matrix which may be rank-deficient. */
/* > */
/* > Several right hand side vectors b and solution vectors x can be */
/* > handled in a single call;
they are stored as the columns of the */
/* > M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix */
/* > X. */
/* > */
/* > The effective rank of A is determined by treating as zero those */
/* > singular values which are less than RCOND times the largest singular */
/* > value. */
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
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, the first fla_min(m,n) rows of A are overwritten with */
/* > its right singular vectors, stored rowwise. */
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
/* > On entry, the M-by-NRHS right hand side matrix B. */
/* > On exit, B is overwritten by the N-by-NRHS solution */
/* > matrix X. If m >= n and RANK = n, the residual */
/* > sum-of-squares for the solution in the i-th column is given */
/* > by the sum of squares of elements n+1:m in that column. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,fla_max(M,N)). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is REAL array, dimension (fla_min(M,N)) */
/* > The singular values of A in decreasing order. */
/* > The condition number of A in the 2-norm = S(1)/S(fla_min(m,n)). */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > RCOND is used to determine the effective rank of A. */
/* > Singular values S(i) <= RCOND*S(1) are treated as zero. */
/* > If RCOND < 0, machine precision is used instead. */
/* > \endverbatim */
/* > */
/* > \param[out] RANK */
/* > \verbatim */
/* > RANK is INTEGER */
/* > The effective rank of A, i.e., the number of singular values */
/* > which are greater than RCOND*S(1). */
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
/* > The dimension of the array WORK. LWORK >= 1, and also: */
/* > LWORK >= 3*fla_min(M,N) + fla_max( 2*fla_min(M,N), fla_max(M,N), NRHS ) */
/* > For good performance, LWORK should generally be larger. */
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
/* > > 0: the algorithm for computing the SVD failed to converge;
 */
/* > if INFO = i, i off-diagonal elements of an intermediate */
/* > bidiagonal form did not converge to zero. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup gelss */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sgelss_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *nrhs, real *a, aocl_int_t *lda, real *b,
             aocl_int_t *ldb, real *s, real *rcond, aocl_int_t *rank, real *work, aocl_int_t *lwork,
             aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t rank_64 = *rank;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgelss(&m_64, &n_64, &nrhs_64, a, &lda_64, b, &ldb_64, s, rcond, &rank_64, work,
                       &lwork_64, &info_64);

    *rank = (aocl_int_t)rank_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sgelss(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs, real *a,
                        aocl_int64_t *lda, real *b, aocl_int64_t *ldb, real *s, real *rcond,
                        aocl_int64_t *rank, real *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgelss inputs: m %" FLA_IS ", n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS ", rank %" FLA_IS "",
                      *m, *n, *nrhs, *lda, *ldb, *rank);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    real r__1;
    /* Local variables */
    aocl_int64_t i__, bl, ie, il, mm;
    real dum[1], eps, thr, anrm, bnrm;
    aocl_int64_t itau, lwork_sgebrd__, lwork_sgeqrf__, lwork_sorgbr__, lwork_sormbr__,
        lwork_sormlq__, iascl, ibscl, lwork_sormqr__, chunk;
    real sfmin;
    aocl_int64_t minmn, maxmn;
    aocl_int64_t itaup, itauq;
    aocl_int64_t mnthr, iwork;
    aocl_int64_t bdspac;
    real bignum;
    aocl_int64_t ldwork;
    aocl_int64_t minwrk, maxwrk;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --s;
    --work;
    /* Function Body */
    *info = 0;
    minmn = fla_min(*m, *n);
    maxmn = fla_max(*m, *n);
    lquery = *lwork == -1;
    mnthr = 0;
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
    else if(*ldb < fla_max(1, maxmn))
    {
        *info = -7;
    }
    /* Compute workspace */
    /* (Note: Comments in the code beginning "Workspace:" describe the */
    /* minimal amount of workspace needed at that point in the code, */
    /* as well as the preferred amount for good performance. */
    /* NB refers to the optimal block size for the immediately */
    /* following subroutine, as returned by ILAENV.) */
    if(*info == 0)
    {
        minwrk = 1;
        maxwrk = 1;
        if(minmn > 0)
        {
            mm = *m;
            mnthr = aocl_lapack_ilaenv(&c__6, "SGELSS", " ", m, n, nrhs, &c_n1);
            if(*m >= *n && *m >= mnthr)
            {
                /* Path 1a - overdetermined, with many more rows than */
                /* columns */
                /* Compute space needed for SGEQRF */
                aocl_lapack_sgeqrf(m, n, &a[a_offset], lda, dum, dum, &c_n1, info);
                lwork_sgeqrf__ = (integer)dum[0];
                /* Compute space needed for SORMQR */
                aocl_lapack_sormqr("L", "T", m, nrhs, n, &a[a_offset], lda, dum, &b[b_offset], ldb,
                                   dum, &c_n1, info);
                lwork_sormqr__ = (integer)dum[0];
                mm = *n;
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + lwork_sgeqrf__; // , expr subst
                maxwrk = fla_max(i__1, i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n + lwork_sormqr__; // , expr subst
                maxwrk = fla_max(i__1, i__2);
            }
            if(*m >= *n)
            {
                /* Path 1 - overdetermined or exactly determined */
                /* Compute workspace needed for SBDSQR */
                /* Computing MAX */
                i__1 = 1;
                i__2 = *n * 5; // , expr subst
                bdspac = fla_max(i__1, i__2);
                /* Compute space needed for SGEBRD */
                aocl_lapack_sgebrd(&mm, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
                                   info);
                lwork_sgebrd__ = (integer)dum[0];
                /* Compute space needed for SORMBR */
                aocl_lapack_sormbr("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, dum,
                                   &b[b_offset], ldb, dum, &c_n1, info);
                lwork_sormbr__ = (integer)dum[0];
                /* Compute space needed for SORGBR */
                aocl_lapack_sorgbr("P", n, n, n, &a[a_offset], lda, dum, dum, &c_n1, info);
                lwork_sorgbr__ = (integer)dum[0];
                /* Compute total workspace needed */
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n * 3 + lwork_sgebrd__; // , expr subst
                maxwrk = fla_max(i__1, i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n * 3 + lwork_sormbr__; // , expr subst
                maxwrk = fla_max(i__1, i__2);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n * 3 + lwork_sorgbr__; // , expr subst
                maxwrk = fla_max(i__1, i__2);
                maxwrk = fla_max(maxwrk, bdspac);
                /* Computing MAX */
                i__1 = maxwrk;
                i__2 = *n * *nrhs; // , expr subst
                maxwrk = fla_max(i__1, i__2);
                /* Computing MAX */
                i__1 = *n * 3 + mm;
                i__2 = *n * 3 + *nrhs;
                i__1 = fla_max(i__1, i__2); // ; expr subst
                minwrk = fla_max(i__1, bdspac);
                maxwrk = fla_max(minwrk, maxwrk);
            }
            if(*n > *m)
            {
                /* Compute workspace needed for SBDSQR */
                /* Computing MAX */
                i__1 = 1;
                i__2 = *m * 5; // , expr subst
                bdspac = fla_max(i__1, i__2);
                /* Computing MAX */
                i__1 = *m * 3 + *nrhs;
                i__2 = *m * 3 + *n;
                i__1 = fla_max(i__1, i__2); // ; expr subst
                minwrk = fla_max(i__1, bdspac);
                if(*n >= mnthr)
                {
                    /* Path 2a - underdetermined, with many more columns */
                    /* than rows */
                    /* Compute space needed for SGEBRD */
                    aocl_lapack_sgebrd(m, m, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
                                       info);
                    lwork_sgebrd__ = (integer)dum[0];
                    /* Compute space needed for SORMBR */
                    aocl_lapack_sormbr("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, dum,
                                       &b[b_offset], ldb, dum, &c_n1, info);
                    lwork_sormbr__ = (integer)dum[0];
                    /* Compute space needed for SORGBR */
                    aocl_lapack_sorgbr("P", m, m, m, &a[a_offset], lda, dum, dum, &c_n1, info);
                    lwork_sorgbr__ = (integer)dum[0];
                    /* Compute space needed for SORMLQ */
                    aocl_lapack_sormlq("L", "T", n, nrhs, m, &a[a_offset], lda, dum, &b[b_offset],
                                       ldb, dum, &c_n1, info);
                    lwork_sormlq__ = (integer)dum[0];
                    /* Compute total workspace needed */
                    maxwrk = *m + *m * aocl_lapack_ilaenv(&c__1, "SGELQF", " ", m, n, &c_n1, &c_n1);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * *m + (*m << 2) + lwork_sgebrd__; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * *m + (*m << 2) + lwork_sormbr__; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * *m + (*m << 2) + lwork_sorgbr__; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * *m + *m + bdspac; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    if(*nrhs > 1)
                    {
                        /* Computing MAX */
                        i__1 = maxwrk;
                        i__2 = *m * *m + *m + *m * *nrhs; // , expr subst
                        maxwrk = fla_max(i__1, i__2);
                    }
                    else
                    {
                        /* Computing MAX */
                        i__1 = maxwrk;
                        i__2 = *m * *m + (*m << 1); // , expr subst
                        maxwrk = fla_max(i__1, i__2);
                    }
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m + lwork_sormlq__; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                }
                else
                {
                    /* Path 2 - underdetermined */
                    /* Compute space needed for SGEBRD */
                    aocl_lapack_sgebrd(m, n, &a[a_offset], lda, &s[1], dum, dum, dum, dum, &c_n1,
                                       info);
                    lwork_sgebrd__ = (integer)dum[0];
                    /* Compute space needed for SORMBR */
                    aocl_lapack_sormbr("Q", "L", "T", m, nrhs, m, &a[a_offset], lda, dum,
                                       &b[b_offset], ldb, dum, &c_n1, info);
                    lwork_sormbr__ = (integer)dum[0];
                    /* Compute space needed for SORGBR */
                    aocl_lapack_sorgbr("P", m, n, m, &a[a_offset], lda, dum, dum, &c_n1, info);
                    lwork_sorgbr__ = (integer)dum[0];
                    maxwrk = *m * 3 + lwork_sgebrd__;
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * 3 + lwork_sormbr__; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *m * 3 + lwork_sorgbr__; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                    maxwrk = fla_max(maxwrk, bdspac);
                    /* Computing MAX */
                    i__1 = maxwrk;
                    i__2 = *n * *nrhs; // , expr subst
                    maxwrk = fla_max(i__1, i__2);
                }
            }
            maxwrk = fla_max(minwrk, maxwrk);
        }
        work[1] = aocl_lapack_sroundup_lwork(&maxwrk);
        if(*lwork < minwrk && !lquery)
        {
            *info = -12;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SGELSS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        *rank = 0;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine parameters */
    eps = slamch_("P");
    sfmin = slamch_("S");
    smlnum = sfmin / eps;
    bignum = 1.f / smlnum;
    /* Scale A if max element outside range [SMLNUM,BIGNUM] */
    anrm = aocl_lapack_slange("M", m, n, &a[a_offset], lda, &work[1]);
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
        aocl_lapack_slaset("F", &i__1, nrhs, &c_b50, &c_b50, &b[b_offset], ldb);
        aocl_lapack_slaset("F", &minmn, &c__1, &c_b50, &c_b50, &s[1], &minmn);
        *rank = 0;
        goto L70;
    }
    /* Scale B if max element outside range [SMLNUM,BIGNUM] */
    bnrm = aocl_lapack_slange("M", m, nrhs, &b[b_offset], ldb, &work[1]);
    ibscl = 0;
    if(bnrm > 0.f && bnrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        aocl_lapack_slascl("G", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[b_offset], ldb, info);
        ibscl = 1;
    }
    else if(bnrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        aocl_lapack_slascl("G", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[b_offset], ldb, info);
        ibscl = 2;
    }
    /* Overdetermined case */
    if(*m >= *n)
    {
        /* Path 1 - overdetermined or exactly determined */
        mm = *m;
        if(*m >= mnthr)
        {
            /* Path 1a - overdetermined, with many more rows than columns */
            mm = *n;
            itau = 1;
            iwork = itau + *n;
            /* Compute A=Q*R */
            /* (Workspace: need 2*N, prefer N+N*NB) */
            i__1 = *lwork - iwork + 1;
            aocl_lapack_sgeqrf(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__1, info);
            /* Multiply B by transpose(Q) */
            /* (Workspace: need N+NRHS, prefer N+NRHS*NB) */
            i__1 = *lwork - iwork + 1;
            aocl_lapack_sormqr("L", "T", m, nrhs, n, &a[a_offset], lda, &work[itau], &b[b_offset],
                               ldb, &work[iwork], &i__1, info);
            /* Zero out below R */
            if(*n > 1)
            {
                i__1 = *n - 1;
                i__2 = *n - 1;
                aocl_lapack_slaset("L", &i__1, &i__2, &c_b50, &c_b50, &a[a_dim1 + 2], lda);
            }
        }
        ie = 1;
        itauq = ie + *n;
        itaup = itauq + *n;
        iwork = itaup + *n;
        /* Bidiagonalize R in A */
        /* (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB) */
        i__1 = *lwork - iwork + 1;
        aocl_lapack_sgebrd(&mm, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq], &work[itaup],
                           &work[iwork], &i__1, info);
        /* Multiply B by transpose of left bidiagonalizing vectors of R */
        /* (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB) */
        i__1 = *lwork - iwork + 1;
        aocl_lapack_sormbr("Q", "L", "T", &mm, nrhs, n, &a[a_offset], lda, &work[itauq],
                           &b[b_offset], ldb, &work[iwork], &i__1, info);
        /* Generate right bidiagonalizing vectors of R in A */
        /* (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB) */
        i__1 = *lwork - iwork + 1;
        aocl_lapack_sorgbr("P", n, n, n, &a[a_offset], lda, &work[itaup], &work[iwork], &i__1,
                           info);
        iwork = ie + *n;
        /* Perform bidiagonal QR iteration */
        /* multiply B by transpose of left singular vectors */
        /* compute right singular vectors in A */
        /* (Workspace: need BDSPAC) */
        aocl_lapack_sbdsqr("U", n, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], lda, dum, &c__1,
                           &b[b_offset], ldb, &work[iwork], info);
        if(*info != 0)
        {
            goto L70;
        }
        /* Multiply B by reciprocals of singular values */
        /* Computing MAX */
        r__1 = *rcond * s[1];
        thr = fla_max(r__1, sfmin);
        if(*rcond < 0.f)
        {
            /* Computing MAX */
            r__1 = eps * s[1];
            thr = fla_max(r__1, sfmin);
        }
        *rank = 0;
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if(s[i__] > thr)
            {
                aocl_lapack_srscl(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                ++(*rank);
            }
            else
            {
                aocl_lapack_slaset("F", &c__1, nrhs, &c_b50, &c_b50, &b[i__ + b_dim1], ldb);
            }
            /* L10: */
        }
        /* Multiply B by right singular vectors */
        /* (Workspace: need N, prefer N*NRHS) */
        if(*lwork >= *ldb * *nrhs && *nrhs > 1)
        {
            aocl_blas_sgemm("T", "N", n, nrhs, n, &c_b83, &a[a_offset], lda, &b[b_offset], ldb,
                            &c_b50, &work[1], ldb);
            aocl_lapack_slacpy("G", n, nrhs, &work[1], ldb, &b[b_offset], ldb);
        }
        else if(*nrhs > 1)
        {
            chunk = *lwork / *n;
            i__1 = *nrhs;
            i__2 = chunk;
            for(i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
            {
                /* Computing MIN */
                i__3 = *nrhs - i__ + 1;
                bl = fla_min(i__3, chunk);
                aocl_blas_sgemm("T", "N", n, &bl, n, &c_b83, &a[a_offset], lda,
                                &b[i__ * b_dim1 + 1], ldb, &c_b50, &work[1], n);
                aocl_lapack_slacpy("G", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb);
                /* L20: */
            }
        }
        else if(*nrhs == 1)
        {
            aocl_blas_sgemv("T", n, n, &c_b83, &a[a_offset], lda, &b[b_offset], &c__1, &c_b50,
                            &work[1], &c__1);
            aocl_blas_scopy(n, &work[1], &c__1, &b[b_offset], &c__1);
        }
    }
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__2 = *m, i__1 = (*m << 1) - 4, i__2 = fla_max(i__2, i__1);
        i__2 = fla_max(i__2, *nrhs);
        i__1 = *n - *m * 3; // ; expr subst
        if(*n >= mnthr && *lwork >= (*m << 2) + *m * *m + fla_max(i__2, i__1))
        {
            /* Path 2a - underdetermined, with many more columns than rows */
            /* and sufficient workspace for an efficient algorithm */
            ldwork = *m;
            /* Computing MAX */
            /* Computing MAX */
            i__3 = *m, i__4 = (*m << 1) - 4, i__3 = fla_max(i__3, i__4);
            i__3 = fla_max(i__3, *nrhs);
            i__4 = *n - *m * 3; // ; expr subst
            i__2 = (*m << 2) + *m * *lda + fla_max(i__3, i__4);
            i__1 = *m * *lda + *m + *m * *nrhs; // , expr subst
            if(*lwork >= fla_max(i__2, i__1))
            {
                ldwork = *lda;
            }
            itau = 1;
            iwork = *m + 1;
            /* Compute A=L*Q */
            /* (Workspace: need 2*M, prefer M+M*NB) */
            i__2 = *lwork - iwork + 1;
            aocl_lapack_sgelqf(m, n, &a[a_offset], lda, &work[itau], &work[iwork], &i__2, info);
            il = iwork;
            /* Copy L to WORK(IL), zeroing out above it */
            aocl_lapack_slacpy("L", m, m, &a[a_offset], lda, &work[il], &ldwork);
            i__2 = *m - 1;
            i__1 = *m - 1;
            aocl_lapack_slaset("U", &i__2, &i__1, &c_b50, &c_b50, &work[il + ldwork], &ldwork);
            ie = il + ldwork * *m;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;
            /* Bidiagonalize L in WORK(IL) */
            /* (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB) */
            i__2 = *lwork - iwork + 1;
            aocl_lapack_sgebrd(m, m, &work[il], &ldwork, &s[1], &work[ie], &work[itauq],
                               &work[itaup], &work[iwork], &i__2, info);
            /* Multiply B by transpose of left bidiagonalizing vectors of L */
            /* (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB) */
            i__2 = *lwork - iwork + 1;
            aocl_lapack_sormbr("Q", "L", "T", m, nrhs, m, &work[il], &ldwork, &work[itauq],
                               &b[b_offset], ldb, &work[iwork], &i__2, info);
            /* Generate right bidiagonalizing vectors of R in WORK(IL) */
            /* (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB) */
            i__2 = *lwork - iwork + 1;
            aocl_lapack_sorgbr("P", m, m, m, &work[il], &ldwork, &work[itaup], &work[iwork], &i__2,
                               info);
            iwork = ie + *m;
            /* Perform bidiagonal QR iteration, */
            /* computing right singular vectors of L in WORK(IL) and */
            /* multiplying B by transpose of left singular vectors */
            /* (Workspace: need M*M+M+BDSPAC) */
            aocl_lapack_sbdsqr("U", m, m, &c__0, nrhs, &s[1], &work[ie], &work[il], &ldwork,
                               &a[a_offset], lda, &b[b_offset], ldb, &work[iwork], info);
            if(*info != 0)
            {
                goto L70;
            }
            /* Multiply B by reciprocals of singular values */
            /* Computing MAX */
            r__1 = *rcond * s[1];
            thr = fla_max(r__1, sfmin);
            if(*rcond < 0.f)
            {
                /* Computing MAX */
                r__1 = eps * s[1];
                thr = fla_max(r__1, sfmin);
            }
            *rank = 0;
            i__2 = *m;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                if(s[i__] > thr)
                {
                    aocl_lapack_srscl(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                    ++(*rank);
                }
                else
                {
                    aocl_lapack_slaset("F", &c__1, nrhs, &c_b50, &c_b50, &b[i__ + b_dim1], ldb);
                }
                /* L30: */
            }
            iwork = ie;
            /* Multiply B by right singular vectors of L in WORK(IL) */
            /* (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS) */
            if(*lwork >= *ldb * *nrhs + iwork - 1 && *nrhs > 1)
            {
                aocl_blas_sgemm("T", "N", m, nrhs, m, &c_b83, &work[il], &ldwork, &b[b_offset], ldb,
                                &c_b50, &work[iwork], ldb);
                aocl_lapack_slacpy("G", m, nrhs, &work[iwork], ldb, &b[b_offset], ldb);
            }
            else if(*nrhs > 1)
            {
                chunk = (*lwork - iwork + 1) / *m;
                i__2 = *nrhs;
                i__1 = chunk;
                for(i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
                {
                    /* Computing MIN */
                    i__3 = *nrhs - i__ + 1;
                    bl = fla_min(i__3, chunk);
                    aocl_blas_sgemm("T", "N", m, &bl, m, &c_b83, &work[il], &ldwork,
                                    &b[i__ * b_dim1 + 1], ldb, &c_b50, &work[iwork], m);
                    aocl_lapack_slacpy("G", m, &bl, &work[iwork], m, &b[i__ * b_dim1 + 1], ldb);
                    /* L40: */
                }
            }
            else if(*nrhs == 1)
            {
                aocl_blas_sgemv("T", m, m, &c_b83, &work[il], &ldwork, &b[b_dim1 + 1], &c__1,
                                &c_b50, &work[iwork], &c__1);
                aocl_blas_scopy(m, &work[iwork], &c__1, &b[b_dim1 + 1], &c__1);
            }
            /* Zero out below first M rows of B */
            i__1 = *n - *m;
            aocl_lapack_slaset("F", &i__1, nrhs, &c_b50, &c_b50, &b[*m + 1 + b_dim1], ldb);
            iwork = itau + *m;
            /* Multiply transpose(Q) by B */
            /* (Workspace: need M+NRHS, prefer M+NRHS*NB) */
            i__1 = *lwork - iwork + 1;
            aocl_lapack_sormlq("L", "T", n, nrhs, m, &a[a_offset], lda, &work[itau], &b[b_offset],
                               ldb, &work[iwork], &i__1, info);
        }
        else
        {
            /* Path 2 - remaining underdetermined cases */
            ie = 1;
            itauq = ie + *m;
            itaup = itauq + *m;
            iwork = itaup + *m;
            /* Bidiagonalize A */
            /* (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB) */
            i__1 = *lwork - iwork + 1;
            aocl_lapack_sgebrd(m, n, &a[a_offset], lda, &s[1], &work[ie], &work[itauq],
                               &work[itaup], &work[iwork], &i__1, info);
            /* Multiply B by transpose of left bidiagonalizing vectors */
            /* (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB) */
            i__1 = *lwork - iwork + 1;
            aocl_lapack_sormbr("Q", "L", "T", m, nrhs, n, &a[a_offset], lda, &work[itauq],
                               &b[b_offset], ldb, &work[iwork], &i__1, info);
            /* Generate right bidiagonalizing vectors in A */
            /* (Workspace: need 4*M, prefer 3*M+M*NB) */
            i__1 = *lwork - iwork + 1;
            aocl_lapack_sorgbr("P", m, n, m, &a[a_offset], lda, &work[itaup], &work[iwork], &i__1,
                               info);
            iwork = ie + *m;
            /* Perform bidiagonal QR iteration, */
            /* computing right singular vectors of A in A and */
            /* multiplying B by transpose of left singular vectors */
            /* (Workspace: need BDSPAC) */
            aocl_lapack_sbdsqr("L", m, n, &c__0, nrhs, &s[1], &work[ie], &a[a_offset], lda, dum,
                               &c__1, &b[b_offset], ldb, &work[iwork], info);
            if(*info != 0)
            {
                goto L70;
            }
            /* Multiply B by reciprocals of singular values */
            /* Computing MAX */
            r__1 = *rcond * s[1];
            thr = fla_max(r__1, sfmin);
            if(*rcond < 0.f)
            {
                /* Computing MAX */
                r__1 = eps * s[1];
                thr = fla_max(r__1, sfmin);
            }
            *rank = 0;
            i__1 = *m;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                if(s[i__] > thr)
                {
                    aocl_lapack_srscl(nrhs, &s[i__], &b[i__ + b_dim1], ldb);
                    ++(*rank);
                }
                else
                {
                    aocl_lapack_slaset("F", &c__1, nrhs, &c_b50, &c_b50, &b[i__ + b_dim1], ldb);
                }
                /* L50: */
            }
            /* Multiply B by right singular vectors of A */
            /* (Workspace: need N, prefer N*NRHS) */
            if(*lwork >= *ldb * *nrhs && *nrhs > 1)
            {
                aocl_blas_sgemm("T", "N", n, nrhs, m, &c_b83, &a[a_offset], lda, &b[b_offset], ldb,
                                &c_b50, &work[1], ldb);
                aocl_lapack_slacpy("F", n, nrhs, &work[1], ldb, &b[b_offset], ldb);
            }
            else if(*nrhs > 1)
            {
                chunk = *lwork / *n;
                i__1 = *nrhs;
                i__2 = chunk;
                for(i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
                {
                    /* Computing MIN */
                    i__3 = *nrhs - i__ + 1;
                    bl = fla_min(i__3, chunk);
                    aocl_blas_sgemm("T", "N", n, &bl, m, &c_b83, &a[a_offset], lda,
                                    &b[i__ * b_dim1 + 1], ldb, &c_b50, &work[1], n);
                    aocl_lapack_slacpy("F", n, &bl, &work[1], n, &b[i__ * b_dim1 + 1], ldb);
                    /* L60: */
                }
            }
            else if(*nrhs == 1)
            {
                aocl_blas_sgemv("T", m, n, &c_b83, &a[a_offset], lda, &b[b_offset], &c__1, &c_b50,
                                &work[1], &c__1);
                aocl_blas_scopy(n, &work[1], &c__1, &b[b_offset], &c__1);
            }
        }
    }
    /* Undo scaling */
    if(iascl == 1)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[b_offset], ldb, info);
        aocl_lapack_slascl("G", &c__0, &c__0, &smlnum, &anrm, &minmn, &c__1, &s[1], &minmn, info);
    }
    else if(iascl == 2)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[b_offset], ldb, info);
        aocl_lapack_slascl("G", &c__0, &c__0, &bignum, &anrm, &minmn, &c__1, &s[1], &minmn, info);
    }
    if(ibscl == 1)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[b_offset], ldb, info);
    }
    else if(ibscl == 2)
    {
        aocl_lapack_slascl("G", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[b_offset], ldb, info);
    }
L70:
    work[1] = aocl_lapack_sroundup_lwork(&maxwrk);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGELSS */
}
/* sgelss_ */
