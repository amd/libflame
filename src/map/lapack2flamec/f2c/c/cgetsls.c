/* ./cgetsls.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b1 = {{0.f}, {0.f}};
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c_n2 = -2;
static aocl_int64_t c__0 = 0;
/* > \brief \b CGETSLS */
/* Definition: */
/* =========== */
/* SUBROUTINE CGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, */
/* $ WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDA, LDB, LWORK, M, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), B( LDB, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGETSLS solves overdetermined or underdetermined scomplex linear systems */
/* > involving an M-by-N matrix A, using a tall skinny QR or short wide LQ */
/* > factorization of A. It is assumed that A has full rank. */
/* > */
/* > */
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
/* > 3. If TRANS = 'C' and m >= n: find the minimum norm solution of */
/* > an undetermined system A**T * X = B. */
/* > */
/* > 4. If TRANS = 'C' and m < n: find the least squares solution of */
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
/* > = 'C': the linear system involves A**H. */
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
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, */
/* > A is overwritten by details of its QR or LQ */
/* > factorization as returned by CGEQR or CGELQ. */
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
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the matrix B of right hand side vectors, stored */
/* > columnwise;
B is M-by-NRHS if TRANS = 'N', or N-by-NRHS */
/* > if TRANS = 'C'. */
/* > On exit, if INFO = 0, B is overwritten by the solution */
/* > vectors, stored columnwise: */
/* > if TRANS = 'N' and m >= n, rows 1 to n of B contain the least */
/* > squares solution vectors. */
/* > if TRANS = 'N' and m < n, rows 1 to N of B contain the */
/* > minimum norm solution vectors;
 */
/* > if TRANS = 'C' and m >= n, rows 1 to M of B contain the */
/* > minimum norm solution vectors;
 */
/* > if TRANS = 'C' and m < n, rows 1 to M of B contain the */
/* > least squares solution vectors. */
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
/* > (workspace) COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) contains optimal (or either minimal */
/* > or optimal, if query was assumed) LWORK. */
/* > See LWORK for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If LWORK = -1 or -2, then a workspace query is assumed. */
/* > If LWORK = -1, the routine calculates optimal size of WORK for the */
/* > optimal performance and returns this value in WORK(1). */
/* > If LWORK = -2, the routine calculates minimal size of WORK and */
/* > returns this value in WORK(1). */
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
/* > \ingroup getsls */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgetsls_(char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *nrhs, scomplex *a,
              aocl_int_t *lda, scomplex *b, aocl_int_t *ldb, scomplex *work, aocl_int_t *lwork,
              aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetsls(trans, &m_64, &n_64, &nrhs_64, a, &lda_64, b, &ldb_64, work, &lwork_64,
                        &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgetsls(char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *nrhs,
                         scomplex *a, aocl_int64_t *lda, scomplex *b, aocl_int64_t *ldb,
                         scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cgetsls inputs: trans %c, m %lld, n %lld, nrhs %lld, lda %lld, ldb %lld, lwork %lld",
             *trans, *m, *n, *nrhs, *lda, *ldb, *lwork);
#else
    snprintf(buffer, 256, "cgetsls inputs: trans %c, m %d, n %d, nrhs %d, lda %d, ldb %d, lwork %d",
             *trans, *m, *n, *nrhs, *lda, *ldb, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    real r__1;
    /* Local variables */
    aocl_int64_t i__, j;
    scomplex tq[5];
    aocl_int64_t lw1, lw2;
    real dum[1];
    aocl_int64_t lwm, lwo;
    real anrm, bnrm;
    logical tran;
    aocl_int64_t brow, tszm, tszo, info2, iascl, ibscl;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t maxmn;
    scomplex workq[1];
    extern real slamch_(char *);
    aocl_int64_t scllen;
    real bignum, smlnum;
    aocl_int64_t wsizem, wsizeo;
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
    maxmn = fla_max(*m, *n);
    tran = lsame_(trans, "C", 1, 1);
    lquery = *lwork == -1 || *lwork == -2;
    if(!(lsame_(trans, "N", 1, 1) || lsame_(trans, "C", 1, 1)))
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
    }
    if(*info == 0)
    {
        /* Determine the optimum and minimum LWORK */
        if(*m >= *n)
        {
            aocl_lapack_cgeqr(m, n, &a[a_offset], lda, tq, &c_n1, workq, &c_n1, &info2);
            tszo = (integer)tq[0].r;
            lwo = (integer)workq[0].r;
            aocl_lapack_cgemqr("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszo, &b[b_offset],
                               ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwo;
            i__2 = (integer)workq[0].r; // , expr subst
            lwo = fla_max(i__1, i__2);
            aocl_lapack_cgeqr(m, n, &a[a_offset], lda, tq, &c_n2, workq, &c_n2, &info2);
            tszm = (integer)tq[0].r;
            lwm = (integer)workq[0].r;
            aocl_lapack_cgemqr("L", trans, m, nrhs, n, &a[a_offset], lda, tq, &tszm, &b[b_offset],
                               ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwm;
            i__2 = (integer)workq[0].r; // , expr subst
            lwm = fla_max(i__1, i__2);
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        }
        else
        {
            aocl_lapack_cgelq(m, n, &a[a_offset], lda, tq, &c_n1, workq, &c_n1, &info2);
            tszo = (integer)tq[0].r;
            lwo = (integer)workq[0].r;
            aocl_lapack_cgemlq("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszo, &b[b_offset],
                               ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwo;
            i__2 = (integer)workq[0].r; // , expr subst
            lwo = fla_max(i__1, i__2);
            aocl_lapack_cgelq(m, n, &a[a_offset], lda, tq, &c_n2, workq, &c_n2, &info2);
            tszm = (integer)tq[0].r;
            lwm = (integer)workq[0].r;
            aocl_lapack_cgemlq("L", trans, n, nrhs, m, &a[a_offset], lda, tq, &tszm, &b[b_offset],
                               ldb, workq, &c_n1, &info2);
            /* Computing MAX */
            i__1 = lwm;
            i__2 = (integer)workq[0].r; // , expr subst
            lwm = fla_max(i__1, i__2);
            wsizeo = tszo + lwo;
            wsizem = tszm + lwm;
        }
        if(*lwork < wsizem && !lquery)
        {
            *info = -10;
        }
        r__1 = aocl_lapack_sroundup_lwork(&wsizeo);
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGETSLS", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(lquery)
    {
        if(*lwork == -2)
        {
            r__1 = aocl_lapack_sroundup_lwork(&wsizem);
            work[1].r = r__1;
            work[1].i = 0.f; // , expr subst
        }
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*lwork < wsizeo)
    {
        lw1 = tszm;
        lw2 = lwm;
    }
    else
    {
        lw1 = tszo;
        lw2 = lwo;
    }
    /* Quick return if possible */
    /* Computing MIN */
    i__1 = fla_min(*m, *n);
    if(fla_min(i__1, *nrhs) == 0)
    {
        i__1 = fla_max(*m, *n);
        aocl_lapack_claset("FULL", &i__1, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Get machine parameters */
    smlnum = slamch_("S") / slamch_("P");
    bignum = 1.f / smlnum;
    /* Scale A, B if max element outside range [SMLNUM,BIGNUM] */
    anrm = aocl_lapack_clange("M", m, n, &a[a_offset], lda, dum);
    iascl = 0;
    if(anrm > 0.f && anrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        aocl_lapack_clascl("G", &c__0, &c__0, &anrm, &smlnum, m, n, &a[a_offset], lda, info);
        iascl = 1;
    }
    else if(anrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        aocl_lapack_clascl("G", &c__0, &c__0, &anrm, &bignum, m, n, &a[a_offset], lda, info);
        iascl = 2;
    }
    else if(anrm == 0.f)
    {
        /* Matrix all zero. Return zero solution. */
        aocl_lapack_claset("F", &maxmn, nrhs, &c_b1, &c_b1, &b[b_offset], ldb);
        goto L50;
    }
    brow = *m;
    if(tran)
    {
        brow = *n;
    }
    bnrm = aocl_lapack_clange("M", &brow, nrhs, &b[b_offset], ldb, dum);
    ibscl = 0;
    if(bnrm > 0.f && bnrm < smlnum)
    {
        /* Scale matrix norm up to SMLNUM */
        aocl_lapack_clascl("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 1;
    }
    else if(bnrm > bignum)
    {
        /* Scale matrix norm down to BIGNUM */
        aocl_lapack_clascl("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &b[b_offset], ldb, info);
        ibscl = 2;
    }
    if(*m >= *n)
    {
        /* compute QR factorization of A */
        aocl_lapack_cgeqr(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, info);
        if(!tran)
        {
            /* Least-Squares Problem min || A * X - B || */
            /* B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS) */
            aocl_lapack_cgemqr("L", "C", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &lw1,
                               &b[b_offset], ldb, &work[1], &lw2, info);
            /* B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */
            aocl_lapack_ctrtrs("U", "N", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
            scllen = *n;
        }
        else
        {
            /* Overdetermined system of equations A**T * X = B */
            /* B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS) */
            aocl_lapack_ctrtrs("U", "C", "N", n, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
            /* B(N+1:M,1:NRHS) = CZERO */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *m;
                for(i__ = *n + 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    b[i__3].r = 0.f;
                    b[i__3].i = 0.f; // , expr subst
                    /* L10: */
                }
                /* L20: */
            }
            /* B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */
            aocl_lapack_cgemqr("L", "N", m, nrhs, n, &a[a_offset], lda, &work[lw2 + 1], &lw1,
                               &b[b_offset], ldb, &work[1], &lw2, info);
            scllen = *m;
        }
    }
    else
    {
        /* Compute LQ factorization of A */
        aocl_lapack_cgelq(m, n, &a[a_offset], lda, &work[lw2 + 1], &lw1, &work[1], &lw2, info);
        /* workspace at least M, optimally M*NB. */
        if(!tran)
        {
            /* underdetermined system of equations A * X = B */
            /* B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */
            aocl_lapack_ctrtrs("L", "N", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
            /* B(M+1:N,1:NRHS) = 0 */
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *n;
                for(i__ = *m + 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    b[i__3].r = 0.f;
                    b[i__3].i = 0.f; // , expr subst
                    /* L30: */
                }
                /* L40: */
            }
            /* B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS) */
            aocl_lapack_cgemlq("L", "C", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &lw1,
                               &b[b_offset], ldb, &work[1], &lw2, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            scllen = *n;
        }
        else
        {
            /* overdetermined system min || A**T * X - B || */
            /* B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */
            aocl_lapack_cgemlq("L", "N", n, nrhs, m, &a[a_offset], lda, &work[lw2 + 1], &lw1,
                               &b[b_offset], ldb, &work[1], &lw2, info);
            /* workspace at least NRHS, optimally NRHS*NB */
            /* B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS) */
            aocl_lapack_ctrtrs("L", "C", "N", m, nrhs, &a[a_offset], lda, &b[b_offset], ldb, info);
            if(*info > 0)
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
            scllen = *m;
        }
    }
    /* Undo scaling */
    if(iascl == 1)
    {
        aocl_lapack_clascl("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
    else if(iascl == 2)
    {
        aocl_lapack_clascl("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
    if(ibscl == 1)
    {
        aocl_lapack_clascl("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
    else if(ibscl == 2)
    {
        aocl_lapack_clascl("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &b[b_offset], ldb,
                           info);
    }
L50:
    i__1 = tszo + lwo;
    r__1 = aocl_lapack_sroundup_lwork(&i__1);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGETSLS */
}
/* cgetsls_ */
