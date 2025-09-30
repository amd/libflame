/* ./sgemlq.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SGEMLQ */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T, */
/* $ TSIZE, C, LDC, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), T( * ), C(LDC, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEMLQ overwrites the general real M-by-N matrix C with */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'T': Q**T * C C * Q**T */
/* > where Q is a real orthogonal matrix defined as the product */
/* > of blocked elementary reflectors computed by short wide LQ */
/* > factorization (SGELQ) */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] SIDE */
/* > \verbatim */
/* > SIDE is CHARACTER*1 */
/* > = 'L': apply Q or Q**T from the Left;
 */
/* > = 'R': apply Q or Q**T from the Right. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > = 'N': No transpose, apply Q;
 */
/* > = 'T': Transpose, apply Q**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix A. M >=0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix C. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines */
/* > the matrix Q. */
/* > If SIDE = 'L', M >= K >= 0;
 */
/* > if SIDE = 'R', N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension */
/* > (LDA,M) if SIDE = 'L', */
/* > (LDA,N) if SIDE = 'R' */
/* > Part of the data structure to represent Q as returned by DGELQ. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,K). */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* > T is REAL array, dimension (MAX(5,TSIZE)). */
/* > Part of the data structure to represent Q as returned by SGELQ. */
/* > \endverbatim */
/* > */
/* > \param[in] TSIZE */
/* > \verbatim */
/* > TSIZE is INTEGER */
/* > The dimension of the array T. TSIZE >= 5. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL array, dimension (LDC,N) */
/* > On entry, the M-by-N matrix C. */
/* > On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* > LDC is INTEGER */
/* > The leading dimension of the array C. LDC >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > (workspace) REAL array, dimension (MAX(1,LWORK)) */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. */
/* > If LWORK = -1, then a workspace query is assumed. The routine */
/* > only calculates the size of the WORK array, returns this */
/* > value as WORK(1), and no error message related to WORK */
/* > is issued by XERBLA. */
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
/* > \par Further Details */
/* ==================== */
/* > */
/* > \verbatim */
/* > */
/* > These details are particular for this LAPACK implementation. Users should not */
/* > take them for granted. These details may change in the future, and are not likely */
/* > true for another LAPACK implementation. These details are relevant if one wants */
/* > to try to understand the code. They are not part of the interface. */
/* > */
/* > In this version, */
/* > */
/* > T(2): row block size (MB) */
/* > T(3): column block size (NB) */
/* > T(6:TSIZE): data structure needed for Q, computed by */
/* > SLASWLQ or SGELQT */
/* > */
/* > Depending on the matrix dimensions M and N, and row and column */
/* > block sizes MB and NB returned by ILAENV, SGELQ will use either */
/* > SLASWLQ (if the matrix is wide-and-short) or SGELQT to compute */
/* > the LQ factorization. */
/* > This version of SGEMLQ will use either SLAMSWLQ or SGEMLQT to */
/* > multiply matrix Q by another matrix. */
/* > Further Details in SLAMSWLQ or SGEMLQT. */
/* > \endverbatim */
/* > */
/* > \ingroup gemlq */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void sgemlq_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, real *a,
             aocl_int_t *lda, real *t, aocl_int_t *tsize, real *c__, aocl_int_t *ldc, real *work,
             aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgemlq(side, trans, m, n, k, a, lda, t, tsize, c__, ldc, work, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t tsize_64 = *tsize;
    aocl_int64_t ldc_64 = *ldc;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgemlq(side, trans, &m_64, &n_64, &k_64, a, &lda_64, t, &tsize_64, c__, &ldc_64,
                       work, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_sgemlq(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k,
                        real *a, aocl_int64_t *lda, real *t, aocl_int64_t *tsize, real *c__,
                        aocl_int64_t *ldc, real *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgemlq inputs: side %c ,trans %c ,m %" FLA_IS ",n %" FLA_IS ",k %" FLA_IS
                      ",lda %" FLA_IS ",tsize %" FLA_IS ",ldc %" FLA_IS ",lwork %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *lda, *tsize, *ldc, *lwork);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, c_dim1, c_offset, i__1;
    /* Local variables */
    aocl_int64_t mb, nb, mn, lw;
    logical left, tran;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    logical right;
    logical notran, lquery;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --t;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    lquery = *lwork == -1;
    notran = lsame_(trans, "N", 1, 1);
    tran = lsame_(trans, "T", 1, 1);
    left = lsame_(side, "L", 1, 1);
    right = lsame_(side, "R", 1, 1);
    mb = (integer)t[2];
    nb = (integer)t[3];
    if(left)
    {
        lw = *n * mb;
        mn = *m;
    }
    else
    {
        lw = *m * mb;
        mn = *n;
    }
    *info = 0;
    if(!left && !right)
    {
        *info = -1;
    }
    else if(!tran && !notran)
    {
        *info = -2;
    }
    else if(*m < 0)
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*k < 0 || *k > mn)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, *k))
    {
        *info = -7;
    }
    else if(*tsize < 5)
    {
        *info = -9;
    }
    else if(*ldc < fla_max(1, *m))
    {
        *info = -11;
    }
    else if(*lwork < fla_max(1, lw) && !lquery)
    {
        *info = -13;
    }
    if(*info == 0)
    {
        work[1] = aocl_lapack_sroundup_lwork(&lw);
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SGEMLQ", &i__1, (ftnlen)6);
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
    if(fla_min(i__1, *k) == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Computing MAX */
    i__1 = fla_max(*m, *n);
    if(left && *m <= *k || right && *n <= *k || nb <= *k || nb >= fla_max(i__1, *k))
    {
        aocl_lapack_sgemlqt(side, trans, m, n, k, &mb, &a[a_offset], lda, &t[6], &mb,
                            &c__[c_offset], ldc, &work[1], info);
    }
    else
    {
        aocl_lapack_slamswlq(side, trans, m, n, k, &mb, &nb, &a[a_offset], lda, &t[6], &mb,
                             &c__[c_offset], ldc, &work[1], lwork, info);
    }
    work[1] = aocl_lapack_sroundup_lwork(&lw);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGEMLQ */
}
/* sgemlq_ */
