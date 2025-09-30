/* slamswlq.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__0 = 0;
/* > \brief \b SLAMSWLQ */
/* Definition: */
/* =========== */
/* SUBROUTINE SLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, */
/* $ LDT, C, LDC, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER SIDE, TRANS */
/* INTEGER INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE A( LDA, * ), WORK( * ), C(LDC, * ), */
/* $ T( LDT, * ) */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAMSWLQ overwrites the general real M-by-N matrix C with */
/* > */
/* > */
/* > SIDE = 'L' SIDE = 'R' */
/* > TRANS = 'N': Q * C C * Q */
/* > TRANS = 'T': Q**T * C C * Q**T */
/* > where Q is a real orthogonal matrix defined as the product of blocked */
/* > elementary reflectors computed by short wide LQ */
/* > factorization (SLASWLQ) */
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
/* > The number of rows of the matrix C. M >=0. */
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
/* > M >= K >= 0;
 */
/* > */
/* > \endverbatim */
/* > \param[in] MB */
/* > \verbatim */
/* > MB is INTEGER */
/* > The row block size to be used in the blocked LQ. */
/* > M >= MB >= 1 */
/* > \endverbatim */
/* > */
/* > \param[in] NB */
/* > \verbatim */
/* > NB is INTEGER */
/* > The column block size to be used in the blocked LQ. */
/* > NB > M. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension */
/* > (LDA,M) if SIDE = 'L', */
/* > (LDA,N) if SIDE = 'R' */
/* > The i-th row must contain the vector which defines the blocked */
/* > elementary reflector H(i), for i = 1,2,...,k, as returned by */
/* > SLASWLQ in the first k rows of its array argument A. */
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
/* > T is REAL array, dimension */
/* > ( M * Number of blocks(CEIL(N-K/NB-K)), */
/* > The blocked upper triangular block reflectors stored in compact form */
/* > as a sequence of upper triangular blocks. See below */
/* > for further details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* > LDT is INTEGER */
/* > The leading dimension of the array T. LDT >= MB. */
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
/* > If SIDE = 'L', LWORK >= fla_max(1,NB) * MB;
 */
/* > if SIDE = 'R', LWORK >= fla_max(1,M) * MB. */
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
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > Short-Wide LQ (SWLQ) performs LQ by a sequence of orthogonal transformations, */
/* > representing Q as a product of other orthogonal matrices */
/* > Q = Q(1) * Q(2) * . . . * Q(k) */
/* > where each Q(i) zeros out upper diagonal entries of a block of NB rows of A: */
/* > Q(1) zeros out the upper diagonal entries of rows 1:NB of A */
/* > Q(2) zeros out the bottom MB-N rows of rows [1:M,NB+1:2*NB-M] of A */
/* > Q(3) zeros out the bottom MB-N rows of rows [1:M,2*NB-M+1:3*NB-2*M] of A */
/* > . . . */
/* > */
/* > Q(1) is computed by GELQT, which represents Q(1) by Householder vectors */
/* > stored under the diagonal of rows 1:MB of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,1:N). */
/* > For more information see Further Details in GELQT. */
/* > */
/* > Q(i) for i>1 is computed by TPLQT, which represents Q(i) by Householder vectors */
/* > stored in columns [(i-1)*(NB-M)+M+1:i*(NB-M)+M] of A, and by upper triangular */
/* > block reflectors, stored in array T(1:LDT,(i-1)*M+1:i*M). */
/* > The last Q(k) may use fewer rows. */
/* > For more information see Further Details in TPLQT. */
/* > */
/* > For more details of the overall algorithm, see the description of */
/* > Sequential TSQR in Section 2.2 of [1]. */
/* > */
/* > [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations, */
/* > J. Demmel, L. Grigori, M. Hoemmen, J. Langou, */
/* > SIAM J. Sci. Comput, vol. 34, no. 1, 2012 */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void slamswlq_(char *side, char *trans, aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, aocl_int_t *mb,
               aocl_int_t *nb, real *a, aocl_int_t *lda, real *t, aocl_int_t *ldt, real *c__,
               aocl_int_t *ldc, real *work, aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_slamswlq(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c__, ldc, work, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t mb_64 = *mb;
    aocl_int64_t nb_64 = *nb;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldt_64 = *ldt;
    aocl_int64_t ldc_64 = *ldc;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_slamswlq(side, trans, &m_64, &n_64, &k_64, &mb_64, &nb_64, a, &lda_64, t, &ldt_64,
                         c__, &ldc_64, work, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_slamswlq(char *side, char *trans, aocl_int64_t *m, aocl_int64_t *n,
                          aocl_int64_t *k, aocl_int64_t *mb, aocl_int64_t *nb, real *a,
                          aocl_int64_t *lda, real *t, aocl_int64_t *ldt, real *c__,
                          aocl_int64_t *ldc, real *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slamswlq inputs: side %c, trans %c, m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS
                      ", mb %" FLA_IS ", nb %" FLA_IS ", lda %" FLA_IS ", ldt %" FLA_IS
                      ", ldc %" FLA_IS ", lwork %" FLA_IS "",
                      *side, *trans, *m, *n, *k, *mb, *nb, *lda, *ldt, *ldc, *lwork);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, c_dim1, c_offset, t_dim1, t_offset, i__1, i__2, i__3;
    /* Local variables */
    aocl_int64_t i__, ii, kk, lw, ctr;
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    lquery = *lwork < 0;
    notran = lsame_(trans, "N", 1, 1);
    tran = lsame_(trans, "T", 1, 1);
    left = lsame_(side, "L", 1, 1);
    right = lsame_(side, "R", 1, 1);
    if(left)
    {
        lw = *n * *mb;
    }
    else
    {
        lw = *m * *mb;
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
    else if(*k < 0)
    {
        *info = -5;
    }
    else if(*m < *k)
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*k < *mb || *mb < 1)
    {
        *info = -6;
    }
    else if(*lda < fla_max(1, *k))
    {
        *info = -9;
    }
    else if(*ldt < fla_max(1, *mb))
    {
        *info = -11;
    }
    else if(*ldc < fla_max(1, *m))
    {
        *info = -13;
    }
    else if(*lwork < fla_max(1, lw) && !lquery)
    {
        *info = -15;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SLAMSWLQ", &i__1, (ftnlen)8);
        work[1] = (real)lw;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        work[1] = (real)lw;
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
    if(*nb <= *k || *nb >= fla_max(i__1, *k))
    {
        aocl_lapack_sgemlqt(side, trans, m, n, k, mb, &a[a_offset], lda, &t[t_offset], ldt,
                            &c__[c_offset], ldc, &work[1], info);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(left && tran)
    {
        /* Multiply Q to the last block of C */
        kk = (*m - *k) % (*nb - *k);
        ctr = (*m - *k) / (*nb - *k);
        if(kk > 0)
        {
            ii = *m - kk + 1;
            aocl_lapack_stpmlqt("L", "T", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[ii + c_dim1], ldc, &work[1], info);
        }
        else
        {
            ii = *m + 1;
        }
        i__1 = *nb + 1;
        i__2 = -(*nb - *k);
        for(i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            /* Multiply Q to the current block of C (1:M,I:I+NB) */
            --ctr;
            i__3 = *nb - *k;
            aocl_lapack_stpmlqt("L", "T", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[i__ + c_dim1], ldc, &work[1], info);
        }
        /* Multiply Q to the first block of C (1:M,1:NB) */
        aocl_lapack_sgemlqt("L", "T", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt,
                            &c__[c_dim1 + 1], ldc, &work[1], info);
    }
    else if(left && notran)
    {
        /* Multiply Q to the first block of C */
        kk = (*m - *k) % (*nb - *k);
        ii = *m - kk + 1;
        ctr = 1;
        aocl_lapack_sgemlqt("L", "N", nb, n, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt,
                            &c__[c_dim1 + 1], ldc, &work[1], info);
        i__2 = ii - *nb + *k;
        i__1 = *nb - *k;
        for(i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
        {
            /* Multiply Q to the current block of C (I:I+NB,1:N) */
            i__3 = *nb - *k;
            aocl_lapack_stpmlqt("L", "N", &i__3, n, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[i__ + c_dim1], ldc, &work[1], info);
            ++ctr;
        }
        if(ii <= *m)
        {
            /* Multiply Q to the last block of C */
            aocl_lapack_stpmlqt("L", "N", &kk, n, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[ii + c_dim1], ldc, &work[1], info);
        }
    }
    else if(right && notran)
    {
        /* Multiply Q to the last block of C */
        kk = (*n - *k) % (*nb - *k);
        ctr = (*n - *k) / (*nb - *k);
        if(kk > 0)
        {
            ii = *n - kk + 1;
            aocl_lapack_stpmlqt("R", "N", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[ii * c_dim1 + 1], ldc, &work[1], info);
        }
        else
        {
            ii = *n + 1;
        }
        i__1 = *nb + 1;
        i__2 = -(*nb - *k);
        for(i__ = ii - (*nb - *k); i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            /* Multiply Q to the current block of C (1:M,I:I+MB) */
            --ctr;
            i__3 = *nb - *k;
            aocl_lapack_stpmlqt("R", "N", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[i__ * c_dim1 + 1], ldc, &work[1], info);
        }
        /* Multiply Q to the first block of C (1:M,1:MB) */
        aocl_lapack_sgemlqt("R", "N", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt,
                            &c__[c_dim1 + 1], ldc, &work[1], info);
    }
    else if(right && tran)
    {
        /* Multiply Q to the first block of C */
        kk = (*n - *k) % (*nb - *k);
        ii = *n - kk + 1;
        ctr = 1;
        aocl_lapack_sgemlqt("R", "T", m, nb, k, mb, &a[a_dim1 + 1], lda, &t[t_offset], ldt,
                            &c__[c_dim1 + 1], ldc, &work[1], info);
        i__2 = ii - *nb + *k;
        i__1 = *nb - *k;
        for(i__ = *nb + 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1)
        {
            /* Multiply Q to the current block of C (1:M,I:I+MB) */
            i__3 = *nb - *k;
            aocl_lapack_stpmlqt("R", "T", m, &i__3, k, &c__0, mb, &a[i__ * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[i__ * c_dim1 + 1], ldc, &work[1], info);
            ++ctr;
        }
        if(ii <= *n)
        {
            /* Multiply Q to the last block of C */
            aocl_lapack_stpmlqt("R", "T", m, &kk, k, &c__0, mb, &a[ii * a_dim1 + 1], lda,
                                &t[(ctr * *k + 1) * t_dim1 + 1], ldt, &c__[c_dim1 + 1], ldc,
                                &c__[ii * c_dim1 + 1], ldc, &work[1], info);
        }
    }
    work[1] = (real)lw;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLAMSWLQ */
}
/* slamswlq_ */
