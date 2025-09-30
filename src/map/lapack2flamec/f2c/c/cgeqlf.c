/* ./cgeqlf.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;
static aocl_int64_t c__3 = 3;
static aocl_int64_t c__2 = 2;
/* > \brief \b CGEQLF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGEQLF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeqlf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeqlf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeqlf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, LWORK, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEQLF computes a QL factorization of a scomplex M-by-N matrix A: */
/* > A = Q * L. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the M-by-N matrix A. */
/* > On exit, */
/* > if m >= n, the lower triangle of the subarray */
/* > A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L;
 */
/* > if m <= n, the elements on and below the (n-m)-th */
/* > superdiagonal contain the M-by-N lower trapezoidal matrix L;
 */
/* > the remaining elements, with the array TAU, represent the */
/* > unitary matrix Q as a product of elementary reflectors */
/* > (see Further Details). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (fla_min(M,N)) */
/* > The scalar factors of the elementary reflectors (see Further */
/* > Details). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* > On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK is INTEGER */
/* > The dimension of the array WORK. LWORK >= fla_max(1,N). */
/* > For optimum performance LWORK >= N*NB, where NB is */
/* > the optimal blocksize. */
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
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup geqlf */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > The matrix Q is represented as a product of elementary reflectors */
/* > */
/* > Q = H(k) . . . H(2) H(1), where k = fla_min(m,n). */
/* > */
/* > Each H(i) has the form */
/* > */
/* > H(i) = I - tau * v * v**H */
/* > */
/* > where tau is a scomplex scalar, and v is a scomplex vector with */
/* > v(m-k+i+1:m) = 0 and v(m-k+i) = 1;
v(1:m-k+i-1) is stored on exit in */
/* > A(1:m-k+i-1,n-k+i), and tau in TAU(i). */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgeqlf_(aocl_int_t *m, aocl_int_t *n, scomplex *a, aocl_int_t *lda, scomplex *tau, scomplex *work,
             aocl_int_t *lwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgeqlf(m, n, a, lda, tau, work, lwork, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t lwork_64 = *lwork;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgeqlf(&m_64, &n_64, a, &lda_64, tau, work, &lwork_64, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgeqlf(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                        scomplex *tau, scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cgeqlf inputs: m %lld, n %lld, lda %lld, lwork %lld", *m, *n, *lda,
             *lwork);
#else
    snprintf(buffer, 256, "cgeqlf inputs: m %d, n %d, lda %d, lwork %d", *m, *n, *lda, *lwork);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1;
    /* Local variables */
    aocl_int64_t i__, k, ib, nb, ki, kk, mu, nu, nx, iws, nbmin, iinfo;
    aocl_int64_t ldwork, lwkopt;
    logical lquery;
    /* -- LAPACK computational routine -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    nb = aocl_lapack_ilaenv(&c__1, "CGEQLF", " ", m, n, &c_n1, &c_n1);
    if(*info == 0)
    {
        k = fla_min(*m, *n);
        if(k == 0)
        {
            lwkopt = 1;
        }
        else
        {
            lwkopt = *n * nb;
        }
        r__1 = aocl_lapack_sroundup_lwork(&lwkopt);
        work[1].r = r__1;
        work[1].i = 0.f; // , expr subst
        if(*lwork < fla_max(1, *n) && !lquery)
        {
            *info = -7;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGEQLF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    else if(lquery)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(k == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    nbmin = 2;
    nx = 1;
    iws = *n;
    if(nb > 1 && nb < k)
    {
        /* Determine when to cross over from blocked to unblocked code. */
        /* Computing MAX */
        i__1 = 0;
        i__2 = aocl_lapack_ilaenv(&c__3, "CGEQLF", " ", m, n, &c_n1, &c_n1); // , expr subst
        nx = fla_max(i__1, i__2);
        if(nx < k)
        {
            /* Determine if workspace is large enough for blocked code. */
            ldwork = *n;
            iws = ldwork * nb;
            if(*lwork < iws)
            {
                /* Not enough workspace to use optimal NB: reduce NB and */
                /* determine the minimum value of NB. */
                nb = *lwork / ldwork;
                /* Computing MAX */
                i__1 = 2;
                i__2 = aocl_lapack_ilaenv(&c__2, "CGEQLF", " ", m, n, &c_n1, &c_n1); // , expr subst
                nbmin = fla_max(i__1, i__2);
            }
        }
    }
    if(nb >= nbmin && nb < k && nx < k)
    {
        /* Use blocked code initially. */
        /* The last kk columns are handled by the block method. */
        ki = (k - nx - 1) / nb * nb;
        /* Computing MIN */
        i__1 = k;
        i__2 = ki + nb; // , expr subst
        kk = fla_min(i__1, i__2);
        i__1 = k - kk + 1;
        i__2 = -nb;
        for(i__ = k - kk + ki + 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
        {
            /* Computing MIN */
            i__3 = k - i__ + 1;
            ib = fla_min(i__3, nb);
            /* Compute the QL factorization of the current block */
            /* A(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1) */
            i__3 = *m - k + i__ + ib - 1;
            aocl_lapack_cgeql2(&i__3, &ib, &a[(*n - k + i__) * a_dim1 + 1], lda, &tau[i__],
                               &work[1], &iinfo);
            if(*n - k + i__ > 1)
            {
                /* Form the triangular factor of the block reflector */
                /* H = H(i+ib-1) . . . H(i+1) H(i) */
                i__3 = *m - k + i__ + ib - 1;
                aocl_lapack_clarft("Backward", "Columnwise", &i__3, &ib,
                                   &a[(*n - k + i__) * a_dim1 + 1], lda, &tau[i__], &work[1],
                                   &ldwork);
                /* Apply H**H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left */
                i__3 = *m - k + i__ + ib - 1;
                i__4 = *n - k + i__ - 1;
                aocl_lapack_clarfb("Left", "Conjugate transpose", "Backward",
                                   "Columnwi"
                                   "se",
                                   &i__3, &i__4, &ib, &a[(*n - k + i__) * a_dim1 + 1], lda,
                                   &work[1], &ldwork, &a[a_offset], lda, &work[ib + 1], &ldwork);
            }
            /* L10: */
        }
        mu = *m - k + i__ + nb - 1;
        nu = *n - k + i__ + nb - 1;
    }
    else
    {
        mu = *m;
        nu = *n;
    }
    /* Use unblocked code to factor the last or only block */
    if(mu > 0 && nu > 0)
    {
        aocl_lapack_cgeql2(&mu, &nu, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
    }
    r__1 = aocl_lapack_sroundup_lwork(&iws);
    work[1].r = r__1;
    work[1].i = 0.f; // , expr subst
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGEQLF */
}
/* cgeqlf_ */
