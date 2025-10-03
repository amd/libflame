/* ../netlib/cung2r.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CUNG2R */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CUNG2R + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cung2r.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cung2r.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cung2r.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CUNG2R( M, N, K, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CUNG2R generates an m by n scomplex matrix Q with orthonormal columns, */
/* > which is defined as the first n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* > Q = H(1) H(2) . . . H(k) */
/* > */
/* > as returned by CGEQRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] M */
/* > \verbatim */
/* > M is INTEGER */
/* > The number of rows of the matrix Q. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of columns of the matrix Q. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines the */
/* > matrix Q. N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the i-th column must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by CGEQRF in the first k columns of its array */
/* > argument A. */
/* > On exit, the m by n matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The first dimension of the array A. LDA >= fla_max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* > TAU is COMPLEX array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by CGEQRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cung2r_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, scomplex *a, aocl_int_t *lda, scomplex *tau,
             scomplex *work, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cung2r(m, n, k, a, lda, tau, work, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cung2r(&m_64, &n_64, &k_64, a, &lda_64, tau, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cung2r(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, scomplex *a,
                        aocl_int64_t *lda, scomplex *tau, scomplex *work, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cung2r inputs: m %lld, n %lld, k %lld, lda %lld", *m, *n, *k, *lda);
#else
    snprintf(buffer, 256, "cung2r inputs: m %d, n %d, k %d, lda %d", *m, *n, *k, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    scomplex q__1;
    /* Local variables */
    aocl_int64_t i__, j, l;
    /* -- LAPACK computational routine (version 3.4.0) -- */
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
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0 || *n > *m)
    {
        *info = -2;
    }
    else if(*k < 0 || *k > *n)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CUNG2R", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Initialise columns k+1:n to columns of the unit matrix */
    i__1 = *n;
    for(j = *k + 1; j <= i__1; ++j)
    {
        i__2 = *m;
        for(l = 1; l <= i__2; ++l)
        {
            i__3 = l + j * a_dim1;
            a[i__3].real = 0.f;
            a[i__3].imag = 0.f; // , expr subst
            /* L10: */
        }
        i__2 = j + j * a_dim1;
        a[i__2].real = 1.f;
        a[i__2].imag = 0.f; // , expr subst
        /* L20: */
    }
    for(i__ = *k; i__ >= 1; --i__)
    {
        /* Apply H(i) to A(i:m,i:n) from the left */
        if(i__ < *n)
        {
            i__1 = i__ + i__ * a_dim1;
            a[i__1].real = 1.f;
            a[i__1].imag = 0.f; // , expr subst
            i__1 = *m - i__ + 1;
            i__2 = *n - i__;
            aocl_lapack_clarf("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[i__],
                              &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
        }
        if(i__ < *m)
        {
            i__1 = *m - i__;
            i__2 = i__;
            q__1.real = -tau[i__2].real;
            q__1.imag = -tau[i__2].imag; // , expr subst
            aocl_blas_cscal(&i__1, &q__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
        }
        i__1 = i__ + i__ * a_dim1;
        i__2 = i__;
        q__1.real = 1.f - tau[i__2].real;
        q__1.imag = 0.f - tau[i__2].imag; // , expr subst
        a[i__1].real = q__1.real;
        a[i__1].imag = q__1.imag; // , expr subst
        /* Set A(1:i-1,i) to zero */
        i__1 = i__ - 1;
        for(l = 1; l <= i__1; ++l)
        {
            i__2 = l + i__ * a_dim1;
            a[i__2].real = 0.f;
            a[i__2].imag = 0.f; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CUNG2R */
}
/* cung2r_ */
