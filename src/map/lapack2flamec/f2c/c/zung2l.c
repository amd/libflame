/* ../netlib/zung2l.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b ZUNG2L generates all or part of the unitary matrix Q from a QL factorization
 * determined by cgeq lf (unblocked algorithm). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZUNG2L + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zung2l.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zung2l.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zung2l.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNG2L generates an m by n scomplex matrix Q with orthonormal columns, */
/* > which is defined as the last n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* > Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by ZGEQLF. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the (n-k+i)-th column must contain the vector which */
/* > defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* > returned by ZGEQLF in the last k columns of its array */
/* > argument A. */
/* > On exit, the m-by-n matrix Q. */
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
/* > TAU is COMPLEX*16 array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by ZGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N) */
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
/* > \date September 2012 */
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zung2l_(aocl_int_t *m, aocl_int_t *n, aocl_int_t *k, dcomplex *a, aocl_int_t *lda,
             dcomplex *tau, dcomplex *work, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zung2l(m, n, k, a, lda, tau, work, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t k_64 = *k;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zung2l(&m_64, &n_64, &k_64, a, &lda_64, tau, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zung2l(aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a,
                        aocl_int64_t *lda, dcomplex *tau, dcomplex *work,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zung2l inputs: m %" FLA_IS ", n %" FLA_IS ", k %" FLA_IS ", lda %" FLA_IS "",
                      *m, *n, *k, *lda);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    dcomplex z__1;
    /* Local variables */
    aocl_int64_t i__, j, l, ii;
    /* -- LAPACK computational routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
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
        aocl_blas_xerbla("ZUNG2L", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Initialise columns 1:n-k to columns of the unit matrix */
    i__1 = *n - *k;
    for(j = 1; j <= i__1; ++j)
    {
        i__2 = *m;
        for(l = 1; l <= i__2; ++l)
        {
            i__3 = l + j * a_dim1;
            a[i__3].real = 0.;
            a[i__3].imag = 0.; // , expr subst
            /* L10: */
        }
        i__2 = *m - *n + j + j * a_dim1;
        a[i__2].real = 1.;
        a[i__2].imag = 0.; // , expr subst
        /* L20: */
    }
    i__1 = *k;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        ii = *n - *k + i__;
        /* Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */
        i__2 = *m - *n + ii + ii * a_dim1;
        a[i__2].real = 1.;
        a[i__2].imag = 0.; // , expr subst
        i__2 = *m - *n + ii;
        i__3 = ii - 1;
        aocl_lapack_zlarf("Left", &i__2, &i__3, &a[ii * a_dim1 + 1], &c__1, &tau[i__], &a[a_offset],
                          lda, &work[1]);
        i__2 = *m - *n + ii - 1;
        i__3 = i__;
        z__1.real = -tau[i__3].real;
        z__1.imag = -tau[i__3].imag; // , expr subst
        aocl_blas_zscal(&i__2, &z__1, &a[ii * a_dim1 + 1], &c__1);
        i__2 = *m - *n + ii + ii * a_dim1;
        i__3 = i__;
        z__1.real = 1. - tau[i__3].real;
        z__1.imag = 0. - tau[i__3].imag; // , expr subst
        a[i__2].real = z__1.real;
        a[i__2].imag = z__1.imag; // , expr subst
        /* Set A(m-k+i+1:m,n-k+i) to zero */
        i__2 = *m;
        for(l = *m - *n + ii + 1; l <= i__2; ++l)
        {
            i__3 = l + ii * a_dim1;
            a[i__3].real = 0.;
            a[i__3].imag = 0.; // , expr subst
            /* L30: */
        }
        /* L40: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZUNG2L */
}
/* zung2l_ */
