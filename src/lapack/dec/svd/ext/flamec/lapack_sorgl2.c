/* ../netlib/sorgl2.f -- translated by f2c (version 20000121). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h" /* > \brief \b SORGL2 */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SORGL2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgl2. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgl2. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgl2. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SORGL2( M, N, K, A, LDA, TAU, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, K, LDA, M, N */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), TAU( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORGL2 generates an m by n real matrix Q with orthonormal rows, */
/* > which is defined as the first m rows of a product of k elementary */
/* > reflectors of order n */
/* > */
/* > Q = H(k) . . . H(2) H(1) */
/* > */
/* > as returned by SGELQF. */
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
/* > The number of columns of the matrix Q. N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* > K is INTEGER */
/* > The number of elementary reflectors whose product defines the */
/* > matrix Q. M >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > On entry, the i-th row must contain the vector which defines */
/* > the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* > by SGELQF in the first k rows of its array argument A. */
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
/* > TAU is REAL array, dimension (K) */
/* > TAU(i) must contain the scalar factor of the elementary */
/* > reflector H(i), as returned by SGELQF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (M) */
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
/* > \ingroup realOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
int lapack_sorgl2(integer *m, integer *n, integer *k, real *a, integer *lda, real *tau, real *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;
    /* Local variables */
    integer i__, j, l;
    extern /* Subroutine */
    void sscal_(integer *, real *, real *, integer *), slarf_(char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *), xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* -- LAPACK computational routine -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input arguments */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < *m)
    {
        *info = -2;
    }
    else if (*k < 0 || *k > *m)
    {
        *info = -3;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -5;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SORGL2", &i__1, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible */
    if (*m <= 0)
    {
        return 0;
    }
    if (*k < *m)
    {
        /* Initialise rows k+1:m to rows of the unit matrix */
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (l = *k + 1;
                    l <= i__2;
                    ++l)
            {
                a[l + j * a_dim1] = 0.f;
                /* L10: */
            }
            if (j > *k && j <= *m)
            {
                a[j + j * a_dim1] = 1.f;
            }
            /* L20: */
        }
    }
    for (i__ = *k;
            i__ >= 1;
            --i__)
    {
        /* Apply H(i) to A(i:m,i:n) from the right */
        if (i__ < *n)
        {
            if (i__ < *m)
            {
                a[i__ + i__ * a_dim1] = 1.f;
                i__1 = *m - i__;
                i__2 = *n - i__ + 1;
                slarf_("Right", &i__1, &i__2, &a[i__ + i__ * a_dim1], lda, & tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1]);
            }
            i__1 = *n - i__;
            r__1 = -tau[i__];
            sscal_(&i__1, &r__1, &a[i__ + (i__ + 1) * a_dim1], lda);
        }
        a[i__ + i__ * a_dim1] = 1.f - tau[i__];
        /* Set A(i,1:i-1) to zero */
        i__1 = i__ - 1;
        for (l = 1;
                l <= i__1;
                ++l)
        {
            a[i__ + l * a_dim1] = 0.f;
            /* L30: */
        }
        /* L40: */
    }
    return 0;
    /* End of SORGL2 */
}
/* lapack_sorgl2 */
