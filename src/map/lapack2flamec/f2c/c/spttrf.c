/* ../netlib/spttrf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SPTTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPTTRF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spttrf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spttrf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spttrf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPTTRF( N, D, E, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* REAL D( * ), E( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPTTRF computes the L*D*L**T factorization of a real symmetric */
/* > positive definite tridiagonal matrix A. The factorization may also */
/* > be regarded as having the form A = U**T*D*U. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* > D is REAL array, dimension (N) */
/* > On entry, the n diagonal elements of the tridiagonal matrix */
/* > A. On exit, the n diagonal elements of the diagonal matrix */
/* > D from the L*D*L**T factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] E */
/* > \verbatim */
/* > E is REAL array, dimension (N-1) */
/* > On entry, the (n-1) subdiagonal elements of the tridiagonal */
/* > matrix A. On exit, the (n-1) subdiagonal elements of the */
/* > unit bidiagonal factor L from the L*D*L**T factorization of A. */
/* > E can also be regarded as the superdiagonal of the unit */
/* > bidiagonal factor U from the U**T*D*U factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -k, the k-th argument had an illegal value */
/* > > 0: if INFO = k, the leading minor of order k is not */
/* > positive definite;
if k < N, the factorization could not */
/* > be completed, while if k = N, the factorization was */
/* > completed, but D(N) <= 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void spttrf_(integer *n, real *d__, real *e, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spttrf inputs: n %" FLA_IS "", *n);
    /* System generated locals */
    integer i__1;
    /* Local variables */
    integer i__, i4;
    real ei;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    /* Test the input parameters. */
    /* Parameter adjustments */
    --e;
    --d__;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
        i__1 = -(*info);
        xerbla_("SPTTRF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Compute the L*D*L**T (or U**T*D*U) factorization of A. */
    i4 = (*n - 1) % 4;
    i__1 = i4;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        if(d__[i__] <= 0.f)
        {
            *info = i__;
            goto L30;
        }
        ei = e[i__];
        e[i__] = ei / d__[i__];
        d__[i__ + 1] -= e[i__] * ei;
        /* L10: */
    }
    i__1 = *n - 4;
    for(i__ = i4 + 1; i__ <= i__1; i__ += 4)
    {
        /* Drop out of the loop if d(i) <= 0: the matrix is not positive */
        /* definite. */
        if(d__[i__] <= 0.f)
        {
            *info = i__;
            goto L30;
        }
        /* Solve for e(i) and d(i+1). */
        ei = e[i__];
        e[i__] = ei / d__[i__];
        d__[i__ + 1] -= e[i__] * ei;
        if(d__[i__ + 1] <= 0.f)
        {
            *info = i__ + 1;
            goto L30;
        }
        /* Solve for e(i+1) and d(i+2). */
        ei = e[i__ + 1];
        e[i__ + 1] = ei / d__[i__ + 1];
        d__[i__ + 2] -= e[i__ + 1] * ei;
        if(d__[i__ + 2] <= 0.f)
        {
            *info = i__ + 2;
            goto L30;
        }
        /* Solve for e(i+2) and d(i+3). */
        ei = e[i__ + 2];
        e[i__ + 2] = ei / d__[i__ + 2];
        d__[i__ + 3] -= e[i__ + 2] * ei;
        if(d__[i__ + 3] <= 0.f)
        {
            *info = i__ + 3;
            goto L30;
        }
        /* Solve for e(i+3) and d(i+4). */
        ei = e[i__ + 3];
        e[i__ + 3] = ei / d__[i__ + 3];
        d__[i__ + 4] -= e[i__ + 3] * ei;
        /* L20: */
    }
    /* Check d(n) for positive definiteness. */
    if(d__[*n] <= 0.f)
    {
        *info = *n;
    }
L30:
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPTTRF */
}
/* spttrf_ */
