/* ../netlib/spoequb.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SPOEQUB */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SPOEQUB + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spoequb
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spoequb
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spoequb
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, N */
/* REAL AMAX, SCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), S( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPOEQU computes row and column scalings intended to equilibrate a */
/* > symmetric positive definite matrix A and reduce its condition number */
/* > (with respect to the two-norm). S contains the scale factors, */
/* > S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/* > elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal. This */
/* > choice of S puts the condition number of B within a factor N of the */
/* > smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The N-by-N symmetric positive definite matrix whose scaling */
/* > factors are to be computed. Only the diagonal elements of A */
/* > are referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is REAL array, dimension (N) */
/* > If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* > SCOND is REAL */
/* > If INFO = 0, S contains the ratio of the smallest S(i) to */
/* > the largest S(i). If SCOND >= 0.1 and AMAX is neither too */
/* > large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* > AMAX is REAL */
/* > Absolute value of largest matrix element. If AMAX is very */
/* > close to overflow or very close to underflow, the matrix */
/* > should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, the i-th diagonal element is nonpositive. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realPOcomputational */
/* ===================================================================== */
/* Subroutine */
void spoequb_(integer *n, real *a, integer *lda, real *s, real *scond, real *amax, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spoequb inputs: n %" FLA_IS ",lda %" FLA_IS "", *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double log(doublereal), pow_ri(real *, integer *), sqrt(doublereal);
    /* Local variables */
    integer i__;
    real tmp, base, smin;
    extern real slamch_(char *);
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
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Positive definite only performs 1 pass of equilibration. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    /* Function Body */
    *info = 0;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -3;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SPOEQUB", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible. */
    if(*n == 0)
    {
        *scond = 1.f;
        *amax = 0.f;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    base = slamch_("B");
    tmp = -.5f / log(base);
    /* Find the minimum and maximum diagonal elements. */
    s[1] = a[a_dim1 + 1];
    smin = s[1];
    *amax = s[1];
    i__1 = *n;
    for(i__ = 2; i__ <= i__1; ++i__)
    {
        s[i__] = a[i__ + i__ * a_dim1];
        /* Computing MIN */
        r__1 = smin;
        r__2 = s[i__]; // , expr subst
        smin = fla_min(r__1, r__2);
        /* Computing MAX */
        r__1 = *amax;
        r__2 = s[i__]; // , expr subst
        *amax = fla_max(r__1, r__2);
        /* L10: */
    }
    if(smin <= 0.f)
    {
        /* Find the first non-positive diagonal element and return. */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if(s[i__] <= 0.f)
            {
                *info = i__;
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* L20: */
        }
    }
    else
    {
        /* Set the scale factors to the reciprocals */
        /* of the diagonal elements. */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            i__2 = (integer)(tmp * log(s[i__]));
            s[i__] = pow_ri(&base, &i__2);
            /* L30: */
        }
        /* Compute SCOND = fla_min(S(I)) / fla_max(S(I)). */
        *scond = sqrt(smin) / sqrt(*amax);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SPOEQUB */
}
/* spoequb_ */
