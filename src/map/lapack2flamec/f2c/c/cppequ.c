/* ../netlib/cppequ.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CPPEQU */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPPEQU + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppequ.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppequ.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppequ.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* REAL AMAX, SCOND */
/* .. */
/* .. Array Arguments .. */
/* REAL S( * ) */
/* COMPLEX AP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPPEQU computes row and column scalings intended to equilibrate a */
/* > Hermitian positive definite matrix A in packed storage and reduce */
/* > its condition number (with respect to the two-norm). S contains the */
/* > scale factors, S(i)=1/sqrt(A(i,i)), chosen so that the scaled matrix */
/* > B with elements B(i,j)=S(i)*A(i,j)*S(j) has ones on the diagonal. */
/* > This choice of S puts the condition number of B within a factor N of */
/* > the smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* > AP is COMPLEX array, dimension (N*(N+1)/2) */
/* > The upper or lower triangle of the Hermitian matrix A, packed */
/* > columnwise in a linear array. The j-th column of A is stored */
/* > in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
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
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void cppequ_(char *uplo, integer *n, complex *ap, real *s, real *scond, real *amax, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "cppequ inputs: uplo %c, n %lld", *uplo, *n);
#else
    snprintf(buffer, 256, "cppequ inputs: uplo %c, n %d", *uplo, *n);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, jj;
    real smin;
    extern logical lsame_(char *, char *, integer, integer);
    logical upper;
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
    /* Parameter adjustments */
    --s;
    --ap;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CPPEQU", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        *scond = 1.f;
        *amax = 0.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Initialize SMIN and AMAX. */
    s[1] = ap[1].r;
    smin = s[1];
    *amax = s[1];
    if(upper)
    {
        /* UPLO = 'U': Upper triangle of A is stored. */
        /* Find the minimum and maximum diagonal elements. */
        jj = 1;
        i__1 = *n;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            jj += i__;
            i__2 = jj;
            s[i__] = ap[i__2].r;
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
    }
    else
    {
        /* UPLO = 'L': Lower triangle of A is stored. */
        /* Find the minimum and maximum diagonal elements. */
        jj = 1;
        i__1 = *n;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            jj = jj + *n - i__ + 2;
            i__2 = jj;
            s[i__] = ap[i__2].r;
            /* Computing MIN */
            r__1 = smin;
            r__2 = s[i__]; // , expr subst
            smin = fla_min(r__1, r__2);
            /* Computing MAX */
            r__1 = *amax;
            r__2 = s[i__]; // , expr subst
            *amax = fla_max(r__1, r__2);
            /* L20: */
        }
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
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
            /* L30: */
        }
    }
    else
    {
        /* Set the scale factors to the reciprocals */
        /* of the diagonal elements. */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            s[i__] = 1.f / sqrt(s[i__]);
            /* L40: */
        }
        /* Compute SCOND = fla_min(S(I)) / fla_max(S(I)) */
        *scond = sqrt(smin) / sqrt(*amax);
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CPPEQU */
}
/* cppequ_ */
