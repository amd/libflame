/* ../netlib/sgbequ.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SGBEQU */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGBEQU + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbequ.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbequ.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbequ.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, */
/* AMAX, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, KL, KU, LDAB, M, N */
/* REAL AMAX, COLCND, ROWCND */
/* .. */
/* .. Array Arguments .. */
/* REAL AB( LDAB, * ), C( * ), R( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGBEQU computes row and column scalings intended to equilibrate an */
/* > M-by-N band matrix A and reduce its condition number. R returns the */
/* > row scale factors and C the column scale factors, chosen to try to */
/* > make the largest element in each row and column of the matrix B with */
/* > elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1. */
/* > */
/* > R(i) and C(j) are restricted to be between SMLNUM = smallest safe */
/* > number and BIGNUM = largest safe number. Use of these scaling */
/* > factors is not guaranteed to reduce the condition number of A but */
/* > works well in practice. */
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
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is REAL array, dimension (LDAB,N) */
/* > The band matrix A, stored in rows 1 to KL+KU+1. The j-th */
/* > column of A is stored in the j-th column of the array AB as */
/* > follows: */
/* > AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=fla_min(m,j+kl). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* > R is REAL array, dimension (M) */
/* > If INFO = 0, or INFO > M, R contains the row scale factors */
/* > for A. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is REAL array, dimension (N) */
/* > If INFO = 0, C contains the column scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] ROWCND */
/* > \verbatim */
/* > ROWCND is REAL */
/* > If INFO = 0 or INFO > M, ROWCND contains the ratio of the */
/* > smallest R(i) to the largest R(i). If ROWCND >= 0.1 and */
/* > AMAX is neither too large nor too small, it is not worth */
/* > scaling by R. */
/* > \endverbatim */
/* > */
/* > \param[out] COLCND */
/* > \verbatim */
/* > COLCND is REAL */
/* > If INFO = 0, COLCND contains the ratio of the smallest */
/* > C(i) to the largest C(i). If COLCND >= 0.1, it is not */
/* > worth scaling by C. */
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
/* > > 0: if INFO = i, and i is */
/* > <= M: the i-th row of A is exactly zero */
/* > > M: the (i-M)-th column of A is exactly zero */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup realGBcomputational */
/* ===================================================================== */
/* Subroutine */
void sgbequ_(integer *m, integer *n, integer *kl, integer *ku, real *ab, integer *ldab, real *r__,
             real *c__, real *rowcnd, real *colcnd, real *amax, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgbequ inputs: m %" FLA_IS ", n %" FLA_IS ", kl %" FLA_IS ", ku %" FLA_IS
                      ", ldab %" FLA_IS "",
                      *m, *n, *kl, *ku, *ldab);
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2, r__3;
    /* Local variables */
    integer i__, j, kd;
    real rcmin, rcmax;
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum, smlnum;
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
    /* Test the input parameters */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --r__;
    --c__;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*kl < 0)
    {
        *info = -3;
    }
    else if(*ku < 0)
    {
        *info = -4;
    }
    else if(*ldab < *kl + *ku + 1)
    {
        *info = -6;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGBEQU", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        *rowcnd = 1.f;
        *colcnd = 1.f;
        *amax = 0.f;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine constants. */
    smlnum = slamch_("S");
    bignum = 1.f / smlnum;
    /* Compute row scale factors. */
    i__1 = *m;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        r__[i__] = 0.f;
        /* L10: */
    }
    /* Find the maximum element in each row. */
    kd = *ku + 1;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        /* Computing MAX */
        i__2 = j - *ku;
        /* Computing MIN */
        i__4 = j + *kl;
        i__3 = fla_min(i__4, *m);
        for(i__ = fla_max(i__2, 1); i__ <= i__3; ++i__)
        {
            /* Computing MAX */
            r__2 = r__[i__];
            r__3 = (r__1 = ab[kd + i__ - j + j * ab_dim1], f2c_abs(r__1)); // , expr subst
            r__[i__] = fla_max(r__2, r__3);
            /* L20: */
        }
        /* L30: */
    }
    /* Find the maximum and minimum scale factors. */
    rcmin = bignum;
    rcmax = 0.f;
    i__1 = *m;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        /* Computing MAX */
        r__1 = rcmax;
        r__2 = r__[i__]; // , expr subst
        rcmax = fla_max(r__1, r__2);
        /* Computing MIN */
        r__1 = rcmin;
        r__2 = r__[i__]; // , expr subst
        rcmin = fla_min(r__1, r__2);
        /* L40: */
    }
    *amax = rcmax;
    if(rcmin == 0.f)
    {
        /* Find the first zero scale factor and return an error code. */
        i__1 = *m;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if(r__[i__] == 0.f)
            {
                *info = i__;
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* L50: */
        }
    }
    else
    {
        /* Invert the scale factors. */
        i__1 = *m;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Computing MIN */
            /* Computing MAX */
            r__2 = r__[i__];
            r__1 = fla_max(r__2, smlnum);
            r__[i__] = 1.f / fla_min(r__1, bignum);
            /* L60: */
        }
        /* Compute ROWCND = fla_min(R(I)) / fla_max(R(I)) */
        *rowcnd = fla_max(rcmin, smlnum) / fla_min(rcmax, bignum);
    }
    /* Compute column scale factors */
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        c__[j] = 0.f;
        /* L70: */
    }
    /* Find the maximum element in each column, */
    /* assuming the row scaling computed above. */
    kd = *ku + 1;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        /* Computing MAX */
        i__3 = j - *ku;
        /* Computing MIN */
        i__4 = j + *kl;
        i__2 = fla_min(i__4, *m);
        for(i__ = fla_max(i__3, 1); i__ <= i__2; ++i__)
        {
            /* Computing MAX */
            r__2 = c__[j];
            r__3
                = (r__1 = ab[kd + i__ - j + j * ab_dim1], f2c_abs(r__1)) * r__[i__]; // , expr subst
            c__[j] = fla_max(r__2, r__3);
            /* L80: */
        }
        /* L90: */
    }
    /* Find the maximum and minimum scale factors. */
    rcmin = bignum;
    rcmax = 0.f;
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        /* Computing MIN */
        r__1 = rcmin;
        r__2 = c__[j]; // , expr subst
        rcmin = fla_min(r__1, r__2);
        /* Computing MAX */
        r__1 = rcmax;
        r__2 = c__[j]; // , expr subst
        rcmax = fla_max(r__1, r__2);
        /* L100: */
    }
    if(rcmin == 0.f)
    {
        /* Find the first zero scale factor and return an error code. */
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            if(c__[j] == 0.f)
            {
                *info = *m + j;
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* L110: */
        }
    }
    else
    {
        /* Invert the scale factors. */
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            /* Computing MIN */
            /* Computing MAX */
            r__2 = c__[j];
            r__1 = fla_max(r__2, smlnum);
            c__[j] = 1.f / fla_min(r__1, bignum);
            /* L120: */
        }
        /* Compute COLCND = fla_min(C(J)) / fla_max(C(J)) */
        *colcnd = fla_max(rcmin, smlnum) / fla_min(rcmax, bignum);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGBEQU */
}
/* sgbequ_ */
