/* ../netlib/ssyequb.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b SSYEQUB */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYEQUB + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyequb
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyequb
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyequb
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, LDA, N */
/* REAL AMAX, SCOND */
/* CHARACTER UPLO */
/* .. */
/* .. Array Arguments .. */
/* REAL A( LDA, * ), S( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYEQUB computes row and column scalings intended to equilibrate a */
/* > symmetric matrix A (with respect to the Euclidean norm) and reduce */
/* > its condition number. The scale factors S are computed by the BIN */
/* > algorithm (see references) so that the scaled matrix B with elements */
/* > B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of */
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
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The N-by-N symmetric matrix whose scaling factors are to be */
/* > computed. */
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
/* > Largest absolute value of any matrix element. If AMAX is */
/* > very close to overflow or very close to underflow, the */
/* > matrix should be scaled. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (2*N) */
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
/* > \date November 2017 */
/* > \ingroup realSYcomputational */
/* > \par References: */
/* ================ */
/* > */
/* > Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n */
/* > Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n */
/* > DOI 10.1023/B:NUMA.0000016606.32820.69 \n */
/* > Tech report version: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.3.1679 */
/* > */
/* ===================================================================== */
/* Subroutine */
void ssyequb_(char *uplo, integer *n, real *a, integer *lda, real *s, real *scond, real *amax,
              real *work, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ssyequb inputs: uplo %c ,n %" FLA_IS ",lda %" FLA_IS "", *uplo, *n, *lda);
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1, r__2, r__3;
    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_ri(real *, integer *);
    /* Local variables */
    real d__;
    integer i__, j;
    real t, u, c0, c1, c2, si;
    logical up;
    real avg, std, tol, base;
    integer iter;
    real smin, smax, scale;
    extern logical lsame_(char *, char *, integer, integer);
    real sumsq;
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    extern /* Subroutine */
        void
        slassq_(integer *, real *, integer *, real *, real *);
    real smlnum;
    /* -- LAPACK computational routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2017 */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --s;
    --work;
    /* Function Body */
    *info = 0;
    if(!(lsame_(uplo, "U", 1, 1) || lsame_(uplo, "L", 1, 1)))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYEQUB", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    up = lsame_(uplo, "U", 1, 1);
    *amax = 0.f;
    /* Quick return if possible. */
    if(*n == 0)
    {
        *scond = 1.f;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        s[i__] = 0.f;
    }
    *amax = 0.f;
    if(up)
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = j - 1;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                /* Computing MAX */
                r__2 = s[i__];
                r__3 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)); // , expr subst
                s[i__] = fla_max(r__2, r__3);
                /* Computing MAX */
                r__2 = s[j];
                r__3 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)); // , expr subst
                s[j] = fla_max(r__2, r__3);
                /* Computing MAX */
                r__2 = *amax;
                r__3 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)); // , expr subst
                *amax = fla_max(r__2, r__3);
            }
            /* Computing MAX */
            r__2 = s[j];
            r__3 = (r__1 = a[j + j * a_dim1], f2c_abs(r__1)); // , expr subst
            s[j] = fla_max(r__2, r__3);
            /* Computing MAX */
            r__2 = *amax;
            r__3 = (r__1 = a[j + j * a_dim1], f2c_abs(r__1)); // , expr subst
            *amax = fla_max(r__2, r__3);
        }
    }
    else
    {
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            r__2 = s[j];
            r__3 = (r__1 = a[j + j * a_dim1], f2c_abs(r__1)); // , expr subst
            s[j] = fla_max(r__2, r__3);
            /* Computing MAX */
            r__2 = *amax;
            r__3 = (r__1 = a[j + j * a_dim1], f2c_abs(r__1)); // , expr subst
            *amax = fla_max(r__2, r__3);
            i__2 = *n;
            for(i__ = j + 1; i__ <= i__2; ++i__)
            {
                /* Computing MAX */
                r__2 = s[i__];
                r__3 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)); // , expr subst
                s[i__] = fla_max(r__2, r__3);
                /* Computing MAX */
                r__2 = s[j];
                r__3 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)); // , expr subst
                s[j] = fla_max(r__2, r__3);
                /* Computing MAX */
                r__2 = *amax;
                r__3 = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)); // , expr subst
                *amax = fla_max(r__2, r__3);
            }
        }
    }
    i__1 = *n;
    for(j = 1; j <= i__1; ++j)
    {
        s[j] = 1.f / s[j];
    }
    tol = 1.f / sqrt(*n * 2.f);
    for(iter = 1; iter <= 100; ++iter)
    {
        scale = 0.f;
        sumsq = 0.f;
        /* beta = |A|s */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            work[i__] = 0.f;
        }
        if(up)
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = j - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    work[i__] += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)) * s[j];
                    work[j] += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)) * s[i__];
                }
                work[j] += (r__1 = a[j + j * a_dim1], f2c_abs(r__1)) * s[j];
            }
        }
        else
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                work[j] += (r__1 = a[j + j * a_dim1], f2c_abs(r__1)) * s[j];
                i__2 = *n;
                for(i__ = j + 1; i__ <= i__2; ++i__)
                {
                    work[i__] += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)) * s[j];
                    work[j] += (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1)) * s[i__];
                }
            }
        }
        /* avg = s^T beta / n */
        avg = 0.f;
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            avg += s[i__] * work[i__];
        }
        avg /= *n;
        std = 0.f;
        i__1 = *n << 1;
        for(i__ = *n + 1; i__ <= i__1; ++i__)
        {
            work[i__] = s[i__ - *n] * work[i__ - *n] - avg;
        }
        slassq_(n, &work[*n + 1], &c__1, &scale, &sumsq);
        std = scale * sqrt(sumsq / *n);
        if(std < tol * avg)
        {
            goto L999;
        }
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            t = (r__1 = a[i__ + i__ * a_dim1], f2c_abs(r__1));
            si = s[i__];
            c2 = (*n - 1) * t;
            c1 = (*n - 2) * (work[i__] - t * si);
            c0 = -(t * si) * si + work[i__] * 2 * si - *n * avg;
            d__ = c1 * c1 - c0 * 4 * c2;
            if(d__ <= 0.f)
            {
                *info = -1;
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            si = c0 * -2 / (c1 + sqrt(d__));
            d__ = si - s[i__];
            u = 0.f;
            if(up)
            {
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    t = (r__1 = a[j + i__ * a_dim1], f2c_abs(r__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
                i__2 = *n;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    t = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
            }
            else
            {
                i__2 = i__;
                for(j = 1; j <= i__2; ++j)
                {
                    t = (r__1 = a[i__ + j * a_dim1], f2c_abs(r__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
                i__2 = *n;
                for(j = i__ + 1; j <= i__2; ++j)
                {
                    t = (r__1 = a[j + i__ * a_dim1], f2c_abs(r__1));
                    u += s[j] * t;
                    work[j] += d__ * t;
                }
            }
            avg += (u + work[i__]) * d__ / *n;
            s[i__] = si;
        }
    }
L999:
    smlnum = slamch_("SAFEMIN");
    bignum = 1.f / smlnum;
    smin = bignum;
    smax = 0.f;
    t = 1.f / sqrt(avg);
    base = slamch_("B");
    u = 1.f / log(base);
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = (integer)(u * log(s[i__] * t));
        s[i__] = pow_ri(&base, &i__2);
        /* Computing MIN */
        r__1 = smin;
        r__2 = s[i__]; // , expr subst
        smin = fla_min(r__1, r__2);
        /* Computing MAX */
        r__1 = smax;
        r__2 = s[i__]; // , expr subst
        smax = fla_max(r__1, r__2);
    }
    *scond = fla_max(smin, smlnum) / fla_min(smax, bignum);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* ssyequb_ */
