/* ../netlib/slarfg.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLARFG generates an elementary reflector (Householder matrix). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLARFG + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfg.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfg.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfg.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* REAL ALPHA, TAU */
/* .. */
/* .. Array Arguments .. */
/* REAL X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARFG generates a real elementary reflector H of order n, such */
/* > that */
/* > */
/* > H * ( alpha ) = ( beta ), H**T * H = I. */
/* > ( x ) ( 0 ) */
/* > */
/* > where alpha and beta are scalars, and x is an (n-1)-element real */
/* > vector. H is represented in the form */
/* > */
/* > H = I - tau * ( 1 ) * ( 1 v**T ) , */
/* > ( v ) */
/* > */
/* > where tau is a real scalar and v is a real (n-1)-element */
/* > vector. */
/* > */
/* > If the elements of x are all zero, then tau = 0 and H is taken to be */
/* > the unit matrix. */
/* > */
/* > Otherwise 1 <= tau <= 2. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the elementary reflector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ALPHA */
/* > \verbatim */
/* > ALPHA is REAL */
/* > On entry, the value alpha. */
/* > On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is REAL array, dimension */
/* > (1+(N-2)*f2c_abs(INCX)) */
/* > On entry, the vector x. */
/* > On exit, it is overwritten with the vector v. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* > TAU is REAL */
/* > The value tau. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup realOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
void slarfg_(integer *n, real *alpha, real *x, integer *incx, real *tau)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slarfg inputs: n %" FLA_IS ", incx %" FLA_IS "", *n, *incx);
    /* System generated locals */
    integer i__1;
    real r__1;
    /* Builtin functions */
    double r_sign(real *, real *);
    /* Local variables */
    integer j, knt;
    real beta;
    extern real snrm2_(integer *, real *, integer *);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    real xnorm;
    extern real slapy2_(real *, real *), slamch_(char *);
    real safmin, rsafmn;
#if FLA_ENABLE_AMD_OPT
    void fla_sscal(integer * n, real * alpha, real * x, integer * incx);
#endif

    /* -- LAPACK auxiliary routine (version 3.8.0) -- */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --x;
    /* Function Body */
    if(*n <= 1)
    {
        *tau = 0.f;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    i__1 = *n - 1;
    xnorm = snrm2_(&i__1, &x[1], incx);
    if(xnorm == 0.f)
    {
        /* H = I */
        *tau = 0.f;
    }
    else
    {
        /* general case */
        r__1 = slapy2_(alpha, &xnorm);
        beta = -r_sign(&r__1, alpha);
        safmin = slamch_("S") / slamch_("E");
        knt = 0;
        if(f2c_abs(beta) < safmin)
        {
            /* XNORM, BETA may be inaccurate;
            scale X and recompute them */
            rsafmn = 1.f / safmin;
        L10:
            ++knt;
            i__1 = *n - 1;
#if FLA_ENABLE_AMD_OPT
            /* Inline SSCAL for small sizes */
            fla_sscal(&i__1, &rsafmn, &x[1], incx);
#else
            sscal_(&i__1, &rsafmn, &x[1], incx);
#endif
            beta *= rsafmn;
            *alpha *= rsafmn;
            if(f2c_abs(beta) < safmin && knt < 20)
            {
                goto L10;
            }
            /* New BETA is at most 1, at least SAFMIN */
            i__1 = *n - 1;
            xnorm = snrm2_(&i__1, &x[1], incx);
            r__1 = slapy2_(alpha, &xnorm);
            beta = -r_sign(&r__1, alpha);
        }
        *tau = (beta - *alpha) / beta;
        i__1 = *n - 1;
        r__1 = 1.f / (*alpha - beta);
#if FLA_ENABLE_AMD_OPT
        /* Inline SSCAL for small sizes */
        fla_sscal(&i__1, &r__1, &x[1], incx);
#else
        sscal_(&i__1, &r__1, &x[1], incx);
#endif
        /* If ALPHA is subnormal, it may lose relative accuracy */
        i__1 = knt;
        for(j = 1; j <= i__1; ++j)
        {
            beta *= safmin;
            /* L20: */
        }
        *alpha = beta;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SLARFG */
}
/* slarfg_ */
