/* ../netlib/dlarfg.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
    Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_f2c.h" /* > \brief \b DLARFG generates an elementary reflector (Householder matrix). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARFG + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* DOUBLE PRECISION ALPHA, TAU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARFG generates a real elementary reflector H of order n, such */
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
/* > ALPHA is DOUBLE PRECISION */
/* > On entry, the value alpha. */
/* > On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is DOUBLE PRECISION array, dimension */
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
/* > TAU is DOUBLE PRECISION */
/* > The value tau. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
void dlarfg_(integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *tau)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlarfg inputs: n %" FLA_IS ", incx %" FLA_IS "", *n, *incx);
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    /* Local variables */
    integer j, knt;
    doublereal beta;
    extern /* Subroutine */
        void
        dscal_(integer *, doublereal *, doublereal *, integer *);
    doublereal xnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *);
#if FLA_ENABLE_AMD_OPT
    extern int fla_dscal(integer * n, doublereal * da, doublereal * dx, integer * incx);
    extern doublereal fla_dnrm2_blas_kernel(integer *, doublereal *, integer *);
#else
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
#endif
    static TLS_CLASS_SPEC integer r_once = 1;
    static TLS_CLASS_SPEC doublereal safmin, rsafmn;
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
        *tau = 0.;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }

    i__1 = *n - 1;
#if FLA_ENABLE_AMD_OPT
    xnorm = fla_dnrm2_blas_kernel(&i__1, &x[1], incx);
#else
    xnorm = dnrm2_(&i__1, &x[1], incx);
#endif
    if(xnorm == 0.)
    {
        /* H = I */
        *tau = 0.;
    }
    else
    {
        /* general case */
        d__1 = dlapy2_(alpha, &xnorm);
        beta = -d_sign(&d__1, alpha);
        if(r_once)
        {
            safmin = dlamch_("S") / dlamch_("E");
            rsafmn = 1. / safmin;
            r_once = 0;
        }
        knt = 0;
        if(f2c_abs(beta) < safmin)
        {
            /* XNORM, BETA may be inaccurate;
            scale X and recompute them */
        L10:
            ++knt;
            i__1 = *n - 1;
#if FLA_ENABLE_AMD_OPT
        /* Inline DSCAL for small sizes */
            fla_dscal(&i__1, &rsafmn, &x[1], incx);
#else
            dscal_(&i__1, &rsafmn, &x[1], incx);
#endif
            beta *= rsafmn;
            *alpha *= rsafmn;
            if(f2c_abs(beta) < safmin && knt < 20)
            {
                goto L10;
            }
            /* New BETA is at most 1, at least SAFMIN */
            i__1 = *n - 1;
#if FLA_ENABLE_AMD_OPT
            xnorm = fla_dnrm2_blas_kernel(&i__1, &x[1], incx);
#else
            xnorm = dnrm2_(&i__1, &x[1], incx);
#endif
            d__1 = dlapy2_(alpha, &xnorm);
            beta = -d_sign(&d__1, alpha);
        }
        *tau = (beta - *alpha) / beta;
        i__1 = *n - 1;
        d__1 = 1. / (*alpha - beta);

#if FLA_ENABLE_AMD_OPT
        /* Inline DSCAL for small sizes */
        fla_dscal(&i__1, &d__1, &x[1], incx);
#else
        dscal_(&i__1, &d__1, &x[1], incx);
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
    /* End of DLARFG */
}
/* dlarfg_ */
