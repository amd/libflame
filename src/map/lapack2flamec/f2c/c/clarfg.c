/* ../netlib/clarfg.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b5 = {1.f, 0.f};
/* > \brief \b CLARFG generates an elementary reflector (Householder matrix). */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARFG + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfg.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfg.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfg.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* COMPLEX ALPHA, TAU */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFG generates a scomplex elementary reflector H of order n, such */
/* > that */
/* > */
/* > H**H * ( alpha ) = ( beta ), H**H * H = I. */
/* > ( x ) ( 0 ) */
/* > */
/* > where alpha and beta are scalars, with beta real, and x is an */
/* > (n-1)-element scomplex vector. H is represented in the form */
/* > */
/* > H = I - tau * ( 1 ) * ( 1 v**H ) , */
/* > ( v ) */
/* > */
/* > where tau is a scomplex scalar and v is a scomplex (n-1)-element */
/* > vector. Note that H is not hermitian. */
/* > */
/* > If the elements of x are all zero and alpha is real, then tau = 0 */
/* > and H is taken to be the unit matrix. */
/* > */
/* > Otherwise 1 <= real(tau) <= 2 and f2c_abs(tau-1) <= 1 . */
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
/* > ALPHA is COMPLEX */
/* > On entry, the value alpha. */
/* > On exit, it is overwritten with the value beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension */
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
/* > TAU is COMPLEX */
/* > The value tau. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2017 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clarfg_(aocl_int_t *n, scomplex *alpha, scomplex *x, aocl_int_t *incx, scomplex *tau)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clarfg(n, alpha, x, incx, tau);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;

    aocl_lapack_clarfg(&n_64, alpha, x, &incx_64, tau);
#endif
}

void aocl_lapack_clarfg(aocl_int64_t *n, scomplex *alpha, scomplex *x, aocl_int64_t *incx,
                        scomplex *tau)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarfg inputs: n %lld, incx %lld", *n, *incx);
#else
    snprintf(buffer, 256, "clarfg inputs: n %d, incx %d", *n, *incx);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t i__1;
    real r__1, r__2;
    scomplex q__1, q__2;
    /* Builtin functions */
    double r_sign(real *, real *);
    /* Local variables */
    aocl_int64_t j, knt;
    real beta;
    real alphi, alphr, xnorm;
    extern real slapy3_(real *, real *, real *);
    extern /* Complex */
        void
        cladiv_f2c_(scomplex *, scomplex *, scomplex *);
    extern real slamch_(char *);
    real safmin, rsafmn;
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
    if(*n <= 0)
    {
        tau->real = 0.f, tau->imag = 0.f;
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    i__1 = *n - 1;
    xnorm = aocl_blas_scnrm2(&i__1, &x[1], incx);
    alphr = alpha->real;
    alphi = alpha->imag;
    if(xnorm == 0.f && alphi == 0.f)
    {
        /* H = I */
        tau->real = 0.f, tau->imag = 0.f;
    }
    else
    {
        /* general case */
        r__1 = slapy3_(&alphr, &alphi, &xnorm);
        beta = -r_sign(&r__1, &alphr);
        safmin = slamch_("S") / slamch_("E");
        rsafmn = 1.f / safmin;
        knt = 0;
        if(f2c_abs(beta) < safmin)
        {
            /* XNORM, BETA may be inaccurate;
            scale X and recompute them */
        L10:
            ++knt;
            i__1 = *n - 1;
            aocl_blas_csscal(&i__1, &rsafmn, &x[1], incx);
            beta *= rsafmn;
            alphi *= rsafmn;
            alphr *= rsafmn;
            if(f2c_abs(beta) < safmin && knt < 20)
            {
                goto L10;
            }
            /* New BETA is at most 1, at least SAFMIN */
            i__1 = *n - 1;
            xnorm = aocl_blas_scnrm2(&i__1, &x[1], incx);
            q__1.real = alphr;
            q__1.imag = alphi; // , expr subst
            alpha->real = q__1.real, alpha->imag = q__1.imag;
            r__1 = slapy3_(&alphr, &alphi, &xnorm);
            beta = -r_sign(&r__1, &alphr);
        }
        r__1 = (beta - alphr) / beta;
        r__2 = -alphi / beta;
        q__1.real = r__1;
        q__1.imag = r__2; // , expr subst
        tau->real = q__1.real, tau->imag = q__1.imag;
        q__2.real = alpha->real - beta;
        q__2.imag = alpha->imag; // , expr subst
        cladiv_f2c_(&q__1, &c_b5, &q__2);
        alpha->real = q__1.real, alpha->imag = q__1.imag;
        i__1 = *n - 1;
        aocl_blas_cscal(&i__1, alpha, &x[1], incx);
        /* If ALPHA is subnormal, it may lose relative accuracy */
        i__1 = knt;
        for(j = 1; j <= i__1; ++j)
        {
            beta *= safmin;
            /* L20: */
        }
        alpha->real = beta, alpha->imag = 0.f;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLARFG */
}
/* clarfg_ */
