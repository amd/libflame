/* ./clarfgp.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b6 = {1.f, 0.f};
/* > \brief \b CLARFGP generates an elementary reflector (Householder matrix) with non-negative
 * beta. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLARFGP + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfgp
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfgp
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfgp
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARFGP( N, ALPHA, X, INCX, TAU ) */
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
/* > CLARFGP generates a scomplex elementary reflector H of order n, such */
/* > that */
/* > */
/* > H**H * ( alpha ) = ( beta ), H**H * H = I. */
/* > ( x ) ( 0 ) */
/* > */
/* > where alpha and beta are scalars, beta is real and non-negative, and */
/* > x is an (n-1)-element scomplex vector. H is represented in the form */
/* > */
/* > H = I - tau * ( 1 ) * ( 1 v**H ) , */
/* > ( v ) */
/* > */
/* > where tau is a scomplex scalar and v is a scomplex (n-1)-element */
/* > vector. Note that H is not hermitian. */
/* > */
/* > If the elements of x are all zero and alpha is real, then tau = 0 */
/* > and H is taken to be the unit matrix. */
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
/* > (1+(N-2)*abs(INCX)) */
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
/* > \ingroup larfgp */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clarfgp_(aocl_int_t *n, scomplex *alpha, scomplex *x, aocl_int_t *incx, scomplex *tau)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clarfgp(n, alpha, x, incx, tau);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;

    aocl_lapack_clarfgp(&n_64, alpha, x, &incx_64, tau);
#endif
}

void aocl_lapack_clarfgp(aocl_int64_t *n, scomplex *alpha, scomplex *x, aocl_int64_t *incx,
                         scomplex *tau)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clarfgp inputs: n %lld, incx %lld", *n, *incx);
#else
    snprintf(buffer, 256, "clarfgp inputs: n %d, incx %d", *n, *incx);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t i__1, i__2;
    real r__1, r__2;
    scomplex q__1, q__2;
    /* Builtin functions */
    double r_imag(scomplex *), c_abs(scomplex *), r_sign(real *, real *);
    /* Local variables */
    aocl_int64_t j;
    scomplex savealpha;
    real eps;
    aocl_int64_t knt;
    real beta;
    real alphi, alphr, xnorm;
    extern real slapy2_(real *, real *),
        slapy3_(real *, real *, real *);
    extern /* Complex */
        void
        cladiv_f2c_(scomplex *, scomplex *, scomplex *);
    extern real slamch_(char *);
    real bignum, smlnum;
    /* -- LAPACK auxiliary routine -- */
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
    eps = slamch_("Precision");
    i__1 = *n - 1;
    xnorm = aocl_blas_scnrm2(&i__1, &x[1], incx);
    alphr = alpha->real;
    alphi = r_imag(alpha);
    if(xnorm <= eps * c_abs(alpha))
    {
        /* H = [1-alpha/abs(alpha) 0;
        0 I], sign chosen so ALPHA >= 0. */
        if(alphi == 0.f)
        {
            if(alphr >= 0.f)
            {
                /* When TAU.eq.ZERO, the vector is special-cased to be */
                /* all zeros in the application routines. We do not need */
                /* to clear it. */
                tau->real = 0.f, tau->imag = 0.f;
            }
            else
            {
                /* However, the application routines rely on explicit */
                /* zero checks when TAU.ne.ZERO, and we must clear X. */
                tau->real = 2.f, tau->imag = 0.f;
                i__1 = *n - 1;
                for(j = 1; j <= i__1; ++j)
                {
                    i__2 = (j - 1) * *incx + 1;
                    x[i__2].real = 0.f;
                    x[i__2].imag = 0.f; // , expr subst
                }
                q__1.real = -alpha->real;
                q__1.imag = -alpha->imag; // , expr subst
                alpha->real = q__1.real, alpha->imag = q__1.imag;
            }
        }
        else
        {
            /* Only "reflecting" the diagonal entry to be real and non-negative. */
            xnorm = slapy2_(&alphr, &alphi);
            r__1 = 1.f - alphr / xnorm;
            r__2 = -alphi / xnorm;
            q__1.real = r__1;
            q__1.imag = r__2; // , expr subst
            tau->real = q__1.real, tau->imag = q__1.imag;
            i__1 = *n - 1;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = (j - 1) * *incx + 1;
                x[i__2].real = 0.f;
                x[i__2].imag = 0.f; // , expr subst
            }
            alpha->real = xnorm, alpha->imag = 0.f;
        }
    }
    else
    {
        /* general case */
        r__1 = slapy3_(&alphr, &alphi, &xnorm);
        beta = r_sign(&r__1, &alphr);
        smlnum = slamch_("S") / slamch_("E");
        bignum = 1.f / smlnum;
        knt = 0;
        if(f2c_abs(beta) < smlnum)
        {
        /* XNORM, BETA may be inaccurate;
        scale X and recompute them */
        L10:
            ++knt;
            i__1 = *n - 1;
            aocl_blas_csscal(&i__1, &bignum, &x[1], incx);
            beta *= bignum;
            alphi *= bignum;
            alphr *= bignum;
            if(f2c_abs(beta) < smlnum && knt < 20)
            {
                goto L10;
            }
            /* New BETA is at most 1, at least SMLNUM */
            i__1 = *n - 1;
            xnorm = aocl_blas_scnrm2(&i__1, &x[1], incx);
            q__1.real = alphr;
            q__1.imag = alphi; // , expr subst
            alpha->real = q__1.real, alpha->imag = q__1.imag;
            r__1 = slapy3_(&alphr, &alphi, &xnorm);
            beta = r_sign(&r__1, &alphr);
        }
        savealpha.real = alpha->real;
        savealpha.imag = alpha->imag; // , expr subst
        q__1.real = alpha->real + beta;
        q__1.imag = alpha->imag; // , expr subst
        alpha->real = q__1.real, alpha->imag = q__1.imag;
        if(beta < 0.f)
        {
            beta = -beta;
            q__2.real = -alpha->real;
            q__2.imag = -alpha->imag; // , expr subst
            q__1.real = q__2.real / beta;
            q__1.imag = q__2.imag / beta; // , expr subst
            tau->real = q__1.real, tau->imag = q__1.imag;
        }
        else
        {
            alphr = alphi * (alphi / alpha->real);
            alphr += xnorm * (xnorm / alpha->real);
            r__1 = alphr / beta;
            r__2 = -alphi / beta;
            q__1.real = r__1;
            q__1.imag = r__2; // , expr subst
            tau->real = q__1.real, tau->imag = q__1.imag;
            r__1 = -alphr;
            q__1.real = r__1;
            q__1.imag = alphi; // , expr subst
            alpha->real = q__1.real, alpha->imag = q__1.imag;
        }
        cladiv_f2c_(&q__1, &c_b6, alpha);
        alpha->real = q__1.real, alpha->imag = q__1.imag;
        if(c_abs(tau) <= smlnum)
        {
            /* In the case where the computed TAU ends up being a denormalized number, */
            /* it loses relative accuracy. This is a BIG problem. Solution: flush TAU */
            /* to ZERO (or TWO or whatever makes a nonnegative real number for BETA). */
            /* (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.) */
            /* (Thanks Pat. Thanks MathWorks.) */
            alphr = savealpha.real;
            alphi = r_imag(&savealpha);
            if(alphi == 0.f)
            {
                if(alphr >= 0.f)
                {
                    tau->real = 0.f, tau->imag = 0.f;
                }
                else
                {
                    tau->real = 2.f, tau->imag = 0.f;
                    i__1 = *n - 1;
                    for(j = 1; j <= i__1; ++j)
                    {
                        i__2 = (j - 1) * *incx + 1;
                        x[i__2].real = 0.f;
                        x[i__2].imag = 0.f; // , expr subst
                    }
                    q__1.real = -savealpha.real;
                    q__1.imag = -savealpha.imag; // , expr subst
                    beta = q__1.real;
                }
            }
            else
            {
                xnorm = slapy2_(&alphr, &alphi);
                r__1 = 1.f - alphr / xnorm;
                r__2 = -alphi / xnorm;
                q__1.real = r__1;
                q__1.imag = r__2; // , expr subst
                tau->real = q__1.real, tau->imag = q__1.imag;
                i__1 = *n - 1;
                for(j = 1; j <= i__1; ++j)
                {
                    i__2 = (j - 1) * *incx + 1;
                    x[i__2].real = 0.f;
                    x[i__2].imag = 0.f; // , expr subst
                }
                beta = xnorm;
            }
        }
        else
        {
            /* This is the general case. */
            i__1 = *n - 1;
            aocl_blas_cscal(&i__1, alpha, &x[1], incx);
        }
        /* If BETA is subnormal, it may lose relative accuracy */
        i__1 = knt;
        for(j = 1; j <= i__1; ++j)
        {
            beta *= smlnum;
            /* L20: */
        }
        alpha->real = beta, alpha->imag = 0.f;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CLARFGP */
}
/* clarfgp_ */
