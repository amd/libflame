/* ../netlib/dlasv2.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b3 = 2.;
static doublereal c_b4 = 1.;
/* > \brief \b DLASV2 computes the singular value decomposition of a 2-by-2 triangular matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLASV2 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasv2.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasv2.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasv2.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASV2 computes the singular value decomposition of a 2-by-2 */
/* > triangular matrix */
/* > [ F G ] */
/* > [ 0 H ]. */
/* > On return, f2c_dabs(SSMAX) is the larger singular value, f2c_dabs(SSMIN) is the */
/* > smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and */
/* > right singular vectors for f2c_dabs(SSMAX), giving the decomposition */
/* > */
/* > [ CSL SNL ] [ F G ] [ CSR -SNR ] = [ SSMAX 0 ] */
/* > [-SNL CSL ] [ 0 H ] [ SNR CSR ] [ 0 SSMIN ]. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] F */
/* > \verbatim */
/* > F is DOUBLE PRECISION */
/* > The (1,1) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* > G is DOUBLE PRECISION */
/* > The (1,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] H */
/* > \verbatim */
/* > H is DOUBLE PRECISION */
/* > The (2,2) element of the 2-by-2 matrix. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMIN */
/* > \verbatim */
/* > SSMIN is DOUBLE PRECISION */
/* > f2c_dabs(SSMIN) is the smaller singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] SSMAX */
/* > \verbatim */
/* > SSMAX is DOUBLE PRECISION */
/* > f2c_dabs(SSMAX) is the larger singular value. */
/* > \endverbatim */
/* > */
/* > \param[out] SNL */
/* > \verbatim */
/* > SNL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] CSL */
/* > \verbatim */
/* > CSL is DOUBLE PRECISION */
/* > The vector (CSL, SNL) is a unit left singular vector for the */
/* > singular value f2c_dabs(SSMAX). */
/* > \endverbatim */
/* > */
/* > \param[out] SNR */
/* > \verbatim */
/* > SNR is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] CSR */
/* > \verbatim */
/* > CSR is DOUBLE PRECISION */
/* > The vector (CSR, SNR) is a unit right singular vector for the */
/* > singular value f2c_dabs(SSMAX). */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup auxOTHERauxiliary */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Any input parameter may be aliased with any output parameter. */
/* > */
/* > Barring over/underflow and assuming a guard digit in subtraction, all */
/* > output quantities are correct to within a few units in the last */
/* > place (ulps). */
/* > */
/* > In IEEE arithmetic, the code works correctly if one matrix element is */
/* > infinite. */
/* > */
/* > Overflow will not occur unless the largest singular value itself */
/* > overflows or is within a few ulps of overflow. (On machines with */
/* > partial overflow, like the Cray, overflow may occur if the largest */
/* > singular value is within a factor of 2 of overflow.) */
/* > */
/* > Underflow is harmless if underflow is gradual. Otherwise, results */
/* > may correspond to a matrix modified by perturbations of size near */
/* > the underflow threshold. */
/* > \endverbatim */
/* > */
/* ===================================================================== */
/* Subroutine */
void dlasv2_(doublereal *f, doublereal *g, doublereal *h__, doublereal *ssmin, doublereal *ssmax,
             doublereal *snr, doublereal *csr, doublereal *snl, doublereal *csl)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlasv2 inputs: f %lf, g %lf, h %lf", *f, *g, *h__);
    /* System generated locals */
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    doublereal a, d__, l, m, r__, s, t, fa, ga, ha, ft, gt, ht, mm, tt, clt, crt, slt, srt;
    integer pmax;
    doublereal temp;
    logical swap;
    doublereal tsign;
    extern doublereal dlamch_(char *);
    logical gasmal;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    ft = *f;
    fa = f2c_dabs(ft);
    ht = *h__;
    ha = f2c_dabs(*h__);
    /* PMAX points to the maximum absolute element of matrix */
    /* PMAX = 1 if F largest in absolute values */
    /* PMAX = 2 if G largest in absolute values */
    /* PMAX = 3 if H largest in absolute values */
    pmax = 1;
    swap = ha > fa;
    if(swap)
    {
        pmax = 3;
        temp = ft;
        ft = ht;
        ht = temp;
        temp = fa;
        fa = ha;
        ha = temp;
        /* Now FA .ge. HA */
    }
    gt = *g;
    ga = f2c_dabs(gt);
    if(ga == 0.)
    {
        /* Diagonal matrix */
        *ssmin = ha;
        *ssmax = fa;
        clt = 1.;
        crt = 1.;
        slt = 0.;
        srt = 0.;
    }
    else
    {
        gasmal = TRUE_;
        if(ga > fa)
        {
            pmax = 2;
            if(fa / ga < dlamch_("EPS"))
            {
                /* Case of very large GA */
                gasmal = FALSE_;
                *ssmax = ga;
                if(ha > 1.)
                {
                    *ssmin = fa / (ga / ha);
                }
                else
                {
                    *ssmin = fa / ga * ha;
                }
                clt = 1.;
                slt = ht / gt;
                srt = 1.;
                crt = ft / gt;
            }
        }
        if(gasmal)
        {
            /* Normal case */
            d__ = fa - ha;
            if(d__ == fa)
            {
                /* Copes with infinite F or H */
                l = 1.;
            }
            else
            {
                l = d__ / fa;
            }
            /* Note that 0 .le. L .le. 1 */
            m = gt / ft;
            /* Note that f2c_dabs(M) .le. 1/macheps */
            t = 2. - l;
            /* Note that T .ge. 1 */
            mm = m * m;
            tt = t * t;
            s = sqrt(tt + mm);
            /* Note that 1 .le. S .le. 1 + 1/macheps */
            if(l == 0.)
            {
                r__ = f2c_dabs(m);
            }
            else
            {
                r__ = sqrt(l * l + mm);
            }
            /* Note that 0 .le. R .le. 1 + 1/macheps */
            a = (s + r__) * .5;
            /* Note that 1 .le. A .le. 1 + f2c_dabs(M) */
            *ssmin = ha / a;
            *ssmax = fa * a;
            if(mm == 0.)
            {
                /* Note that M is very tiny */
                if(l == 0.)
                {
                    t = d_sign(&c_b3, &ft) * d_sign(&c_b4, &gt);
                }
                else
                {
                    t = gt / d_sign(&d__, &ft) + m / t;
                }
            }
            else
            {
                t = (m / (s + t) + m / (r__ + l)) * (a + 1.);
            }
            l = sqrt(t * t + 4.);
            crt = 2. / l;
            srt = t / l;
            clt = (crt + srt * m) / a;
            slt = ht / ft * srt / a;
        }
    }
    if(swap)
    {
        *csl = srt;
        *snl = crt;
        *csr = slt;
        *snr = clt;
    }
    else
    {
        *csl = clt;
        *snl = slt;
        *csr = crt;
        *snr = srt;
    }
    /* Correct signs of SSMAX and SSMIN */
    if(pmax == 1)
    {
        tsign = d_sign(&c_b4, csr) * d_sign(&c_b4, csl) * d_sign(&c_b4, f);
    }
    if(pmax == 2)
    {
        tsign = d_sign(&c_b4, snr) * d_sign(&c_b4, csl) * d_sign(&c_b4, g);
    }
    if(pmax == 3)
    {
        tsign = d_sign(&c_b4, snr) * d_sign(&c_b4, snl) * d_sign(&c_b4, h__);
    }
    *ssmax = d_sign(ssmax, &tsign);
    d__1 = tsign * d_sign(&c_b4, f) * d_sign(&c_b4, h__);
    *ssmin = d_sign(ssmin, &d__1);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLASV2 */
}
/* dlasv2_ */
