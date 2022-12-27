/* ./clartg.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
/**
 * Modifications Copyright (C) 2014-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
#include "FLA_f2c.h" /* > \brief \b CLARTG generates a plane rotation with real cosine and scomplex sine. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/  */
/* > \htmlonly */
/* > Download CLARTG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clartg. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clartg. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clartg. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLARTG( F, G, C, S, R ) */
/* .. Scalar Arguments .. */
/* REAL(wp) C */
/* COMPLEX(wp) F, G, R, S */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARTG generates a plane rotation so that */
/* > */
/* > [ C S ] . [ F ] = [ R ] */
/* > [ -conjg(S) C ] [ G ] [ 0 ] */
/* > */
/* > where C is real and C**2 + |S|**2 = 1. */
/* > */
/* > The mathematical formulas used for C and S are */
/* > */
/* > sgn(x) = {
x / |x|, x != 0 */
/* > {
1, x = 0 */
/* > */
/* > R = sgn(F) * sqrt(|F|**2 + |G|**2) */
/* > */
/* > C = |F| / sqrt(|F|**2 + |G|**2) */
/* > */
/* > S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2) */
/* > */
/* > Special conditions: */
/* > If G=0, then C=1 and S=0. */
/* > If F=0, then C=0 and S is chosen so that R is real. */
/* > */
/* > When F and G are real, the formulas simplify to C = F/R and */
/* > S = G/R, and the returned values of C, S, and R should be */
/* > identical to those returned by SLARTG. */
/* > */
/* > The algorithm used to compute these quantities incorporates scaling */
/* > to avoid overflow or underflow in computing the square root of the */
/* > sum of squares. */
/* > */
/* > This is the same routine CROTG fom BLAS1, except that */
/* > F and G are unchanged on return. */
/* > */
/* > Below, wp=>sp stands for single precision from LA_CONSTANTS module. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] F */
/* > \verbatim */
/* > F is COMPLEX(wp) */
/* > The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* > G is COMPLEX(wp) */
/* > The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* > C is REAL(wp) */
/* > The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* > S is COMPLEX(wp) */
/* > The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* > R is COMPLEX(wp) */
/* > The nonzero component of the rotated vector. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Weslley Pereira, University of Colorado Denver, USA */
/* > \date December 2021 */
/* > \ingroup lartg */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Based on the algorithm from */
/* > */
/* > Anderson E. (2017) */
/* > Algorithm 978: Safe Scaling in the Level 1 BLAS */
/* > ACM Trans Math Softw 44:1--28 */
/* > https://doi.org/10.1145/3061665 */
/* > */
/* > \endverbatim */
/* Subroutine */
void clartg_(scomplex *f, scomplex *g, real *c__, scomplex *s, scomplex *r__)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    double log(doublereal), pow_ri(real *, integer *), r_imag(complex *), c_abs(complex *), sqrt(doublereal);
    /* Local variables */
    real d__, u, v, w, f1, f2, g1, g2, h2;
    scomplex fs, gs, f__t, g__t;
    real rtmin, rtmax, safmin, safmax;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 00:33:35 2/21/25 */
    /* ...Switches: */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Constants .. */
    safmin = 1.1754943508222875e-38f;
    safmax = 8.5070591730234616e37f;
    rtmin = sqrt(safmin);
    /* .. */
    /* .. Executable Statements .. */
    safmin = slamch_("S");
    eps = slamch_("E");
    r__1 = slamch_("B");
    i__1 = (integer) (log(safmin / eps) / log(slamch_("B")) / 2.f);
    safmn2 = pow_ri(&r__1, &i__1);
    safmx2 = 1.f / safmn2;
    /* Computing MAX */
    /* Computing MAX */
    r__7 = (r__1 = f->r, f2c_abs(r__1));
    r__8 = (r__2 = f->i, f2c_abs(r__2)); // , expr subst
    /* Computing MAX */
    r__9 = (r__3 = g->r, f2c_abs(r__3));
    r__10 = (r__4 = g->i, f2c_abs(r__4)); // , expr subst
    r__5 = fla_max(r__7,r__8);
    r__6 = fla_max(r__9,r__10); // , expr subst
    scale = fla_max(r__5,r__6);
    fs.r = f->r;
    fs.i = f->i; // , expr subst
    gs.r = g->r;
    gs.i = g->i; // , expr subst
    count = 0;
    if (scale >= safmx2)
    {
        *c__ = 1.f;
        s->real = 0.f, s->imag = 0.f;
        r__->real = f__t.real, r__->imag = f__t.imag;
    }
    else if(f__t.real == 0.f && f__t.imag == 0.f)
    {
        r__1 = c_abs(g);
        if (g->r == 0.f && g->i == 0.f || sisnan_(&r__1))
        {
            *cs = 1.f;
            sn->r = 0.f, sn->i = 0.f;
            r__->r = f->r, r__->i = f->i;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return 0;
        }
        else if(r_imag(&g__t) == 0.f)
        {
            goto L20;
        }
    }
    /* Computing 2nd power */
    r__1 = fs.r;
    /* Computing 2nd power */
    r__2 = fs.i;
    f2 = r__1 * r__1 + r__2 * r__2;
    /* Computing 2nd power */
    r__1 = gs.r;
    /* Computing 2nd power */
    r__2 = gs.i;
    g2 = r__1 * r__1 + r__2 * r__2;
    if (f2 <= fla_max(g2,1.f) * safmin)
    {
        /* This is a rare case: F is very small. */
        if (f->r == 0.f && f->i == 0.f)
        {
            *cs = 0.f;
            r__2 = g->r;
            r__3 = g->i;
            r__1 = slapy2_(&r__2, &r__3);
            r__->r = r__1, r__->i = 0.f;
            /* Do complex/real division explicitly with two real divisions */
            r__1 = gs.r;
            r__2 = gs.i;
            d__ = slapy2_(&r__1, &r__2);
            r__1 = gs.r / d__;
            r__2 = -gs.i / d__;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            sn->r = q__1.r, sn->i = q__1.i;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return 0;
        }
        r__1 = fs.r;
        r__2 = fs.i;
        f2s = slapy2_(&r__1, &r__2);
        /* G2 and G2S are accurate */
        /* G2 is at least SAFMIN, and G2S is at least SAFMN2 */
        g2s = sqrt(g2);
        /* Error in CS from underflow in F2S is at most */
        /* UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS */
        /* If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN, */
        /* and so CS .lt. sqrt(SAFMIN) */
        /* If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN */
        /* and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS) */
        /* Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S */
        *cs = f2s / g2s;
        /* Make sure f2c_abs(FF) = 1 */
        /* Do complex/real division explicitly with 2 real divisions */
        /* Computing MAX */
        r__3 = (r__1 = f->r, f2c_abs(r__1));
        r__4 = (r__2 = f->i, f2c_abs(r__2)); // , expr subst
        if (fla_max(r__3,r__4) > 1.f)
        {
            r__1 = f->r;
            r__2 = f->i;
            d__ = slapy2_(&r__1, &r__2);
            r__1 = f->r / d__;
            r__2 = f->i / d__;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            ff.r = q__1.r;
            ff.i = q__1.i; // , expr subst
        }
        else
        {
            dr = safmx2 * f->r;
            di = safmx2 * f->i;
            d__ = slapy2_(&dr, &di);
            r__1 = dr / d__;
            r__2 = di / d__;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            ff.r = q__1.r;
            ff.i = q__1.i; // , expr subst
        }
        r__1 = gs.r / g2s;
        r__2 = -gs.i / g2s;
        q__2.r = r__1;
        q__2.i = r__2; // , expr subst
        q__1.r = ff.r * q__2.r - ff.i * q__2.i;
        q__1.i = ff.r * q__2.i + ff.i * q__2.r; // , expr subst
        sn->r = q__1.r, sn->i = q__1.i;
        q__2.r = *cs * f->r;
        q__2.i = *cs * f->i; // , expr subst
        q__3.r = sn->r * g->r - sn->i * g->i;
        q__3.i = sn->r * g->i + sn->i * g->r; // , expr subst
        q__1.r = q__2.r + q__3.r;
        q__1.i = q__2.i + q__3.i; // , expr subst
        r__->r = q__1.r, r__->i = q__1.i;
    }
    else
    {
        /* This is the most common case. */
        /* Neither F2 nor F2/G2 are less than SAFMIN */
        /* F2S cannot overflow, and it is accurate */
        f2s = sqrt(g2 / f2 + 1.f);
        /* Do the F2S(real)*FS(complex) multiply with two real multiplies */
        r__1 = f2s * fs.r;
        r__2 = f2s * fs.i;
        q__1.r = r__1;
        q__1.i = r__2; // , expr subst
        r__->r = q__1.r, r__->i = q__1.i;
        *cs = 1.f / f2s;
        d__ = f2 + g2;
        /* Do complex/real division explicitly with two real divisions */
        r__1 = r__->r / d__;
        r__2 = r__->i / d__;
        q__1.r = r__1;
        q__1.i = r__2; // , expr subst
        sn->r = q__1.r, sn->i = q__1.i;
        q__2.r = gs.r;
        q__2.i = -gs.i;
        q__1.r = sn->r * q__2.r - sn->i * q__2.i;
        q__1.i = sn->r * q__2.i + sn->i * q__2.r; // , expr subst
        sn->r = q__1.r, sn->i = q__1.i;
        if (count != 0)
        {
            /* Use unscaled algorithm */
            /* Computing 2nd power */
            r__1 = f__t.real;
            /* Computing 2nd power */
            r__2 = r_imag(&f__t);
            f2 = r__1 * r__1 + r__2 * r__2;
            /* Computing 2nd power */
            r__1 = g__t.real;
            /* Computing 2nd power */
            r__2 = r_imag(&g__t);
            g2 = r__1 * r__1 + r__2 * r__2;
            h2 = f2 + g2;
            /* safmin <= f2 <= h2 <= safmax */
            if(f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                q__1.real = f__t.real / *c__;
                q__1.imag = f__t.imag / *c__; // , expr subst
                r__->real = q__1.real, r__->imag = q__1.imag;
                rtmax *= 2;
                if(f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    r_cnjg(&q__2, &g__t);
                    r__1 = sqrt(f2 * h2);
                    q__3.real = f__t.real / r__1;
                    q__3.imag = f__t.imag / r__1; // , expr subst
                    q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                    q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                    s->real = q__1.real, s->imag = q__1.imag;
                }
                else
                {
                    r_cnjg(&q__2, &g__t);
                    q__3.real = r__->real / h2;
                    q__3.imag = r__->imag / h2; // , expr subst
                    q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                    q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                    s->real = q__1.real, s->imag = q__1.imag;
                }
            }
            else
            {
                /* f2/h2 <= safmin may be subnormal, and h2/f2 may overflow. */
                /* Moreover, */
                /* safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= sa */
                /* sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax). */
                /* Also, */
                /* g2 >> f2, which means that h2 = g2. */
                d__ = sqrt(f2 * h2);
                *c__ = f2 / d__;
                if(*c__ >= safmin)
                {
                    q__1.real = f__t.real / *c__;
                    q__1.imag = f__t.imag / *c__; // , expr subst
                    r__->real = q__1.real, r__->imag = q__1.imag;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    r__1 = h2 / d__;
                    q__1.real = r__1 * f__t.real;
                    q__1.imag = r__1 * f__t.imag; // , expr subst
                    r__->real = q__1.real, r__->imag = q__1.imag;
                }
                r_cnjg(&q__2, &g__t);
                q__3.real = f__t.real / d__;
                q__3.imag = f__t.imag / d__; // , expr subst
                q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
            }
        }
        else
        {
            /* Use scaled algorithm */
            /* Computing MIN */
            /* Computing MAX */
            r__3 = fla_max(safmin, f1);
            r__1 = safmax;
            r__2 = fla_max(r__3, g1); // , expr subst
            u = fla_min(r__1, r__2);
            q__1.real = g__t.real / u;
            q__1.imag = g__t.imag / u; // , expr subst
            gs.real = q__1.real;
            gs.imag = q__1.imag; // , expr subst
            /* Computing 2nd power */
            r__1 = gs.real;
            /* Computing 2nd power */
            r__2 = r_imag(&gs);
            g2 = r__1 * r__1 + r__2 * r__2;
            if(f1 / u < rtmin)
            {
                /* f is not well-scaled when scaled by g1. */
                /* Use a different scaling for f. */
                /* Computing MIN */
                r__1 = safmax;
                r__2 = fla_max(safmin, f1); // , expr subst
                v = fla_min(r__1, r__2);
                w = v / u;
                q__1.real = f__t.real / v;
                q__1.imag = f__t.imag / v; // , expr subst
                fs.real = q__1.real;
                fs.imag = q__1.imag; // , expr subst
                /* Computing 2nd power */
                r__1 = fs.real;
                /* Computing 2nd power */
                r__2 = r_imag(&fs);
                f2 = r__1 * r__1 + r__2 * r__2;
                /* Computing 2nd power */
                r__1 = w;
                h2 = f2 * (r__1 * r__1) + g2;
            }
            else
            {
                /* Otherwise use the same scaling for f and g. */
                w = 1.f;
                q__1.real = f__t.real / u;
                q__1.imag = f__t.imag / u; // , expr subst
                fs.real = q__1.real;
                fs.imag = q__1.imag; // , expr subst
                /* Computing 2nd power */
                r__1 = fs.real;
                /* Computing 2nd power */
                r__2 = r_imag(&fs);
                f2 = r__1 * r__1 + r__2 * r__2;
                h2 = f2 + g2;
            }
            /* safmin <= f2 <= h2 <= safmax */
            if(f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                q__1.real = fs.real / *c__;
                q__1.imag = fs.imag / *c__; // , expr subst
                r__->real = q__1.real, r__->imag = q__1.imag;
                rtmax *= 2;
                if(f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    r_cnjg(&q__2, &gs);
                    r__1 = sqrt(f2 * h2);
                    q__3.real = fs.real / r__1;
                    q__3.imag = fs.imag / r__1; // , expr subst
                    q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                    q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                    s->real = q__1.real, s->imag = q__1.imag;
                }
                else
                {
                    r_cnjg(&q__2, &gs);
                    q__3.real = r__->real / h2;
                    q__3.imag = r__->imag / h2; // , expr subst
                    q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                    q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                    s->real = q__1.real, s->imag = q__1.imag;
                }
            }
            else
            {
                /* f2/h2 <= safmin may be subnormal, and h2/f2 may overflow. */
                /* Moreover, */
                /* safmin <= f2*f2 * safmax < f2 * h2 < h2*h2 * safmin <= sa */
                /* sqrt(safmin) <= sqrt(f2 * h2) <= sqrt(safmax). */
                /* Also, */
                /* g2 >> f2, which means that h2 = g2. */
                d__ = sqrt(f2 * h2);
                *c__ = f2 / d__;
                if(*c__ >= safmin)
                {
                    q__1.real = fs.real / *c__;
                    q__1.imag = fs.imag / *c__; // , expr subst
                    r__->real = q__1.real, r__->imag = q__1.imag;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    r__1 = h2 / d__;
                    q__1.real = r__1 * fs.real;
                    q__1.imag = r__1 * fs.imag; // , expr subst
                    r__->real = q__1.real, r__->imag = q__1.imag;
                }
                r_cnjg(&q__2, &gs);
                q__3.real = fs.real / d__;
                q__3.imag = fs.imag / d__; // , expr subst
                q__1.real = q__2.real * q__3.real - q__2.imag * q__3.imag;
                q__1.imag = q__2.real * q__3.imag + q__2.imag * q__3.real; // , expr subst
                s->real = q__1.real, s->imag = q__1.imag;
            }
            /* Rescale c and r */
            *c__ *= w;
            q__1.real = u * r__->real;
            q__1.imag = u * r__->imag; // , expr subst
            r__->real = q__1.real, r__->imag = q__1.imag;
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of CLARTG */
}
/* clartg_ */