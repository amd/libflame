/* zlartg.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Subroutine */
/* > \brief \b ZLARTG generates a plane rotation with real cosine and scomplex sine. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARTG( F, G, C, S, R ) */

/*       .. Scalar Arguments .. */
/*       REAL(wp)              C */
/*       COMPLEX(wp)           F, G, R, S */
/*       .. */

/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARTG generates a plane rotation so that */
/* > */
/* >    [  C         S  ] . [ F ]  =  [ R ] */
/* >    [ -conjg(S)  C  ]   [ G ]     [ 0 ] */
/* > */
/* > where C is real and C**2 + |S|**2 = 1. */
/* > */
/* > The mathematical formulas used for C and S are */
/* > */
/* >    sgn(x) = {  x / |x|,   x != 0 */
/* >             {  1,         x  = 0 */
/* > */
/* >    R = sgn(F) * sqrt(|F|**2 + |G|**2) */
/* > */
/* >    C = |F| / sqrt(|F|**2 + |G|**2) */
/* > */
/* >    S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2) */
/* > */
/* > Special conditions: */
/* >    If G=0, then C=1 and S=0. */
/* >    If F=0, then C=0 and S is chosen so that R is real. */
/* > */
/* > When F and G are real, the formulas simplify to C = F/R and */
/* > S = G/R, and the returned values of C, S, and R should be */
/* > identical to those returned by DLARTG. */
/* > */
/* > The algorithm used to compute these quantities incorporates scaling */
/* > to avoid overflow or underflow in computing the square root of the */
/* > sum of squares. */
/* > */
/* > This is the same routine ZROTG fom BLAS1, except that */
/* > F and G are unchanged on return. */
/* > */
/* > Below, wp=>dp stands for double precision from LA_CONSTANTS module. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] F */
/* > \verbatim */
/* >          F is COMPLEX(wp) */
/* >          The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is COMPLEX(wp) */
/* >          The second component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[out] C */
/* > \verbatim */
/* >          C is REAL(wp) */
/* >          The cosine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is COMPLEX(wp) */
/* >          The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is COMPLEX(wp) */
/* >          The nonzero component of the rotated vector. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Weslley Pereira, University of Colorado Denver, USA */

/* > \date December 2021 */

/* > \ingroup OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Based on the algorithm from */
/* > */
/* >  Anderson E. (2017) */
/* >  Algorithm 978: Safe Scaling in the Level 1 BLAS */
/* >  ACM Trans Math Softw 44:1--28 */
/* >  https://doi.org/10.1145/3061665 */
/* > */
/* > \endverbatim */
void zlartg_(dcomplex *f, dcomplex *g, doublereal *c__, dcomplex *s,
             dcomplex *r__)
{
    AOCL_DTL_TRACE_ENTRY_INDENT
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    double sqrt(doublereal), d_imag(dcomplex *);
    void d_cnjg(dcomplex *, dcomplex *),
        z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    doublereal d__, u, v, w, f1, f2, g1, g2, h2;
    doublereal d__1, d__2, d__3, d__4;
    dcomplex fs, gs;
    doublereal rtmin, rtmax, safmin, safmax;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 04:17:29 1/20/23 */
    /*  -- LAPACK auxiliary routine -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     February 2021 */
    /* .. Constants .. */
    safmin = 2.2250738585072014e-308;
    safmax = 4.4942328371557898e307;
    rtmin = sqrt(safmin);
    /* .. */
    /* .. Executable Statements .. */
    if(g->real == 0. && g->imag == 0.)
    {
        *c__ = 1.;
        s->real = 0., s->imag = 0.;
        r__->real = f->real, r__->imag = f->imag;
    }
    else if(f->real == 0. && f->imag == 0.)
    {
        *c__ = 0.;
        if(g->real == 0.)
        {
            d__2 = (d__1 = d_imag(g), f2c_dabs(d__1));
            r__->real = d__2, r__->imag = 0.;
            d_cnjg(&z__2, g);
            z_div(&z__1, &z__2, r__);
            s->real = z__1.real, s->imag = z__1.imag;
        }
        else if(d_imag(g) == 0.)
        {
            d__2 = (d__1 = g->real, f2c_dabs(d__1));
            r__->real = d__2, r__->imag = 0.;
            d_cnjg(&z__2, g);
            z_div(&z__1, &z__2, r__);
            s->real = z__1.real, s->imag = z__1.imag;
        }
        else
        {
            /* Computing MAX */
            d__3 = (d__1 = g->real, f2c_dabs(d__1));
            d__4 = (d__2 = d_imag(g), f2c_dabs(d__2)); // , expr subst
            g1 = fla_max(d__3, d__4);
            rtmax = sqrt(safmax / 2);
            if(g1 > rtmin && g1 < rtmax)
            {
                /* Use unscaled algorithm */
                /* The following two lines can be replaced by `d = f2c_dabs( g )`. */
                /* This algorithm do not use the intrinsic scomplex abs. */
                /* Computing 2nd power */
                d__1 = g->real;
                /* Computing 2nd power */
                d__2 = d_imag(g);
                g2 = d__1 * d__1 + d__2 * d__2;
                d__ = sqrt(g2);
                d_cnjg(&z__2, g);
                z__1.real = z__2.real / d__;
                z__1.imag = z__2.imag / d__; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                r__->real = d__, r__->imag = 0.;
            }
            else
            {
                /* Use scaled algorithm */
                /* Computing MIN */
                d__1 = safmax;
                d__2 = fla_max(safmin, g1); // , expr subst
                u = fla_min(d__1, d__2);
                z__1.real = g->real / u;
                z__1.imag = g->imag / u; // , expr subst
                gs.real = z__1.real;
                gs.imag = z__1.imag; // , expr subst
                /* The following two lines can be replaced by `d = f2c_dabs( gs )`. */
                /* This algorithm do not use the intrinsic scomplex abs. */
                /* Computing 2nd power */
                d__1 = gs.real;
                /* Computing 2nd power */
                d__2 = d_imag(&gs);
                g2 = d__1 * d__1 + d__2 * d__2;
                d__ = sqrt(g2);
                d_cnjg(&z__2, &gs);
                z__1.real = z__2.real / d__;
                z__1.imag = z__2.imag / d__; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
                d__1 = d__ * u;
                r__->real = d__1, r__->imag = 0.;
            }
        }
    }
    else
    {
        /* Computing MAX */
        d__3 = (d__1 = f->real, f2c_dabs(d__1));
        d__4 = (d__2 = d_imag(f), f2c_dabs(d__2)); // , expr subst
        f1 = fla_max(d__3, d__4);
        /* Computing MAX */
        d__3 = (d__1 = g->real, f2c_dabs(d__1));
        d__4 = (d__2 = d_imag(g), f2c_dabs(d__2)); // , expr subst
        g1 = fla_max(d__3, d__4);
        rtmax = sqrt(safmax / 4);
        if(f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax)
        {
            /* Use unscaled algorithm */
            /* Computing 2nd power */
            d__1 = f->real;
            /* Computing 2nd power */
            d__2 = d_imag(f);
            f2 = d__1 * d__1 + d__2 * d__2;
            /* Computing 2nd power */
            d__1 = g->real;
            /* Computing 2nd power */
            d__2 = d_imag(g);
            g2 = d__1 * d__1 + d__2 * d__2;
            h2 = f2 + g2;
            /* safmin <= f2 <= h2 <= safmax */
            if(f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                z__1.real = f->real / *c__;
                z__1.imag = f->imag / *c__; // , expr subst
                r__->real = z__1.real, r__->imag = z__1.imag;
                rtmax *= 2;
                if(f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    d_cnjg(&z__2, g);
                    d__1 = sqrt(f2 * h2);
                    z__3.real = f->real / d__1;
                    z__3.imag = f->imag / d__1; // , expr subst
                    z__1.real = z__2.real * z__3.real - z__2.imag * z__3.imag;
                    z__1.imag = z__2.real * z__3.imag + z__2.imag * z__3.real; // , expr subst
                    s->real = z__1.real, s->imag = z__1.imag;
                }
                else
                {
                    d_cnjg(&z__2, g);
                    z__3.real = r__->real / h2;
                    z__3.imag = r__->imag / h2; // , expr subst
                    z__1.real = z__2.real * z__3.real - z__2.imag * z__3.imag;
                    z__1.imag = z__2.real * z__3.imag + z__2.imag * z__3.real; // , expr subst
                    s->real = z__1.real, s->imag = z__1.imag;
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
                    z__1.real = f->real / *c__;
                    z__1.imag = f->imag / *c__; // , expr subst
                    r__->real = z__1.real, r__->imag = z__1.imag;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    d__1 = h2 / d__;
                    z__1.real = d__1 * f->real;
                    z__1.imag = d__1 * f->imag; // , expr subst
                    r__->real = z__1.real, r__->imag = z__1.imag;
                }
                d_cnjg(&z__2, g);
                z__3.real = f->real / d__;
                z__3.imag = f->imag / d__; // , expr subst
                z__1.real = z__2.real * z__3.real - z__2.imag * z__3.imag;
                z__1.imag = z__2.real * z__3.imag + z__2.imag * z__3.real; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
            }
        }
        else
        {
            /* Use scaled algorithm */
            /* Computing MIN */
            /* Computing MAX */
            d__3 = fla_max(safmin, f1);
            d__1 = safmax;
            d__2 = fla_max(d__3, g1); // , expr subst
            u = fla_min(d__1, d__2);
            z__1.real = g->real / u;
            z__1.imag = g->imag / u; // , expr subst
            gs.real = z__1.real;
            gs.imag = z__1.imag; // , expr subst
            /* Computing 2nd power */
            d__1 = gs.real;
            /* Computing 2nd power */
            d__2 = d_imag(&gs);
            g2 = d__1 * d__1 + d__2 * d__2;
            if(f1 / u < rtmin)
            {
                /* f is not well-scaled when scaled by g1. */
                /* Use a different scaling for f. */
                /* Computing MIN */
                d__1 = safmax;
                d__2 = fla_max(safmin, f1); // , expr subst
                v = fla_min(d__1, d__2);
                w = v / u;
                z__1.real = f->real / v;
                z__1.imag = f->imag / v; // , expr subst
                fs.real = z__1.real;
                fs.imag = z__1.imag; // , expr subst
                /* Computing 2nd power */
                d__1 = fs.real;
                /* Computing 2nd power */
                d__2 = d_imag(&fs);
                f2 = d__1 * d__1 + d__2 * d__2;
                /* Computing 2nd power */
                d__1 = w;
                h2 = f2 * (d__1 * d__1) + g2;
            }
            else
            {
                /* Otherwise use the same scaling for f and g. */
                w = 1.;
                z__1.real = f->real / u;
                z__1.imag = f->imag / u; // , expr subst
                fs.real = z__1.real;
                fs.imag = z__1.imag; // , expr subst
                /* Computing 2nd power */
                d__1 = fs.real;
                /* Computing 2nd power */
                d__2 = d_imag(&fs);
                f2 = d__1 * d__1 + d__2 * d__2;
                h2 = f2 + g2;
            }
            /* safmin <= f2 <= h2 <= safmax */
            if(f2 >= h2 * safmin)
            {
                /* safmin <= f2/h2 <= 1, and h2/f2 is finite */
                *c__ = sqrt(f2 / h2);
                z__1.real = fs.real / *c__;
                z__1.imag = fs.imag / *c__; // , expr subst
                r__->real = z__1.real, r__->imag = z__1.imag;
                rtmax *= 2;
                if(f2 > rtmin && h2 < rtmax)
                {
                    /* safmin <= sqrt( f2*h2 ) <= safmax */
                    d_cnjg(&z__2, &gs);
                    d__1 = sqrt(f2 * h2);
                    z__3.real = fs.real / d__1;
                    z__3.imag = fs.imag / d__1; // , expr subst
                    z__1.real = z__2.real * z__3.real - z__2.imag * z__3.imag;
                    z__1.imag = z__2.real * z__3.imag + z__2.imag * z__3.real; // , expr subst
                    s->real = z__1.real, s->imag = z__1.imag;
                }
                else
                {
                    d_cnjg(&z__2, &gs);
                    z__3.real = r__->real / h2;
                    z__3.imag = r__->imag / h2; // , expr subst
                    z__1.real = z__2.real * z__3.real - z__2.imag * z__3.imag;
                    z__1.imag = z__2.real * z__3.imag + z__2.imag * z__3.real; // , expr subst
                    s->real = z__1.real, s->imag = z__1.imag;
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
                    z__1.real = fs.real / *c__;
                    z__1.imag = fs.imag / *c__; // , expr subst
                    r__->real = z__1.real, r__->imag = z__1.imag;
                }
                else
                {
                    /* f2 / sqrt(f2 * h2) < safmin, then */
                    /* sqrt(safmin) <= f2 * sqrt(safmax) <= h2 / sqrt(f2 * h2 */
                    d__1 = h2 / d__;
                    z__1.real = d__1 * fs.real;
                    z__1.imag = d__1 * fs.imag; // , expr subst
                    r__->real = z__1.real, r__->imag = z__1.imag;
                }
                d_cnjg(&z__2, &gs);
                z__3.real = fs.real / d__;
                z__3.imag = fs.imag / d__; // , expr subst
                z__1.real = z__2.real * z__3.real - z__2.imag * z__3.imag;
                z__1.imag = z__2.real * z__3.imag + z__2.imag * z__3.real; // , expr subst
                s->real = z__1.real, s->imag = z__1.imag;
            }
            /* Rescale c and r */
            *c__ *= w;
            z__1.real = u * r__->real;
            z__1.imag = u * r__->imag; // , expr subst
            r__->real = z__1.real, r__->imag = z__1.imag;
        }
    }
    AOCL_DTL_TRACE_EXIT_INDENT
    return;
}
/* zlartg_ */
