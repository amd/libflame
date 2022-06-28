/* ../netlib/dlartg.f -- translated by f2c (version 20100827). You must link the resulting object file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLARTG generates a plane rotation with real cosine and real sine. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at  */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARTG + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartg. f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartg. f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartg. f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARTG( F, G, CS, SN, R ) */
/* .. Scalar Arguments .. */
/* DOUBLE PRECISION CS, F, G, R, SN */
/* .. */
/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARTG generates a plane rotation so that */
/* > */
/* >    [  C  S  ]  .  [ F ]  =  [ R ] */
/* >    [ -S  C  ]     [ G ]     [ 0 ] */
/* > */
/* > where C**2 + S**2 = 1. */
/* > */
/* > The mathematical formulas used for C and S are */
/* >    R = sign(F) * sqrt(F**2 + G**2) */
/* >    C = F / R */
/* >    S = G / R */
/* > Hence C >= 0. The algorithm used to compute these quantities */
/* > incorporates scaling to avoid overflow or underflow in computing the */
/* > square root of the sum of squares. */
/* > */
/* > This version is discontinuous in R at F = 0 but it returns the same */
/* > C and S as ZLARTG for scomplex inputs (F,0) and (G,0). */
/* > */
/* > This is a more accurate version of the BLAS1 routine DROTG, */
/* > with the following other differences: */
/* >    F and G are unchanged on return. */
/* >    If G=0, then C=1 and S=0. */
/* >    If F=0 and (G .ne. 0), then C=0 and S=sign(1,G) without doing any */
/* >       floating point operations (saves work in DBDSQR when */
/* >       there are zeros on the diagonal). */
/* > */
/* > Below, wp=>dp stands for double precision from LA_CONSTANTS module. */
/* > \endverbatim */
/*  Arguments: */
/*  ========== */
/* > \param[in] F */
/* > \verbatim */
/* >          F is REAL(wp) */
/* >          The first component of vector to be rotated. */
/* > \endverbatim */
/* > */
/* > \param[in] G */
/* > \verbatim */
/* >          G is REAL(wp) */
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
/* >          S is REAL(wp) */
/* >          The sine of the rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] R */
/* > \verbatim */
/* >          R is REAL(wp) */
/* >          The nonzero component of the rotated vector. */
/* > \endverbatim */
/*  Authors: */
/*  ======== */
/* > \author Edward Anderson, Lockheed Martin */
/* > \date July 2016 */
/* > \ingroup lartg */
/* > \par Contributors: */
/*  ================== */
/* > */
/* > Weslley Pereira, University of Colorado Denver, USA */
/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Anderson E. (2017) */
/* >  Algorithm 978: Safe Scaling in the Level 1 BLAS */
/* >  ACM Trans Math Softw 44:1--28 */
/* >  https://doi.org/10.1145/3061665 */
/* > */
/* > \endverbatim */
static doublereal c_b2 = 1.;
/* Subroutine */
void dlartg_(doublereal *f, doublereal *g, doublereal *c__, doublereal *s, doublereal *r__)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlartg inputs : f %lf, g %lf", *f, *g);
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    doublereal d__, u, f1, g1, fs, gs, rtmin, rtmax, safmin, safmax;
    doublereal d__1, d__2, d__3, f__t, g__t;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 05:19:29 1/25/23 */
    /*  -- LAPACK auxiliary routine -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     February 2021 */
    /* .. */
    /* .. Constants .. */
    safmin = 2.2250738585072014e-308;
    safmax = 4.4942328371557898e307;
    rtmin = sqrt(safmin);
    rtmax = sqrt(safmax / 2);
    /* .. */
    /* .. Executable Statements .. */
    f__t = *f;
    g__t = *g;
    f1 = f2c_dabs(f__t);
    g1 = f2c_dabs(g__t);
    if(g__t == 0.)
    {
        *c__ = 1.;
        *s = 0.;
        *r__ = f__t;
    }
    else if(f__t == 0.)
    {
        *c__ = 0.;
        *s = d_sign(&c_b2, &g__t);
        *r__ = g1;
    }
    else if(f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax)
    {
        d__ = sqrt(f__t * f__t + g__t * g__t);
        *c__ = f1 / d__;
        *r__ = d_sign(&d__, &f__t);
        *s = g__t / *r__;
    }
    else
    {
        /* Computing MIN */
        /* Computing MAX */
        d__1 = f2c_dabs(f1);
        d__2 = f2c_dabs(g1); // , expr subst
        scale = max(d__1,d__2);
        if (scale >= safmx2)
        {
            count = 0;
L10:
            ++count;
            f1 *= safmn2;
            g1 *= safmn2;
            /* Computing MAX */
            d__1 = f2c_dabs(f1);
            d__2 = f2c_dabs(g1); // , expr subst
            scale = max(d__1,d__2);
            if (scale >= safmx2)
            {
                goto L10;
            }
            /* Computing 2nd power */
            d__1 = f1;
            /* Computing 2nd power */
            d__2 = g1;
            *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
            *cs = f1 / *r__;
            *sn = g1 / *r__;
            i__1 = count;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                *r__ *= safmx2;
                /* L20: */
            }
        }
        else if (scale <= safmn2)
        {
            count = 0;
L30:
            ++count;
            f1 *= safmx2;
            g1 *= safmx2;
            /* Computing MAX */
            d__1 = f2c_dabs(f1);
            d__2 = f2c_dabs(g1); // , expr subst
            scale = max(d__1,d__2);
            if (scale <= safmn2)
            {
                goto L30;
            }
            /* Computing 2nd power */
            d__1 = f1;
            /* Computing 2nd power */
            d__2 = g1;
            *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
            *cs = f1 / *r__;
            *sn = g1 / *r__;
            i__1 = count;
            for (i__ = 1;
                    i__ <= i__1;
                    ++i__)
            {
                *r__ *= safmn2;
                /* L40: */
            }
        }
        else
        {
            /* Computing 2nd power */
            d__1 = f1;
            /* Computing 2nd power */
            d__2 = g1;
            *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
            *cs = f1 / *r__;
            *sn = g1 / *r__;
        }
        if (f2c_dabs(*f) > f2c_dabs(*g) && *cs < 0.)
        {
            *cs = -(*cs);
            *sn = -(*sn);
            *r__ = -(*r__);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return 0;
    /* End of DLARTG */
}
/* dlartg_ */
