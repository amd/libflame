/* ./zrscl.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZDRSCL multiplies a vector by the reciprocal of a real scalar. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZDRSCL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zdrscl.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zdrscl.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zdrscl.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZRSCL( N, A, X, INCX ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* COMPLEX*16 A */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZRSCL multiplies an n-element complex vector x by the complex scalar */
/* > 1/a. This is done without overflow or underflow as long as */
/* > the final result x/a does not overflow or underflow. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of components of the vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX*16 */
/* > The scalar a which is used to divide each component of x. */
/* > A must not be 0, or the subroutine will divide by zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension */
/* > (1+(N-1)*f2c_dabs(INCX)) */
/* > The n-element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of the vector SX. */
/* > > 0: SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i), 1< i<= n */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
void zrscl_(integer *n, doublecomplex *a, doublecomplex *x, integer *incx)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zrscl inputs: n %" FLA_IS ",incx %" FLA_IS "", *n, *incx);
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    doublereal ai, ar, ui, ov, ur, absi, absr;
    extern /* Subroutine */
        void
        zscal_(integer *, doublecomplex *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *);
    doublereal safmin, safmax;
    extern /* Subroutine */
        void
        zdscal_(integer *, doublereal *, doublecomplex *, integer *),
        zdrscl_(integer *, doublereal *, doublecomplex *, integer *);
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
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Quick return if possible */
    /* Parameter adjustments */
    --x;
    /* Function Body */
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Get machine parameters */
    safmin = dlamch_("S");
    safmax = 1. / safmin;
    ov = dlamch_("O");
    /* Initialize constants related to A. */
    ar = a->r;
    ai = d_imag(a);
    absr = f2c_dabs(ar);
    absi = f2c_dabs(ai);
    if(ai == 0.)
    {
        /* If alpha is real, then we can use csrscl */
        zdrscl_(n, &ar, &x[1], incx);
    }
    else if(ar == 0.)
    {
        /* If alpha has a zero real part, then we follow the same rules as if */
        /* alpha were real. */
        if(absi > safmax)
        {
            zdscal_(n, &safmin, &x[1], incx);
            d__1 = -safmax / ai;
            z__1.r = 0.;
            z__1.i = d__1; // , expr subst
            zscal_(n, &z__1, &x[1], incx);
        }
        else if(absi < safmin)
        {
            d__1 = -safmin / ai;
            z__1.r = 0.;
            z__1.i = d__1; // , expr subst
            zscal_(n, &z__1, &x[1], incx);
            zdscal_(n, &safmax, &x[1], incx);
        }
        else
        {
            d__1 = -1. / ai;
            z__1.r = 0.;
            z__1.i = d__1; // , expr subst
            zscal_(n, &z__1, &x[1], incx);
        }
    }
    else
    {
        /* The following numbers can be computed. */
        /* They are the inverse of the real and imaginary parts of 1/alpha. */
        /* Note that a and b are always different from zero. */
        /* NaNs are only possible if either: */
        /* 1. alphaR or alphaI is NaN. */
        /* 2. alphaR and alphaI are both infinite, in which case it makes sense */
        /* to propagate a NaN. */
        ur = ar + ai * (ai / ar);
        ui = ai + ar * (ar / ai);
        if(f2c_dabs(ur) < safmin || f2c_dabs(ui) < safmin)
        {
            /* This means that both alphaR and alphaI are very small. */
            d__1 = safmin / ur;
            d__2 = -safmin / ui;
            z__1.r = d__1;
            z__1.i = d__2; // , expr subst
            zscal_(n, &z__1, &x[1], incx);
            zdscal_(n, &safmax, &x[1], incx);
        }
        else if(f2c_dabs(ur) > safmax || f2c_dabs(ui) > safmax)
        {
            if(absr > ov || absi > ov)
            {
                /* This means that a and b are both Inf. No need for scaling. */
                d__1 = 1. / ur;
                d__2 = -1. / ui;
                z__1.r = d__1;
                z__1.i = d__2; // , expr subst
                zscal_(n, &z__1, &x[1], incx);
            }
            else
            {
                zdscal_(n, &safmin, &x[1], incx);
                if(f2c_dabs(ur) > ov || f2c_dabs(ui) > ov)
                {
                    /* Infs were generated. We do proper scaling to avoid them. */
                    if(absr >= absi)
                    {
                        /* f2c_dabs( UR ) <= f2c_dabs( UI ) */
                        ur = safmin * ar + safmin * (ai * (ai / ar));
                        ui = safmin * ai + ar * (safmin * ar / ai);
                    }
                    else
                    {
                        /* f2c_dabs( UR ) > f2c_dabs( UI ) */
                        ur = safmin * ar + ai * (safmin * ai / ar);
                        ui = safmin * ai + safmin * (ar * (ar / ai));
                    }
                    d__1 = 1. / ur;
                    d__2 = -1. / ui;
                    z__1.r = d__1;
                    z__1.i = d__2; // , expr subst
                    zscal_(n, &z__1, &x[1], incx);
                }
                else
                {
                    d__1 = safmax / ur;
                    d__2 = -safmax / ui;
                    z__1.r = d__1;
                    z__1.i = d__2; // , expr subst
                    zscal_(n, &z__1, &x[1], incx);
                }
            }
        }
        else
        {
            d__1 = 1. / ur;
            d__2 = -1. / ui;
            z__1.r = d__1;
            z__1.i = d__2; // , expr subst
            zscal_(n, &z__1, &x[1], incx);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZRSCL */
}
/* zrscl_ */
