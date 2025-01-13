/* ./crscl.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CRSCL multiplies a vector by the reciprocal of a real scalar. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CRSCL + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/crscl.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/crscl.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/crscl.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CRSCL( N, A, X, INCX ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* COMPLEX A */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CRSCL multiplies an n-element complex vector x by the complex scalar */
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
/* > A is COMPLEX */
/* > The scalar a which is used to divide each component of x. */
/* > A must not be 0, or the subroutine will divide by zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension */
/* > (1+(N-1)*abs(INCX)) */
/* > The n-element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of the vector X. */
/* > > 0: X(1) = X(1) and X(1+(i-1)*INCX) = x(i), 1< i<= n */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
void crscl_(integer *n, complex *a, complex *x, integer *incx)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("crscl inputs: n %" FLA_IS ",incx %" FLA_IS "", *n, *incx);

    /* System generated locals */
    real r__1, r__2;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    real ai, ar, ui, ov, ur, absi, absr;
    extern /* Subroutine */
        void
        cscal_(integer *, complex *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        csscal_(integer *, real *, complex *, integer *);
    real safmin, safmax;
    extern /* Subroutine */
        void
        csrscl_(integer *, real *, complex *, integer *);
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
    safmin = slamch_("S");
    safmax = 1.f / safmin;
    ov = slamch_("O");
    /* Initialize constants related to A. */
    ar = a->r;
    ai = r_imag(a);
    absr = f2c_abs(ar);
    absi = f2c_abs(ai);
    if(ai == 0.f)
    {
        /* If alpha is real, then we can use csrscl */
        csrscl_(n, &ar, &x[1], incx);
    }
    else if(ar == 0.f)
    {
        /* If alpha has a zero real part, then we follow the same rules as if */
        /* alpha were real. */
        if(absi > safmax)
        {
            csscal_(n, &safmin, &x[1], incx);
            r__1 = -safmax / ai;
            q__1.r = 0.f;
            q__1.i = r__1; // , expr subst
            cscal_(n, &q__1, &x[1], incx);
        }
        else if(absi < safmin)
        {
            r__1 = -safmin / ai;
            q__1.r = 0.f;
            q__1.i = r__1; // , expr subst
            cscal_(n, &q__1, &x[1], incx);
            csscal_(n, &safmax, &x[1], incx);
        }
        else
        {
            r__1 = -1.f / ai;
            q__1.r = 0.f;
            q__1.i = r__1; // , expr subst
            cscal_(n, &q__1, &x[1], incx);
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
        if(f2c_abs(ur) < safmin || f2c_abs(ui) < safmin)
        {
            /* This means that both alphaR and alphaI are very small. */
            r__1 = safmin / ur;
            r__2 = -safmin / ui;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            cscal_(n, &q__1, &x[1], incx);
            csscal_(n, &safmax, &x[1], incx);
        }
        else if(f2c_abs(ur) > safmax || f2c_abs(ui) > safmax)
        {
            if(absr > ov || absi > ov)
            {
                /* This means that a and b are both Inf. No need for scaling. */
                r__1 = 1.f / ur;
                r__2 = -1.f / ui;
                q__1.r = r__1;
                q__1.i = r__2; // , expr subst
                cscal_(n, &q__1, &x[1], incx);
            }
            else
            {
                csscal_(n, &safmin, &x[1], incx);
                if(f2c_abs(ur) > ov || f2c_abs(ui) > ov)
                {
                    /* Infs were generated. We do proper scaling to avoid them. */
                    if(absr >= absi)
                    {
                        /* ABS( UR ) <= ABS( UI ) */
                        ur = safmin * ar + safmin * (ai * (ai / ar));
                        ui = safmin * ai + ar * (safmin * ar / ai);
                    }
                    else
                    {
                        /* ABS( UR ) > ABS( UI ) */
                        ur = safmin * ar + ai * (safmin * ai / ar);
                        ui = safmin * ai + safmin * (ar * (ar / ai));
                    }
                    r__1 = 1.f / ur;
                    r__2 = -1.f / ui;
                    q__1.r = r__1;
                    q__1.i = r__2; // , expr subst
                    cscal_(n, &q__1, &x[1], incx);
                }
                else
                {
                    r__1 = safmax / ur;
                    r__2 = -safmax / ui;
                    q__1.r = r__1;
                    q__1.i = r__2; // , expr subst
                    cscal_(n, &q__1, &x[1], incx);
                }
            }
        }
        else
        {
            r__1 = 1.f / ur;
            r__2 = -1.f / ui;
            q__1.r = r__1;
            q__1.i = r__2; // , expr subst
            cscal_(n, &q__1, &x[1], incx);
        }
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of CRSCL */
}
/* crscl_ */
