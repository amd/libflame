/* ../netlib/dlasq3.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLASQ3 checks for deflation, computes a shift and calls dqds. Used by sbdsqr. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLASQ3 + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq3.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasq3.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasq3.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, */
/* ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, */
/* DN2, G, TAU ) */
/* .. Scalar Arguments .. */
/* LOGICAL IEEE */
/* INTEGER I0, ITER, N0, NDIV, NFAIL, PP */
/* DOUBLE PRECISION DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, */
/* $ QMAX, SIGMA, TAU */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION Z( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds. */
/* > In case of failure it changes shifts, and tries again until output */
/* > is positive. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] I0 */
/* > \verbatim */
/* > I0 is INTEGER */
/* > First index. */
/* > \endverbatim */
/* > */
/* > \param[in,out] N0 */
/* > \verbatim */
/* > N0 is INTEGER */
/* > Last index. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* > Z is DOUBLE PRECISION array, dimension ( 4*N0 ) */
/* > Z holds the qd array. */
/* > \endverbatim */
/* > */
/* > \param[in,out] PP */
/* > \verbatim */
/* > PP is INTEGER */
/* > PP=0 for ping, PP=1 for pong. */
/* > PP=2 indicates that flipping was applied to the Z array */
/* > and that the initial tests for deflation should not be */
/* > performed. */
/* > \endverbatim */
/* > */
/* > \param[out] DMIN */
/* > \verbatim */
/* > DMIN is DOUBLE PRECISION */
/* > Minimum value of d. */
/* > \endverbatim */
/* > */
/* > \param[out] SIGMA */
/* > \verbatim */
/* > SIGMA is DOUBLE PRECISION */
/* > Sum of shifts used in current segment. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DESIG */
/* > \verbatim */
/* > DESIG is DOUBLE PRECISION */
/* > Lower order part of SIGMA */
/* > \endverbatim */
/* > */
/* > \param[in] QMAX */
/* > \verbatim */
/* > QMAX is DOUBLE PRECISION */
/* > Maximum value of q. */
/* > \endverbatim */
/* > */
/* > \param[in,out] NFAIL */
/* > \verbatim */
/* > NFAIL is INTEGER */
/* > Increment NFAIL by 1 each time the shift was too big. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ITER */
/* > \verbatim */
/* > ITER is INTEGER */
/* > Increment ITER by 1 for each iteration. */
/* > \endverbatim */
/* > */
/* > \param[in,out] NDIV */
/* > \verbatim */
/* > NDIV is INTEGER */
/* > Increment NDIV by 1 for each division. */
/* > \endverbatim */
/* > */
/* > \param[in] IEEE */
/* > \verbatim */
/* > IEEE is LOGICAL */
/* > Flag for IEEE or non IEEE arithmetic (passed to DLASQ5). */
/* > \endverbatim */
/* > */
/* > \param[in,out] TTYPE */
/* > \verbatim */
/* > TTYPE is INTEGER */
/* > Shift type. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DMIN1 */
/* > \verbatim */
/* > DMIN1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DMIN2 */
/* > \verbatim */
/* > DMIN2 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DN */
/* > \verbatim */
/* > DN is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DN1 */
/* > \verbatim */
/* > DN1 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] DN2 */
/* > \verbatim */
/* > DN2 is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] G */
/* > \verbatim */
/* > G is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in,out] TAU */
/* > \verbatim */
/* > TAU is DOUBLE PRECISION */
/* > */
/* > These are passed as arguments in order to save their values */
/* > between calls to DLASQ3. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2016 */
/* > \ingroup auxOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dlasq3_(aocl_int_t *i0, aocl_int_t *n0, doublereal *z__, aocl_int_t *pp, doublereal *dmin__,
             doublereal *sigma, doublereal *desig, doublereal *qmax, aocl_int_t *nfail,
             aocl_int_t *iter, aocl_int_t *ndiv, logical *ieee, aocl_int_t *ttype,
             doublereal *dmin1, doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2,
             doublereal *g, doublereal *tau)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlasq3(i0, n0, z__, pp, dmin__, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype,
                       dmin1, dmin2, dn, dn1, dn2, g, tau);
#else
    aocl_int64_t i0_64 = *i0;
    aocl_int64_t n0_64 = *n0;
    aocl_int64_t pp_64 = *pp;
    aocl_int64_t nfail_64 = *nfail;
    aocl_int64_t iter_64 = *iter;
    aocl_int64_t ndiv_64 = *ndiv;
    aocl_int64_t ttype_64 = *ttype;

    aocl_lapack_dlasq3(&i0_64, &n0_64, z__, &pp_64, dmin__, sigma, desig, qmax, &nfail_64, &iter_64,
                       &ndiv_64, ieee, &ttype_64, dmin1, dmin2, dn, dn1, dn2, g, tau);

    *n0 = (aocl_int_t)n0_64;
    *pp = (aocl_int_t)pp_64;
    *nfail = (aocl_int_t)nfail_64;
    *iter = (aocl_int_t)iter_64;
    *ndiv = (aocl_int_t)ndiv_64;
    *ttype = (aocl_int_t)ttype_64;
#endif
}

void aocl_lapack_dlasq3(aocl_int64_t *i0, aocl_int64_t *n0, doublereal *z__, aocl_int64_t *pp,
                        doublereal *dmin__, doublereal *sigma, doublereal *desig, doublereal *qmax,
                        aocl_int64_t *nfail, aocl_int64_t *iter, aocl_int64_t *ndiv, logical *ieee,
                        aocl_int64_t *ttype, doublereal *dmin1, doublereal *dmin2, doublereal *dn,
                        doublereal *dn1, doublereal *dn2, doublereal *g, doublereal *tau)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlasq3 inputs: i0 %" FLA_IS ", n0 %" FLA_IS ", pp %" FLA_IS
                      ", nfail %" FLA_IS ", iter %" FLA_IS ", ndiv %" FLA_IS ", ttype %" FLA_IS "",
                      *i0, *n0, *pp, *nfail, *iter, *ndiv, *ttype);
    /* System generated locals */
    aocl_int64_t i__1;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    doublereal s, t;
    aocl_int64_t j4, nn;
    doublereal eps, tol;
    aocl_int64_t n0in, ipn4;
    doublereal tol2, temp;
    extern doublereal dlamch_(char *);
    extern logical disnan_(doublereal *);
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Function .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --z__;
    /* Function Body */
    n0in = *n0;
    eps = dlamch_("Precision");
    tol = eps * 100.;
    /* Computing 2nd power */
    d__1 = tol;
    tol2 = d__1 * d__1;
    /* Check for deflation. */
L10:
    if(*n0 < *i0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*n0 == *i0)
    {
        goto L20;
    }
    nn = (*n0 << 2) + *pp;
    if(*n0 == *i0 + 1)
    {
        goto L40;
    }
    /* Check whether E(N0-1) is negligible, 1 eigenvalue. */
    if(z__[nn - 5] > tol2 * (*sigma + z__[nn - 3]) && z__[nn - (*pp << 1) - 4] > tol2 * z__[nn - 7])
    {
        goto L30;
    }
L20:
    z__[(*n0 << 2) - 3] = z__[(*n0 << 2) + *pp - 3] + *sigma;
    --(*n0);
    goto L10;
    /* Check whether E(N0-2) is negligible, 2 eigenvalues. */
L30:
    if(z__[nn - 9] > tol2 * *sigma && z__[nn - (*pp << 1) - 8] > tol2 * z__[nn - 11])
    {
        goto L50;
    }
L40:
    if(z__[nn - 3] > z__[nn - 7])
    {
        s = z__[nn - 3];
        z__[nn - 3] = z__[nn - 7];
        z__[nn - 7] = s;
    }
    t = (z__[nn - 7] - z__[nn - 3] + z__[nn - 5]) * .5;
    if(z__[nn - 5] > z__[nn - 3] * tol2 && t != 0.)
    {
        s = z__[nn - 3] * (z__[nn - 5] / t);
        if(s <= t)
        {
            s = z__[nn - 3] * (z__[nn - 5] / (t * (sqrt(s / t + 1.) + 1.)));
        }
        else
        {
            s = z__[nn - 3] * (z__[nn - 5] / (t + sqrt(t) * sqrt(t + s)));
        }
        t = z__[nn - 7] + (s + z__[nn - 5]);
        z__[nn - 3] *= z__[nn - 7] / t;
        z__[nn - 7] = t;
    }
    z__[(*n0 << 2) - 7] = z__[nn - 7] + *sigma;
    z__[(*n0 << 2) - 3] = z__[nn - 3] + *sigma;
    *n0 += -2;
    goto L10;
L50:
    if(*pp == 2)
    {
        *pp = 0;
    }
    /* Reverse the qd-array, if warranted. */
    if(*dmin__ <= 0. || *n0 < n0in)
    {
        if(z__[(*i0 << 2) + *pp - 3] * 1.5 < z__[(*n0 << 2) + *pp - 3])
        {
            ipn4 = *i0 + *n0 << 2;
            i__1 = *i0 + *n0 - 1 << 1;
            for(j4 = *i0 << 2; j4 <= i__1; j4 += 4)
            {
                temp = z__[j4 - 3];
                z__[j4 - 3] = z__[ipn4 - j4 - 3];
                z__[ipn4 - j4 - 3] = temp;
                temp = z__[j4 - 2];
                z__[j4 - 2] = z__[ipn4 - j4 - 2];
                z__[ipn4 - j4 - 2] = temp;
                temp = z__[j4 - 1];
                z__[j4 - 1] = z__[ipn4 - j4 - 5];
                z__[ipn4 - j4 - 5] = temp;
                temp = z__[j4];
                z__[j4] = z__[ipn4 - j4 - 4];
                z__[ipn4 - j4 - 4] = temp;
                /* L60: */
            }
            if(*n0 - *i0 <= 4)
            {
                z__[(*n0 << 2) + *pp - 1] = z__[(*i0 << 2) + *pp - 1];
                z__[(*n0 << 2) - *pp] = z__[(*i0 << 2) - *pp];
            }
            /* Computing MIN */
            d__1 = *dmin2;
            d__2 = z__[(*n0 << 2) + *pp - 1]; // , expr subst
            *dmin2 = fla_min(d__1, d__2);
            /* Computing MIN */
            d__1 = z__[(*n0 << 2) + *pp - 1], d__2 = z__[(*i0 << 2) + *pp - 1];
            d__1 = fla_min(d__1, d__2);
            d__2 = z__[(*i0 << 2) + *pp + 3]; // ; expr subst
            z__[(*n0 << 2) + *pp - 1] = fla_min(d__1, d__2);
            /* Computing MIN */
            d__1 = z__[(*n0 << 2) - *pp], d__2 = z__[(*i0 << 2) - *pp];
            d__1 = fla_min(d__1, d__2);
            d__2 = z__[(*i0 << 2) - *pp + 4]; // ; expr subst
            z__[(*n0 << 2) - *pp] = fla_min(d__1, d__2);
            /* Computing MAX */
            d__1 = *qmax, d__2 = z__[(*i0 << 2) + *pp - 3];
            d__1 = fla_max(d__1, d__2);
            d__2 = z__[(*i0 << 2) + *pp + 1]; // ; expr subst
            *qmax = fla_max(d__1, d__2);
            *dmin__ = -0.;
        }
    }
    /* Choose a shift. */
    aocl_lapack_dlasq4(i0, n0, &z__[1], pp, &n0in, dmin__, dmin1, dmin2, dn, dn1, dn2, tau, ttype,
                       g);
    /* Call dqds until DMIN > 0. */
L70:
    aocl_lapack_dlasq5(i0, n0, &z__[1], pp, tau, sigma, dmin__, dmin1, dmin2, dn, dn1, dn2, ieee,
                       &eps);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    /* Check status. */
    if(*dmin__ >= 0. && *dmin1 >= 0.)
    {
        /* Success. */
        goto L90;
    }
    else if(*dmin__ < 0. && *dmin1 > 0. && z__[(*n0 - 1 << 2) - *pp] < tol * (*sigma + *dn1)
            && f2c_abs(*dn) < tol * *sigma)
    {
        /* Convergence hidden by negative DN. */
        z__[(*n0 - 1 << 2) - *pp + 2] = 0.;
        *dmin__ = 0.;
        goto L90;
    }
    else if(*dmin__ < 0.)
    {
        /* TAU too big. Select new TAU and try again. */
        ++(*nfail);
        if(*ttype < -22)
        {
            /* Failed twice. Play it safe. */
            *tau = 0.;
        }
        else if(*dmin1 > 0.)
        {
            /* Late failure. Gives excellent shift. */
            *tau = (*tau + *dmin__) * (1. - eps * 2.);
            *ttype += -11;
        }
        else
        {
            /* Early failure. Divide by 4. */
            *tau *= .25;
            *ttype += -12;
        }
        goto L70;
    }
    else if(disnan_(dmin__))
    {
        /* NaN. */
        if(*tau == 0.)
        {
            goto L80;
        }
        else
        {
            *tau = 0.;
            goto L70;
        }
    }
    else
    {
        /* Possible underflow. Play it safe. */
        goto L80;
    }
    /* Risk of underflow. */
L80:
    aocl_lapack_dlasq6(i0, n0, &z__[1], pp, dmin__, dmin1, dmin2, dn, dn1, dn2);
    *ndiv += *n0 - *i0 + 2;
    ++(*iter);
    *tau = 0.;
L90:
    if(*tau < *sigma)
    {
        *desig += *tau;
        t = *sigma + *desig;
        *desig -= t - *sigma;
    }
    else
    {
        t = *sigma + *tau;
        *desig = *sigma - (t - *tau) + *desig;
    }
    *sigma = t;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLASQ3 */
}
/* dlasq3_ */
