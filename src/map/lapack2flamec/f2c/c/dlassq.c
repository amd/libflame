/* dlassq.f -- translated by f2c (version 20160102). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h"
/* > \brief \b DLASSQ updates a sum of squares represented in scaled form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLASSQ + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlassq.f90">
 */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlassq.f90">
 */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlassq.f90">
 */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* DOUBLE PRECISION    SCALE, SUMSQ */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION    X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASSQ returns the values scl and smsq such that */
/* > */
/* > ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, */
/* > */
/* > where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is */
/* > assumed to be non-negative and scl returns the value */
/* > */
/* > scl = fla_max( scale, f2c_abs( x( i ) ) ). */
/* > */
/* > scale and sumsq must be supplied in SCALE and SUMSQ and */
/* > scl and smsq are overwritten on SCALE and SUMSQ respectively. */
/* > */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of elements to be used from the vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is DOUBLE PRECISION array, dimension (1+(N-1)*f2c_dabs(INCX)) */
/* > The vector for which a scaled sum of squares is computed. */
/* > x( i ) = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of the vector x. */
/* > If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n */
/* > If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n */
/* > If INCX = 0, x isn't a vector so there is no need to call */
/* > this subroutine. If you call it anyway, it will count x(1) */
/* > in the vector norm N times. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SCALE */
/* > \verbatim */
/* > SCALE is DOUBLE PRECISION */
/* > On entry, the value scale in the equation above. */
/* > On exit, SCALE is overwritten with scl , the scaling factor */
/* > for the sum of squares. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SUMSQ */
/* > \verbatim */
/* > SUMSQ is DOUBLE PRECISION */
/* > On entry, the value sumsq in the equation above. */
/* > On exit, SUMSQ is overwritten with smsq , the basic sum of */
/* > squares from which scl has been factored out. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Edward Anderson, Lockheed Martin */
/* > \par Contributors: */
/* ================== */
/* > */
/* > Weslley Pereira, University of Colorado Denver, USA */
/* > Nick Papior, Technical University of Denmark, DK */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > Anderson E. (2017) */
/* > Algorithm 978: Safe Scaling in the Level 1 BLAS */
/* > ACM Trans Math Softw 44:1--28 */
/* > https://doi.org/10.1145/3061665 */
/* > */
/* > Blue, James L. (1978) */
/* > A Portable Fortran Program to Find the Euclidean Norm of a Vector */
/* > ACM Trans Math Softw 4:15--23 */
/* > https://doi.org/10.1145/355769.355771 */
/* > */
/* > \endverbatim */
/* > \ingroup lassq */
/* ===================================================================== */

/* la_isnan__ */
/* Subroutine */
void dlassq_(integer *n, doublereal *x, integer *incx, doublereal *scl, doublereal *sumsq)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlassq inputs: n %" FLA_IS ", incx %" FLA_IS "", *n, *incx);
    /* System generated locals */
    integer i__1;
    doublereal r__1, r__2;
    /* Builtin functions */
    double pow_ri(doublereal *, doublereal *), sqrt(doublereal);
    integer i__;
    doublereal ax;
    integer ix;
    doublereal abig, amed, sbig, tbig, asml, ymin, ssml, tsml, ymax;
    logical notbig;
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 09:17:33 8/30/21 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* .. Local Scalars .. */
    /* Parameter adjustments */
    --x;
    /* Function Body */
    tsml = 1.4916681462400413E-154;
    tbig = 1.9979190722022350E+146;
    ssml = 4.4989137945431964E+161;
    sbig = 1.1113793747425387E-162;

    /* .. */
    /* Quick return if possible */
    if(*scl != *scl || *sumsq != *sumsq)
    {
        return;
    }
    if(*sumsq == 0.)
    {
        *scl = 1.;
    }
    if(*scl == 0.)
    {
        *scl = 1.;
        *sumsq = 0.;
    }
    if(*n <= 0)
    {
        return;
    }
    /* Compute the sum of squares in 3 accumulators: */
    /* abig -- sums of squares scaled down to avoid overflow */
    /* asml -- sums of squares scaled up to avoid underflow */
    /* amed -- sums of squares that do not require scaling */
    /* The thresholds and multipliers are */
    /* dtbig -- values bigger than this are scaled down by dsbig */
    /* dtsml -- values smaller than this are scaled up by dssml */
    notbig = TRUE_;
    asml = 0.;
    amed = 0.;
    abig = 0.;
    ix = 1;
    if(*incx < 0)
    {
        ix = 1 - (*n - 1) * *incx;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        ax = f2c_dabs(x[ix]);
        if(ax > tbig)
        {
            /* Computing 2nd power */
            r__1 = ax * sbig;
            abig += r__1 * r__1;
            notbig = FALSE_;
        }
        else if(ax < tsml)
        {
            if(notbig)
            {
                /* Computing 2nd power */
                r__1 = ax * ssml;
                asml += r__1 * r__1;
            }
        }
        else
        {
            /* Computing 2nd power */
            r__1 = ax;
            amed += r__1 * r__1;
        }
        ix += *incx;
    }
    /* Put the existing sum of squares into one of the accumulators */
    if(*sumsq > 0.)
    {
        ax = *scl * sqrt(*sumsq);
        if(ax > tbig)
        {
            if(*scl > 1.0)
            {
                *scl *= sbig;
                abig += (*scl) * (*scl) * (*sumsq);
            }
            else
            {
                abig += (*scl) * (*scl) * (sbig * (sbig * (*sumsq)));
            }
        }
        else if(ax < tsml)
        {
            if(notbig)
            {
                if(*scl < 1.0)
                {
                    *scl *= ssml;
                    asml += (*scl) * (*scl) * (*sumsq);
                }
                else
                {
                    asml += (*scl) * (*scl) * (ssml * (ssml * (*sumsq)));
                }
            }
        }
        else
        {
            /* Computing 2nd power */
            r__1 = *scl;
            amed += (r__1 * r__1) * *sumsq;
        }
    }
    /* Combine abig and amed or amed and asml if more than one */
    /* accumulator was used. */
    if(abig > 0.)
    {
        if(amed > 0. || amed != amed)
        {
            abig += amed * sbig * sbig;
        }
        *scl = 1. / sbig;
        *sumsq = abig;
    }
    else if(asml > 0.)
    {
        /* Combine amed and asml if asml > 0. */
        if(amed > 0. || amed != amed)
        {
            amed = sqrt(amed);
            asml = sqrt(asml) / ssml;
            if(asml > amed)
            {
                ymin = amed;
                ymax = asml;
            }
            else
            {
                ymin = asml;
                ymax = amed;
            }
            *scl = 1.;
            /* Computing 2nd power */
            r__1 = ymax;
            /* Computing 2nd power */
            r__2 = ymin / ymax;
            *sumsq = r__1 * r__1 * (r__2 * r__2 + 1.);
        }
        else
        {
            *scl = 1. / ssml;
            *sumsq = asml;
        }
    }
    else
    {
        /* Otherwise all values are mid-range or zero */
        *scl = 1.;
        *sumsq = amed;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* dlassq_ */
