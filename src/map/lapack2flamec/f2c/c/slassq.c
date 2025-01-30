/* ./slassq.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SLASSQ updates a sum of squares represented in scaled form. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SLASSQ + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slassq.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slassq.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slassq.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, N */
/* REAL SCALE, SUMSQ */
/* .. */
/* .. Array Arguments .. */
/* REAL X( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASSQ returns the values scale_out and sumsq_out such that */
/* > */
/* > (scale_out**2)*sumsq_out = x( 1 )**2 +...+ x( n )**2 + (scale**2)*sumsq, */
/* > */
/* > where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is */
/* > assumed to be non-negative. */
/* > */
/* > scale and sumsq must be supplied in SCALE and SUMSQ and */
/* > scale_out and sumsq_out are overwritten on SCALE and SUMSQ respectively. */
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
/* > X is REAL array, dimension (1+(N-1)*abs(INCX)) */
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
/* > SCALE is REAL */
/* > On entry, the value scale in the equation above. */
/* > On exit, SCALE is overwritten by scale_out, the scaling factor */
/* > for the sum of squares. */
/* > \endverbatim */
/* > */
/* > \param[in,out] SUMSQ */
/* > \verbatim */
/* > SUMSQ is REAL */
/* > On entry, the value sumsq in the equation above. */
/* > On exit, SUMSQ is overwritten by sumsq_out, the basic sum of */
/* > squares from which scale_out has been factored out. */
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
/* Subroutine */
void slassq_(integer *n, real *x, integer *incx, real *scale, real *sumsq)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slassq inputs: n %" FLA_IS ", incx %" FLA_IS "", *n, *incx);
    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__;
    real ax;
    integer ix;
    real  abig, amed, sbig, tbig, asml, ymin, ymax, tsml, ssml;
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
    sbig = 1.32348898E-23;
    ssml = 3.77789319E+22;
    tsml = 1.08420217E-19;
    tbig = 4.50359963E+15;
    /* .. */
    /* Quick return if possible */
    if(*scale != *scale || *sumsq != *sumsq)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(*sumsq == 0.f)
    {
        *scale = 1.f;
    }
    if(*scale == 0.f)
    {
        *scale = 1.f;
        *sumsq = 0.f;
    }
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Compute the sum of squares in 3 accumulators: */
    /* abig -- sums of squares scaled down to avoid overflow */
    /* asml -- sums of squares scaled up to avoid underflow */
    /* amed -- sums of squares that do not require scaling */
    /* The thresholds and multipliers are */
    /* tbig -- values bigger than this are scaled down by sbig */
    /* tsml -- values smaller than this are scaled up by ssml */
    notbig = TRUE_;
    asml = 0.f;
    amed = 0.f;
    abig = 0.f;
    ix = 1;
    if(*incx < 0)
    {
        ix = 1 - (*n - 1) * *incx;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        ax = (r__1 = x[ix], f2c_abs(r__1));
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
    if(*sumsq > 0.f)
    {
        ax = *scale * sqrt(*sumsq);
        if(ax > tbig)
        {
            if(*scale > 1.f)
            {
                *scale *= sbig;
                abig += *scale * (*scale * *sumsq);
            }
            else
            {
                /* sumsq > tbig^2 => (sbig * (sbig * sumsq)) is representable */
                abig += *scale * (*scale * (sbig * (sbig * *sumsq)));
            }
        }
        else if(ax < tsml)
        {
            if(notbig)
            {
                if(*scale < 1.f)
                {
                    *scale *= ssml;
                    asml += *scale * (*scale * *sumsq);
                }
                else
                {
                    /* sumsq < tsml^2 => (ssml * (ssml * sumsq)) is representa */
                    asml += *scale * (*scale * (ssml * (ssml * *sumsq)));
                }
            }
        }
        else
        {
            amed += *scale * (*scale * *sumsq);
        }
    }
    /* Combine abig and amed or amed and asml if more than one */
    /* accumulator was used. */
    if(abig > 0.f)
    {
        if(amed > 0.f || amed != amed)
        {
            abig += amed * sbig * sbig;
        }
        *scale = 1.f / sbig;
        *sumsq = abig;
    }
    else if(asml > 0.f)
    {
        /* Combine amed and asml if asml > 0. */
        if(amed > 0.f || amed != amed)
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
            *scale = 1.f;
            /* Computing 2nd power */
            r__1 = ymax;
            /* Computing 2nd power */
            r__2 = ymin / ymax;
            *sumsq = r__1 * r__1 * (1.f + r__2 * r__2);
        }
        else
        {
            *scale = 1.f / ssml;
            *sumsq = asml;
        }
    }
    else
    {
        /* Otherwise all values are mid-range or zero */
        *scale = 1.f;
        *sumsq = amed;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* slassq_ */
