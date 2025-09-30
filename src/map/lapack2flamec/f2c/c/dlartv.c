/* ../netlib/dlartv.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b DLARTV applies a vector of plane rotations with real cosines and real sines to the elements of a pair of vectors. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DLARTV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartv.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartv.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartv.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DLARTV( N, X, INCX, Y, INCY, C, S, INCC ) */
/* .. Scalar Arguments .. */
/* INTEGER INCC, INCX, INCY, N */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION C( * ), S( * ), X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARTV applies a vector of real plane rotations to elements of the */
/* > real vectors x and y. For i = 1,2,...,n */
/* > */
/* > ( x(i) ) := ( c(i) s(i) ) ( x(i) ) */
/* > ( y(i) ) ( -s(i) c(i) ) ( y(i) ) */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of plane rotations to be applied. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is DOUBLE PRECISION array, */
/* > dimension (1+(N-1)*INCX) */
/* > The vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between elements of X. INCX > 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is DOUBLE PRECISION array, */
/* > dimension (1+(N-1)*INCY) */
/* > The vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between elements of Y. INCY > 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
/* > The cosines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is DOUBLE PRECISION array, dimension (1+(N-1)*INCC) */
/* > The sines of the plane rotations. */
/* > \endverbatim */
/* > */
/* > \param[in] INCC */
/* > \verbatim */
/* > INCC is INTEGER */
/* > The increment between elements of C and S. INCC > 0. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup doubleOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dlartv_(aocl_int_t *n, doublereal *x, aocl_int_t *incx, doublereal *y, aocl_int_t *incy,
             doublereal *c__, doublereal *s, aocl_int_t *incc)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dlartv(n, x, incx, y, incy, c__, s, incc);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;
    aocl_int64_t incy_64 = *incy;
    aocl_int64_t incc_64 = *incc;

    aocl_lapack_dlartv(&n_64, x, &incx_64, y, &incy_64, c__, s, &incc_64);
#endif
}

void aocl_lapack_dlartv(aocl_int64_t *n, doublereal *x, aocl_int64_t *incx, doublereal *y,
                        aocl_int64_t *incy, doublereal *c__, doublereal *s, aocl_int64_t *incc)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlartv inputs: n %" FLA_IS ", incx %" FLA_IS ", incy %" FLA_IS
                      ", incc %" FLA_IS "",
                      *n, *incx, *incy, *incc);
    /* System generated locals */
    aocl_int64_t i__1;
    /* Local variables */
    aocl_int64_t i__, ic, ix, iy;
    doublereal xi, yi;
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Local Scalars .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    --s;
    --c__;
    --y;
    --x;
    /* Function Body */
    ix = 1;
    iy = 1;
    ic = 1;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        xi = x[ix];
        yi = y[iy];
        x[ix] = c__[ic] * xi + s[ic] * yi;
        y[iy] = c__[ic] * yi - s[ic] * xi;
        ix += *incx;
        iy += *incy;
        ic += *incc;
        /* L10: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DLARTV */
}
/* dlartv_ */
