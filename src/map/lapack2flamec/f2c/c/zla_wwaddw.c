/* ../netlib/zla_wwaddw.f -- translated by f2c (version 20100827). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b ZLA_WWADDW adds a vector into a doubled-single vector. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZLA_WWADDW + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_wwa
 * ddw.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_wwa
 * ddw.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_wwa
 * ddw.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZLA_WWADDW( N, X, Y, W ) */
/* .. Scalar Arguments .. */
/* INTEGER N */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 X( * ), Y( * ), W( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLA_WWADDW adds a vector W into a doubled-single vector (X, Y). */
/* > */
/* > This works for all extant IBM's hex and binary floating point */
/* > arithmetics, but not for decimal. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The length of vectors X, Y, and W. */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (N) */
/* > The first part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX*16 array, dimension (N) */
/* > The second part of the doubled-single accumulation vector. */
/* > \endverbatim */
/* > */
/* > \param[in] W */
/* > \verbatim */
/* > W is COMPLEX*16 array, dimension (N) */
/* > The vector to be added. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zla_wwaddw_(aocl_int_t *n, dcomplex *x, dcomplex *y, dcomplex *w)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zla_wwaddw(n, x, y, w);
#else
    aocl_int64_t n_64 = *n;

    aocl_lapack_zla_wwaddw(&n_64, x, y, w);
#endif
}

void aocl_lapack_zla_wwaddw(aocl_int64_t *n, dcomplex *x, dcomplex *y, dcomplex *w)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zla_wwaddw inputs: n %" FLA_IS "", *n);
    /* System generated locals */
    aocl_int64_t i__1, i__2, i__3, i__4, i__5;
    dcomplex z__1, z__2, z__3;
    /* Local variables */
    aocl_int64_t i__;
    dcomplex s;
    /* -- LAPACK computational routine (version 3.4.2) -- */
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
    --w;
    --y;
    --x;
    /* Function Body */
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        i__3 = i__;
        z__1.real = x[i__2].real + w[i__3].real;
        z__1.imag = x[i__2].imag + w[i__3].imag; // , expr subst
        s.real = z__1.real;
        s.imag = z__1.imag; // , expr subst
        z__2.real = s.real + s.real;
        z__2.imag = s.imag + s.imag; // , expr subst
        z__1.real = z__2.real - s.real;
        z__1.imag = z__2.imag - s.imag; // , expr subst
        s.real = z__1.real;
        s.imag = z__1.imag; // , expr subst
        i__2 = i__;
        i__3 = i__;
        z__3.real = x[i__3].real - s.real;
        z__3.imag = x[i__3].imag - s.imag; // , expr subst
        i__4 = i__;
        z__2.real = z__3.real + w[i__4].real;
        z__2.imag = z__3.imag + w[i__4].imag; // , expr subst
        i__5 = i__;
        z__1.real = z__2.real + y[i__5].real;
        z__1.imag = z__2.imag + y[i__5].imag; // , expr subst
        y[i__2].real = z__1.real;
        y[i__2].imag = z__1.imag; // , expr subst
        i__2 = i__;
        x[i__2].real = s.real;
        x[i__2].imag = s.imag; // , expr subst
        /* L10: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* zla_wwaddw__ */
