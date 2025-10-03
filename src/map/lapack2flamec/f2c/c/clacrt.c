/* ../netlib/clacrt.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CLACRT performs a linear transformation of a pair of scomplex vectors. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CLACRT + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clacrt.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clacrt.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clacrt.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CLACRT( N, CX, INCX, CY, INCY, C, S ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* COMPLEX C, S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX CX( * ), CY( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLACRT performs the operation */
/* > */
/* > ( c s )( x ) ==> ( x ) */
/* > ( -s c )( y ) ( y ) */
/* > */
/* > where c and s are scomplex and the vectors x and y are scomplex. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of elements in the vectors CX and CY. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CX */
/* > \verbatim */
/* > CX is COMPLEX array, dimension (N) */
/* > On input, the vector x. */
/* > On output, CX is overwritten with c*x + s*y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of CX. INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* > CY is COMPLEX array, dimension (N) */
/* > On input, the vector y. */
/* > On output, CY is overwritten with -s*x + c*y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between successive values of CY. INCY <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is COMPLEX */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is COMPLEX */
/* > C and S define the matrix */
/* > [ C S ]. */
/* > [ -S C ] */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexOTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void clacrt_(aocl_int_t *n, scomplex *cx, aocl_int_t *incx, scomplex *cy, aocl_int_t *incy,
             scomplex *c__, scomplex *s)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_clacrt(n, cx, incx, cy, incy, c__, s);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;
    aocl_int64_t incy_64 = *incy;

    aocl_lapack_clacrt(&n_64, cx, &incx_64, cy, &incy_64, c__, s);
#endif
}

void aocl_lapack_clacrt(aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx, scomplex *cy,
                        aocl_int64_t *incy, scomplex *c__, scomplex *s)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "clacrt inputs: n %lld, incx %lld, incy %lld", *n, *incx, *incy);
#else
    snprintf(buffer, 256, "clacrt inputs: n %d, incx %d, incy %d", *n, *incx, *incy);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t i__1, i__2, i__3, i__4;
    scomplex q__1, q__2, q__3;
    /* Local variables */
    aocl_int64_t i__, ix, iy;
    scomplex ctemp;
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
    --cy;
    --cx;
    /* Function Body */
    if(*n <= 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(*incx == 1 && *incy == 1)
    {
        goto L20;
    }
    /* Code for unequal increments or equal increments not equal to 1 */
    ix = 1;
    iy = 1;
    if(*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if(*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = ix;
        q__2.real = c__->real * cx[i__2].real - c__->imag * cx[i__2].imag;
        q__2.imag = c__->real * cx[i__2].imag + c__->imag * cx[i__2].real; // , expr subst
        i__3 = iy;
        q__3.real = s->real * cy[i__3].real - s->imag * cy[i__3].imag;
        q__3.imag = s->real * cy[i__3].imag + s->imag * cy[i__3].real; // , expr subst
        q__1.real = q__2.real + q__3.real;
        q__1.imag = q__2.imag + q__3.imag; // , expr subst
        ctemp.real = q__1.real;
        ctemp.imag = q__1.imag; // , expr subst
        i__2 = iy;
        i__3 = iy;
        q__2.real = c__->real * cy[i__3].real - c__->imag * cy[i__3].imag;
        q__2.imag = c__->real * cy[i__3].imag + c__->imag * cy[i__3].real; // , expr subst
        i__4 = ix;
        q__3.real = s->real * cx[i__4].real - s->imag * cx[i__4].imag;
        q__3.imag = s->real * cx[i__4].imag + s->imag * cx[i__4].real; // , expr subst
        q__1.real = q__2.real - q__3.real;
        q__1.imag = q__2.imag - q__3.imag; // , expr subst
        cy[i__2].real = q__1.real;
        cy[i__2].imag = q__1.imag; // , expr subst
        i__2 = ix;
        cx[i__2].real = ctemp.real;
        cx[i__2].imag = ctemp.imag; // , expr subst
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        q__2.real = c__->real * cx[i__2].real - c__->imag * cx[i__2].imag;
        q__2.imag = c__->real * cx[i__2].imag + c__->imag * cx[i__2].real; // , expr subst
        i__3 = i__;
        q__3.real = s->real * cy[i__3].real - s->imag * cy[i__3].imag;
        q__3.imag = s->real * cy[i__3].imag + s->imag * cy[i__3].real; // , expr subst
        q__1.real = q__2.real + q__3.real;
        q__1.imag = q__2.imag + q__3.imag; // , expr subst
        ctemp.real = q__1.real;
        ctemp.imag = q__1.imag; // , expr subst
        i__2 = i__;
        i__3 = i__;
        q__2.real = c__->real * cy[i__3].real - c__->imag * cy[i__3].imag;
        q__2.imag = c__->real * cy[i__3].imag + c__->imag * cy[i__3].real; // , expr subst
        i__4 = i__;
        q__3.real = s->real * cx[i__4].real - s->imag * cx[i__4].imag;
        q__3.imag = s->real * cx[i__4].imag + s->imag * cx[i__4].real; // , expr subst
        q__1.real = q__2.real - q__3.real;
        q__1.imag = q__2.imag - q__3.imag; // , expr subst
        cy[i__2].real = q__1.real;
        cy[i__2].imag = q__1.imag; // , expr subst
        i__2 = i__;
        cx[i__2].real = ctemp.real;
        cx[i__2].imag = ctemp.imag; // , expr subst
        /* L30: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* clacrt_ */
