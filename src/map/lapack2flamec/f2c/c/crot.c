/* ../netlib/crot.f -- translated by f2c (version 20100827). You must link the resulting object file
 with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CROT applies a plane rotation with real cosine and scomplex sine to a pair of scomplex vectors. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CROT + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/crot.f"
 * > */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/crot.f"
 * > */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/crot.f"
 * > */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CROT( N, CX, INCX, CY, INCY, C, S ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* REAL C */
/* COMPLEX S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX CX( * ), CY( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CROT applies a plane rotation, where the cos (C) is real and the */
/* > sin (S) is scomplex, and the vectors CX and CY are scomplex. */
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
/* > On input, the vector X. */
/* > On output, CX is overwritten with C*X + S*Y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > The increment between successive values of CY. INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* > CY is COMPLEX array, dimension (N) */
/* > On input, the vector Y. */
/* > On output, CY is overwritten with -CONJG(S)*X + C*Y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > The increment between successive values of CY. INCX <> 0. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* > C is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is COMPLEX */
/* > C and S define a rotation */
/* > [ C S ] */
/* > [ -conjg(S) C ] */
/* > where C*C + S*CONJG(S) = 1.0. */
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
void crot_(aocl_int_t *n, scomplex *cx, aocl_int_t *incx, scomplex *cy, aocl_int_t *incy, real *c__,
           scomplex *s)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_crot(n, cx, incx, cy, incy, c__, s);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;
    aocl_int64_t incy_64 = *incy;

    aocl_lapack_crot(&n_64, cx, &incx_64, cy, &incy_64, c__, s);
#endif
}

void aocl_lapack_crot(aocl_int64_t *n, scomplex *cx, aocl_int64_t *incx, scomplex *cy,
                      aocl_int64_t *incy, real *c__, scomplex *s)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "crot inputs: n %lld, incx %lld, incy %lld", *n, *incx, *incy);
#else
    snprintf(buffer, 256, "crot inputs: n %d, incx %d, incy %d", *n, *incx, *incy);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t i__1;
    scomplex q__1, q__2, q__3;
    /* Local variables */
    aocl_int64_t i__, ix, iy;
    scomplex stemp;
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
    /* .. Intrinsic Functions .. */
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
    real sr = s->real;
    real si = s->imag;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        q__2.real = *c__ * cx[ix].real;
        q__2.imag = *c__ * cx[ix].imag; // , expr subst
        q__3.real = sr * cy[iy].real - si * cy[iy].imag;
        q__3.imag = sr * cy[iy].imag + si * cy[iy].real; // , expr subst
        q__1.real = q__2.real + q__3.real;
        q__1.imag = q__2.imag + q__3.imag; // , expr subst
        stemp.real = q__1.real;
        stemp.imag = q__1.imag; // , expr subst
        q__2.real = *c__ * cy[iy].real;
        q__2.imag = *c__ * cy[iy].imag; // , expr subst
        q__3.real = sr * cx[ix].real + si * cx[ix].imag;
        q__3.imag = sr * cx[ix].imag - si * cx[ix].real; // , expr subst
        q__1.real = q__2.real - q__3.real;
        q__1.imag = q__2.imag - q__3.imag; // , expr subst
        cy[iy].real = q__1.real;
        cy[iy].imag = q__1.imag; // , expr subst
        cx[ix].real = stemp.real;
        cx[ix].imag = stemp.imag; // , expr subst
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* Code for both increments equal to 1 */
L20:
    sr = s->real;
    si = s->imag;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        q__2.real = *c__ * cx[i__].real;
        q__2.imag = *c__ * cx[i__].imag; // , expr subst
        q__3.real = sr * cy[i__].real - si * cy[i__].imag;
        q__3.imag = sr * cy[i__].imag + si * cy[i__].real; // , expr subst
        q__1.real = q__2.real + q__3.real;
        q__1.imag = q__2.imag + q__3.imag; // , expr subst
        stemp.real = q__1.real;
        stemp.imag = q__1.imag; // , expr subst
        q__2.real = *c__ * cy[i__].real;
        q__2.imag = *c__ * cy[i__].imag; // , expr subst
        q__3.real = sr * cx[i__].real + si * cx[i__].imag;
        q__3.imag = sr * cx[i__].imag - si * cx[i__].real; // , expr subst
        q__1.real = q__2.real - q__3.real;
        q__1.imag = q__2.imag - q__3.imag; // , expr subst
        cy[i__].real = q__1.real;
        cy[i__].imag = q__1.imag; // , expr subst
        cx[i__].real = stemp.real;
        cx[i__].imag = stemp.imag; // , expr subst
        /* L30: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* crot_ */