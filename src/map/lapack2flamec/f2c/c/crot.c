/* ../netlib/crot.f -- translated by f2c (version 20100827). You must link the resulting object file
 with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors. */
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
/* > sin (S) is complex, and the vectors CX and CY are complex. */
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
void crot_(integer *n, complex *cx, integer *incx, complex *cy, integer *incy, real *c__,
           complex *s)
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
    integer i__1;
    complex q__1, q__2, q__3;
    /* Local variables */
    integer i__, ix, iy;
    complex stemp;
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
    real sr = s->r;
    real si = s->i;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        q__2.r = *c__ * cx[ix].r;
        q__2.i = *c__ * cx[ix].i; // , expr subst
        q__3.r = sr * cy[iy].r - si * cy[iy].i;
        q__3.i = sr * cy[iy].i + si * cy[iy].r; // , expr subst
        q__1.r = q__2.r + q__3.r;
        q__1.i = q__2.i + q__3.i; // , expr subst
        stemp.r = q__1.r;
        stemp.i = q__1.i; // , expr subst
        q__2.r = *c__ * cy[iy].r;
        q__2.i = *c__ * cy[iy].i; // , expr subst
        q__3.r = sr * cx[ix].r + si * cx[ix].i;
        q__3.i = sr * cx[ix].i - si * cx[ix].r; // , expr subst
        q__1.r = q__2.r - q__3.r;
        q__1.i = q__2.i - q__3.i; // , expr subst
        cy[iy].r = q__1.r;
        cy[iy].i = q__1.i; // , expr subst
        cx[ix].r = stemp.r;
        cx[ix].i = stemp.i; // , expr subst
        ix += *incx;
        iy += *incy;
        /* L10: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* Code for both increments equal to 1 */
L20:
    sr = s->r;
    si = s->i;
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        q__2.r = *c__ * cx[i__].r;
        q__2.i = *c__ * cx[i__].i; // , expr subst
        q__3.r = sr * cy[i__].r - si * cy[i__].i;
        q__3.i = sr * cy[i__].i + si * cy[i__].r; // , expr subst
        q__1.r = q__2.r + q__3.r;
        q__1.i = q__2.i + q__3.i; // , expr subst
        stemp.r = q__1.r;
        stemp.i = q__1.i; // , expr subst
        q__2.r = *c__ * cy[i__].r;
        q__2.i = *c__ * cy[i__].i; // , expr subst
        q__3.r = sr * cx[i__].r + si * cx[i__].i;
        q__3.i = sr * cx[i__].i - si * cx[i__].r; // , expr subst
        q__1.r = q__2.r - q__3.r;
        q__1.i = q__2.i - q__3.i; // , expr subst
        cy[i__].r = q__1.r;
        cy[i__].i = q__1.i; // , expr subst
        cx[i__].r = stemp.r;
        cx[i__].i = stemp.i; // , expr subst
        /* L30: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
}
/* crot_ */