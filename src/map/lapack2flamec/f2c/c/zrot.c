/* ../netlib/zrot.f -- translated by f2c (version 20100827). You must link the resulting object file
 with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
    Copyright (c) 2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLA_f2c.h" /* > \brief \b ZROT applies a plane rotation with real cosine and scomplex sine to a pair of scomplex vectors. */
#ifdef FLA_ENABLE_AMD_OPT
#include "immintrin.h"
#endif

/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZROT + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zrot.f"
 * > */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zrot.f"
 * > */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zrot.f"
 * > */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S ) */
/* .. Scalar Arguments .. */
/* INTEGER INCX, INCY, N */
/* DOUBLE PRECISION C */
/* COMPLEX*16 S */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX*16 CX( * ), CY( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZROT applies a plane rotation, where the cos (C) is real and the */
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
/* > CX is COMPLEX*16 array, dimension (N) */
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
/* > CY is COMPLEX*16 array, dimension (N) */
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
/* > C is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] S */
/* > \verbatim */
/* > S is COMPLEX*16 */
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
/* > \ingroup complex16OTHERauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zrot_(aocl_int_t *n, dcomplex *cx, aocl_int_t *incx, dcomplex *cy, aocl_int_t *incy,
           doublereal *c__, dcomplex *s)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zrot(n, cx, incx, cy, incy, c__, s);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t incx_64 = *incx;
    aocl_int64_t incy_64 = *incy;

    aocl_lapack_zrot(&n_64, cx, &incx_64, cy, &incy_64, c__, s);
#endif
}

void aocl_lapack_zrot(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy,
                      aocl_int64_t *incy, doublereal *c__, dcomplex *s)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zrot inputs: n %" FLA_IS ", incx %" FLA_IS ", incy %" FLA_IS "", *n, *incx,
                      *incy);
    extern fla_context fla_global_context;
    extern void fla_zrot(aocl_int64_t * n, dcomplex * cx, aocl_int64_t * incx,
                         dcomplex * cy, aocl_int64_t * incy, doublereal * c__,
                         dcomplex * s);
    extern void fla_zrot_native(aocl_int64_t * n, dcomplex * cx, aocl_int64_t * incx,
                                dcomplex * cy, aocl_int64_t * incy, doublereal * c__,
                                dcomplex * s);

    /* Initialize global context data */
    aocl_fla_init();

#ifdef FLA_ENABLE_AMD_OPT
    if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_zrot(n, cx, incx, cy, incy, c__, s);
    }
    else
    {
        fla_zrot_native(n, cx, incx, cy, incy, c__, s);
    }
#else
    fla_zrot_native(n, cx, incx, cy, incy, c__, s);
#endif

    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

void fla_zrot_native(aocl_int64_t *n, dcomplex *cx, aocl_int64_t *incx, dcomplex *cy,
                     aocl_int64_t *incy, doublereal *c__, dcomplex *s)
{
    /* System generated locals */
    aocl_int64_t i__1;
    dcomplex z__1, z__2, z__3;
    /* Local variables */
    aocl_int64_t i__, ix, iy;
    doublereal lc, sr, si;
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
        return;
    }
    lc = *c__;
    sr = s->real;
    si = s->imag;

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
    if(*incx != *incy)
    {
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            z__2.real = lc * cx[ix].real;
            z__2.imag = lc * cx[ix].imag; // , expr subst
            z__3.real = sr * cy[iy].real - si * cy[iy].imag;
            z__3.imag = sr * cy[iy].imag + si * cy[iy].real; // , expr subst
            z__1.real = z__2.real + z__3.real;
            z__1.imag = z__2.imag + z__3.imag; // , expr subst

            z__2.real = lc * cy[iy].real;
            z__2.imag = lc * cy[iy].imag; // , expr subst
            z__3.real = sr * cx[ix].real + si * cx[ix].imag;
            z__3.imag = sr * cx[ix].imag - si * cx[ix].real; // , expr subst

            cy[iy].real = z__2.real - z__3.real;
            cy[iy].imag = z__2.imag - z__3.imag; // , expr subst
            cx[ix].real = z__1.real;
            cx[ix].imag = z__1.imag; // , expr subst
            ix += *incx;
            iy += *incy;
        }
    }
    else
    {
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            z__2.real = lc * cx[ix].real;
            z__2.imag = lc * cx[ix].imag; // , expr subst
            z__3.real = sr * cy[ix].real - si * cy[ix].imag;
            z__3.imag = sr * cy[ix].imag + si * cy[ix].real; // , expr subst
            z__1.real = z__2.real + z__3.real;
            z__1.imag = z__2.imag + z__3.imag; // , expr subst

            z__2.real = lc * cy[ix].real;
            z__2.imag = lc * cy[ix].imag; // , expr subst
            z__3.real = sr * cx[ix].real + si * cx[ix].imag;
            z__3.imag = sr * cx[ix].imag - si * cx[ix].real; // , expr subst

            cy[ix].real = z__2.real - z__3.real;
            cy[ix].imag = z__2.imag - z__3.imag; // , expr subst
            cx[ix].real = z__1.real;
            cx[ix].imag = z__1.imag; // , expr subst
            ix += *incx;
        }
    }
    return;
    /* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        z__2.real = lc * cx[i__].real;
        z__2.imag = lc * cx[i__].imag; // , expr subst
        z__3.real = sr * cy[i__].real - si * cy[i__].imag;
        z__3.imag = sr * cy[i__].imag + si * cy[i__].real; // , expr subst
        z__1.real = z__2.real + z__3.real;
        z__1.imag = z__2.imag + z__3.imag; // , expr subst

        z__2.real = lc * cy[i__].real;
        z__2.imag = lc * cy[i__].imag; // , expr subst
        z__3.real = sr * cx[i__].real + si * cx[i__].imag;
        z__3.imag = sr * cx[i__].imag - si * cx[i__].real; // , expr subst

        cy[i__].real = z__2.real - z__3.real;
        cy[i__].imag = z__2.imag - z__3.imag; // , expr subst
        cx[i__].real = z__1.real;
        cx[i__].imag = z__1.imag; // , expr subst
    }
    return;
}
/* zrot_ */
