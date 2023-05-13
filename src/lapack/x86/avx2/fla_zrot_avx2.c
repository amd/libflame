/******************************************************************************
 * * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_zrot_avx2.c
 *  *  @brief Plane rotations in AVX2.
 *   *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT

/* Application of 2x2 Plane Rotation on two vectors */
int fla_zrot_avx2(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s)
{
    /* System generated locals */
    integer i__1;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    integer i__, ix, iy;
    doublereal lc, sr, si, msi;

    __m256d cmm, srmm, simm, sinm;
    __m256d sirmm, srimm, msirmm, msrimm;
    __m256d xmm0, ymm0, xmm1, ymm1;
    __m256d xrmm0, yrmm0, ximm0, yimm0;
    __m256d xrmm1, yrmm1, ximm1, yimm1;
    __m256d oxm0, oym0, oxm1, oym1;

    __m128d cm, srm, sim, sin;
    __m128d sirm, srim, msirm, msrim;
    __m128d xmm, ymm;
    __m128d xrmm, yrmm, ximm, yimm;
    __m128d oxm, oym;
    __m128d hxmm0, hxmm1, hymm0, hymm1;
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

    if (*n <= 0)
    {
        return 0;
    }
    lc  = *c__;
    sr  = s->r;
    si  = s->i;
    msi = -si;

    cmm  = _mm256_broadcast_sd((double const *) &lc);
    srmm = _mm256_broadcast_sd((double const *) &sr);
    simm = _mm256_broadcast_sd((double const *) &si);
    sinm = _mm256_broadcast_sd((double const *) &msi);

    sirmm = _mm256_shuffle_pd(srmm, simm, 0xA);  
    srimm = _mm256_shuffle_pd(simm, srmm, 0x5);  
    msirmm = _mm256_shuffle_pd(srmm, sinm, 0xA); 
    msrimm = _mm256_shuffle_pd(sinm, srmm, 0x5); 

    cm  = _mm_loaddup_pd ((double const *) &lc);
    srm = _mm_loaddup_pd ((double const *) &sr);
    sim = _mm_loaddup_pd ((double const *) &si);
    sin = _mm_loaddup_pd ((double const *) &msi);

    sirm = _mm_shuffle_pd(srm, sim, 0x2); 
    srim = _mm_shuffle_pd(sim, srm, 0x1);
    msirm = _mm_shuffle_pd(srm, sin, 0x2);
    msrim = _mm_shuffle_pd(sin, srm, 0x1);


    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }
    /* Code for unequal increments or equal increments not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }

    i__1 = *n;
    if (*incx != *incy)
    {
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            z__2.r = lc * cx[ix].r;
            z__2.i = lc * cx[ix].i; // , expr subst
            z__3.r = sr * cy[iy].r - si * cy[iy].i;
            z__3.i = sr * cy[iy].i + si * cy[iy].r; // , expr subst
            z__1.r = z__2.r + z__3.r;
            z__1.i = z__2.i + z__3.i; // , expr subst

            z__2.r = lc * cy[iy].r;
            z__2.i = lc * cy[iy].i; // , expr subst
            z__3.r = sr * cx[ix].r + si * cx[ix].i;
            z__3.i = sr * cx[ix].i - si * cx[ix].r; // , expr subst

            cy[iy].r = z__2.r - z__3.r;
            cy[iy].i = z__2.i - z__3.i; // , expr subst
            cx[ix].r = z__1.r;
            cx[ix].i = z__1.i; // , expr subst
            ix += *incx;
            iy += *incy;
        }
    }
    else
    {
        for (i__ = 1; i__ <= (i__1 - 1); i__ += 2)
        {
            /* load complex inputs from x & y */
            xmm0   = _mm256_loadu_pd((double const *) &cx[ix]);
            hxmm1 = _mm_loadu_pd((double const *) &cx[ix + *incx]);
            ymm0   = _mm256_loadu_pd((double const *) &cy[ix]);
            hymm1 = _mm_loadu_pd((double const *) &cy[ix + *incx]);

            /* pack the inputs into 256-bit registers */
            xmm0 = _mm256_insertf128_pd(xmm0, hxmm1, 0x1);
            ymm0 = _mm256_insertf128_pd(ymm0, hymm1, 0x1);

            /* shuffle the loaded inputs */
            xrmm0 = _mm256_movedup_pd(xmm0);
            ximm0 = _mm256_unpackhi_pd(xmm0, xmm0);
            yrmm0 = _mm256_movedup_pd(ymm0);
            yimm0 = _mm256_unpackhi_pd(ymm0, ymm0);

            /* compute x outputs */
            oxm0 = _mm256_mul_pd(srimm, yimm0);
            oxm0 = _mm256_fmaddsub_pd(sirmm, yrmm0, oxm0);
            oxm0 = _mm256_fmadd_pd(cmm, xmm0, oxm0);

            /* compute y outputs */
            oym0 = _mm256_mul_pd(msrimm, ximm0);
            oym0 = _mm256_fmaddsub_pd(msirmm, xrmm0, oym0);
            oym0 = _mm256_fmsub_pd(cmm, ymm0, oym0);

            /* extract the results */
            hxmm0 = _mm256_extractf128_pd(oxm0, 0x0);
            hxmm1 = _mm256_extractf128_pd(oxm0, 0x1);
            hymm0 = _mm256_extractf128_pd(oym0, 0x0);
            hymm1 = _mm256_extractf128_pd(oym0, 0x1);

            /* store the results */
            _mm_storeu_pd((double *) &cx[ix], hxmm0);
            _mm_storeu_pd((double *) &cx[ix + *incx], hxmm1);
            _mm_storeu_pd((double *) &cy[ix], hymm0);
            _mm_storeu_pd((double *) &cy[ix + *incx], hymm1);

            ix += 2 * *incx;
        }
        for ( ; i__ <= i__1; ++i__)
        {
            /* load complex inputs from x & y */
            xmm  = _mm_loadu_pd((double const *) &cx[ix]);
            ymm  = _mm_loadu_pd((double const *) &cy[ix]);

            /* shuffle the loaded inputs */
            xrmm = _mm_movedup_pd(xmm);
            ximm = _mm_unpackhi_pd(xmm, xmm);
            yrmm = _mm_movedup_pd(ymm);
            yimm = _mm_unpackhi_pd(ymm, ymm);   

            /* compute x outputs */
            oxm = _mm_mul_pd(srim, yimm);
            oxm = _mm_fmaddsub_pd(sirm, yrmm, oxm);
            oxm = _mm_fmadd_pd(cm, xmm, oxm);

            /* compute y outputs */
            oym = _mm_mul_pd(msrim, ximm);
            oym = _mm_fmaddsub_pd(msirm, xrmm, oym);
            oym = _mm_fmsub_pd(cm, ymm, oym);

            /* store the results */
            _mm_storeu_pd((double *) &cx[ix], oxm);
            _mm_storeu_pd((double *) &cy[ix], oym);

            ix += *incx;
        }
    }
    return 0;
    /* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= (i__1 - 3); i__ += 4)
    {
        /* load complex inputs from x & y */
        xmm0 = _mm256_loadu_pd((double const *) &cx[i__]);
        ymm0 = _mm256_loadu_pd((double const *) &cy[i__]);
        xmm1 = _mm256_loadu_pd((double const *) &cx[i__ + 2]);
        ymm1 = _mm256_loadu_pd((double const *) &cy[i__ + 2]);

        /* shuffle the loaded inputs */
        xrmm0 = _mm256_movedup_pd(xmm0);
        ximm0 = _mm256_unpackhi_pd(xmm0, xmm0);
        yrmm0 = _mm256_movedup_pd(ymm0);
        yimm0 = _mm256_unpackhi_pd(ymm0, ymm0);

        xrmm1 = _mm256_movedup_pd(xmm1);
        ximm1 = _mm256_unpackhi_pd(xmm1, xmm1);
        yrmm1 = _mm256_movedup_pd(ymm1);
        yimm1 = _mm256_unpackhi_pd(ymm1, ymm1);

        /* compute x outputs */
        oxm0 = _mm256_mul_pd(srimm, yimm0);
        oxm0 = _mm256_fmaddsub_pd(sirmm, yrmm0, oxm0);
        oxm0 = _mm256_fmadd_pd(cmm, xmm0, oxm0);

        oxm1 = _mm256_mul_pd(srimm, yimm1);
        oxm1 = _mm256_fmaddsub_pd(sirmm, yrmm1, oxm1);
        oxm1 = _mm256_fmadd_pd(cmm, xmm1, oxm1);

        /* compute y outputs */
        oym0 = _mm256_mul_pd(msrimm, ximm0);
        oym0 = _mm256_fmaddsub_pd(msirmm, xrmm0, oym0);
        oym0 = _mm256_fmsub_pd(cmm, ymm0, oym0);

        oym1 = _mm256_mul_pd(msrimm, ximm1);
        oym1 = _mm256_fmaddsub_pd(msirmm, xrmm1, oym1);
        oym1 = _mm256_fmsub_pd(cmm, ymm1, oym1);

        /* store the results */
        _mm256_storeu_pd((double *) &cx[i__], oxm0);
        _mm256_storeu_pd((double *) &cy[i__], oym0);
        _mm256_storeu_pd((double *) &cx[i__ + 2], oxm1);
        _mm256_storeu_pd((double *) &cy[i__ + 2], oym1);
    }

    for ( ; i__ <= i__1; ++i__)
    {
        /* load complex inputs from x & y */
        xmm  = _mm_loadu_pd((double const *) &cx[i__]);
        ymm  = _mm_loadu_pd((double const *) &cy[i__]);

        /* shuffle the loaded inputs */
        xrmm = _mm_movedup_pd(xmm);
        ximm = _mm_unpackhi_pd(xmm, xmm);
        yrmm = _mm_movedup_pd(ymm);
        yimm = _mm_unpackhi_pd(ymm, ymm);   

        /* compute x outputs */
        oxm = _mm_mul_pd(srim, yimm);
        oxm = _mm_fmaddsub_pd(sirm, yrmm, oxm);
        oxm = _mm_fmadd_pd(cm, xmm, oxm);

        /* compute y outputs */
        oym = _mm_mul_pd(msrim, ximm);
        oym = _mm_fmaddsub_pd(msirm, xrmm, oym);
        oym = _mm_fmsub_pd(cm, ymm, oym);

        /* store the results */
        _mm_storeu_pd((double *) &cx[i__], oxm);
        _mm_storeu_pd((double *) &cy[i__], oym);
    }

    return 0;
}
#endif