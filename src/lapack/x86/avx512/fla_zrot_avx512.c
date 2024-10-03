/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file fla_zrot_avx512.c
 *  *  @brief Plane rotations in AVX512.
 *   *  */

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Application of 2x2 Plane Rotation on two vectors */
int fla_zrot_avx512(integer *n, doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy,
                    doublereal *c__, doublecomplex *s)
{
    /* System generated locals */
    integer i__1;
    doublecomplex z__1, z__2, z__3;
    /* Local variables */
    integer i__, ix, iy;
    integer aix, aiy;
    doublereal lc, sr, si, msi;

    __m512d vd8_cmm, vd8_srmm, vd8_simm, vd8_sinm;
    __m512d vd8_sirmm, vd8_srimm, vd8_msirmm, vd8_msrimm;
    __m512d vd8_xmm0, vd8_ymm0, vd8_xmm1, vd8_ymm1;
    __m512d vd8_xrmm0, vd8_yrmm0, vd8_ximm0, vd8_yimm0;
    __m512d vd8_xrmm1, vd8_yrmm1, vd8_ximm1, vd8_yimm1;
    __m512d vd8_oxm0, vd8_oym0, vd8_oxm1, vd8_oym1;

    __m256d vd4_cmm, vd4_srmm, vd4_simm, vd4_sinm;
    __m256d vd4_sirmm, vd4_srimm, vd4_msirmm, vd4_msrimm;
    __m256d vd4_xmm0, vd4_ymm0;
    __m256d vd4_xrmm0, vd4_yrmm0, vd4_ximm0, vd4_yimm0;
    __m256d vd4_oxm0, vd4_oym0;

    __m128d vd2_cm, vd2_srm, vd2_sim, vd2_sin;
    __m128d vd2_sirm, vd2_srim, vd2_msirm, vd2_msrim;
    __m128d vd2_xmm, vd2_ymm;
    __m128d vd2_xrmm, vd2_yrmm, vd2_ximm, vd2_yimm;
    __m128d vd2_oxm, vd2_oym;
    __m128d vd2_hxmm0, vd2_hxmm1, vd2_hymm0, vd2_hymm1;
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
        return 0;
    }
    lc = *c__;
    sr = s->r;
    si = s->i;
    msi = -si;

    vd4_cmm = _mm256_broadcast_sd((double const *)&lc);
    vd4_srmm = _mm256_broadcast_sd((double const *)&sr);
    vd4_simm = _mm256_broadcast_sd((double const *)&si);
    vd4_sinm = _mm256_broadcast_sd((double const *)&msi);

    vd4_sirmm = _mm256_shuffle_pd(vd4_srmm, vd4_simm, 0xA);
    vd4_srimm = _mm256_shuffle_pd(vd4_simm, vd4_srmm, 0x5);
    vd4_msirmm = _mm256_shuffle_pd(vd4_srmm, vd4_sinm, 0xA);
    vd4_msrimm = _mm256_shuffle_pd(vd4_sinm, vd4_srmm, 0x5);

    vd2_cm = _mm_loaddup_pd((double const *)&lc);
    vd2_srm = _mm_loaddup_pd((double const *)&sr);
    vd2_sim = _mm_loaddup_pd((double const *)&si);
    vd2_sin = _mm_loaddup_pd((double const *)&msi);

    vd2_sirm = _mm_shuffle_pd(vd2_srm, vd2_sim, 0x2);
    vd2_srim = _mm_shuffle_pd(vd2_sim, vd2_srm, 0x1);
    vd2_msirm = _mm_shuffle_pd(vd2_srm, vd2_sin, 0x2);
    vd2_msrim = _mm_shuffle_pd(vd2_sin, vd2_srm, 0x1);

    vd8_cmm = _mm512_broadcastsd_pd(vd2_cm);
    vd8_srmm = _mm512_broadcastsd_pd(vd2_srm);
    vd8_simm = _mm512_broadcastsd_pd(vd2_sim);
    vd8_sinm = _mm512_broadcastsd_pd(vd2_sin);

    vd8_sirmm = _mm512_shuffle_pd(vd8_srmm, vd8_simm, 0xA);
    vd8_srimm = _mm512_shuffle_pd(vd8_simm, vd8_srmm, 0x5);
    vd8_msirmm = _mm512_shuffle_pd(vd8_srmm, vd8_sinm, 0xA);
    vd8_msrimm = _mm512_shuffle_pd(vd8_sinm, vd8_srmm, 0x5);

    aix = f2c_abs(*incx);
    aiy = f2c_abs(*incy);

    if(aix == 1 && aiy == 1)
    {
        goto L20;
    }

    /* Code for unequal increments or equal increments not equal to 1 */
    ix = 1;
    iy = 1;

    i__1 = *n;
    if(*incx != *incy)
    {
        for(i__ = 1; i__ <= i__1; ++i__)
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
            ix += aix;
            iy += aiy;
        }
    }
    else
    {
        for(i__ = 1; i__ <= (i__1 - 1); i__ += 2)
        {
            /* load complex inputs from x & y */
            vd4_xmm0 = _mm256_loadu_pd((double const *)&cx[ix]);
            vd2_hxmm1 = _mm_loadu_pd((double const *)&cx[ix + aix]);
            vd4_ymm0 = _mm256_loadu_pd((double const *)&cy[ix]);
            vd2_hymm1 = _mm_loadu_pd((double const *)&cy[ix + aix]);

            /* pack the inputs into 256-bit registers */
            vd4_xmm0 = _mm256_insertf128_pd(vd4_xmm0, vd2_hxmm1, 0x1);
            vd4_ymm0 = _mm256_insertf128_pd(vd4_ymm0, vd2_hymm1, 0x1);

            /* shuffle the loaded inputs */
            vd4_xrmm0 = _mm256_movedup_pd(vd4_xmm0);
            vd4_ximm0 = _mm256_unpackhi_pd(vd4_xmm0, vd4_xmm0);
            vd4_yrmm0 = _mm256_movedup_pd(vd4_ymm0);
            vd4_yimm0 = _mm256_unpackhi_pd(vd4_ymm0, vd4_ymm0);

            /* compute x outputs */
            vd4_oxm0 = _mm256_mul_pd(vd4_srimm, vd4_yimm0);
            vd4_oxm0 = _mm256_fmaddsub_pd(vd4_sirmm, vd4_yrmm0, vd4_oxm0);
            vd4_oxm0 = _mm256_fmadd_pd(vd4_cmm, vd4_xmm0, vd4_oxm0);

            /* compute y outputs */
            vd4_oym0 = _mm256_mul_pd(vd4_msrimm, vd4_ximm0);
            vd4_oym0 = _mm256_fmaddsub_pd(vd4_msirmm, vd4_xrmm0, vd4_oym0);
            vd4_oym0 = _mm256_fmsub_pd(vd4_cmm, vd4_ymm0, vd4_oym0);

            /* extract the results */
            vd2_hxmm0 = _mm256_extractf128_pd(vd4_oxm0, 0x0);
            vd2_hxmm1 = _mm256_extractf128_pd(vd4_oxm0, 0x1);
            vd2_hymm0 = _mm256_extractf128_pd(vd4_oym0, 0x0);
            vd2_hymm1 = _mm256_extractf128_pd(vd4_oym0, 0x1);

            /* store the results */
            _mm_storeu_pd((double *)&cx[ix], vd2_hxmm0);
            _mm_storeu_pd((double *)&cx[ix + aix], vd2_hxmm1);
            _mm_storeu_pd((double *)&cy[ix], vd2_hymm0);
            _mm_storeu_pd((double *)&cy[ix + aix], vd2_hymm1);

            ix += 2 * aix;
        }
        for(; i__ <= i__1; ++i__)
        {
            /* load complex inputs from x & y */
            vd2_xmm = _mm_loadu_pd((double const *)&cx[ix]);
            vd2_ymm = _mm_loadu_pd((double const *)&cy[ix]);

            /* shuffle the loaded inputs */
            vd2_xrmm = _mm_movedup_pd(vd2_xmm);
            vd2_ximm = _mm_unpackhi_pd(vd2_xmm, vd2_xmm);
            vd2_yrmm = _mm_movedup_pd(vd2_ymm);
            vd2_yimm = _mm_unpackhi_pd(vd2_ymm, vd2_ymm);

            /* compute x outputs */
            vd2_oxm = _mm_mul_pd(vd2_srim, vd2_yimm);
            vd2_oxm = _mm_fmaddsub_pd(vd2_sirm, vd2_yrmm, vd2_oxm);
            vd2_oxm = _mm_fmadd_pd(vd2_cm, vd2_xmm, vd2_oxm);

            /* compute y outputs */
            vd2_oym = _mm_mul_pd(vd2_msrim, vd2_ximm);
            vd2_oym = _mm_fmaddsub_pd(vd2_msirm, vd2_xrmm, vd2_oym);
            vd2_oym = _mm_fmsub_pd(vd2_cm, vd2_ymm, vd2_oym);

            /* store the results */
            _mm_storeu_pd((double *)&cx[ix], vd2_oxm);
            _mm_storeu_pd((double *)&cy[ix], vd2_oym);

            ix += aix;
        }
    }
    return 0;
    /* Code for both increments equal to 1 */
L20:
    i__1 = *n;
    for(i__ = 1; i__ <= (i__1 - 7); i__ += 8)
    {
        /* load complex inputs from x & y */
        vd8_xmm0 = _mm512_loadu_pd((double const *)&cx[i__]);
        vd8_ymm0 = _mm512_loadu_pd((double const *)&cy[i__]);
        vd8_xmm1 = _mm512_loadu_pd((double const *)&cx[i__ + 4]);
        vd8_ymm1 = _mm512_loadu_pd((double const *)&cy[i__ + 4]);

        /* shuffle the loaded inputs */
        vd8_xrmm0 = _mm512_movedup_pd(vd8_xmm0);
        vd8_ximm0 = _mm512_unpackhi_pd(vd8_xmm0, vd8_xmm0);
        vd8_yrmm0 = _mm512_movedup_pd(vd8_ymm0);
        vd8_yimm0 = _mm512_unpackhi_pd(vd8_ymm0, vd8_ymm0);

        vd8_xrmm1 = _mm512_movedup_pd(vd8_xmm1);
        vd8_ximm1 = _mm512_unpackhi_pd(vd8_xmm1, vd8_xmm1);
        vd8_yrmm1 = _mm512_movedup_pd(vd8_ymm1);
        vd8_yimm1 = _mm512_unpackhi_pd(vd8_ymm1, vd8_ymm1);

        /* compute x outputs */
        vd8_oxm0 = _mm512_mul_pd(vd8_srimm, vd8_yimm0);
        vd8_oxm0 = _mm512_fmaddsub_pd(vd8_sirmm, vd8_yrmm0, vd8_oxm0);
        vd8_oxm0 = _mm512_fmadd_pd(vd8_cmm, vd8_xmm0, vd8_oxm0);

        vd8_oxm1 = _mm512_mul_pd(vd8_srimm, vd8_yimm1);
        vd8_oxm1 = _mm512_fmaddsub_pd(vd8_sirmm, vd8_yrmm1, vd8_oxm1);
        vd8_oxm1 = _mm512_fmadd_pd(vd8_cmm, vd8_xmm1, vd8_oxm1);

        /* compute y outputs */
        vd8_oym0 = _mm512_mul_pd(vd8_msrimm, vd8_ximm0);
        vd8_oym0 = _mm512_fmaddsub_pd(vd8_msirmm, vd8_xrmm0, vd8_oym0);
        vd8_oym0 = _mm512_fmsub_pd(vd8_cmm, vd8_ymm0, vd8_oym0);

        vd8_oym1 = _mm512_mul_pd(vd8_msrimm, vd8_ximm1);
        vd8_oym1 = _mm512_fmaddsub_pd(vd8_msirmm, vd8_xrmm1, vd8_oym1);
        vd8_oym1 = _mm512_fmsub_pd(vd8_cmm, vd8_ymm1, vd8_oym1);

        /* store the results */
        _mm512_storeu_pd((double *)&cx[i__], vd8_oxm0);
        _mm512_storeu_pd((double *)&cy[i__], vd8_oym0);
        _mm512_storeu_pd((double *)&cx[i__ + 4], vd8_oxm1);
        _mm512_storeu_pd((double *)&cy[i__ + 4], vd8_oym1);
    }
    for(; i__ <= (i__1 - 3); i__ += 4)
    {
        /* load complex inputs from x & y */
        vd8_xmm0 = _mm512_loadu_pd((double const *)&cx[i__]);
        vd8_ymm0 = _mm512_loadu_pd((double const *)&cy[i__]);

        /* shuffle the loaded inputs */
        vd8_xrmm0 = _mm512_movedup_pd(vd8_xmm0);
        vd8_ximm0 = _mm512_unpackhi_pd(vd8_xmm0, vd8_xmm0);
        vd8_yrmm0 = _mm512_movedup_pd(vd8_ymm0);
        vd8_yimm0 = _mm512_unpackhi_pd(vd8_ymm0, vd8_ymm0);

        /* compute x outputs */
        vd8_oxm0 = _mm512_mul_pd(vd8_srimm, vd8_yimm0);
        vd8_oxm0 = _mm512_fmaddsub_pd(vd8_sirmm, vd8_yrmm0, vd8_oxm0);
        vd8_oxm0 = _mm512_fmadd_pd(vd8_cmm, vd8_xmm0, vd8_oxm0);

        /* compute y outputs */
        vd8_oym0 = _mm512_mul_pd(vd8_msrimm, vd8_ximm0);
        vd8_oym0 = _mm512_fmaddsub_pd(vd8_msirmm, vd8_xrmm0, vd8_oym0);
        vd8_oym0 = _mm512_fmsub_pd(vd8_cmm, vd8_ymm0, vd8_oym0);

        /* store the results */
        _mm512_storeu_pd((double *)&cx[i__], vd8_oxm0);
        _mm512_storeu_pd((double *)&cy[i__], vd8_oym0);
    }

    for(; i__ <= (i__1 - 1); i__ += 2)
    {
        /* load complex inputs from x & y */
        vd4_xmm0 = _mm256_loadu_pd((double const *)&cx[i__]);
        vd4_ymm0 = _mm256_loadu_pd((double const *)&cy[i__]);

        /* shuffle the loaded inputs */
        vd4_xrmm0 = _mm256_movedup_pd(vd4_xmm0);
        vd4_ximm0 = _mm256_unpackhi_pd(vd4_xmm0, vd4_xmm0);
        vd4_yrmm0 = _mm256_movedup_pd(vd4_ymm0);
        vd4_yimm0 = _mm256_unpackhi_pd(vd4_ymm0, vd4_ymm0);

        /* compute x outputs */
        vd4_oxm0 = _mm256_mul_pd(vd4_srimm, vd4_yimm0);
        vd4_oxm0 = _mm256_fmaddsub_pd(vd4_sirmm, vd4_yrmm0, vd4_oxm0);
        vd4_oxm0 = _mm256_fmadd_pd(vd4_cmm, vd4_xmm0, vd4_oxm0);

        /* compute y outputs */
        vd4_oym0 = _mm256_mul_pd(vd4_msrimm, vd4_ximm0);
        vd4_oym0 = _mm256_fmaddsub_pd(vd4_msirmm, vd4_xrmm0, vd4_oym0);
        vd4_oym0 = _mm256_fmsub_pd(vd4_cmm, vd4_ymm0, vd4_oym0);

        /* store the results */
        _mm256_storeu_pd((double *)&cx[i__], vd4_oxm0);
        _mm256_storeu_pd((double *)&cy[i__], vd4_oym0);
    }

    for(; i__ <= i__1; ++i__)
    {
        /* load complex inputs from x & y */
        vd2_xmm = _mm_loadu_pd((double const *)&cx[i__]);
        vd2_ymm = _mm_loadu_pd((double const *)&cy[i__]);

        /* shuffle the loaded inputs */
        vd2_xrmm = _mm_movedup_pd(vd2_xmm);
        vd2_ximm = _mm_unpackhi_pd(vd2_xmm, vd2_xmm);
        vd2_yrmm = _mm_movedup_pd(vd2_ymm);
        vd2_yimm = _mm_unpackhi_pd(vd2_ymm, vd2_ymm);

        /* compute x outputs */
        vd2_oxm = _mm_mul_pd(vd2_srim, vd2_yimm);
        vd2_oxm = _mm_fmaddsub_pd(vd2_sirm, vd2_yrmm, vd2_oxm);
        vd2_oxm = _mm_fmadd_pd(vd2_cm, vd2_xmm, vd2_oxm);

        /* compute y outputs */
        vd2_oym = _mm_mul_pd(vd2_msrim, vd2_ximm);
        vd2_oym = _mm_fmaddsub_pd(vd2_msirm, vd2_xrmm, vd2_oym);
        vd2_oym = _mm_fmsub_pd(vd2_cm, vd2_ymm, vd2_oym);

        /* store the results */
        _mm_storeu_pd((double *)&cx[i__], vd2_oxm);
        _mm_storeu_pd((double *)&cy[i__], vd2_oym);
    }

    return 0;
}
#endif
