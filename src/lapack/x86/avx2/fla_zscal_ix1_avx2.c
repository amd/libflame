/******************************************************************************
 * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_zscal_ix1_avx2.c
 *  @brief ZSCAL scales a vector by a scalar constant using AVX2 intrinsics.
 *         The vector elements are assumed to be contiguosly stored in memory.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
int fla_zscal_ix1_avx2(integer *n, doublecomplex *alpha, doublecomplex *x)
{
    /* Local variables */
    integer i__1, i;
    __m256d alpha256_real, alpha256_img, ar256, ai256;
    __m128d alpha128_real, alpha128_img, ar128, ai128;
    __m256d xv256[4], yv256[4], sv256[4], pv256[4];
    __m128d xv128, yv128, sv128, pv128;

    i__1 = *n;
    if(i__1 <= 0)
    {
        return 0;
    }

    /* load scale factor in 256 bit register */
    alpha256_real = _mm256_broadcast_sd((double const *)&alpha->r);
    alpha256_img = _mm256_broadcast_sd((double const *)&alpha->i);
    ai256 = _mm256_shuffle_pd(alpha256_real, alpha256_img, 0xA);
    ar256 = _mm256_shuffle_pd(alpha256_img, alpha256_real, 0x5);

    /* load scale factor in 128 bit register */
    alpha128_real = _mm_loaddup_pd((double const *)&alpha->r);
    alpha128_img = _mm_loaddup_pd((double const *)&alpha->i);
    ai128 = _mm_shuffle_pd(alpha128_real, alpha128_img, 0x2);
    ar128 = _mm_shuffle_pd(alpha128_img, alpha128_real, 0x1);

    /* Code for increments equal to 1 only */
    for(i = 0; i < (i__1 - 7); i += 8)
    {
        /* load complex inputs */
        xv256[0] = _mm256_loadu_pd((double const *)&x[i]);
        xv256[1] = _mm256_loadu_pd((double const *)&x[i + 2]);
        xv256[2] = _mm256_loadu_pd((double const *)&x[i + 4]);
        xv256[3] = _mm256_loadu_pd((double const *)&x[i + 6]);

        /* shuffle the loaded inputs */
        yv256[0] = _mm256_movedup_pd(xv256[0]);
        yv256[1] = _mm256_movedup_pd(xv256[1]);
        yv256[2] = _mm256_movedup_pd(xv256[2]);
        yv256[3] = _mm256_movedup_pd(xv256[3]);
        sv256[0] = _mm256_unpackhi_pd(xv256[0], xv256[0]);
        sv256[1] = _mm256_unpackhi_pd(xv256[1], xv256[1]);
        sv256[2] = _mm256_unpackhi_pd(xv256[2], xv256[2]);
        sv256[3] = _mm256_unpackhi_pd(xv256[3], xv256[3]);

        /* performs the scaling */
        pv256[0] = _mm256_mul_pd(ar256, sv256[0]);
        pv256[1] = _mm256_mul_pd(ar256, sv256[1]);
        pv256[2] = _mm256_mul_pd(ar256, sv256[2]);
        pv256[3] = _mm256_mul_pd(ar256, sv256[3]);
        pv256[0] = _mm256_fmaddsub_pd(ai256, yv256[0], pv256[0]);
        pv256[1] = _mm256_fmaddsub_pd(ai256, yv256[1], pv256[1]);
        pv256[2] = _mm256_fmaddsub_pd(ai256, yv256[2], pv256[2]);
        pv256[3] = _mm256_fmaddsub_pd(ai256, yv256[3], pv256[3]);

        /* store the results */
        _mm256_storeu_pd((double *)&x[i], pv256[0]);
        _mm256_storeu_pd((double *)&x[i + 2], pv256[1]);
        _mm256_storeu_pd((double *)&x[i + 4], pv256[2]);
        _mm256_storeu_pd((double *)&x[i + 6], pv256[3]);
    }

    for(; i < (i__1 - 3); i += 4)
    {
        /* load complex inputs */
        xv256[0] = _mm256_loadu_pd((double const *)&x[i]);
        xv256[1] = _mm256_loadu_pd((double const *)&x[i + 2]);

        /* shuffle the loaded inputs */
        yv256[0] = _mm256_movedup_pd(xv256[0]);
        yv256[1] = _mm256_movedup_pd(xv256[1]);
        sv256[0] = _mm256_unpackhi_pd(xv256[0], xv256[0]);
        sv256[1] = _mm256_unpackhi_pd(xv256[1], xv256[1]);

        /* performs the scaling */
        pv256[0] = _mm256_mul_pd(ar256, sv256[0]);
        pv256[1] = _mm256_mul_pd(ar256, sv256[1]);
        pv256[0] = _mm256_fmaddsub_pd(ai256, yv256[0], pv256[0]);
        pv256[1] = _mm256_fmaddsub_pd(ai256, yv256[1], pv256[1]);

        /* store the results */
        _mm256_storeu_pd((double *)&x[i], pv256[0]);
        _mm256_storeu_pd((double *)&x[i + 2], pv256[1]);
    }

    /* remainder iterations */
    for(; i < i__1; ++i)
    {
        /* load inputs */
        xv128 = _mm_loadu_pd((double const *)&x[i]);

        /* shuffle inputs */
        yv128 = _mm_movedup_pd(xv128);
        sv128 = _mm_unpackhi_pd(xv128, xv128);

        /* performs scaling */
        pv128 = _mm_mul_pd(ar128, sv128);
        pv128 = _mm_fmaddsub_pd(ai128, yv128, pv128);

        /* store result */
        _mm_storeu_pd((double *)&x[i], pv128);
    }
    return 0;
}
#endif
