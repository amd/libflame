/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dcopy_scal_avx512.c
 *  @brief y = alpha * x operation for double precision vectors
 *  */

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/**
 * @brief Computes y = alpha * x for double precision vectors using AVX512.
 * @param[in] n Length of vectors x and y
 * @param[in] da Scalar multiplier alpha
 * @param[in] dx Input vector x
 * @param[out] dy Output vector y
 * @return 0 on success
 * @note Strides of both vectors is assumed to be 1.
 */
int fla_dcopy_scal_ix1_avx512(aocl_int64_t n, doublereal da, const doublereal *dx, doublereal *dy)
{
    aocl_int64_t i;
    __m512d alpha512, xv512, pv512;
    __m256d alpha256, xv256, pv256;
    __m128d alpha128, xv128, pv128;

    /* Function Body */
    if(n <= 0)
    {
        return 0;
    }

    /* Load scaling factor alpha*/
    alpha128 = _mm_set_pd1(da);
    alpha256 = _mm256_broadcastsd_pd(alpha128);
    alpha512 = _mm512_broadcastsd_pd(alpha128);

    /* Process 8 elements at a time */
    for(i = 0; i < (n - 7); i += 8)
    {
        /* Load the input values */
        xv512 = _mm512_loadu_pd((double const *)&dx[i]);

        /* perform alpha * x  */
        pv512 = _mm512_mul_pd(alpha512, xv512);

        /* Store the output */
        _mm512_storeu_pd((double *)&dy[i], pv512);
    }

    /* Process 4 elements at a time */
    for(; i < (n - 3); i += 4)
    {
        /* Load the input values */
        xv256 = _mm256_loadu_pd((double const *)&dx[i]);

        /* perform alpha * x  */
        pv256 = _mm256_mul_pd(alpha256, xv256);

        /* Store the output */
        _mm256_storeu_pd((double *)&dy[i], pv256);
    }

    /* Process 2 elements at a time */
    for(; i < (n - 1); i += 2)
    {
        xv128 = _mm_loadu_pd((double const *)&dx[i]);
        pv128 = _mm_mul_pd(alpha128, xv128);
        _mm_storeu_pd((double *)&dy[i], pv128);
    }

    /* last iteration */
    if(i < n)
    {
        xv128 = _mm_load1_pd((double const *)&dx[i]);
        pv128 = _mm_mul_pd(alpha128, xv128);
        _mm_storel_pd((double *)&dy[i], pv128);
    }

    return 0;
}

#endif