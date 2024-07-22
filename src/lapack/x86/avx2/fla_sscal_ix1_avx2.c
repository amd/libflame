/******************************************************************************
 * Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_sscal_ix1_avx2.c
 *  @brief SSCAL scales a vector by a constant in AVX2.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

int fla_sscal_ix1_avx2(integer *n, real *sa, real *sx)
{
    /* Local variables */
    integer i__1, i;
    __m256 xv256[2], pv256[2], alpha256;
    __m128 xv128, pv128, alpha128;

    i__1 = *n;

    if(i__1 <= 0)
    {
        return 0;
    }

    /* Load scaling factor alpha*/
    alpha128 = _mm_set_ps1(*sa);
    alpha256 = _mm256_broadcastss_ps(alpha128);

    // Process 16 elements at a time
    for(i = 0; i < (i__1 - 15); i += 16)
    {
        /* load complex inputs */
        xv256[0] = _mm256_loadu_ps((float const *)&sx[i]);
        xv256[1] = _mm256_loadu_ps((float const *)&sx[i + 8]);

        /* performs the scaling */
        pv256[0] = _mm256_mul_ps(alpha256, xv256[0]);
        pv256[1] = _mm256_mul_ps(alpha256, xv256[1]);

        /* store the results */
        _mm256_storeu_ps((float *)&sx[i], pv256[0]);
        _mm256_storeu_ps((float *)&sx[i + 8], pv256[1]);
    }

    // Process 8 elements at a time
    for(; i < (i__1 - 7); i += 8)
    {
        /* load complex inputs */
        xv256[0] = _mm256_loadu_ps((float const *)&sx[i]);

        /* performs the scaling */
        pv256[0] = _mm256_mul_ps(alpha256, xv256[0]);

        /* store the results */
        _mm256_storeu_ps((float *)&sx[i], pv256[0]);
    }

    // Process 4 elements at a time
    for(; i < (i__1 - 3); i += 4)
    {
        /* load complex inputs */
        xv128 = _mm_loadu_ps((float const *)&sx[i]);

        /* performs the scaling */
        pv128 = _mm_mul_ps(alpha128, xv128);

        /* store the results */
        _mm_storeu_ps((float *)&sx[i], pv128);
    }

    /* remainder iterations */
    for(; i < i__1; i += 1)
    {
        /* load complex inputs */
        xv128 = _mm_load_ss((float const *)&sx[i]);

        /* performs the scaling */
        pv128 = _mm_mul_ss(alpha128, xv128);

        /* store the results */
        _mm_store_ss((float *)&sx[i], pv128);
    }
    return 0;
}

#endif
