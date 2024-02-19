/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *   Portions of this file consist of AI-generated content
 * *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
doublereal fla_get_max_abs_element_vector_avx2(integer m, doublereal *a, integer a_dim)
{
    __m256d vec1, vec2;
    __m256d vd4_zero = _mm256_set1_pd(-0.0f);
    __m256d vmax = {0};
    integer i__, lasti;
    doublereal max_value, temp;
    max_value = 0.0f;
    lasti = m - (m % 8);

    /* Get the maxmium elements upto size lasti = m - (m % 8)*/
    for(i__ = 1; i__ <= lasti; i__ += 8)
    {
        vec1 = _mm256_loadu_pd(&a[i__ + a_dim]);
        vec2 = _mm256_loadu_pd(&a[i__ + 4 + a_dim]);

        vec1 = _mm256_andnot_pd(vd4_zero, vec1);
        vec2 = _mm256_andnot_pd(vd4_zero, vec2);

        vec1 = _mm256_max_pd(vec1, vec2);
        vmax = _mm256_max_pd(vec1, vmax);
    }

    /* Permute and find the maximum element*/
    vec1 = _mm256_permute2f128_pd(vmax, vmax, 1);
    vmax = _mm256_max_pd(vec1, vmax);
    vec1 = _mm256_permute_pd(vmax, 5);
    vmax = _mm256_max_pd(vec1, vmax);
    vec1 = _mm256_permute2f128_pd(vmax, vmax, 1);
    vmax = _mm256_max_pd(vec1, vmax);

    max_value = vmax[0];

    /* Find maximum from reminder loop*/
    for(i__ = lasti + 1; i__ <= m; ++i__)
    {
        temp = f2c_abs(a[i__ + a_dim]);
        if(max_value < temp || temp != temp)
        {
            max_value = temp;
        }
    }

    return max_value;
}
#endif
