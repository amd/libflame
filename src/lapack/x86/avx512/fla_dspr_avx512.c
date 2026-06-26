/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

/**
 * @brief Symmetric rank-1 update (lower packed storage) using AVX-512.
 * @details Computes: A := alpha*x*x^T + A
 * @param[in] len Length of vector x and dimension of matrix A
 * @param[in] alpha Scalar multiplier
 * @param[in] x Input vector
 * @param[in,out] ap Symmetric matrix A in lower packed storage format
 */
void fla_dspr_lower_avx512(aocl_int64_t len, doublereal alpha, const doublereal *x, doublereal *ap)
{
    aocl_int64_t i, j;
    doublereal *ap_col = ap;

    for(j = 0; j < len; j++)
    {
        doublereal temp = alpha * x[j];

        if(temp == 0.0)
        {
            ap_col += (len - j);
            continue;
        }

        aocl_int64_t col_len = len - j;
        const doublereal *xi = &x[j];

        i = 0;

        /* Process 8 doubles per iteration */

        if(col_len >= 8)
        {
            __m512d v_temp = _mm512_set1_pd(temp);
            aocl_int64_t ilast = col_len - 8;
            for(; i <= ilast; i += 8)
            {
                __m512d v_x0 = _mm512_loadu_pd(&xi[i]);
                __m512d v_a0 = _mm512_loadu_pd(&ap_col[i]);
                v_a0 = _mm512_fmadd_pd(v_x0, v_temp, v_a0);
                _mm512_storeu_pd(&ap_col[i], v_a0);
            }
        }

        if(col_len - i > 3)
        {
            __m256d v_temp256 = _mm256_set1_pd(temp);
            __m256d v_x = _mm256_loadu_pd(&xi[i]);
            __m256d v_a = _mm256_loadu_pd(&ap_col[i]);
            v_a = _mm256_fmadd_pd(v_x, v_temp256, v_a);
            _mm256_storeu_pd(&ap_col[i], v_a);
            i += 4;
        }

        if(col_len - i > 1)
        {
            __m128d v_temp128 = _mm_set1_pd(temp);
            __m128d v_x = _mm_loadu_pd(&xi[i]);
            __m128d v_a = _mm_loadu_pd(&ap_col[i]);
            v_a = _mm_fmadd_pd(v_x, v_temp128, v_a);
            _mm_storeu_pd(&ap_col[i], v_a);
            i += 2;
        }

        if(i < col_len)
        {
            ap_col[i] += xi[i] * temp;
        }

        ap_col += col_len;
    }
}

#endif /* FLA_ENABLE_AMD_OPT */