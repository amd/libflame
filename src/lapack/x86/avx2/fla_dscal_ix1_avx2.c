/******************************************************************************
 * Copyright (C) 2023-24, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dscal_ix1_avx2.c
 *  @brief scales a vector by a constant
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

int fla_dscal_ix1_avx2(integer *n, doublereal *da, doublereal *dx, integer *incx)
{
    /* Parameter adjustments */
    --dx;
    /* Function Body */
    if(*n <= 0)
    {
        return 0;
    }

    integer i, i__1;
    doublereal *d__1;
    i__1 = *n;
    d__1 = da;

    /* Load scaling factor alpha*/
    __m128d alpha128 = _mm_set_pd1(*d__1);
    __m256d alpha256 = _mm256_broadcastsd_pd(alpha128);
    /* Process 8 elements at a time */
    for(i = 1; i <= (i__1 - 7); i += 8)
    {
        /* Load the input values */
        __m256d x0v = _mm256_loadu_pd((double const *)&dx[i]);
        __m256d x1v = _mm256_loadu_pd((double const *)&dx[i + 4]);
        /* perform alpha * x  */
        x0v = _mm256_mul_pd(alpha256, x0v);
        x1v = _mm256_mul_pd(alpha256, x1v);
        /* Store the output */
        _mm256_storeu_pd((double *)&dx[i], x0v);
        _mm256_storeu_pd((double *)&dx[i + 4], x1v);
    }
    /* Remainder iterations */
    if((i__1 - i) >= 3)
    {
        __m256d x0v = _mm256_loadu_pd((double const *)&dx[i]);
        x0v = _mm256_mul_pd(alpha256, x0v);
        _mm256_storeu_pd((double *)&dx[i], x0v);
        i += 4;
    }
    /* last two iterations */
    if(i__1 > i)
    {
        __m128d x0v128 = _mm_loadu_pd((double const *)&dx[i]);
        x0v128 = _mm_mul_pd(alpha128, x0v128);
        _mm_storeu_pd((double *)&dx[i], x0v128);
        i += 2;
    }
    /* last iteration */
    if(i__1 == i)
    {
        __m128d x1v128 = _mm_load1_pd((double const *)&dx[i]);
        x1v128 = _mm_mul_pd(alpha128, x1v128);
        _mm_storel_pd((double *)&dx[i], x1v128);
    }
    return 0;
}
#endif
