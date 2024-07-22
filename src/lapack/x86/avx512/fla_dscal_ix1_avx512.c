/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dscal_ix1_avx512.c
 *  @brief scales a vector by a constant
 *  */

#include "FLAME.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT

int fla_dscal_ix1_avx512(integer *n, doublereal *da, doublereal *dx, integer *incx)
{
    integer i, i__1;
    __m512d alpha512, xv512, pv512;
    __m256d alpha256, xv256, pv256;
    __m128d alpha128, xv128, pv128;

    i__1 = *n;

    /* Function Body */
    if(i__1 <= 0)
    {
        return 0;
    }

    /* Load scaling factor alpha*/
    alpha128 = _mm_set_pd1(*da);
    alpha256 = _mm256_broadcastsd_pd(alpha128);
    alpha512 = _mm512_broadcastsd_pd(alpha128);

    /* Process 8 elements at a time */
    for(i = 0; i < (i__1 - 7); i += 8)
    {
        /* Load the input values */
        xv512 = _mm512_loadu_pd((double const *)&dx[i]);

        /* perform alpha * x  */
        pv512 = _mm512_mul_pd(alpha512, xv512);

        /* Store the output */
        _mm512_storeu_pd((double *)&dx[i], pv512);
    }

    /* Process 4 elements at a time */
    for(; i < (i__1 - 3); i += 4)
    {
        /* Load the input values */
        xv256 = _mm256_loadu_pd((double const *)&dx[i]);

        /* perform alpha * x  */
        pv256 = _mm256_mul_pd(alpha256, xv256);

        /* Store the output */
        _mm256_storeu_pd((double *)&dx[i], pv256);
    }

    /* Process 2 elements at a time */
    for(; i < (i__1 - 1); i += 2)
    {
        xv128 = _mm_loadu_pd((double const *)&dx[i]);
        pv128 = _mm_mul_pd(alpha128, xv128);
        _mm_storeu_pd((double *)&dx[i], pv128);
    }

    /* last iteration */
    if(i < i__1)
    {
        xv128 = _mm_load1_pd((double const *)&dx[i]);
        pv128 = _mm_mul_pd(alpha128, xv128);
        _mm_storel_pd((double *)&dx[i], pv128);
    }

    return 0;
}

#endif