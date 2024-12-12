/******************************************************************************
 * Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dnrm2blas_avx2.c
 *  @brief norm2 for AVX2 kernel.
 *  */

#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
    /* dnrm2 for small sizes */
#define CMP256_df(v, t, T) \
    _mm256_or_pd(_mm256_cmp_pd(v, t, _CMP_LE_OS), _mm256_cmp_pd(v, T, _CMP_GE_OS));
// One 256-bit AVX register holds 4 DP elements.
typedef union
{
    __m256d v;
    double d[4] __attribute__((aligned(64)));
} v4df_t;

static inline int bli_horizontal_or_df(__m256d a)
{
    return !_mm256_testz_pd(a, a);
}

doublereal fla_dnrm2_blas_avx2(integer *sd, doublereal *a, integer *incx)
{
    double sumsq = 0;

    double *xt = a;

    // Compute the sum of squares on 3 accumulators to avoid overflow
    // and underflow, depending on the vector element value.
    // Accumulator for small values; using scaling to avoid underflow.
    double sum_sml = 0;
    // Accumulator for medium values; no scaling required.
    double sum_med = 0;
    // Accumulator for big values; using scaling to avoid overflow.
    double sum_big = 0;
    // Constants chosen to minimize roundoff, according to Blue's algorithm.
    static TLS_CLASS_SPEC double thres_sml = 1.4916681462400413e-154;
    static TLS_CLASS_SPEC double thres_big = 1.9979190722022350e+146;
    static TLS_CLASS_SPEC double scale_sml = 4.4989137945431964e+161;
    static TLS_CLASS_SPEC double scale_big = 1.1113793747425387e-162;

    double scale;
    int isbig = FALSE;

    integer i = 0;

    if(*incx == 1)
    {
        // AVX-2 code-section
        // Partial sums used for scaling.
        v4df_t sum_med_vec0, sum_big_vec0, sum_sml_vec0;
        v4df_t sum_med_vec1, sum_big_vec1, sum_sml_vec1;

        // Vectors used for comparisons and getting absolute values.
        v4df_t thres_sml_vec, thres_big_vec, scale_sml_vec, scale_big_vec;
        v4df_t temp, zerov;

        sum_med_vec0.v = _mm256_setzero_pd();
        sum_big_vec0.v = _mm256_setzero_pd();
        sum_sml_vec0.v = _mm256_setzero_pd();
        sum_med_vec1.v = _mm256_setzero_pd();
        sum_big_vec1.v = _mm256_setzero_pd();
        sum_sml_vec1.v = _mm256_setzero_pd();

        // Pre-broadcasting the thresholds and scale factors before entering the loops
        thres_sml_vec.v = _mm256_broadcast_sd(&thres_sml);
        thres_big_vec.v = _mm256_broadcast_sd(&thres_big);
        scale_sml_vec.v = _mm256_broadcast_sd(&scale_sml);
        scale_big_vec.v = _mm256_broadcast_sd(&scale_big);

        // This is used to convert the values in a vector to their absolute value
        temp.v = _mm256_set1_pd(-0.0);

        // Vectors used for loading from memory and setting masks
        v4df_t x0v, x1v, mask_vec0, mask_vec1;

        for(; (i + 8) <= *sd; i = i + 8)
        {
            x0v.v = _mm256_loadu_pd(xt);
            x1v.v = _mm256_loadu_pd(xt + 4);

            // Getting the abs of the vector elements.
            x0v.v = _mm256_andnot_pd(temp.v, x0v.v);
            x1v.v = _mm256_andnot_pd(temp.v, x1v.v);

            // Mask vectors which indicate whether
            // xi <= thres_sml or xi >= thres_big.
            mask_vec0.v = CMP256_df(x0v.v, thres_sml_vec.v, thres_big_vec.v);
            mask_vec1.v = CMP256_df(x1v.v, thres_sml_vec.v, thres_big_vec.v);

            if(!bli_horizontal_or_df(mask_vec0.v))
            {
                // Scaling is not necessary; only medium values.
                sum_med_vec0.v = _mm256_fmadd_pd(x0v.v, x0v.v, sum_med_vec0.v);
            }
            else
            {
                // Mask vector which indicate whether xi > thres_big.
                mask_vec0.v = _mm256_cmp_pd(x0v.v, thres_big_vec.v, _CMP_GT_OQ);
                zerov.v = _mm256_setzero_pd();

                if(bli_horizontal_or_df(mask_vec0.v))
                {
                    isbig = TRUE;

                    // Fill sum_med vector without scaling.
                    zerov.v = _mm256_blendv_pd(x0v.v, zerov.v, mask_vec0.v);
                    sum_med_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_med_vec0.v);

                    // Fill sum_big vector using scaling.
                    zerov.v = _mm256_setzero_pd();
                    zerov.v = _mm256_blendv_pd(zerov.v, scale_big_vec.v, mask_vec0.v);
                    zerov.v = _mm256_mul_pd(x0v.v, zerov.v);
                    sum_big_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_big_vec0.v);
                }
                else
                {
                    // Mask vector which indicates whether xi > thres_small.
                    mask_vec0.v = _mm256_cmp_pd(x0v.v, thres_sml_vec.v, _CMP_LT_OQ);
                    // Fill sum_med vector without scaling.
                    zerov.v = _mm256_blendv_pd(x0v.v, zerov.v, mask_vec0.v);
                    sum_med_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_med_vec0.v);

                    // Accumulate small values only if there have not been any big values so far.
                    if(!isbig)
                    {
                        // Fill sum_sml vector using scaling.
                        zerov.v = _mm256_setzero_pd();
                        zerov.v = _mm256_blendv_pd(zerov.v, scale_sml_vec.v, mask_vec0.v);
                        zerov.v = _mm256_mul_pd(x0v.v, zerov.v);
                        sum_sml_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_sml_vec0.v);
                    }
                }
            }

            if(!bli_horizontal_or_df(mask_vec1.v))
            {
                // Scaling is not necessary; only medium values.
                sum_med_vec1.v = _mm256_fmadd_pd(x1v.v, x1v.v, sum_med_vec1.v);
            }
            else
            {
                // Mask vector which indicate whether xi > thres_big.
                mask_vec1.v = _mm256_cmp_pd(x1v.v, thres_big_vec.v, _CMP_GT_OQ);

                zerov.v = _mm256_setzero_pd();

                if(bli_horizontal_or_df(mask_vec1.v))
                {
                    isbig = TRUE;

                    // Fill sum_med vector without scaling.
                    zerov.v = _mm256_blendv_pd(x1v.v, zerov.v, mask_vec1.v);
                    sum_med_vec1.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_med_vec1.v);

                    // Fill sum_big vector using scaling.
                    zerov.v = _mm256_setzero_pd();
                    zerov.v = _mm256_blendv_pd(zerov.v, scale_big_vec.v, mask_vec1.v);
                    zerov.v = _mm256_mul_pd(x1v.v, zerov.v);
                    sum_big_vec1.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_big_vec1.v);
                }
                else
                {
                    // Mask vector which indicates whether xi > thres_small.
                    mask_vec1.v = _mm256_cmp_pd(x1v.v, thres_sml_vec.v, _CMP_LT_OQ);
                    // Fill sum_med vector without scaling.
                    zerov.v = _mm256_blendv_pd(x1v.v, zerov.v, mask_vec1.v);
                    sum_med_vec1.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_med_vec1.v);

                    // Accumulate small values only if there have not been any big values so far.
                    if(!isbig)
                    {
                        // Fill sum_sml vector using scaling.
                        zerov.v = _mm256_setzero_pd();
                        zerov.v = _mm256_blendv_pd(zerov.v, scale_sml_vec.v, mask_vec1.v);
                        zerov.v = _mm256_mul_pd(x1v.v, zerov.v);
                        sum_sml_vec1.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_sml_vec1.v);
                    }
                }
            }

            xt += 8;
        }

        for(; (i + 4) <= *sd; i = i + 4)
        {
            x0v.v = _mm256_loadu_pd(xt);

            // Getting the abs of the vector elements.
            x0v.v = _mm256_andnot_pd(temp.v, x0v.v);

            // Mask vectors which indicate whether
            // xi<=thres_sml or xi>=thres_big.
            mask_vec0.v = CMP256_df(x0v.v, thres_sml_vec.v, thres_big_vec.v);

            if(!bli_horizontal_or_df(mask_vec0.v))
            {
                // Scaling is not necessary; only medium values.
                sum_med_vec0.v = _mm256_fmadd_pd(x0v.v, x0v.v, sum_med_vec0.v);
            }
            else
            {
                // Mask vector which indicate whether xi > thres_big.
                mask_vec0.v = _mm256_cmp_pd(x0v.v, thres_big_vec.v, _CMP_GT_OQ);
                zerov.v = _mm256_setzero_pd();

                if(bli_horizontal_or_df(mask_vec0.v))
                {
                    isbig = TRUE;

                    // Fill sum_med vector without scaling.
                    zerov.v = _mm256_blendv_pd(x0v.v, zerov.v, mask_vec0.v);
                    sum_med_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_med_vec0.v);

                    // Fill sum_big vector using scaling.
                    zerov.v = _mm256_setzero_pd();
                    zerov.v = _mm256_blendv_pd(zerov.v, scale_big_vec.v, mask_vec0.v);
                    zerov.v = _mm256_mul_pd(x0v.v, zerov.v);
                    sum_big_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_big_vec0.v);
                }
                else
                {
                    // Mask vector which indicates whether xi > thres_small.
                    mask_vec0.v = _mm256_cmp_pd(x0v.v, thres_sml_vec.v, _CMP_LT_OQ);
                    // Fill sum_med vector without scaling.
                    zerov.v = _mm256_blendv_pd(x0v.v, zerov.v, mask_vec0.v);
                    sum_med_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_med_vec0.v);

                    // Accumulate small values only if there have not been any big values so far.
                    if(!isbig)
                    {
                        // Fill sum_sml vector using scaling.
                        zerov.v = _mm256_setzero_pd();
                        zerov.v = _mm256_blendv_pd(zerov.v, scale_sml_vec.v, mask_vec0.v);
                        zerov.v = _mm256_mul_pd(x0v.v, zerov.v);
                        sum_sml_vec0.v = _mm256_fmadd_pd(zerov.v, zerov.v, sum_sml_vec0.v);
                    }
                }
            }
            xt += 4;
        }

        sum_sml_vec0.v = _mm256_add_pd(sum_sml_vec0.v, sum_sml_vec1.v);
        sum_med_vec0.v = _mm256_add_pd(sum_med_vec0.v, sum_med_vec1.v);
        sum_big_vec0.v = _mm256_add_pd(sum_big_vec0.v, sum_big_vec1.v);

        sum_sml += sum_sml_vec0.v[0] + sum_sml_vec0.v[1] + sum_sml_vec0.v[2] + sum_sml_vec0.v[3];
        sum_med += sum_med_vec0.v[0] + sum_med_vec0.v[1] + sum_med_vec0.v[2] + sum_med_vec0.v[3];
        sum_big += sum_big_vec0.v[0] + sum_big_vec0.v[1] + sum_big_vec0.v[2] + sum_big_vec0.v[3];
    }

    // Dealing with fringe cases
    double abs_chi;
    for(; i < *sd; i += 1)
    {
        abs_chi = fabs(*xt);
        
        if((abs_chi <= thres_big) && (abs_chi >= thres_sml))
        {
            sum_med += abs_chi * abs_chi;
        }
        // Case where there could be an overflow. Scaling is required.
        else if(abs_chi > thres_big)
        {
            sum_big += (abs_chi * scale_big) * (abs_chi * scale_big);
            isbig = TRUE;
        }
        // Case where there could be an underflow. Scaling is required.
        else if((!isbig) && (abs_chi < thres_sml))
        {
            sum_sml += (abs_chi * scale_sml) * (abs_chi * scale_sml);
        }

        xt += *incx;
    }

    // Combine accumulators.
    if(isbig)
    {
        // Combine sum_big and sum_med if sum_med > 0.
        if(sum_med > 0.0)
        {
            sum_big += (sum_med * scale_big) * scale_big;
        }
        scale = 1.0 / scale_big;
        sumsq = sum_big;
    }
    else if(sum_sml > 0.0)
    {
        // Combine sum_med and sum_sml if sum_sml>0.
        if(sum_med > 0.0)
        {
            sum_med = sqrt(sum_med);
            sum_sml = sqrt(sum_sml) / scale_sml;
            double ymin, ymax;
            if(sum_sml > sum_med)
            {
                ymin = sum_med;
                ymax = sum_sml;
            }
            else
            {
                ymin = sum_sml;
                ymax = sum_med;
            }
            scale = 1.0;
            sumsq = ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax));
        }
        else
        {
            scale = 1.0 / scale_sml;
            sumsq = sum_sml;
        }
    }
    else
    {
        // If all values are mid-range:
        scale = 1.0;
        sumsq = sum_med;
    }

    return scale * sqrt(sumsq);
}
#endif