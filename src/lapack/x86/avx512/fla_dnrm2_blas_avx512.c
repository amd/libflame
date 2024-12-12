/******************************************************************************
 * Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dnrm2_blas_avx512.c
 *  @brief norm2 for avx512 kernel.
 *  */

#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT
// Union data structure to access AVX registers
// One 512-bit AVX register holds 8 DP elements.
typedef union
{
    __m512d v;
    double d[8] __attribute__((aligned(64)));
} v8df_t;
/* dnrm2 for small sizes */
doublereal fla_dnrm2_blas_avx512(integer *sd, doublereal *a, integer *incx)
{
    // Local variables and pointers used for the computation
    double sumsq = 0;

    // Local pointer alias to the input vector
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
    static TLS_CLASS_SPEC double thresh_sml = 1.4916681462400413e-154;
    static TLS_CLASS_SPEC double thresh_big = 1.9979190722022350e+146;
    static TLS_CLASS_SPEC double scale_sml = 4.4989137945431964e+161;
    static TLS_CLASS_SPEC double scale_big = 1.1113793747425387e-162;

    // Scaling factor to be set and used in the final accumulation
    double scale;

    // Boolean to check if any value > thresh_big has been encountered
    int isbig = FALSE;

    // Iterator
    integer i = 0;

    // In case of unit-strided input
    if(*incx == 1)
    {
        // AVX-512 code-section
        // Declaring registers for loading, accumulation, thresholds and scale factors
        v8df_t x_vec[4], sum_sml_vec[4], sum_med_vec[4], sum_big_vec[4], temp[4];
        v8df_t thresh_sml_vec, thresh_big_vec, scale_sml_vec, scale_big_vec;
        v8df_t zero_reg;

        // Masks to be used in computation
        __mmask8 k_mask[8];

        // Containers to hold the results of operations on mask registers
        // Bitwise operations on 8-bit mask registers would return an
        // unsigned char as its result(0 or 1)
        unsigned char truth_val[4];

        // Setting the thresholds and scaling factors
        thresh_sml_vec.v = _mm512_set1_pd(thresh_sml);
        thresh_big_vec.v = _mm512_set1_pd(thresh_big);
        scale_sml_vec.v = _mm512_set1_pd(scale_sml);
        scale_big_vec.v = _mm512_set1_pd(scale_big);

        // Resetting the accumulators
        sum_sml_vec[0].v = _mm512_setzero_pd();
        sum_sml_vec[1].v = _mm512_setzero_pd();
        sum_sml_vec[2].v = _mm512_setzero_pd();
        sum_sml_vec[3].v = _mm512_setzero_pd();

        sum_med_vec[0].v = _mm512_setzero_pd();
        sum_med_vec[1].v = _mm512_setzero_pd();
        sum_med_vec[2].v = _mm512_setzero_pd();
        sum_med_vec[3].v = _mm512_setzero_pd();

        sum_big_vec[0].v = _mm512_setzero_pd();
        sum_big_vec[1].v = _mm512_setzero_pd();
        sum_big_vec[2].v = _mm512_setzero_pd();
        sum_big_vec[3].v = _mm512_setzero_pd();

        zero_reg.v = _mm512_setzero_pd();

        // Computing in blocks of 32
        for(; (i + 32) <= *sd; i = i + 32)
        {
            // Set temp[0..3] to zero
            temp[0].v = _mm512_setzero_pd();
            temp[1].v = _mm512_setzero_pd();
            temp[2].v = _mm512_setzero_pd();
            temp[3].v = _mm512_setzero_pd();

            // Loading the vectors
            x_vec[0].v = _mm512_loadu_pd(xt);
            x_vec[1].v = _mm512_loadu_pd(xt + 8);
            x_vec[2].v = _mm512_loadu_pd(xt + 16);
            x_vec[3].v = _mm512_loadu_pd(xt + 24);

            // Getting the absoulte values of elements in the vectors
            x_vec[0].v = _mm512_abs_pd(x_vec[0].v);
            x_vec[1].v = _mm512_abs_pd(x_vec[1].v);
            x_vec[2].v = _mm512_abs_pd(x_vec[2].v);
            x_vec[3].v = _mm512_abs_pd(x_vec[3].v);

            // Setting the masks by comparing with thresh_sml_vec.v
            // That is, k_mask[0][i] = 1 if x_vec[0].v[i] > thresh_sml_vec.v
            //          k_mask[1][i] = 1 if x_vec[1].v[i] > thresh_sml_vec.v
            //          k_mask[2][i] = 1 if x_vec[2].v[i] > thresh_sml_vec.v
            //          k_mask[3][i] = 1 if x_vec[3].v[i] > thresh_sml_vec.v
            k_mask[0] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_sml_vec.v, _CMP_GT_OS);
            k_mask[1] = _mm512_cmp_pd_mask(x_vec[1].v, thresh_sml_vec.v, _CMP_GT_OS);
            k_mask[2] = _mm512_cmp_pd_mask(x_vec[2].v, thresh_sml_vec.v, _CMP_GT_OS);
            k_mask[3] = _mm512_cmp_pd_mask(x_vec[3].v, thresh_sml_vec.v, _CMP_GT_OS);

            // Setting the masks by comparing with thresh_big_vec.v
            // That is, k_mask[4][i] = 1 if x_vec[0].v[i] < thresh_big_vec.v
            //          k_mask[5][i] = 1 if x_vec[1].v[i] < thresh_big_vec.v
            //          k_mask[6][i] = 1 if x_vec[2].v[i] < thresh_big_vec.v
            //          k_mask[7][i] = 1 if x_vec[3].v[i] < thresh_big_vec.v
            k_mask[4] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_big_vec.v, _CMP_LT_OS);
            k_mask[5] = _mm512_cmp_pd_mask(x_vec[1].v, thresh_big_vec.v, _CMP_LT_OS);
            k_mask[6] = _mm512_cmp_pd_mask(x_vec[2].v, thresh_big_vec.v, _CMP_LT_OS);
            k_mask[7] = _mm512_cmp_pd_mask(x_vec[3].v, thresh_big_vec.v, _CMP_LT_OS);

            // Setting the masks to filter only the elements within the thresholds
            // k_mask[0 ... 3] contain masks for elements > thresh_sml
            // k_mask[4 ... 7] contain masks for elements < thresh_big
            // Thus, AND operation on these would give elements within these thresholds
            k_mask[4] = _kand_mask8(k_mask[0], k_mask[4]);
            k_mask[5] = _kand_mask8(k_mask[1], k_mask[5]);
            k_mask[6] = _kand_mask8(k_mask[2], k_mask[6]);
            k_mask[7] = _kand_mask8(k_mask[3], k_mask[7]);

            // Setting booleans to check for underflow/overflow handling
            // In case of having values outside threshold, the associated
            // bit in k_mask[4 ... 7] is 0.
            // Thus, truth_val[0] = 0 if x_vec[0].v has elements outside thresholds
            //       truth_val[1] = 0 if x_vec[1].v has elements outside thresholds
            //       truth_val[2] = 0 if x_vec[2].v has elements outside thresholds
            //       truth_val[3] = 0 if x_vec[3].v has elements outside thresholds
            truth_val[0] = _kortestc_mask8_u8(k_mask[4], k_mask[4]);
            truth_val[1] = _kortestc_mask8_u8(k_mask[5], k_mask[5]);
            truth_val[2] = _kortestc_mask8_u8(k_mask[6], k_mask[6]);
            truth_val[3] = _kortestc_mask8_u8(k_mask[7], k_mask[7]);

            // Computing using masked fmadds, that carries over values from
            // accumulator register if the mask bit is 0
            sum_med_vec[0].v
                = _mm512_mask3_fmadd_pd(x_vec[0].v, x_vec[0].v, sum_med_vec[0].v, k_mask[4]);
            sum_med_vec[1].v
                = _mm512_mask3_fmadd_pd(x_vec[1].v, x_vec[1].v, sum_med_vec[1].v, k_mask[5]);
            sum_med_vec[2].v
                = _mm512_mask3_fmadd_pd(x_vec[2].v, x_vec[2].v, sum_med_vec[2].v, k_mask[6]);
            sum_med_vec[3].v
                = _mm512_mask3_fmadd_pd(x_vec[3].v, x_vec[3].v, sum_med_vec[3].v, k_mask[7]);

            // In case of having elements outside the threshold
            if(!(truth_val[0] && truth_val[1] && truth_val[2] && truth_val[3]))
            {
                // Acquiring the masks for numbers greater than thresh_big
                // k_mask[4 ... 7] contain masks for elements within the thresholds
                // k_mask[0 ... 3] contain masks for elements > thresh_sml. This would
                // include both elements < thresh_big and >= thresh_big
                // XOR on these will produce masks for elements >= thresh_big
                // That is, k_mask[4][i] = 1 if x_vec[0].v[i] >= thresh_big_vec.v
                //          k_mask[5][i] = 1 if x_vec[1].v[i] >= thresh_big_vec.v
                //          k_mask[6][i] = 1 if x_vec[2].v[i] >= thresh_big_vec.v
                //          k_mask[7][i] = 1 if x_vec[3].v[i] >= thresh_big_vec.v
                k_mask[4] = _kxor_mask8(k_mask[0], k_mask[4]);
                k_mask[5] = _kxor_mask8(k_mask[1], k_mask[5]);
                k_mask[6] = _kxor_mask8(k_mask[2], k_mask[6]);
                k_mask[7] = _kxor_mask8(k_mask[3], k_mask[7]);

                // Inverting k_mask[0 ... 3], to obtain masks for elements <= thresh_sml
                // That is, k_mask[0][i] = 1 if x_vec[0].v[i] <= thresh_sml_vec.v
                //          k_mask[1][i] = 1 if x_vec[1].v[i] <= thresh_sml_vec.v
                //          k_mask[2][i] = 1 if x_vec[2].v[i] <= thresh_sml_vec.v
                //          k_mask[3][i] = 1 if x_vec[3].v[i] <= thresh_sml_vec.v
                k_mask[0] = _knot_mask8(k_mask[0]);
                k_mask[1] = _knot_mask8(k_mask[1]);
                k_mask[2] = _knot_mask8(k_mask[2]);
                k_mask[3] = _knot_mask8(k_mask[3]);

                // Checking whether we have values greater than thresh_big
                // The truth_val is set to 0 if any bit in the mask is 1
                // Thus, truth_val[2] = 0 if x_vec[0].v or x_vec[1].v has elements >=
                // thresh_big_vec.v
                //       truth_val[3] = 0 if x_vec[2].v or x_vec[3].v has elements >=
                //       thresh_big_vec.v
                truth_val[2] = _kortestz_mask8_u8(k_mask[4], k_mask[5]);
                truth_val[3] = _kortestz_mask8_u8(k_mask[6], k_mask[7]);

                // In case of having values greater than thresh_big
                if(!(truth_val[2] && truth_val[3]))
                {
                    // Set isbig to true
                    isbig = TRUE;

                    // Computing by breaking it into masked muls and fmadds
                    // This computation involves only the elements that
                    // are greater than thresh_big

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[4], scale_big_vec.v, x_vec[0].v);
                    temp[1].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[5], scale_big_vec.v, x_vec[1].v);
                    temp[2].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[6], scale_big_vec.v, x_vec[2].v);
                    temp[3].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[7], scale_big_vec.v, x_vec[3].v);

                    // Square and add the elements to the accumulators
                    sum_big_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_big_vec[0].v);
                    sum_big_vec[1].v = _mm512_fmadd_pd(temp[1].v, temp[1].v, sum_big_vec[1].v);
                    sum_big_vec[2].v = _mm512_fmadd_pd(temp[2].v, temp[2].v, sum_big_vec[2].v);
                    sum_big_vec[3].v = _mm512_fmadd_pd(temp[3].v, temp[3].v, sum_big_vec[3].v);
                }
                else if(!isbig)
                {
                    // Computing by breaking it into muls and adds
                    // This computation involves only the elements that
                    // are lesser than thresh_sml, if needed

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[0], scale_sml_vec.v, x_vec[0].v);
                    temp[1].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[1], scale_sml_vec.v, x_vec[1].v);
                    temp[2].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[2], scale_sml_vec.v, x_vec[2].v);
                    temp[3].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[3], scale_sml_vec.v, x_vec[3].v);

                    // Square and add the elements to the accumulators
                    sum_sml_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_sml_vec[0].v);
                    sum_sml_vec[1].v = _mm512_fmadd_pd(temp[1].v, temp[1].v, sum_sml_vec[1].v);
                    sum_sml_vec[2].v = _mm512_fmadd_pd(temp[2].v, temp[2].v, sum_sml_vec[2].v);
                    sum_sml_vec[3].v = _mm512_fmadd_pd(temp[3].v, temp[3].v, sum_sml_vec[3].v);
                }
            }

            // Updating the pointer for the next iteration
            xt += 32;
        }

        // Computing in blocks of 16
        for(; (i + 16) <= *sd; i = i + 16)
        {
            // Set temp[0..1] to zero
            temp[0].v = _mm512_setzero_pd();
            temp[1].v = _mm512_setzero_pd();

            // Loading the vectors
            x_vec[0].v = _mm512_loadu_pd(xt);
            x_vec[1].v = _mm512_loadu_pd(xt + 8);

            // Getting the absoulte values of elements in the vectors
            x_vec[0].v = _mm512_abs_pd(x_vec[0].v);
            x_vec[1].v = _mm512_abs_pd(x_vec[1].v);

            // Setting the masks by comparing with thresh_sml_vec.v
            // That is, k_mask[0][i] = 1 if x_vec[0].v[i] > thresh_sml_vec.v
            //          k_mask[1][i] = 1 if x_vec[1].v[i] > thresh_sml_vec.v
            k_mask[0] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_sml_vec.v, _CMP_GT_OS);
            k_mask[1] = _mm512_cmp_pd_mask(x_vec[1].v, thresh_sml_vec.v, _CMP_GT_OS);

            // Setting the masks by comparing with thresh_big_vec.v
            // That is, k_mask[4][i] = 1 if x_vec[0].v[i] < thresh_big_vec.v
            //          k_mask[5][i] = 1 if x_vec[1].v[i] < thresh_big_vec.v
            k_mask[4] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_big_vec.v, _CMP_LT_OS);
            k_mask[5] = _mm512_cmp_pd_mask(x_vec[1].v, thresh_big_vec.v, _CMP_LT_OS);

            // Setting the masks to filter only the elements within the thresholds
            // k_mask[0 ... 1] contain masks for elements > thresh_sml
            // k_mask[4 ... 5] contain masks for elements < thresh_big
            // Thus, AND operation on these would give elements within these thresholds
            k_mask[4] = _kand_mask8(k_mask[0], k_mask[4]);
            k_mask[5] = _kand_mask8(k_mask[1], k_mask[5]);

            // Setting booleans to check for underflow/overflow handling
            // In case of having values outside threshold, the associated
            // bit in k_mask[4 ... 7] is 0.
            // Thus, truth_val[0] = 0 if x_vec[0].v has elements outside thresholds
            //       truth_val[1] = 0 if x_vec[1].v has elements outside thresholds
            truth_val[0] = _kortestc_mask8_u8(k_mask[4], k_mask[4]);
            truth_val[1] = _kortestc_mask8_u8(k_mask[5], k_mask[5]);

            // Computing using masked fmadds, that carries over values from
            // accumulator register if the mask bit is 0
            sum_med_vec[0].v
                = _mm512_mask3_fmadd_pd(x_vec[0].v, x_vec[0].v, sum_med_vec[0].v, k_mask[4]);
            sum_med_vec[1].v
                = _mm512_mask3_fmadd_pd(x_vec[1].v, x_vec[1].v, sum_med_vec[1].v, k_mask[5]);

            // In case of having elements outside the threshold
            if(!(truth_val[0] && truth_val[1]))
            {
                // Acquiring the masks for numbers greater than thresh_big
                // k_mask[4 ... 5] contain masks for elements within the thresholds
                // k_mask[0 ... 1] contain masks for elements > thresh_sml. This would
                // include both elements < thresh_big and >= thresh_big
                // XOR on these will produce masks for elements >= thresh_big
                // That is, k_mask[4][i] = 1 if x_vec[0].v[i] >= thresh_big_vec.v
                //          k_mask[5][i] = 1 if x_vec[1].v[i] >= thresh_big_vec.v
                k_mask[4] = _kxor_mask8(k_mask[0], k_mask[4]);
                k_mask[5] = _kxor_mask8(k_mask[1], k_mask[5]);

                // Inverting k_mask[0 ... 1], to obtain masks for elements <= thresh_sml
                // That is, k_mask[0][i] = 1 if x_vec[0].v[i] <= thresh_sml_vec.v
                //          k_mask[1][i] = 1 if x_vec[1].v[i] <= thresh_sml_vec.v
                k_mask[0] = _knot_mask8(k_mask[0]);
                k_mask[1] = _knot_mask8(k_mask[1]);

                // Checking whether we have values greater than thresh_big
                // The truth_val is set to 0 if any bit in the mask is 1
                // Thus, truth_val[2] = 0 if x_vec[0].v or x_vec[1].v has elements >=
                // thresh_big_vec.v
                truth_val[2] = _kortestz_mask8_u8(k_mask[4], k_mask[5]);

                // In case of having values greater than thresh_big
                if(!truth_val[2])
                {
                    // Set isbig to true
                    isbig = TRUE;

                    // Computing by breaking it into masked muls and fmadds
                    // This computation involves only the elements that
                    // are greater than thresh_big

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[4], scale_big_vec.v, x_vec[0].v);
                    temp[1].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[5], scale_big_vec.v, x_vec[1].v);

                    // Square and add the elements to the accumulators
                    sum_big_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_big_vec[0].v);
                    sum_big_vec[1].v = _mm512_fmadd_pd(temp[1].v, temp[1].v, sum_big_vec[1].v);
                }
                else if(!isbig)
                {
                    // Computing by breaking it into muls and adds
                    // This computation involves only the elements that
                    // are lesser than thresh_sml, if needed

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[0], scale_sml_vec.v, x_vec[0].v);
                    temp[1].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[1], scale_sml_vec.v, x_vec[1].v);

                    // Square and add the elements to the accumulators
                    sum_sml_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_sml_vec[0].v);
                    sum_sml_vec[1].v = _mm512_fmadd_pd(temp[1].v, temp[1].v, sum_sml_vec[1].v);
                }
            }

            // Updating the pointer for the next iteration
            xt += 16;
        }
        for(; (i + 8) <= *sd; i = i + 8)
        {
            // Set temp[0].v to zero
            temp[0].v = _mm512_setzero_pd();

            // Loading the vectors
            x_vec[0].v = _mm512_loadu_pd(xt);

            // Getting the absoulte values of elements in the vectors
            x_vec[0].v = _mm512_abs_pd(x_vec[0].v);

            // Setting the masks by comparing with thresh_sml_vec.v
            // That is, k_mask[0][i] = 1 if x_vec[0].v[i] > thresh_sml_vec.v
            k_mask[0] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_sml_vec.v, _CMP_GT_OS);

            // Setting the masks by comparing with thresh_big_vec.v
            // That is, k_mask[4][i] = 1 if x_vec[0].v[i] < thresh_big_vec.v
            k_mask[4] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_big_vec.v, _CMP_LT_OS);

            // Setting the masks to filter only the elements within the thresholds
            // k_mask[0] contain masks for elements > thresh_sml
            // k_mask[4] contain masks for elements < thresh_big
            // Thus, AND operation on these would give elements within these thresholds
            k_mask[4] = _kand_mask8(k_mask[0], k_mask[4]);

            // Setting booleans to check for underflow/overflow handling
            // In case of having values outside threshold, the associated
            // bit in k_mask[4] is 0.
            // Thus, truth_val[0] = 0 if x_vec[0].v has elements outside thresholds
            truth_val[0] = _kortestc_mask8_u8(k_mask[4], k_mask[4]);

            // Computing using masked fmadds, that carries over values from
            // accumulator register if the mask bit is 0
            sum_med_vec[0].v
                = _mm512_mask3_fmadd_pd(x_vec[0].v, x_vec[0].v, sum_med_vec[0].v, k_mask[4]);

            // In case of having elements outside the threshold
            if(!truth_val[0])
            {
                // Acquiring the masks for numbers greater than thresh_big
                // k_mask[4 ... 5] contain masks for elements within the thresholds
                // k_mask[0 ... 1] contain masks for elements > thresh_sml. This would
                // include both elements < thresh_big and >= thresh_big
                // XOR on these will produce masks for elements >= thresh_big
                // That is, k_mask[4][i] = 1 if x_vec[0].v[i] >= thresh_big_vec.v
                //          k_mask[5][i] = 1 if x_vec[1].v[i] >= thresh_big_vec.v
                k_mask[4] = _kxor_mask8(k_mask[0], k_mask[4]);

                // Inverting k_mask[0 ... 1], to obtain masks for elements <= thresh_sml
                // That is, k_mask[0][i] = 1 if x_vec[0].v[i] <= thresh_sml_vec.v
                //          k_mask[1][i] = 1 if x_vec[1].v[i] <= thresh_sml_vec.v
                k_mask[0] = _knot_mask8(k_mask[0]);

                // Checking whether we have values greater than thresh_big
                // The truth_val is set to 0 if any bit in the mask is 1
                // Thus, truth_val[2] = 0 if x_vec[0].v or x_vec[1].v has elements >=
                // thresh_big_vec.v
                truth_val[2] = _kortestz_mask8_u8(k_mask[4], k_mask[4]);

                // In case of having values greater than thresh_big
                if(!truth_val[2])
                {
                    // Set isbig to true
                    isbig = TRUE;

                    // Computing by breaking it into masked muls and fmadds
                    // This computation involves only the elements that
                    // are greater than thresh_big

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[4], scale_big_vec.v, x_vec[0].v);

                    // Square and add the elements to the accumulators
                    sum_big_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_big_vec[0].v);
                }
                else if(!isbig)
                {
                    // Computing by breaking it into muls and adds
                    // This computation involves only the elements that
                    // are lesser than thresh_sml, if needed

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[0], scale_sml_vec.v, x_vec[0].v);

                    // Square and add the elements to the accumulators
                    sum_sml_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_sml_vec[0].v);
                }
            }

            // Updating the pointer for the next iteration
            xt += 8;
        }
        if(i < *sd)
        {
            // Set temp[0].v to zero
            temp[0].v = _mm512_setzero_pd();

            // Setting the mask to load
            k_mask[0] = (1 << (*sd - i)) - 1;

            // Loading the vectors
            x_vec[0].v = _mm512_maskz_loadu_pd(k_mask[0], xt);

            // Getting the absoulte values of elements in the vectors
            x_vec[0].v = _mm512_abs_pd(x_vec[0].v);

            // Setting the masks by comparing with thresh_sml_vec.v
            // That is, k_mask[0][i] = 1 if x_vec[0].v[i] > thresh_sml_vec.v
            k_mask[0] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_sml_vec.v, _CMP_GT_OS);

            // Setting the masks by comparing with thresh_big_vec.v
            // That is, k_mask[4][i] = 1 if x_vec[0].v[i] < thresh_big_vec.v
            k_mask[4] = _mm512_cmp_pd_mask(x_vec[0].v, thresh_big_vec.v, _CMP_LT_OS);

            // Setting the masks to filter only the elements within the thresholds
            // k_mask[0] contain masks for elements > thresh_sml
            // k_mask[4] contain masks for elements < thresh_big
            // Thus, AND operation on these would give elements within these thresholds
            k_mask[4] = _kand_mask8(k_mask[0], k_mask[4]);

            // Setting booleans to check for underflow/overflow handling
            // In case of having values outside threshold, the associated
            // bit in k_mask[4] is 0.
            // Thus, truth_val[0] = 0 if x_vec[0].v has elements outside thresholds
            truth_val[0] = _kortestc_mask8_u8(k_mask[4], k_mask[4]);

            // Computing using masked fmadds, that carries over values from
            // accumulator register if the mask bit is 0
            sum_med_vec[0].v
                = _mm512_mask3_fmadd_pd(x_vec[0].v, x_vec[0].v, sum_med_vec[0].v, k_mask[4]);

            // In case of having elements outside the threshold
            if(!truth_val[0])
            {
                // Acquiring the masks for numbers greater than thresh_big
                // k_mask[4 ... 5] contain masks for elements within the thresholds
                // k_mask[0 ... 1] contain masks for elements > thresh_sml. This would
                // include both elements < thresh_big and >= thresh_big
                // XOR on these will produce masks for elements >= thresh_big
                // That is, k_mask[4][i] = 1 if x_vec[0].v[i] >= thresh_big_vec.v
                //          k_mask[5][i] = 1 if x_vec[1].v[i] >= thresh_big_vec.v
                k_mask[4] = _kxor_mask8(k_mask[0], k_mask[4]);

                // Inverting k_mask[0 ... 1], to obtain masks for elements <= thresh_sml
                // That is, k_mask[0][i] = 1 if x_vec[0].v[i] <= thresh_sml_vec.v
                //          k_mask[1][i] = 1 if x_vec[1].v[i] <= thresh_sml_vec.v
                k_mask[0] = _knot_mask8(k_mask[0]);

                // Checking whether we have values greater than thresh_big
                // The truth_val is set to 0 if any bit in the mask is 1
                // Thus, truth_val[2] = 0 if x_vec[0].v or x_vec[1].v has elements >=
                // thresh_big_vec.v
                truth_val[2] = _kortestz_mask8_u8(k_mask[4], k_mask[4]);

                // In case of having values greater than thresh_big
                if(!truth_val[2])
                {
                    // Set isbig to true
                    isbig = TRUE;

                    // Computing by breaking it into masked muls and fmadds
                    // This computation involves only the elements that
                    // are greater than thresh_big

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[4], scale_big_vec.v, x_vec[0].v);

                    // Square and add the elements to the accumulators
                    sum_big_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_big_vec[0].v);
                }
                else if(!isbig)
                {
                    // Computing by breaking it into muls and adds
                    // This computation involves only the elements that
                    // are lesser than thresh_sml, if needed

                    // Scale the required elements in x_vec[0..3] by scale_smal
                    temp[0].v
                        = _mm512_mask_mul_pd(zero_reg.v, k_mask[0], scale_sml_vec.v, x_vec[0].v);

                    // Square and add the elements to the accumulators
                    sum_sml_vec[0].v = _mm512_fmadd_pd(temp[0].v, temp[0].v, sum_sml_vec[0].v);
                }
            }
        }

        // Reduction step
        // Combining the results of accumulators for each category
        sum_med_vec[0].v = _mm512_add_pd(sum_med_vec[0].v, sum_med_vec[1].v);
        sum_med_vec[2].v = _mm512_add_pd(sum_med_vec[2].v, sum_med_vec[3].v);
        sum_med_vec[0].v = _mm512_add_pd(sum_med_vec[0].v, sum_med_vec[2].v);

        sum_big_vec[0].v = _mm512_add_pd(sum_big_vec[0].v, sum_big_vec[1].v);
        sum_big_vec[2].v = _mm512_add_pd(sum_big_vec[2].v, sum_big_vec[3].v);
        sum_big_vec[0].v = _mm512_add_pd(sum_big_vec[0].v, sum_big_vec[2].v);

        sum_sml_vec[0].v = _mm512_add_pd(sum_sml_vec[0].v, sum_sml_vec[1].v);
        sum_sml_vec[2].v = _mm512_add_pd(sum_sml_vec[2].v, sum_sml_vec[3].v);
        sum_sml_vec[0].v = _mm512_add_pd(sum_sml_vec[0].v, sum_sml_vec[2].v);

        // Final accumulation on the scalars
        sum_sml += sum_sml_vec[0].d[0] + sum_sml_vec[0].d[1] + sum_sml_vec[0].d[2]
                   + sum_sml_vec[0].d[3] + sum_sml_vec[0].d[4] + sum_sml_vec[0].d[5]
                   + sum_sml_vec[0].d[6] + sum_sml_vec[0].d[7];
        sum_med += sum_med_vec[0].d[0] + sum_med_vec[0].d[1] + sum_med_vec[0].d[2]
                   + sum_med_vec[0].d[3] + sum_med_vec[0].d[4] + sum_med_vec[0].d[5]
                   + sum_med_vec[0].d[6] + sum_med_vec[0].d[7];
        sum_big += sum_big_vec[0].d[0] + sum_big_vec[0].d[1] + sum_big_vec[0].d[2]
                   + sum_big_vec[0].d[3] + sum_big_vec[0].d[4] + sum_big_vec[0].d[5]
                   + sum_big_vec[0].d[6] + sum_big_vec[0].d[7];
    }
    // Dealing with non-unit strided inputs
    else
    {
        // Dealing with fringe cases
        double abs_chi;
        for(; i < *sd; i += 1)
        {
            abs_chi = fabs(*xt);
            // Most likely case: medium values, not over/under-flow.
            if((abs_chi <= thresh_big) && (abs_chi >= thresh_sml))
            {
                sum_med += abs_chi * abs_chi;
            }
            // Case where there could be an overflow. Scaling is required.
            else if(abs_chi > thresh_big)
            {
                sum_big += (abs_chi * scale_big) * (abs_chi * scale_big);
                isbig = TRUE;
            }
            // Case where there could be an underflow. Scaling is required.
            else if((!isbig) && (abs_chi < thresh_sml))
            {
                sum_sml += (abs_chi * scale_sml) * (abs_chi * scale_sml);
            }

            xt += *incx;
        }
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