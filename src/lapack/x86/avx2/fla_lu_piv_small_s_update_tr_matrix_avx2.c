/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
void fla_lu_piv_small_s_update_tr_matrix_avx2(integer i_1, integer mi, integer ni, real *acur,
                                              integer lda_t)
{

    integer j, b_off, y_off;
    real p_val;
    __m256 result[4], tempY[4], tempb[4], tempx, p_val8;

    p_val = *acur;
    p_val = 1 / p_val;
    p_val8 = _mm256_broadcast_ss(&p_val);
    for(; i_1 < mi - 7; i_1 += 8)
    {
        tempx = _mm256_mul_ps(_mm256_loadu_ps(&acur[i_1]), p_val8);
        _mm256_storeu_ps(&acur[i_1], tempx);

        for(j = 1; j < ni - 3; j += 4)
        {
            b_off = j * lda_t;
            y_off = i_1 + j * lda_t;
            tempb[0] = _mm256_broadcast_ss(&acur[b_off]);
            tempb[1] = _mm256_broadcast_ss(&acur[b_off] + 1 * lda_t);
            tempb[2] = _mm256_broadcast_ss(&acur[b_off] + 2 * lda_t);
            tempb[3] = _mm256_broadcast_ss(&acur[b_off] + 3 * lda_t);
            tempY[0] = _mm256_loadu_ps(&acur[y_off]);
            tempY[1] = _mm256_loadu_ps(&acur[y_off] + 1 * lda_t);
            tempY[2] = _mm256_loadu_ps(&acur[y_off] + 2 * lda_t);
            tempY[3] = _mm256_loadu_ps(&acur[y_off] + 3 * lda_t);
            // Y := Y - b * x
            result[0] = _mm256_fnmadd_ps(tempb[0], tempx, tempY[0]);
            result[1] = _mm256_fnmadd_ps(tempb[1], tempx, tempY[1]);
            result[2] = _mm256_fnmadd_ps(tempb[2], tempx, tempY[2]);
            result[3] = _mm256_fnmadd_ps(tempb[3], tempx, tempY[3]);
            _mm256_storeu_ps(&acur[y_off], result[0]);
            _mm256_storeu_ps((&acur[y_off] + 1 * lda_t), result[1]);
            _mm256_storeu_ps((&acur[y_off] + 2 * lda_t), result[2]);
            _mm256_storeu_ps((&acur[y_off] + 3 * lda_t), result[3]);
        }
        /* remining inner loop updation */
        for(; j < ni; j++)
        {
            b_off = j * lda_t;
            y_off = i_1 + j * lda_t;
            tempb[0] = _mm256_broadcast_ss(&acur[b_off]);
            tempY[0] = _mm256_loadu_ps(&acur[y_off]);
            result[0] = _mm256_fnmadd_ps(tempb[0], tempx, tempY[0]);
            _mm256_storeu_ps(&acur[y_off], result[0]);
        }
    }
    /* remining outer loop updation */
    for(; i_1 < mi; i_1++)
    {
        acur[i_1] = acur[i_1] * p_val;
        for(j = 1; j < ni; j++)
        {
            acur[i_1 + j * lda_t] = acur[i_1 + j * lda_t] - acur[j * lda_t] * acur[i_1];
        }
    }
}
#endif
