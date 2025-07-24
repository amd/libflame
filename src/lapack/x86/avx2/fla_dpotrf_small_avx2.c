/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
extern int lapack_dpotf2(char *uplo, integer *n, double *a, integer *lda, integer *info);

int fla_dpotrf_small_avx2(char *uplo, integer *n, doublereal *a, integer *lda, integer *info)
{
    *info = 0;
    logical upper = (*uplo == 'U' || *uplo == 'u');
    logical lower = (*uplo == 'L' || *uplo == 'l');
    integer i = 1, i_1 = 1, j = 1, ni = 1, k = 0;
    doublereal *acur = NULL, p_val = 1, sum = 0, temp = 0;
    __m256d p_val_vec_256, acur_vec_256, acur_i1_vec_256, acur_j_vec_256, acur_update_vec_256;
    if(*n == 1)
    {
        /*
         * Check if the diagonal element is non-positive or NaN (not-a-number).
         * This validates that the matrix is positive definite at this pivot.
         * If the element is <= 0 or NaN, the matrix is not positive definite
         * and the Cholesky decomposition cannot proceed.
         */
        if(a[0] <= 0.f || a[0] != a[0])
        {
            *info = 1;
            return 0;
        }
        a[0] = sqrt(a[0]);
        *info = 0;
        return 0;
    }
    if(lower && *n == 2)
    {
        if(a[0] <= 0.f || a[0] != a[0])
        {
            *info = 1;
            return 0;
        }
        a[0] = sqrt(a[0]);
        a[1] = a[1] / a[0];
        p_val = (a[1] * a[1]);
        if(a[*lda + 1] <= p_val || a[*lda + 1] != a[*lda + 1])
        {
            *info = 2;
            return 0;
        }
        a[*lda + 1] = sqrt(a[*lda + 1] - p_val);
        return 0;
    }
    /* Performs Cholesky decomposition for lower triangular matrix when n < 18
       Uses loop unrolling optimization for better performance when j >= 4 */
    if(lower && *n < FLA_DPOTRF_LOWER_SMALL)
    {
        for(j = 0; j < *n; j++)
        {
            for(i = j; i < *n; i++)
            {
                if(i == j)
                {
                    {
                        sum = 0.;
                        if(j < 4)
                        {
                            for(k = 0; k < j; k++)
                            {
                                sum += a[i + (*lda * k)] * a[i + (*lda * k)];
                            }
                        }
                        else
                        {
                            k = 0;
                            for(k = 0; k < (j - 3); k = k + 4)
                            {
                                sum += a[i + (*lda * k)] * a[i + (*lda * k)];
                                sum += a[i + (*lda * (k + 1))] * a[i + (*lda * (k + 1))];
                                sum += a[i + (*lda * (k + 2))] * a[i + (*lda * (k + 2))];
                                sum += a[i + (*lda * (k + 3))] * a[i + (*lda * (k + 3))];
                            }
                            if(k < j)
                            {
                                for(; k < j; k++)
                                {
                                    sum += a[i + (*lda * k)] * a[i + (*lda * k)];
                                }
                            }
                        }
                        temp = (a[(i * *lda) + i] - sum);
                        if(temp <= 0.f || temp != temp)
                        {
                            *info = (i + 1);
                            return 0;
                        }
                        a[(i * *lda) + i] = sqrt(temp);
                    }
                }
                else
                {
                    sum = 0.;
                    if(j < 4)
                    {
                        for(k = 0; k < j; k++)
                        {
                            sum += a[i + (*lda * k)] * a[j + (*lda * k)];
                        }
                    }
                    else
                    {
                        k = 0;
                        for(k = 0; k < (j - 3); k = k + 4)
                        {
                            sum += a[i + (*lda * k)] * a[j + (*lda * k)];
                            sum += a[i + (*lda * (k + 1))] * a[j + (*lda * (k + 1))];
                            sum += a[i + (*lda * (k + 2))] * a[j + (*lda * (k + 2))];
                            sum += a[i + (*lda * (k + 3))] * a[j + (*lda * (k + 3))];
                        }
                        if(k < j)
                        {
                            for(; k < j; k++)
                            {
                                sum += a[i + (*lda * k)] * a[j + (*lda * k)];
                            }
                        }
                    }
                    a[i + (*lda * j)] = (a[i + (*lda * j)] - sum) / a[j + (*lda * j)];
                }
            }
        }
        return 0;
    }
    /* Cholesky decomposition for lower triangular matrix with medium size optimization
       Uses AVX2 vectorization for better performance on matrix with n < 75 */
    if(lower && *n < FLA_DPOTRF_LOWER_MEDIUM)
    {
        for(i = 0; i < *n; i++)
        {
            ni = *n - i;
            acur = &a[i + (*lda * i)];
            if(*acur <= 0.f || *acur != *acur)
            {
                *info = i + 1;
                return 0;
            }
            *acur = sqrt(*acur);
            p_val = 1 / *acur;
            /* Broadcast p_val to all elements of an AVX register*/
            p_val_vec_256 = _mm256_broadcast_sd(&p_val);

            /* Vectorized loop for multipling acur[i_1] by p_val */
            for(i_1 = 1; i_1 < ni - 3; i_1 += 4)
            {
                /* Load 4 elements of acur[i_1] */
                acur_vec_256 = _mm256_loadu_pd(&acur[i_1]);
                /* Multipy by p_val */
                acur_vec_256 = _mm256_mul_pd(acur_vec_256, p_val_vec_256);
                /* Store the result back to acur[i_1] */
                _mm256_storeu_pd(&acur[i_1], acur_vec_256);
            }

            /* Handle remaining elements (if any) with scalar code */
            for(; i_1 < ni; i_1++)
            {
                acur[i_1] = acur[i_1] * p_val;
            }

            /* Vectorized update of acur[i_1 + (j * *lda)] using AVX2 fused multiply-add
               Performs rank-1 update operation for Cholesky factorization */
            for(i_1 = 1; i_1 < ni; i_1++)
            {
                integer i_1xlda = (i_1 * *lda);
                /* Broadcast acur[i_1] */
                acur_i1_vec_256 = _mm256_broadcast_sd(&acur[i_1]);
                for(j = i_1; j < (ni - 3); j += 4)
                {   
                    /* Load 4 elements of acur[j] */
                    acur_j_vec_256 = _mm256_loadu_pd(&acur[j]);
                    /* Load 4 elements of acur[i_1 + (j * *lda)] */
                    acur_update_vec_256 = _mm256_loadu_pd(&acur[j + i_1xlda]);
                    /* Perform fused multiply-subtract */
                    acur_update_vec_256
                        = _mm256_fnmadd_pd(acur_j_vec_256, acur_i1_vec_256, acur_update_vec_256);
                    /* Store the result back */
                    _mm256_storeu_pd(&acur[j + i_1xlda], acur_update_vec_256);
                }

                /* Handle remaining elements (if any) with scalar code */
                for(; j < ni; j++)
                {
                    acur[j + i_1xlda] -= (acur[j] * acur[i_1]);
                }
            }
        }
        return 0;
    }
    else if(upper && *n < FLA_DPOTRF_UPPER_MEDIUM)
    {
        for(i = 0; i < *n; i++)
        {
            ni = *n - i;
            acur = &a[i + (*lda * i)];
            if(*acur <= 0.f || *acur != *acur)
            {
                *info = i + 1;
                return 0;
            }
            *acur = sqrt(*acur);
            p_val = 1 / *acur;

            for(i_1 = 1; i_1 < ni; i_1++)
            {
                acur[i_1 * *lda] = acur[i_1 * *lda] * p_val;
            }
            for(i_1 = 1; i_1 < ni; i_1++)
            {
                integer i_1xlda = (i_1 * *lda);
                for(j = 1; j <= i_1; j++)
                {
                    acur[j + i_1xlda]
                        = acur[j + i_1xlda] - acur[j * *lda] * acur[i_1 * *lda];
                }
            }
        }
        return 0;
    }
    else
    {
        lapack_dpotf2(uplo, n, a, lda, info);
    }
    return 0;
}
#endif