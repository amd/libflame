/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
void fla_dpotri_small_avx2(char *uplo, aocl_int64_t *n, double *A, aocl_int64_t *lda, aocl_int64_t *info)
{
    aocl_int64_t N = *n, i, j, k, LDA = *lda, ilda, jlda;
    char u = toupper(*uplo);
    double sum;
    __m256d vai, vaj, vsum;

    if (u == 'L')
    {
        // Step 1: Invert the triangular Cholesky factor in-place in A
        for (j = 0; j < N; ++j)
        {
            jlda = j * LDA;
            A[j + jlda] = 1.0 / A[j + jlda];
            
            for (i = j + 1; i < N; ++i)
            {
                sum = 0.0;
                for (k = j; k < i; ++k)
                    sum -= A[i + k * LDA] * A[k + jlda];
                A[i + jlda] = sum / A[i + i * LDA];
            }
        }
        
        // Step 2: Compute A⁻¹ = Linvᵀ * Linv
        for (i = 0; i < N; ++i)
        {
            ilda = i * LDA;
            
            for (j = 0; j <= i; ++j)
            {
                jlda = j * LDA;
                k = i;
                sum = 0.0;

                if (k + 3 < N)
                {
                    vsum = _mm256_setzero_pd();
                
                    for (; k + 3 < N; k += 4)
                    {
                        vai = _mm256_loadu_pd(&A[k + ilda]);
                        vaj = _mm256_loadu_pd(&A[k + jlda]);
                        vsum = _mm256_fmadd_pd(vai, vaj, vsum);
                    }
                    
                    sum = vsum[0] + vsum[1] + vsum[2] + vsum[3];
                }
                for (; k < N; ++k)
                {
                    sum += A[k + ilda] * A[k + jlda];
                }
                
                A[i + jlda] = sum;
                A[j + ilda] = sum; // Symmetric
            }
        }
    }
    else if (u == 'U')
    {
        // Step 1: Invert the triangular Cholesky factor in-place in A
        for (j = N - 1; j >= 0; --j)
        {
            jlda = j * LDA;
            A[j + jlda] = 1.0 / A[j + jlda];
            
            for (i = j - 1; i >= 0; --i)
            {
                sum = 0.0;
                for (k = i + 1; k <= j; ++k)
                    sum -= A[i + k * LDA] * A[k + jlda];
                A[i + jlda] = sum / A[i + i * LDA];
            }
        }

        // Step 2: Compute A⁻¹ = Uinv * Uinvᵀ
        for (i = 0; i < N; ++i)
        {   
            for (j = 0; j <= i; ++j)
            {
                k = i;
                sum = 0.0;

                if (k + 3 < N)
                {
                    vsum = _mm256_setzero_pd();
                    for (; k + 3 < N; k += 4)
                    {
                        vai = _mm256_set_pd(A[i + (k+3)*LDA], A[i + (k+2)*LDA], A[i + (k+1)*LDA], A[i + k*LDA]);
                        vaj = _mm256_set_pd(A[j + (k+3)*LDA], A[j + (k+2)*LDA], A[j + (k+1)*LDA], A[j + k*LDA]);
                        vsum = _mm256_fmadd_pd(vai, vaj, vsum);
                    }
                    sum = vsum[0] + vsum[1] + vsum[2] + vsum[3];
                }

                for (; k < N; ++k)
                {
                    sum += A[i + k*LDA] * A[j + k*LDA];
                }

                A[i + j * LDA] = sum;
                A[j + i * LDA] = sum; // Symmetric
            }
        }
    }

    *info = 0;
}
#endif
