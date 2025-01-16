/******************************************************************************
 * Copyright (C) 2023-24, Advanced Micro Devices, Inc. All rights reserved.
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
integer fla_sgetrf_small_avx2(integer *m, integer *n, real *a, integer *lda, integer *ipiv,
                              integer *info)
{
    integer mi, ni;
    integer i, i_1, lda_t;

    real max_val, t_val;
    real *acur, *apiv, *asrc;
    integer p_idx;
    integer min_m_n = fla_min(*m, *n);
    lda_t = *lda;

    for(i = 0; i < min_m_n; i++)
    {
        mi = *m - i;
        ni = *n - i;

        acur = &a[i + lda_t * i];

        /* Find the pivot element */
        max_val = 0;
        p_idx = i;
        for(i_1 = 0; i_1 < mi; i_1++)
        {
            t_val = acur[i_1];
            t_val = (t_val < 0.0) ? -t_val : t_val;
            if(t_val > max_val)
            {
                max_val = t_val;
                p_idx = i + i_1;
            }
        }

        apiv = a + p_idx;
        asrc = a + i;
        ipiv[i] = p_idx + 1;

        /* Swap rows and calculate a column of L */
        if(max_val != 0.0)
        {
            /* Swap entire rows */
            if(p_idx != i)
            {
                for(i_1 = 0; i_1 < *n; i_1++)
                {
                    t_val = apiv[i_1 * lda_t];
                    apiv[i_1 * *lda] = asrc[i_1 * lda_t];
                    asrc[i_1 * *lda] = t_val;
                }
            }

            /* Calculate scalefactors (L)  & update trailing matrix */
            if(mi > 1)
            {
                fla_lu_piv_small_s_update_tr_matrix_avx2(1, mi, ni, acur, *lda);
            }
        }
        else
        {
            *info = (*info == 0) ? p_idx + 1 : *info;
        }
    }

    return *info;
}
#endif
