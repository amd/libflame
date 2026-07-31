/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_sgeqrf_small_avx2.c
 *  @brief QR for small inputs in AVX2.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"
#include "fla_sgeqrf_small_avx2.h"

#if FLA_ENABLE_AMD_OPT

static aocl_int64_t c__1 = 1;
/* QR for small sizes */
int fla_sgeqrf_small_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda, real *tau,
                          real *work)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_SSMALL();

    aocl_int64_t min_m_n;

    /* Adjust pointers */
    a -= (1 + *lda * 1);
    tau--;
    work--;

    min_m_n = fla_min(*m, *n);
    for(i = 1; i <= min_m_n; i++)
    {
        slen = *m - i;
        /* input address */
        real *iptr = &a[i + 1 + i * *lda - 1];
        aocl_int64_t has_outliers = 0;

        if(slen <= 0)
        {
            tau[i] = 0;
        }
        else if(slen < FLA_GEQRF_SLARGE_MIN_LEN)
        {
            FLA_LARF_GEN_SSMALL_COL(i, m, n, tau);
            FLA_LARF_APPLY_SSMALL_COL(i, m, n, a, lda, tau);
        }
        else
        {
            FLA_LARF_GEN_SLARGE_COL(i, m, n, tau);
            FLA_LARF_APPLY_SLARGE_COL(i, m, n, a, lda, tau);
        }
    }
    return 0;
}
#endif
