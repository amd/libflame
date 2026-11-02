/******************************************************************************
 * * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/
/*! @file fla_lapack_x86_common_nofma.c
 *  @brief Common front-end functions for which strictly 
           FMA is not used.
 *  *  */

#include "fla_lapack_x86_common.h"
#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_avx512_kernels.h"

#if FLA_ENABLE_AMD_OPT


/* DLARF1F for small sizes
 * To be used only when vectorized code via avx2/avx512 is enabled
 * */
void fla_dlarf1f_small_left_incv1_simd(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                       doublereal *work)
{
    /* Below one full vector block of rows a scalar loop beats the kernels */
    if(m < 4)
    {
        aocl_int64_t i, j;
        doublereal prod;
        for(j = 0; j < n; j++)
        {
            doublereal *acol = &a_buff[j * ldr];
            doublereal dtmp = acol[0]; /* v(1) = 1 */
            for(i = 1; i < m; i++)
            {
                prod = acol[i] * v[i];
                dtmp += prod;
            }
            work[j] = dtmp;
            dtmp *= ntau;
            acol[0] += dtmp;
            for(i = 1; i < m; i++)
            {
                prod = dtmp * v[i];
                acol[i] += prod;
            }
        }
        return;
    }
    /* Select AVX512 kernel based on preset threshold and ISA support  */
    if(m > FLA_DGEMV_DGER_SIMD_AVX512_THRESH_M && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_dlarf1f_left_apply_incv1_avx512(m, n, a_buff, ldr, v, ntau, work);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dlarf1f_left_apply_incv1_avx2(m, n, a_buff, ldr, v, ntau, work);
    }
    return;
}

/* DLARF1L for small sizes
 * To be used only when vectorized code via avx2/avx512 is enabled
 * */
void fla_dlarf1l_small_left_incv1_simd(aocl_int64_t m, aocl_int64_t n, doublereal *a_buff,
                                       aocl_int64_t ldr, doublereal *v, doublereal ntau,
                                       doublereal *work)
{
    /* Below one full vector block of rows a scalar loop beats the kernels */
    if(m < 4)
    {
        aocl_int64_t i, j, mlast = m - 1;
        doublereal prod;
        for(j = 0; j < n; j++)
        {
            doublereal *acol = &a_buff[j * ldr];
            doublereal dtmp = acol[mlast]; /* v(m) = 1 */
            for(i = 0; i < mlast; i++)
            {
                prod = acol[i] * v[i];
                dtmp += prod;
            }
            work[j] = dtmp;
            dtmp *= ntau;
            acol[mlast] += dtmp;
            for(i = 0; i < mlast; i++)
            {
                prod = dtmp * v[i];
                acol[i] += prod;
            }
        }
        return;
    }
    /* Select AVX512 kernel based on preset threshold and ISA support  */
    if(m > FLA_DGEMV_DGER_SIMD_AVX512_THRESH_M && FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512))
    {
        fla_dlarf1l_left_apply_incv1_avx512(m, n, a_buff, ldr, v, ntau, work);
    }
    else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
    {
        fla_dlarf1l_left_apply_incv1_avx2(m, n, a_buff, ldr, v, ntau, work);
    }
    return;
}

#endif /* FLA_ENABLE_AMD_OPT */