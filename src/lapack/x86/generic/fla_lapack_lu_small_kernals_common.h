/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/

/*! @file ffla_lapack_lu_small_kernals_comon.h
 *  @brief Common front-end functions
 *         for single and double precision
 *         to choose optimized paths
 *  *  */

#include "FLAME.h"

#if FLA_ENABLE_AMD_OPT
/*
 * LU 2x2  with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_GEN_2x2(datatype, i, n, buff_A, ldim_A, buff_p, info) \
    datatype max_val_2 = 0.0, t_val_2 = 0.0;                                   \
    datatype *acur_2, *apiv_2, *asrc_2;                                        \
    integer i_2, p_idx_2 = i, lda2 = *ldim_A;                                  \
    /* #####################     First iteration     ###################### */ \
    acur_2 = &buff_A[i + lda2 * i];                                            \
                                                                               \
    /* Find the pivot element */                                               \
    for(i_2 = 0; i_2 < 2; i_2++)                                               \
    {                                                                          \
        t_val_2 = (datatype)f2c_abs(acur_2[i_2]);                              \
        if(t_val_2 > max_val_2)                                                \
        {                                                                      \
            max_val_2 = t_val_2;                                               \
            p_idx_2 = i + i_2;                                                 \
        }                                                                      \
    }                                                                          \
    apiv_2 = buff_A + p_idx_2;                                                 \
    asrc_2 = buff_A + i;                                                       \
    buff_p[i] = p_idx_2 + 1;                                                   \
    /* Swap rows and calculate a column of L */                                \
    if(max_val_2 != 0.0)                                                       \
    {                                                                          \
        /* Swap entire rows */                                                 \
        if(p_idx_2 != i)                                                       \
        {                                                                      \
            t_val_2 = apiv_2[0];                                               \
            apiv_2[0] = asrc_2[0];                                             \
            asrc_2[0] = t_val_2;                                               \
            t_val_2 = apiv_2[lda2];                                            \
            apiv_2[lda2] = asrc_2[lda2];                                       \
            asrc_2[lda2] = t_val_2;                                            \
            if(*n >= 3)                                                        \
            {                                                                  \
                t_val_2 = apiv_2[2 * lda2];                                    \
                apiv_2[2 * lda2] = asrc_2[2 * lda2];                           \
                asrc_2[2 * lda2] = t_val_2;                                    \
                if(*n == 4)                                                    \
                {                                                              \
                    t_val_2 = apiv_2[3 * lda2];                                \
                    apiv_2[3 * lda2] = asrc_2[3 * lda2];                       \
                    asrc_2[3 * lda2] = t_val_2;                                \
                }                                                              \
            }                                                                  \
        }                                                                      \
        /* Calculate scalefactors (L)  & update trailing matrix */             \
        acur_2[1] = acur_2[1] / *acur_2;                                       \
        acur_2[1 + lda2] = acur_2[1 + lda2] - acur_2[lda2] * acur_2[1];        \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        *info = (*info == 0) ? p_idx_2 + 1 : *info;                            \
    }                                                                          \
    /* #####################     second iteration #########################*/  \
    i = i + 1;                                                                 \
    acur_2 = &buff_A[i + lda2 * i];                                            \
    p_idx_2 = i;                                                               \
    max_val_2 = 0.0;                                                           \
    /* Find the pivot element */                                               \
    t_val_2 = (datatype)f2c_abs(acur_2[0]);                                    \
    if(t_val_2 > 0.0)                                                          \
    {                                                                          \
        max_val_2 = t_val_2;                                                   \
        p_idx_2 = i;                                                           \
    }                                                                          \
                                                                               \
    apiv_2 = buff_A + p_idx_2;                                                 \
    asrc_2 = buff_A + i;                                                       \
    buff_p[i] = p_idx_2 + 1;                                                   \
    /* Swap rows and calculate a column of L */                                \
    if(max_val_2 != 0.0)                                                       \
    {                                                                          \
        /* Swap entire rows */                                                 \
        if(p_idx_2 != i)                                                       \
        {                                                                      \
            t_val_2 = apiv_2[0];                                               \
            apiv_2[0] = asrc_2[0];                                             \
            asrc_2[0] = t_val_2;                                               \
            t_val_2 = apiv_2[lda2];                                            \
            apiv_2[lda2] = asrc_2[lda2];                                       \
            asrc_2[lda2] = t_val_2;                                            \
            if(*n >= 3)                                                        \
            {                                                                  \
                t_val_2 = apiv_2[2 * lda2];                                    \
                apiv_2[2 * lda2] = asrc_2[2 * lda2];                           \
                asrc_2[2 * lda2] = t_val_2;                                    \
                if(*n == 4)                                                    \
                {                                                              \
                    t_val_2 = apiv_2[3 * lda2];                                \
                    apiv_2[3 * lda2] = asrc_2[3 * lda2];                       \
                    asrc_2[3 * lda2] = t_val_2;                                \
                }                                                              \
            }                                                                  \
        }                                                                      \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        *info = (*info == 0) ? p_idx_2 + 1 : *info;                            \
    }
/*
 * LU 3x3 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_GEN_3x3(datatype, i, n, buff_A, ldim_A, buff_p, info)      \
    integer i_3;                                                                    \
    datatype max_val_3, t_val_3 = 0.0;                                              \
    datatype *acur_3, *apiv_3, *asrc_3;                                             \
    integer p_idx_3 = i, lda3 = *ldim_A;                                            \
    /* ################### 1st loop coputation  (i = 0) ####################*/      \
    acur_3 = &buff_A[i + lda3 * i];                                                 \
    /* Find the pivot element */                                                    \
    max_val_3 = 0.0;                                                                \
    p_idx_3 = i;                                                                    \
    for(i_3 = 0; i_3 < 3; i_3++)                                                    \
    {                                                                               \
        t_val_3 = (datatype)f2c_abs(acur_3[i_3]);                                   \
        if(t_val_3 > max_val_3)                                                     \
        {                                                                           \
            max_val_3 = t_val_3;                                                    \
            p_idx_3 = i + i_3;                                                      \
        }                                                                           \
    }                                                                               \
    apiv_3 = buff_A + p_idx_3;                                                      \
    asrc_3 = buff_A + i;                                                            \
    buff_p[i] = p_idx_3 + 1;                                                        \
    /* Swap rows and calculate a column of L */                                     \
    if(max_val_3 != 0.0)                                                            \
    {                                                                               \
        /* Swap entire rows */                                                      \
        if(p_idx_3 != i)                                                            \
        {                                                                           \
            t_val_3 = apiv_3[0];                                                    \
            apiv_3[0] = asrc_3[0];                                                  \
            asrc_3[0] = t_val_3;                                                    \
            t_val_3 = apiv_3[lda3];                                                 \
            apiv_3[lda3] = asrc_3[lda3];                                            \
            asrc_3[lda3] = t_val_3;                                                 \
            t_val_3 = apiv_3[2 * lda3];                                             \
            apiv_3[2 * lda3] = asrc_3[2 * lda3];                                    \
            asrc_3[2 * lda3] = t_val_3;                                             \
            if(*n == 4)                                                             \
            {                                                                       \
                t_val_3 = apiv_3[3 * lda3];                                         \
                apiv_3[3 * lda3] = asrc_3[3 * lda3];                                \
                asrc_3[3 * lda3] = t_val_3;                                         \
            }                                                                       \
        }                                                                           \
        /* Calculate scalefactors (L)  & update trailing matrix */                  \
        t_val_3 = 1 / *acur_3;                                                      \
        acur_3[1] = acur_3[1] * t_val_3;                                            \
        acur_3[1 + lda3] = acur_3[1 + lda3] - acur_3[lda3] * acur_3[1];             \
        acur_3[1 + 2 * lda3] = acur_3[1 + 2 * lda3] - acur_3[2 * lda3] * acur_3[1]; \
        acur_3[2] = acur_3[2] * t_val_3;                                            \
        acur_3[2 + lda3] = acur_3[2 + lda3] - acur_3[lda3] * acur_3[2];             \
        acur_3[2 + 2 * lda3] = acur_3[2 + 2 * lda3] - acur_3[2 * lda3] * acur_3[2]; \
    }                                                                               \
    else                                                                            \
    {                                                                               \
        *info = (*info == 0) ? p_idx_3 + 1 : *info;                                 \
    }                                                                               \
    i = i + 1;                                                                      \
    FLA_LU_PIV_SMALL_GEN_2x2(datatype, i, n, buff_A, ldim_A, buff_p, info);
/*
 * LU 4x4 with partial pivoting for tiny matrices
 */
#define FLA_LU_PIV_SMALL_GEN_4x4(datatype, i, n, buff_A, ldim_A, buff_p, info) \
    integer i_1;                                                               \
    datatype max_val, t_val = 0.0;                                             \
    datatype *acur, *apiv, *asrc;                                              \
    integer p_idx, lda = *ldim_A;                                              \
    /* ########### 1st loop coputation  (i = 0) ####################*/         \
    acur = &buff_A[i + lda * i];                                               \
                                                                               \
    /* Find the pivot element */                                               \
    max_val = 0.0;                                                             \
    p_idx = i;                                                                 \
    for(i_1 = 0; i_1 < 4; i_1++)                                               \
    {                                                                          \
        t_val = (datatype)f2c_abs(acur[i_1]);                                  \
        if(t_val > max_val)                                                    \
        {                                                                      \
            max_val = t_val;                                                   \
            p_idx = i + i_1;                                                   \
        }                                                                      \
    }                                                                          \
                                                                               \
    apiv = buff_A + p_idx;                                                     \
    asrc = buff_A + i;                                                         \
    buff_p[i] = p_idx + 1;                                                     \
                                                                               \
    /* Swap rows and calculate a column of L */                                \
    if(max_val != 0.0)                                                         \
    {                                                                          \
        /* Swap entire rows */                                                 \
        if(p_idx != i)                                                         \
        {                                                                      \
            t_val = apiv[0];                                                   \
            apiv[0] = asrc[0];                                                 \
            asrc[0] = t_val;                                                   \
            t_val = apiv[lda];                                                 \
            apiv[lda] = asrc[lda];                                             \
            asrc[lda] = t_val;                                                 \
            t_val = apiv[2 * lda];                                             \
            apiv[2 * lda] = asrc[2 * lda];                                     \
            asrc[2 * lda] = t_val;                                             \
            t_val = apiv[3 * lda];                                             \
            apiv[3 * lda] = asrc[3 * lda];                                     \
            asrc[3 * lda] = t_val;                                             \
        }                                                                      \
        /* Calculate scalefactors (L)  & update trailing matrix */             \
        t_val = 1 / *acur;                                                     \
        acur[1] = acur[1] * t_val;                                             \
        acur[1 + lda] = acur[1 + lda] - acur[lda] * acur[1];                   \
        acur[1 + 2 * lda] = acur[1 + 2 * lda] - acur[2 * lda] * acur[1];       \
        acur[1 + 3 * lda] = acur[1 + 3 * lda] - acur[3 * lda] * acur[1];       \
        acur[2] = acur[2] * t_val;                                             \
        acur[2 + lda] = acur[2 + lda] - acur[lda] * acur[2];                   \
        acur[2 + 2 * lda] = acur[2 + 2 * lda] - acur[2 * lda] * acur[2];       \
        acur[2 + 3 * lda] = acur[2 + 3 * lda] - acur[3 * lda] * acur[2];       \
        acur[3] = acur[3] * t_val;                                             \
        acur[3 + lda] = acur[3 + lda] - acur[lda] * acur[3];                   \
        acur[3 + 2 * lda] = acur[3 + 2 * lda] - acur[2 * lda] * acur[3];       \
        acur[3 + 3 * lda] = acur[3 + 3 * lda] - acur[3 * lda] * acur[3];       \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        *info = (*info == 0) ? p_idx + 1 : *info;                              \
    }                                                                          \
    i = i + 1;                                                                 \
    FLA_LU_PIV_SMALL_GEN_3x3(datatype, i, n, buff_A, ldim_A, buff_p, info);

#endif