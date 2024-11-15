/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_DGESVD_SMALL_AVX2_DEFS_H
#define FLA_DGESVD_SMALL_AVX2_DEFS_H

/*! @file fla_dgesvd_small_avx2.h
 *  @brief SVD sub-routines for small sizes.
 *  */

#if FLA_ENABLE_AMD_OPT
#include "fla_dgeqrf_small_avx2.h"

doublereal d_sign(doublereal *, doublereal *);

#define FLA_BIDIAGONALIZE_SMALL(nr, nc, ia, ldia, qtau, ptau, dv, ev)    \
    for(i = 1; i <= fla_min(nr, nc); i++)                                \
    {                                                                    \
        slen = nr - i;                                                   \
        /* input address */                                              \
        doublereal *iptr;                                                \
        integer has_outliers = 0;                                        \
                                                                         \
        /* Annihilate elements in current column */                      \
        iptr = (doublereal *)&ia[i + 1 + i * *ldia - 1];                 \
        if(slen == 0)                                                    \
        {                                                                \
            qtau[i] = 0.;                                                \
            beta = 0.;                                                   \
        }                                                                \
        else if(slen < 4)                                                \
        {                                                                \
            /* Generate elementary reflector to annihilate               \
             * elements below diagonal A(i+1:nr,i) */                    \
            FLA_LARF_GEN_DSMALL_COL(i, &nr, &nc, qtau);                  \
            /* Apply the reflector on A(i:nr,i+1:nc) from the left */    \
            FLA_LARF_APPLY_DSMALL_COL(i, &nr, &nc, ia, ldia, qtau);      \
        }                                                                \
        else                                                             \
        {                                                                \
            /* Generate elementary reflector to annihilate               \
             * elements below diagonal A(i+1:nr,i) */                    \
            FLA_LARF_GEN_DLARGE_COL(i, &nr, &nc, qtau);                  \
            /* Apply the reflector on A(i:nr,i+1:nc) from the left */    \
            FLA_LARF_APPLY_DLARGE_COL(i, &nr, &nc, ia, ldia, qtau);      \
        }                                                                \
        dv[i] = *iptr;                                                   \
                                                                         \
        /* Annihilate elements in current row */                         \
        beta = 0.;                                                       \
        rlen = nc - i - 1;                                               \
        tau = ptau;                                                      \
        if(rlen <= 0)                                                    \
        {                                                                \
            tau[i] = 0.;                                                 \
        }                                                                \
        else                                                             \
        {                                                                \
            /* Generate elementary reflector to annihilate               \
             * elements to the right of current row's                    \
             * super diagonalA(i,i+2:nr) */                              \
            FLA_LARF_GEN_DSMALL_ROW(i, &nr, &nc, iptr, ldia, tau);       \
            /* Apply the reflector on A(i+1:nr,i+1:nc) from the right */ \
            FLA_LARF_APPLY_DSMALL_ROW(i, &nr, &nc, iptr, ldia, tau);     \
        }                                                                \
        if(rlen >= 0)                                                    \
            ev[i] = iptr[*ldia];                                         \
    }

#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_DGESVD_SMALL_AVX2_DEFS_H */
