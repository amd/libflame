/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_nn_small1T_avx2_.c
 *  @brief DGESVD Small path (path 1T)
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#ifdef FLA_ENABLE_AMD_OPT

void fla_dgesvd_nn_small1T_avx2(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                                doublereal *work, integer *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    doublereal d__1;
    doublereal *tau, *tauq, *taup;
    doublereal *e;

    integer c__0 = 0;
    integer c__1 = 1;

    integer ie;
    integer itauq, itaup;
    integer rlen, knt;

    /* indices for partitioning work buffer */
    ie = 1;
    itauq = ie + *m;
    itaup = itauq + *m;

    /* parameter adjustments */
    a -= (1 + *lda);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    FLA_BIDIAGONALIZE_SMALL(*m, *m, a, lda, tauq, taup, s, e);

    /* Compute Singular Values */
    lapack_dbdsqr_small("U", m, &c__0, &c__0, &s[1], &e[1], NULL, &c__1, NULL, &c__1, info);

    return;
}
#endif
