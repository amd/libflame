/******************************************************************************
* Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dgesvd_small6T_avx2_.c
 *  @brief DGESVD Small path (path 6T)
 *  without the LQ Factorization.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* SVD for small fat-matrices with LQ factorization
 * already computed
 */
void fla_dgesvd_small6T_avx2(integer *m, integer *n,
                             doublereal *a, integer *lda,
                             doublereal *ql, integer *ldql,
                             doublereal *s,
                             doublereal *u, integer *ldu,
                             doublereal *vt, integer *ldvt,
                             doublereal *work,
                             integer *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    integer iu, ie;
    integer itau, itauq, itaup;
    integer i__1, rlen, knt;
    integer c__1 = 1;

    doublereal *tau, *tauq, *taup;
    doublereal *e, *vtau, *avt;
    doublereal stau, d__1;

    /* indices for partitioning work buffer */
    iu = 1;
    itau = iu + *lda * *m;
    ie = itau + *m;
    itauq = ie + *m;
    itaup = itauq + *m;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    ql -= (1 + *ldql);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    FLA_BIDIAGONALIZE_SMALL(*m, *m, a, lda, tauq, taup, s, e);

    /* Generate Qr (from bidiag) in vt from work[iu] (a here) */
    for (i = 1; i <= *m; i++)
        for (j = 1; j <= *n; j++)
            vt[i + j * *ldvt] = 0.;
    FLA_LARF_VTAPPLY_DSMALL_SQR(m, a, lda, taup, vt, ldvt);

    /* Generate Ql (from bidiag) in u from a */
    FLA_LARF_UAPPLY_DSMALL_SQR(m, a, lda, tauq, u, ldu, taup);

    lapack_dbdsqr_small("U", m, m, m, &s[1], &e[1],
                        &vt[1 + *ldvt], ldvt,
                        &u[1 + *ldu], ldu,
                        info);

    /* Apply HH from LQ factorization (ql) on vt from right */

    tau = &work[itau - 1];
    vtau = tau + *m;
    avt = vtau + *n;
    /* First Iteration corresponding to HH(m) */
    i = *m;
    for (j = i + 1; j <= *n; j++)
    {
        /* - ql[i][j] * tau[i] */
        d__1 = - ql[i + j * *ldql] * tau[i];

        /* vt[1:m, j] = d__1 * vt[1:m, j] */
        for (k = 1; k <= *m; k++)
        {
            vt[k + j * *ldvt] = d__1 * vt[k + i * *ldvt];
        }
    }
    /* vt[m, 1:m] = vt[m, 1:m] * (1 - tau) */
    d__1 = 1 - tau[i];
    for (j = 1; j <= *m; j++)
    {
        vt[j + *m * *ldvt] = vt[j + *m * *ldvt] * d__1;
    }

    /* Second Iteration onwards */
    for (i = *m - 1; i >= 1; i--)
    {
        /* Scale HH vector by tau, store in vtau */
        vtau[1] = - tau[i];
        for (j = 2; j <= (*n - i + 1); j++)
        {
            vtau[j] = vtau[1] * ql[i + (j + i - 1) * *ldql];
        }

        /* avt = Vt * vtau (gemv) */
        for (j = 1; j <= *m; j++)
        {
            avt[j] = 0.;
        }
        for (j = 1; j <= (*n - i + 1); j++) /* for every column of Vt */
        {
            for (k = 1; k <= *m; k++) /* Scale the col and accumulate */
            {
                avt[k] = avt[k] + vtau[j] * vt[k + (j + i - 1) * *ldvt];
            }
        }

        /* Vt = Vt + avt * v' (ger) */
        for (k = 1; k <= *m; k++)
        {
            vt[k + i * *ldvt] = vt[k + i * *ldvt] + avt[k];
        }
        for (j = 2; j <= (*n - i + 1); j++)
        {
            for (k = 1; k <= *m; k++)
            {
                vt[k + (j + i - 1) * *ldvt] = vt[k + (j + i - 1) * *ldvt] + avt[k] * ql[i + (j + i - 1) * *ldql];
            }
        }
    }

    return;
}
#endif
