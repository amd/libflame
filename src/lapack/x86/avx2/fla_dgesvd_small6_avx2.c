/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_small6_avx2_.c
 *  @brief DGESVD Small path (path 6)
 *  without the LQ Factorization.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_x86_common.h"

#if FLA_ENABLE_AMD_OPT

/* SVD for small tall-matrices with QR factorization
 * already computed
 */
void fla_dgesvd_small6_avx2(integer wntus, integer wntvs, integer *m, integer *n, doublereal *a,
                            integer *lda, doublereal *qr, integer *ldqr, doublereal *s,
                            doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
                            doublereal *work, integer *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    integer ie;
    integer itau, itauq, itaup;
    integer rlen, knt;
    integer ni;
    integer tn;
    integer ncvt, nru;
    integer *ldau;
    integer c__1 = 1;

    doublereal *tau, *tauq, *taup;
    doublereal *e, *au;
    doublereal stau, d__1;
    doublereal dum[2];
    doublereal c_zero = 0.;

    /* indices for partitioning work buffer */
    itau = 1;
    ie = itau + *n;
    itauq = ie + *n;
    itaup = itauq + *n;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    qr -= (1 + *ldqr);
    --s;
    --work;

    /* local variables initialization */
    v = &dum[0];
    ncvt = 0;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* QR Factorization */
    fla_dgeqrf_small(m, n, &a[1 + *lda], lda, &work[itau], &work[ie]);

    /* Upper Bidiagonalization */
    if(wntus)
    {
        nru = *n;
        au = u;
        ldau = ldu;
        /* Copy R to U */
        dlacpy_("U", n, n, &a[1 + *lda], lda, &au[1 + *ldau], ldau);
    }
    else
    {
        nru = 0;
        au = a;
        ldau = lda;
    }
    /* Set lower part of U to zero */
    tn = *n - 1;
    dlaset_("L", &tn, &tn, &c_zero, &c_zero, &au[2 + *ldau], ldau);
    FLA_BIDIAGONALIZE_SMALL(*n, *n, au, ldau, tauq, taup, s, e);

    /* Form Vt' in vt from HH vectors in U (right bi-diagonalizing Q) */
    if(wntvs)
    {
        ncvt = *n;
        for(i = 1; i <= *n; i++)
            for(j = 1; j <= *n; j++)
                vt[i + j * *ldvt] = 0.;
        FLA_LARF_VTAPPLY_DSMALL_SQR(n, au, ldau, taup, vt, ldvt);
    }

    /* Form U' in U (left bi-diagonalizing Q) */
    if(wntus)
    {
        FLA_LARF_UAPPLY_DSMALL_SQR(n, au, ldau, tauq, u, ldu, taup);
    }

    /* Compute SVD for bi-diagonal matrix
     * (dbdsqr with no lwork)
     * */
    if(*n == 2)
    {
        /* 2 by 2 block, handle separately */
        doublereal sigmn, sigmx, sinr, cosr, sinl, cosl;
        doublereal scl1, scl2;

        dlasv2_(&s[1], &e[1], &s[2], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);
        scl1 = (sigmx < 0.) ? -1. : 1.;
        scl2 = (sigmn < 0.) ? -1. : 1.;
        s[1] = f2c_abs(sigmx);
        s[2] = f2c_abs(sigmn);
        /* Compute singular vectors, if desired */
        if(ncvt > 0)
        {
            vt[1 + *ldvt] = scl1 * cosr;
            vt[2 + *ldvt] = scl1 * sinr;

            vt[1 + 2 * *ldvt] = scl2 * -sinr;
            vt[2 + 2 * *ldvt] = scl2 * cosr;
        }
        if(nru > 0)
        {
            fla_drot_avx2(&nru, &u[1 + *ldu], &c__1, &u[1 + 2 * *ldu], &c__1, &cosl, &sinl);
        }
    }
    else
    {
        lapack_dbdsqr_small("U", n, &ncvt, &nru, &s[1], &e[1], &vt[1 + *ldvt], ldvt, &u[1 + *ldu],
                            ldu, info);
    }

    /* Compute U by updating U' by applying from the left the Q from QR */
    if(wntus)
    {
        tau = &work[itau - 1];
        /* First Iteration corresponding to HH(n) */
        i = *n;
        for(j = 1; j <= *n; j++)
        {
            /* - u[i][j] * tau[i] */
            d__1 = -u[i + j * *ldu] * tau[i];

            /* u[n+1:m, j] = d__1 * u[n+1:m, j] */
            for(k = *n + 1; k <= *m; k++)
            {
                u[k + j * *ldu] = d__1 * qr[k + *n * *ldqr];
            }
        }
        /* u[m, 1:m] = u[m, 1:m] * (1 - tau) */
        d__1 = 1 - tau[i];
        for(j = 1; j <= *n; j++)
        {
            u[*n + j * *ldu] = u[*n + j * *ldu] * d__1;
        }

        /* Second Iteration onwards */
        beta = 0;
        xnorm = 1.;
        for(i = *n - 1; i >= 1; i--)
        {
            /* incrementing n by i to compensate for decrement
             * by i done in FLA_LARF_APPLY_DLARGE_COL
             */
            ni = *n + i;

            au = &u[-i * *ldu];
            v = &qr[i + i * *ldqr - 1];
            FLA_LARF_APPLY_DLARGE_COL(i, m, &ni, au, ldu, tau);
        }
    }

    return;
}
#endif
