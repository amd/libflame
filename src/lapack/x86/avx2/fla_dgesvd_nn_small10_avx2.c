/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_nx_small10_avx2.c
 *  @brief DGESVD Small path (Path 10)
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

#ifndef FLA_ENABLE_AOCL_BLAS
extern void drot_(integer *, doublereal *, integer *, doublereal *,
                  integer *, doublereal *, doublereal *);
#endif
extern void dlartg_(doublereal *da, doublereal *db, doublereal *c__, doublereal *s, doublereal *r);

void fla_dgesvd_xx_small10_avx2(integer wntus, integer wntvs, integer *m, integer *n, doublereal *a,
                                integer *lda, doublereal *s, doublereal *u, integer *ldu,
                                doublereal *vt, integer *ldvt, doublereal *work, integer *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    doublereal d__1;
    doublereal *tau, *tauq, *taup;
    doublereal *e;
    doublereal stau;

    integer ncvt, nru;
    integer c__1 = 1;

    integer ie;
    integer itauq, itaup;
    integer i__1, rlen, knt;

    /* indices for partitioning work buffer */
    ie = 1;
    itauq = ie + *n;
    itaup = itauq + *n;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    if(*m == 2 && *n == 2)
    {
        /* 2x2 matrix Bi-Diag using Givens */
        doublereal s0, cosr, sinr;

        dlartg_(&a[1 + *lda], &a[2 + *lda], &cosr, &sinr, &s0);
        s[1] = s0;

        /* Update 2nd columns of A */
        dtmp = cosr * a[1 + 2 * *lda] + sinr * a[2 + 2 * *lda];
        a[2 + 2 * *lda] = cosr * a[2 + 2 * *lda] - sinr * a[1 + 2 * *lda];
        a[1 + 2 * *lda] = dtmp;

        /* Update Singular values and vectors */
        s[2] = a[2 + 2 * *lda];
        e[1] = a[1 + 2 * *lda];
        if(wntus)
        {
            u[1 + 1 * *ldu] = cosr;
            u[2 + 1 * *ldu] = -sinr;
            u[1 + 2 * *ldu] = sinr;
            u[2 + 2 * *ldu] = cosr;
        }
    }
    else
    {
        FLA_BIDIAGONALIZE_SMALL(*m, *n, a, lda, tauq, taup, s, e);

        /* Generate Qr (from bidiag) in vt from work[iu] (a here) */
        if(wntvs)
        {
            for(i = 1; i <= *n; i++)
                for(j = 1; j <= *n; j++)
                    vt[i + j * *ldvt] = 0.;
            FLA_LARF_VTAPPLY_DSMALL_SQR(n, a, lda, taup, vt, ldvt);
        }
        /* Generate Ql (from bidiag) in u from a */
        if(wntus)
        {
            /* for all HH vectors from the end */
            for(i = *n; i >= 1; i--)
            {
                /* Update current column */
                stau = -tauq[i];
                for(j = i + 1; j <= *m; j++)
                {
                    u[j + i * *ldu] = stau * a[j + i * *lda];
                }
                u[i + i * *ldu] = 1 + stau;

                /* Update rest of the columns from (i + 1) to n */
                for(k = i + 1; k <= *n; k++)
                {
                    dtmp = 0.;
                    for(j = i + 1; j <= *m; j++)
                    {
                        dtmp = dtmp + a[j + i * *lda] * u[j + k * *ldu];
                    }
                    dtmp = stau * dtmp;

                    for(j = i + 1; j <= *m; j++)
                    {
                        u[j + k * *ldu] = u[j + k * *ldu] + dtmp * a[j + i * *lda];
                    }
                    u[i + k * *ldu] = dtmp;
                }
            }
        }
    }

    /* Compute final Singular Values/Vectors */
    ncvt = 0;
    nru = 0;
    if(wntvs)
    {
        ncvt = *n;
    }
    if(wntus)
    {
        nru = *m;
    }
    if(*m == 2 && *n == 2)
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
            drot_(&nru, &u[1 + *ldu], &c__1, &u[1 + 2 * *ldu], &c__1, &cosl, &sinl);
        }
    }
    else
    {
        lapack_dbdsqr_small("U", n, &ncvt, &nru, &s[1], &e[1], &vt[1 + *ldvt], ldvt, &u[1 + *ldu],
                            ldu, info);
    }
    return;
}
#endif
