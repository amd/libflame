/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_xx_small10_avx2.c
 *  @brief DGESVD Small path (Path 10)
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

extern void dlartg_(doublereal *da, doublereal *db, doublereal *c__, doublereal *s, doublereal *r);

void fla_dgesvd_xx_small10_avx2(aocl_int64_t wntu, aocl_int64_t wntv, aocl_int64_t *m,
                                aocl_int64_t *n, aocl_int64_t *ncu, doublereal *a,
                                aocl_int64_t *lda, doublereal *s, doublereal *u, aocl_int64_t *ldu,
                                doublereal *vt, aocl_int64_t *ldvt, doublereal *work,
                                aocl_int64_t *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    doublereal d__1;
    doublereal *tau, *tauq, *taup;
    doublereal *e;
    doublereal stau;
    doublereal cosu = 0., sinu = 0.;

    aocl_int64_t ncvt, nru;
    aocl_int64_t c__1 = 1;

    aocl_int64_t ie;
    aocl_int64_t itauq, itaup;
    aocl_int64_t rlen, knt;

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
        doublereal s0;

        dlartg_(&a[1 + *lda], &a[2 + *lda], &cosu, &sinu, &s0);
        s[1] = s0;

        /* Update 2nd columns of A */
        dtmp = cosu * a[1 + 2 * *lda] + sinu * a[2 + 2 * *lda];
        a[2 + 2 * *lda] = cosu * a[2 + 2 * *lda] - sinu * a[1 + 2 * *lda];
        a[1 + 2 * *lda] = dtmp;

        /* Update Singular values and vectors */
        s[2] = a[2 + 2 * *lda];
        e[1] = a[1 + 2 * *lda];
    }
    else
    {
        FLA_BIDIAGONALIZE_SMALL(*m, *n, a, lda, tauq, taup, s, e);

        /* Generate Qr (from bidiag) in vt from work[iu] (a here) */
        if(wntv)
        {
            for(i = 1; i <= *n; i++)
                for(j = 1; j <= *n; j++)
                    vt[i + j * *ldvt] = 0.;
            FLA_LARF_VTAPPLY_DSMALL_SQR(n, a, lda, taup, vt, ldvt);
        }
        /* Generate Ql (from bidiag) in u from a */
        if(wntu)
        {
            /* Initialize columns n to ncu of U to eye */
            for(i = *n + 1; i <= *ncu; i++)
            {
                for(j = *n + 1; j <= *ncu; j++)
                {
                    u[i + j * *ldu] = 0.;
                }
                u[i + i * *ldu] = 1.;
            }
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
                for(k = i + 1; k <= *ncu; k++)
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
    if(wntv)
    {
        ncvt = *n;
    }
    if(wntu)
    {
        nru = *m;
    }
    if(*m == 2 && *n == 2)
    {
        /* 2 by 2 block, handle separately */
        doublereal sigmn, sigmx, sinr, cosr, sinl, cosl;

        dlasv2_(&s[1], &e[1], &s[2], &sigmn, &sigmx, &sinr, &cosr, &sinl, &cosl);
        s[1] = f2c_abs(sigmx);
        s[2] = f2c_abs(sigmn);
        /* Compute singular vectors, if desired */
        if(ncvt > 0)
        {
            FLA_COMPUTE_VT_2X2(vt, ldvt, sigmx, sigmn, cosr, sinr);
        }
        if(nru > 0)
        {
            doublereal p0, p1, p2, p3;

            p0 = cosl * cosu;
            p1 = sinl * sinu;
            p2 = sinl * cosu;
            p3 = cosl * sinu;

            u[1 + *ldu] = p0 - p1;
            u[2 + *ldu] = p2 + p3;

            u[1 + 2 * *ldu] = -(p3 + p2);
            u[2 + 2 * *ldu] = p0 - p1;
        }

        /* Normalize singular values and scale corresponding vectors for 2x2 case */
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(1, wntu);
        FLA_NORMALIZE_SINGULAR_VALUE_AND_VECTORS_2X2(2, wntu);
    }
    else
    {
        /* Compute Singular Values and Vectors */
        lapack_dbdsqr_small("U", n, &ncvt, &nru, &s[1], &e[1], &vt[1 + *ldvt], ldvt,
                            &u[1 + *ldu], ldu, info);
    }
    return;
}
#endif
