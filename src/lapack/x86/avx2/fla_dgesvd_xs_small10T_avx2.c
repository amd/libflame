/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_small6_avx2_.c
 *  @brief DGESVD Small path (path 6)
 *  without the LQ Factorization.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* SVD for small fat-matrices
 */
void fla_dgesvd_xs_small10T_avx2(integer *m, integer *n, doublereal *a, integer *lda, doublereal *s,
                                 doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
                                 doublereal *work, integer *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    integer ie;
    integer itauq, itaup;
    integer i__1, rlen, knt;
    integer tm, tn;
    integer c__1 = 1;

    doublereal *tau, *tauq, *taup;
    doublereal *e;
    doublereal *iptr;
    doublereal stau, d__1;

    doublereal *ta, *ts;

    /* indices for partitioning work buffer */
    ie = 1;
    itauq = ie + *m;
    itaup = itauq + *m;

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

    /* Lower Bidiagonalization */
    {
        /* Annihilate first row elements to the right of the diagonal */
        rlen = *n - 1;
        slen = *m - 1;
        iptr = a + 1;
        tau = taup;
        FLA_LARF_GEN_DSMALL_ROW(1, m, n, tau);
        s[1] = beta;
        FLA_LARF_APPLY_DSMALL_ROW(1, m, n, tau);

        /* Upper Bidiagonalize the matrix excluding the first row */
        tm = *m - 1;
        tn = *n;
        ta = a + 1;
        tau = taup + 1;
        ts = s + 1;
        FLA_BIDIAGONALIZE_SMALL(tm, tn, ta, lda, tauq, tau, e, ts);
    }

    /* Generate Qr (from bidiag) in vt */
    xnorm = 1.0;
    for(i = *m; i >= 1; i--)
    {
        /* Update current row */
        for(j = 1; j <= i - 1; j++)
        {
            vt[i + j * *ldvt] = 0.;
        }
        vt[i + i * *ldvt] = 1 - taup[i];
        for(j = i + 1; j <= *n; j++)
        {
            vt[i + j * *ldvt] = -taup[i] * a[i + j * *lda];
        }

        /* Update rows below current row using row-apply */
        v = &a[i + i * *lda - *lda];
        beta = v[*lda];
        slen = *m - i;
        rlen = *n - i;
        ta = &vt[i + i * *ldvt - *ldvt];
        FLA_LARF_VTAPPLY_DSMALL_ROW(i, m, n, taup, ta, ldvt);
    }

    /* Generate Ql (from bidiag) in u from a */
    u[1 + *ldu] = 1.;
    if(*m > 2)
    {
        /* iteration corresponding to (m - 2) HH(m-2) */
        i = *m - 2;
        stau = tauq[i];
        d__1 = a[*m + i * *lda];
        dtmp = -(stau * d__1);

        u[*m - 1 + (*m - 1) * *ldu] = 1.0 - stau; /* 1 - tau */
        u[*m + (*m - 1) * *ldu] = dtmp; /* tau * v2 */
        u[*m - 1 + *m * *ldu] = dtmp; /* tau * v2 */
        u[*m + *m * *ldu] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */
    }
    else if(*m > 1)
    {
        u[1 + *ldu] = 1.0;
        u[2 + *ldu] = 0;
        u[1 + 2 * *ldu] = 0;
        u[2 + 2 * *ldu] = 1.0;
    }
    else
    {
        u[1 + *ldu] = 1.0;
    }
    /* for HH vectors [m-3:1] */
    for(i = *m - 3; i >= 1; i--)
    {
        stau = -tauq[i];
        /* scale col (i + 1) by -tau and dlarf for rest of the columns */
        for(j = i + 2; j <= *m; j++)
        {
            u[j + (i + 1) * *ldu] = stau * a[j + i * *lda];
        }
        /* Columns (i + 2) to m */
        for(j = i + 2; j <= *m; j++)
        {
            /* GEMV part of dlarf excluding zero first row .
               Store the dot product in u.
            */
            dtmp = 0;
            for(k = i + 2; k <= *m; k++)
            {
                dtmp = dtmp + u[k + j * *ldu] * a[k + i * *lda];
            }
            u[i + 1 + j * *ldu] = stau * dtmp;
        }
        u[i + 1 + (i + 1) * *ldu] = 1.0 + stau;

        for(j = i + 2; j <= *m; j++)
        {
            for(k = i + 2; k <= *m; k++)
            {
                u[k + j * *ldu] = u[k + j * *ldu] + a[k + i * *lda] * u[i + 1 + j * *ldu];
            }
        }
    }

    lapack_dbdsqr_small("L", m, n, m, &s[1], &e[1], &vt[1 + *ldvt], ldvt, &u[1 + *ldu], ldu, info);

    return;
}
#endif
