/******************************************************************************
 * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/
#ifndef FLA_DGEQRF_SMALL_AVX2_DEFS_H
#define FLA_DGEQRF_SMALL_AVX2_DEFS_H

/*! @file fla_dgeqrf_small_avx2.h
 *  @brief QR Kernels for small sizes.
 *  */

#if FLA_ENABLE_AMD_OPT


#include "fla_geqrf_small_avx2_kernel.h"

/* Application of Givens Rotation ** T)
 * over rows row & row + 1
 * from the left */
#define FLA_APPLY_GIVENS_LVX(idim, imat, ldi, row, cs, sn)   \
    {                                                        \
        aocl_int64_t im;                                     \
        doublereal tv0, tv1;                                 \
        for(im = 1; im <= *idim; im++)                       \
        {                                                    \
            tv0 = imat[row + 0 + im * *ldi];                 \
            tv1 = imat[row + 1 + im * *ldi];                 \
                                                             \
            imat[row + 0 + im * *ldi] = cs * tv0 + sn * tv1; \
            imat[row + 1 + im * *ldi] = cs * tv1 - sn * tv0; \
        }                                                    \
    }
/* Application of Givens Rotation ** T)
 * over columns col & col + 1
 * from the right */
#define FLA_APPLY_GIVENS_RVX(idim, imat, ldi, col, cs, sn)     \
    {                                                          \
        aocl_int64_t im;                                       \
        doublereal tv0, tv1;                                   \
        for(im = 1; im <= *idim; im++)                         \
        {                                                      \
            tv0 = imat[im + (col + 0) * *ldi];                 \
            tv1 = imat[im + (col + 1) * *ldi];                 \
                                                               \
            imat[im + (col + 0) * *ldi] = cs * tv0 + sn * tv1; \
            imat[im + (col + 1) * *ldi] = cs * tv1 - sn * tv0; \
        }                                                      \
    }

#define FLA_LARF_GEN_DSMALL_ROW(i, m, n, iptr, ldia, tau)              \
    /* Compute norm2 */                                                \
    xnorm = aocl_blas_dnrm2(&rlen, &iptr[2 * *ldia], ldia);            \
    if(xnorm == 0.)                                                    \
    {                                                                  \
        tau[i] = 0.;                                                   \
        beta = iptr[*ldia];                                            \
    }                                                                  \
    else                                                               \
    {                                                                  \
        knt = 0;                                                       \
        v = iptr;                                                      \
        alpha = v[*ldia];                                              \
        d__1 = dlapy2_(&v[*ldia], &xnorm);                             \
        beta = -d_sign(&d__1, &alpha);                                 \
        if(f2c_abs(beta) < safmin)                                     \
        {                                                              \
            for(knt = 0; f2c_abs(beta) < safmin && knt < 20; knt++)    \
            {                                                          \
                aocl_blas_dscal(&rlen, &rsafmin, &v[2 * *ldia], ldia); \
                beta *= rsafmin;                                       \
                alpha *= rsafmin;                                      \
            }                                                          \
            /* New BETA is at most 1, at least SAFMIN */               \
            xnorm = aocl_blas_dnrm2(&rlen, &v[2 * *ldia], ldia);       \
            d__1 = dlapy2_(&alpha, &xnorm);                            \
            beta = -d_sign(&d__1, &alpha);                             \
        }                                                              \
        tau[i] = (beta - alpha) / beta;                                \
        d__1 = 1. / (alpha - beta);                                    \
        aocl_blas_dscal(&rlen, &d__1, &v[2 * *ldia], ldia);            \
        for(j = 1; j <= knt; ++j)                                      \
        {                                                              \
            beta *= safmin;                                            \
        }                                                              \
    }

#define FLA_LARF_APPLY_DSMALL_ROW(i, m, n, iptr, ldia, tau)                  \
    if(xnorm == 0.)                                                          \
    {                                                                        \
        tau[i] = 0.;                                                         \
    }                                                                        \
    else                                                                     \
    {                                                                        \
        /* for every row ac of A(i+1:nr,i+1:nc) */                           \
        ac = iptr;                                                           \
        v[*ldia] = 1;                                                        \
        for(j = 1; j <= slen; j++)                                           \
        {                                                                    \
            dtmp = 0;                                                        \
            /* w = (ac .* v) */                                              \
            for(k = 1; k <= rlen + 1; k++)                                   \
            {                                                                \
                dtmp = dtmp + ac[j + k * *ldia] * v[k * *ldia];              \
            }                                                                \
                                                                             \
            /* (ac .* v) * tau */                                            \
            dtmp = dtmp * tau[i];                                            \
                                                                             \
            /* ac = ac - ac * dtmp */                                        \
            for(k = 1; k <= rlen + 1; k++)                                   \
            {                                                                \
                ac[j + k * *ldia] = ac[j + k * *ldia] - v[k * *ldia] * dtmp; \
            }                                                                \
        }                                                                    \
        v[*ldia] = beta;                                                     \
    }

#define FLA_LARF_UAPPLY_DSMALL_SQR(m, a, lda, tauq, u, ldu, twork)              \
    if(*m > 1)                                                                  \
    {                                                                           \
        /* iteration corresponding to (m - 1) HH(m-1) */                        \
        stau = tauq[*m - 1];                                                    \
        d__1 = a[*m + (*m - 1) * *lda];                                         \
        dtmp = -(stau * d__1);                                                  \
                                                                                \
        u[*m - 1 + (*m - 1) * *ldu] = 1.0 - stau; /* 1 - tau */                 \
        u[*m + (*m - 1) * *ldu] = dtmp; /* tau * v2 */                          \
        u[*m - 1 + *m * *ldu] = dtmp; /* tau * v2 */                            \
        u[*m + *m * *ldu] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */           \
    }                                                                           \
    else                                                                        \
    {                                                                           \
        u[1 + *ldu] = 1.0;                                                      \
    }                                                                           \
    for(i = *m - 2; i >= 1; i--)                                                \
    {                                                                           \
        stau = -tauq[i];                                                        \
                                                                                \
        /* scale col i by -tau and dlarf for rest of the columns */             \
        for(j = i + 1; j <= *m; j++)                                            \
        {                                                                       \
            twork[j] = a[j + i * *lda];                                         \
                                                                                \
            /* GEMV part of dlarf excluding zero first row */                   \
            dtmp = 0;                                                           \
            for(k = i + 1; k <= *m; k++)                                        \
            {                                                                   \
                dtmp = dtmp + u[k + j * *ldu] * a[k + i * *lda];                \
            }                                                                   \
            u[i + j * *ldu] = stau * dtmp;                                      \
        }                                                                       \
        u[i + i * *ldu] = 1.0 + stau;                                           \
                                                                                \
        /* Update all columns except current column*/                           \
        for(j = i + 1; j <= *m; j++)                                            \
        {                                                                       \
            for(k = i + 1; k <= *m; k++)                                        \
            {                                                                   \
                u[k + j * *ldu] = u[k + j * *ldu] + twork[k] * u[i + j * *ldu]; \
            }                                                                   \
        }                                                                       \
        /* Updating the current column */                                       \
        for(j = i + 1; j <= *m; j++)                                            \
        {                                                                       \
            u[j + i * *ldu] = stau * a[j + i * *lda];                           \
        }                                                                       \
    }

#define FLA_LARF_VTAPPLY_DSMALL_SQR(m, a, lda, taup, vt, ldvt)                           \
    if(*m > 2)                                                                           \
    {                                                                                    \
        /* iteration corresponding to (m - 2) HH[m-2] */                                 \
        stau = taup[*m - 2];                                                             \
        d__1 = a[*m - 2 + *m * *lda];                                                    \
        dtmp = -(stau * d__1); /* tau * v2 */                                            \
                                                                                         \
        vt[*m - 1 + (*m - 1) * *ldvt] = 1.0 - stau; /* 1 - tau */                        \
        vt[*m + (*m - 1) * *ldvt] = dtmp; /* tau * v2 */                                 \
        vt[*m - 1 + *m * *ldvt] = dtmp; /* tau * v2 */                                   \
        vt[*m + *m * *ldvt] = 1.0 + (dtmp * d__1); /* 1 - tau * v2^2 */                  \
                                                                                         \
        /* for HH vectors [m-3:1] */                                                     \
        for(i = *m - 3; i >= 1; i--)                                                     \
        {                                                                                \
            stau = -taup[i];                                                             \
                                                                                         \
            /* Compute tmp = tau * v' * in */                                            \
            for(j = i + 2; j <= *m; j++)                                                 \
            {                                                                            \
                /* Scale row i by -tau and dlarf for rest of the rows */                 \
                vt[i + 1 + j * *ldvt] = stau * a[i + j * *lda];                          \
                                                                                         \
                /* DOT part of the dlarf excluding zero first column */                  \
                dtmp = 0.;                                                               \
                for(k = i + 2; k <= *m; k++)                                             \
                {                                                                        \
                    dtmp = dtmp + vt[j + k * *ldvt] * a[i + k * *lda];                   \
                }                                                                        \
                vt[j + (i + 1) * *ldvt] = stau * dtmp;                                   \
            }                                                                            \
            vt[i + 1 + (i + 1) * *ldvt] = 1.0 + stau;                                    \
                                                                                         \
            /* compute (in - v * tmp */                                                  \
            for(j = i + 2; j <= *m; j++)                                                 \
            {                                                                            \
                for(k = i + 2; k <= *m; k++)                                             \
                {                                                                        \
                    vt[j + k * *ldvt]                                                    \
                        = vt[j + k * *ldvt] + a[i + k * *lda] * vt[j + (i + 1) * *ldvt]; \
                }                                                                        \
            }                                                                            \
        }                                                                                \
    }                                                                                    \
    else                                                                                 \
    {                                                                                    \
        for(i = 1; i <= *m; i++)                                                         \
        {                                                                                \
            vt[i + i * *ldvt] = 1.;                                                      \
        }                                                                                \
    }                                                                                    \
    vt[1 + *ldvt] = 1.;

#define FLA_LARF_VTAPPLY_DSMALL_ROW(i, m, n, tau, sv, ldsv)             \
    /* for every row ac of A(i+1:nr,i+1:nc) */                          \
    v[*lda] = 1;                                                        \
    for(j = 1; j <= slen; j++)                                          \
    {                                                                   \
        dtmp = 0;                                                       \
        /* w = (ac .* v) */                                             \
        for(k = 1; k <= rlen + 1; k++)                                  \
        {                                                               \
            dtmp = dtmp + sv[j + k * *ldsv] * v[k * *lda];              \
        }                                                               \
                                                                        \
        /* (ac .* v) * tau */                                           \
        dtmp = dtmp * tau[i];                                           \
                                                                        \
        /* ac = ac - ac * dtmp */                                       \
        for(k = 1; k <= rlen + 1; k++)                                  \
        {                                                               \
            sv[j + k * *ldsv] = sv[j + k * *ldsv] - v[k * *lda] * dtmp; \
        }                                                               \
    }                                                                   \
    v[*lda] = beta;

#endif /* FLA_ENABLE_AMD_OPT */
#endif /* FLA_DGEQRF_SMALL_AVX2_DEFS_H */
