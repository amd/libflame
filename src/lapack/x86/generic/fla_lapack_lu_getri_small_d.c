/******************************************************************************
 * * Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/
/*! @file fla_lapack_lu_getri_small_d.c
 *  @brief getri kernels for small inputs
 *         to choose optimized paths
 *  *  */

#include "FLA_f2c.h"
#include "fla_lapack_lu_small_kernels_d.h"

#if FLA_ENABLE_AMD_OPT
static inline void lapack_getri_small_unblocked_d(aocl_int64_t n, doublereal *restrict a, aocl_int64_t a_dim1,
                                                  aocl_int_t *restrict ipiv, doublereal *restrict work)
{
    aocl_int64_t i, j, k;
    for(i = 1; i <= n; ++i)
    {
        if(a[i + i * a_dim1] == 0.0)
        {
            return;
        }
    }
    if(n == 4)
    {
        LAPACK_GETRI_TRIANGULAR_INVERSION_D_4X4(a, a_dim1);
    }
    else if(n == 5)
    {
        LAPACK_GETRI_TRIANGULAR_INVERSION_D_5X5(a, a_dim1);
    }
    else
    {
        /* Generic triangular inversion */
        for(j = 1; j <= n; ++j)
        {
            doublereal inv_djj = 1.0 / a[j + j * a_dim1];
            a[j + j * a_dim1] = inv_djj;
            if(j > 1)
            {
                for(i = 1; i <= j - 1; ++i)
                {
                    doublereal sum = 0.0;
                    for(k = i; k <= j - 1; ++k)
                    {
                        sum += a[i + k * a_dim1] * a[k + j * a_dim1];
                    }
                    a[i + j * a_dim1] = -inv_djj * sum;
                }
            }
        }
    }
    /* Lower solve phase, macro-based per size */
    if(n == 4)
    {
        LAPACK_GETRI_LOWER_SOLVE_D_4X4(a, a_dim1, work);
    }
    else if(n == 5)
    {
        LAPACK_GETRI_LOWER_SOLVE_D_5X5(a, a_dim1, work);
    }
    else if(n == 6)
    {
        LAPACK_GETRI_LOWER_SOLVE_D_6X6(a, a_dim1, work);
    }
    else
    {
        /* Use generic lower-solve for sizes 7-9 */
        for(j = n; j >= 1; --j)
        {
            for(i = j + 1; i <= n; ++i)
            {
                work[i] = a[i + j * a_dim1];
                a[i + j * a_dim1] = 0.0;
            }
            if(j < n)
            {
                for(i = 1; i <= n; ++i)
                {
                    doublereal acc = 0.0;
                    for(k = j + 1; k <= n; ++k)
                    {
                        acc += a[i + k * a_dim1] * work[k];
                    }
                    a[i + j * a_dim1] -= acc;
                }
            }
        }
    }
    /* Column interchanges (shared) */
    for(j = n - 1; j >= 1; --j)
    {
        aocl_int64_t jp = ipiv[j];
        if(jp != j)
        {
            doublereal *apiv = &a[1 + jp * a_dim1];
            doublereal *asrc = &a[1 + j * a_dim1];
            for(i = 1; i <= n; ++i)
            {
                doublereal t = apiv[i - 1];
                apiv[i - 1] = asrc[i - 1];
                asrc[i - 1] = t;
            }
        }
    }
    work[1] = (doublereal)n;
}

void lapack_getri_small_d(aocl_int64_t *n, doublereal *a, aocl_int64_t *lda, aocl_int_t *ipiv, doublereal *work,
                          aocl_int64_t *info)
{
    aocl_int64_t a_dim1, i__1, i__;
    doublereal t_val, *acur, *apiv, *asrc, *ax, *ay;
    aocl_int64_t i, j, jp;

    a_dim1 = *lda;
    /* check singularity of triangular input matrix  */
    for(*info = 1; *info <= *n; ++(*info))
    {
        if(a[*info + *info * a_dim1] == 0.)
        {
            return;
        }
    }
    *info = 0;
    if(*n == 2)
    {
        LAPACK_GETRI_SMALL_D_2x2(n, a, a_dim1, ipiv, work);
    }
    else if(*n == 3)
    {
        LAPACK_GETRI_SMALL_D_3x3(n, a, a_dim1, ipiv, work);
    }
    else if(*n >= 4 && *n <= 9)
    {
        lapack_getri_small_unblocked_d(*n, a, a_dim1, ipiv, work);
    }
    else
    {
        /* Compute inverse of upper triangular matrix. */
        a = a + (a_dim1 + 1);
        a[0] = 1. / a[0];
        for(j = 1; j < *n; ++j)
        {
            a[j + j * a_dim1] = 1. / a[j + j * a_dim1];
            /* trmv inlined code */
            a[j * a_dim1] = a[j * a_dim1] * a[0];
            for(i = 1; i < j; ++i)
            {
                for(i__ = 0; i__ < i; ++i__)
                {
                    a[i__ + j * a_dim1]
                        = a[i__ + j * a_dim1] + (a[i + j * a_dim1] * a[i__ + i * a_dim1]);
                }
                a[i + j * a_dim1] = a[i + j * a_dim1] * a[i + i * a_dim1];
            }
            /* scal */
            for(i = 0; i < j; ++i)
            {
                a[i + j * a_dim1] = -a[j + j * a_dim1] * a[i + j * a_dim1];
            }
        }
        a = a - (a_dim1 + 1);
        /* Use unblocked code. */
        for(j = *n; j >= 1; --j)
        {
            /* Copy current column of L to WORK and replace with zeros. */
            for(i = j + 1; i <= *n; ++i)
            {
                work[i] = a[i + j * a_dim1];
                a[i + j * a_dim1] = 0.;
            }
            /* Compute current column of inv(A). */
            if(j < *n)
            {
                i__1 = *n - j;
                /* gemv inline */
                acur = &a[j * a_dim1];
                ax = &work[j];
                ay = &a[j * a_dim1];

                for(i = 1; i <= i__1; ++i)
                {
                    for(i__ = 1; i__ <= *n; ++i__)
                    {
                        ay[i__] = ay[i__] - (ax[i] * acur[i__ + i * a_dim1]);
                    }
                }
            }
        }

        /* Apply column interchanges. */
        for(j = *n - 1; j >= 1; --j)
        {
            jp = ipiv[j];
            if(jp != j)
            {
                /*  swap */
                apiv = &a[jp * a_dim1 + 1];
                asrc = &a[j * a_dim1 + 1];
                for(i = 0; i < *n; i++)
                {
                    t_val = apiv[i];
                    apiv[i] = asrc[i];
                    asrc[i] = t_val;
                }
            }
        }
    }
    work[1] = (doublereal)*n;
    return;
}
#endif
