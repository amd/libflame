/******************************************************************************
 * * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 * *******************************************************************************/
/*! @file fla_lapack_lu_getri_small_d.c
 *  @brief getri kernals for small inputs
 *         to choose optimized paths
 *  *  */
#if FLA_ENABLE_AMD_OPT
#include "FLA_f2c.h"
#include "fla_lapack_lu_small_kernals_d.h"

void lapack_getri_small_d(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work,
                          integer *info)
{
    integer a_dim1, i__1, i__;
    doublereal t_val, *acur, *apiv, *asrc, *ax, *ay;
    integer i, j, jp;

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
    if(*n == 3)
    {
        LAPACK_GETRI_SMALL_D_3x3(n, a, a_dim1, ipiv, work);
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