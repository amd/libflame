/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/* > \brief \b validate_gelsd.c                                              */
/* =========== DOCUMENTATION ===========                                     */
/* Definition:                                                               */
/* ===========                                                               */
/* SUBROUTINE validate_gelsd(integer m, integer n, integer NRHS, void *A,    */
/*                           integer lda, void *B, integer ldb, void *S,     */
/*                           void *X, void *rcond, integer *rank,            */
/*                           integer datatype, double* residual);            */
/* > \par Purpose:                                                           */
/* =============                                                             */
/* >                                                                         */
/* > \verbatim                                                               */
/* >                                                                         */
/* > Defines validate function of GELSD() to use in test suite               */
/* >                                                                         */
/* > This API validates GELSD() functionality by performing                  */
/* > the folllowing tests:                                                   */
/* >                                                                         */
/* > TEST 1: If m >= n and RANK = n, the residual sum-of-squares             */
/* >         for the solution in the i-th column is given by the sum of      */
/* >         squares of elements n+1:m in that column                        */
/* >                                                                         */
/* > TEST 2: Check AX = B                                                    */
/* >                                                                         */
/* > TEST 3 (TODO): checks whether X is in the row space of A or A'          */
/* >                                                                         */
/* > Test SVD: Check if Singular values generated are in order               */
/* >                                                                         */
/* > \endverbatim                                                            */
/* ========================================================================= */

#include "test_common.h"

void validate_gelsd(integer m, integer n, integer nrhs, void *A, integer lda, void *B, integer ldb,
                    void *S, void *X, void *rcond, integer *rank, integer datatype,
                    double *residual)
{
    if(m == 0 || n == 0 || nrhs == 0)
        return;

    void *work = NULL;
    double resid1;
    integer ldx;
    ldx = ldb;

    /* Test SVD
       Check the order of Singular values generated */
    resid1 = svd_check_order(datatype, S, m, n, *residual);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_slamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column */
            if((m > n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, residual);
            }
            else
            {
                /* Test 2 */
                /* Compute |B-AX| = 0 */
                norm_a = fla_lapack_slange("1", &m, &n, A, &lda, work);
                norm_b = fla_lapack_slange("1", &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_slange("1", &n, &nrhs, X, &ldx, work);

                /* Compute B-AX */
                sgemm_("N", "N", &m, &nrhs, &n, &s_n_one, A, &lda, X, &ldx, &s_one, B, &ldb);
                norm = fla_lapack_slange("1", &m, &nrhs, B, &ldb, work);

                *residual = (double)(norm / ((norm_a * norm_x + norm_b) * (float)fla_max(m, n) * eps));
            }
            break;
        }
        case DOUBLE:
        {
            double norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_dlamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column*/
            if((m > n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, residual);
            }
            else
            {
                /* Compute |B-AX| = 0 */
                norm_a = fla_lapack_dlange("1", &m, &n, A, &lda, work);
                norm_b = fla_lapack_dlange("1", &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_dlange("1", &n, &nrhs, X, &ldx, work);

                /* Compute B-AX */
                dgemm_("N", "N", &m, &nrhs, &n, &d_n_one, A, &lda, X, &ldx, &d_one, B, &ldb);
                norm = fla_lapack_dlange("1", &m, &nrhs, B, &ldb, work);

                *residual = norm / ((norm_a * norm_x + norm_b) * (double)fla_max(m, n) * eps);
            }
            break;
        }
        case COMPLEX:
        {
            float norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_slamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column*/
            if((m > n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, residual);
            }
            else
            {
                /* Compute |B-AX| = 0 */
                norm_a = fla_lapack_clange("1", &m, &n, A, &lda, work);
                norm_b = fla_lapack_clange("1", &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_clange("1", &n, &nrhs, X, &ldx, work);
                eps = fla_lapack_slamch("E");

                /* Compute B-AX */
                cgemm_("N", "N", &m, &nrhs, &n, &c_n_one, A, &lda, X, &ldx, &c_one, B, &ldb);
                norm = fla_lapack_clange("1", &m, &nrhs, B, &ldb, work);

                *residual = (double)(norm / ((norm_a * norm_x + norm_b) * (float)fla_max(m, n) * eps));
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_dlamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column*/
            if((m > n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, residual);
            }
            else
            {
                /* Compute |B-AX| = 0 */
                norm_a = fla_lapack_zlange("1", &m, &n, A, &lda, work);
                norm_b = fla_lapack_zlange("1", &m, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_zlange("1", &n, &nrhs, X, &ldx, work);
                eps = fla_lapack_dlamch("E");

                /* Compute B-AX */
                zgemm_("N", "N", &m, &nrhs, &n, &z_n_one, A, &lda, X, &ldx, &z_one, B, &ldb);
                norm = fla_lapack_zlange("1", &m, &nrhs, B, &ldb, work);

                *residual = norm / ((norm_a * norm_x + norm_b) * (double)fla_max(m, n) * eps);
            }
            break;
        }
    }
    *residual = fla_max(*residual, resid1);
}