/******************************************************************************
 * Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_potrs.c
 *  @brief Defines validate function of POTRS() to use in test suite.
 *  */

#include "test_common.h"

void validate_potrs(integer n, integer nrhs, void *A, integer lda, void *X, void *B, integer ldb,
                    integer datatype, double *residual, integer *info, char imatrix)
{
    if(n == 0 || nrhs == 0)
        return;
    void *work = NULL;
    integer ldx;
    *info = 0;
    ldx = n;
    char NORM = '1';

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_a, norm_b, norm_x, norm, eps, resid;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            /* Compute AX-B */
            sgemm_("N", "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldx, &s_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid = (norm / (norm_a * norm_x + norm_b)) / ((float)n * eps);

            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
            double norm_a, norm_b, norm_x, norm, eps, resid;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");

            /* Compute AX-B */
            dgemm_("N", "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldx, &d_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid = (norm / (norm_a * norm_x + norm_b)) / ((float)n * eps);

            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
            float norm_a, norm_b, norm_x, norm, eps, resid;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            /* Compute AX-B */
            cgemm_("N", "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldx, &c_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid = (norm / (norm_a * norm_x + norm_b)) / ((float)n * eps);

            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_a, norm_b, norm_x, norm, eps, resid;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");

            /* Compute AX-B */
            zgemm_("N", "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldx, &z_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid = (norm / (norm_a * norm_x + norm_b)) / ((float)n * eps);

            *residual = (double)resid;
            break;
        }
    }
}