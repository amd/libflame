/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_potrs.c
 *  @brief Defines validate function of POTRS() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_potrs(char *tst_api, integer n, integer nrhs, void *A, integer lda, void *X, void *B,
                    integer ldb, integer datatype, double err_thresh, char imatrix)
{
    void *work = NULL;
    integer ldx;
    ldx = n;
    char NORM = '1';
    double residual;

    /* Early return conditions */
    if(n == 0 || nrhs == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_a, norm_b, norm_x, norm, eps;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            /* Compute AX-B */
            sgemm_("N", "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldx, &s_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = (double)((norm / (norm_a * norm_x + norm_b)) / ((float)n * eps));
            break;
        }
        case DOUBLE:
        {
            double norm_a, norm_b, norm_x, norm, eps;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");

            /* Compute AX-B */
            dgemm_("N", "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldx, &d_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = (norm / (norm_a * norm_x + norm_b)) / ((float)n * eps);
            break;
        }
        case COMPLEX:
        {
            float norm_a, norm_b, norm_x, norm, eps;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            /* Compute AX-B */
            cgemm_("N", "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldx, &c_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = (double)((norm / (norm_a * norm_x + norm_b)) / ((float)n * eps));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_a, norm_b, norm_x, norm, eps;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");

            /* Compute AX-B */
            zgemm_("N", "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldx, &z_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = (norm / (norm_a * norm_x + norm_b)) / ((float)n * eps);
            break;
        }
        default:
            residual = err_thresh;
            break;
    }

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}
