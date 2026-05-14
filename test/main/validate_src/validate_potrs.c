/*
    Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_potrs.c
 *  @brief Defines validate function of POTRS() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_potrs(char *tst_api, integer n, integer nrhs, void *A, integer lda, void *A_fact,
                    void *A_save, void *X, void *B, void *B_test, integer ldb, integer datatype,
                    double err_thresh, char imatrix, void *params)
{
    void *work = NULL;
    integer ldx;
    ldx = n;
    char NORM = '1';
    double residual, resid1 = 0., resid2 = 0., resid3 = 0., resid4 = 0.;

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
            float norm_a, norm_b, norm_x, norm;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);

            /* Compute AX-B */
            sgemm_("N", "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldx, &s_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1
                = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        case DOUBLE:
        {
            double norm_a, norm_b, norm_x, norm;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);

            /* Compute AX-B */
            dgemm_("N", "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldx, &d_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1
                = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        case COMPLEX:
        {
            float norm_a, norm_b, norm_x, norm;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);

            /* Compute AX-B */
            cgemm_("N", "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldx, &c_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1
                = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_a, norm_b, norm_x, norm;

            /* Test 1 */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);

            /* Compute AX-B */
            zgemm_("N", "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldx, &z_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1
                = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        default:
            resid1 = err_thresh;
            break;
    }

    /* Test 2: Ensure input matrix A (factored) was not modified by POTRS */
    resid2 = compare_matrix(datatype, "full", n, n, A_save, lda, A_fact, lda);

    /* Test 3: Check padding rows of B not modified */
    resid3 = check_padding(datatype, n, nrhs, B_test, ldb);

    /* Test 4: Check padding rows of A not modified */
    resid4 = check_padding(datatype, n, n, A_fact, lda);

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
}
