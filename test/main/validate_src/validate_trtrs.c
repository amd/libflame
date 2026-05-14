/*
    Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_trtrs.c
 *  @brief Defines validate function of TRTRS() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_trtrs(char *tst_api, integer datatype, char *uplo, char *trans, char *diag, integer n,
                    integer nrhs, void *A, void *A_save, integer lda, void *X, void *B, integer ldb,
                    double err_thresh, char imatrix, void *params)
{
    void *work = NULL;
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

            /* Test 1: Compute residual ||AX - B|| / (||A|| * ||X|| + ||B||) */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            
            /* Compute AX using TRMM: X := A*X */
            strmm_("L", uplo, trans, diag, &n, &nrhs, &s_one, A, &lda, X, &ldb);

            /* Compute AX - B by computing X := X - B */
            matrix_difference(datatype, n, nrhs, X, ldb, B, ldb);

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        case DOUBLE:
        {
            double norm_a, norm_b, norm_x, norm;

            /* Test 1: Compute residual ||AX - B|| / (||A|| * ||X|| + ||B||) */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);

            /* Compute AX using TRMM: X := A*X */
            dtrmm_("L", uplo, trans, diag, &n, &nrhs, &d_one, A, &lda, X, &ldb);

            /* Compute AX - B by computing X := X - B */
            matrix_difference(datatype, n, nrhs, X, ldb, B, ldb);

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        case COMPLEX:
        {
            float norm_a, norm_b, norm_x, norm;

            /* Test 1: Compute residual ||AX - B|| / (||A|| * ||X|| + ||B||) */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);

            /* Compute AX using TRMM: X := A*X */
            ctrmm_("L", uplo, trans, diag, &n, &nrhs, &c_one, A, &lda, X, &ldb);

            /* Compute AX - B by computing X := X - B */
            matrix_difference(datatype, n, nrhs, X, ldb, B, ldb);

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_a, norm_b, norm_x, norm;

            /* Test 1: Compute residual ||AX - B|| / (||A|| * ||X|| + ||B||) */
            compute_matrix_norm(datatype, NORM, n, n, A, lda, &norm_a, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm_b, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);

            /* Compute AX using TRMM: X := A*X */
            ztrmm_("L", uplo, trans, diag, &n, &nrhs, &z_one, A, &lda, X, &ldb);

            /* Compute AX - B by computing X := X - B */
            matrix_difference(datatype, n, nrhs, X, ldb, B, ldb);

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, (norm_a * norm_x + norm_b), n, params);
            break;
        }
        default:
            resid1 = err_thresh;
            break;
    }

    /* Test 2: Check padding rows of A not modified */
    resid2 = check_padding(datatype, n, n, A, lda);

    /* Test 3: Check padding rows of X not modified */
    resid3 = check_padding(datatype, n, nrhs, X, ldb);

    /* Test 4: Ensure input matrix A was not modified by TRTRS */
    resid4 = compare_matrix(datatype, "full", n, n, A_save, lda, A, lda);

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
}
