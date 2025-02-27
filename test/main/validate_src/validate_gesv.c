/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gesv.c
 *  @brief Defines validate function of GESV() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_gesv(char *tst_api, integer n, integer nrhs, void *A, integer lda, void *B,
                   integer ldb, void *X, integer datatype, double err_thresh, char imatrix,
                   void *scal)
{
    void *work = NULL;
    char NORM = '1';
    integer ldx;
    ldx = ldb;
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
            float norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            sgemm_("N", "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldx, &s_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
            {
                sscal_(&n, scal, X, &i_one);
            }
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            residual = norm / (norm_x * n * eps);
            break;
        }
        case DOUBLE:
        {
            double norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            dgemm_("N", "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldx, &d_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
            {
                dscal_(&n, scal, X, &i_one);
            }
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = norm / (norm_x * n * eps);
            break;
        }
        case COMPLEX:
        {
            float norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            cgemm_("N", "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldx, &c_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
            {
                cscal_(&n, scal, X, &i_one);
            }
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = norm / (norm_x * n * eps);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            zgemm_("N", "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldx, &z_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
            {
                zscal_(&n, scal, X, &i_one);
            }
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = norm / (norm_x * n * eps);
            break;
        }
        default:
            residual = err_thresh;
            break;
    }

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}
