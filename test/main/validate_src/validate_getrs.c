/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_getrs.c
 *  @brief Defines validate function of GETRS() to use in test suite.
 *  */

#include "test_common.h"

extern double perf;
extern double time_min;

void validate_getrs(char *tst_api, char *trans, integer n, integer nrhs, void *A, integer lda,
                    void *B, integer ldb, void *X, integer datatype, double err_thresh,
                    char imatrix, void *scal)
{
    void *work = NULL, *Y = NULL;
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

    if(imatrix == 'U')
    {
        create_vector(datatype, &Y, i_one);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            sgemm_(trans, "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldb, &s_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((imatrix == 'O') && (scal != NULL))
            {
                sscal_(&n, scal, X, &i_one);
            }

            if((imatrix == 'U') && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = ((norm / norm_x) / (n * eps));
            break;
        }
        case DOUBLE:
        {
            double norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            dgemm_(trans, "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldb, &d_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((imatrix == 'O') && (scal != NULL))
            {
                dscal_(&n, scal, X, &i_one);
            }

            if((imatrix == 'U') && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = ((norm / norm_x) / (n * eps));
            break;
        }
        case COMPLEX:
        {
            float norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            cgemm_(trans, "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldb, &c_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((imatrix == 'O') && (scal != NULL))
            {
                cscal_(&n, scal, X, &i_one);
            }

            if((imatrix == 'U') && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = ((norm / norm_x) / (n * eps));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_x, norm, eps;

            /* Test 1 */
            /* Compute AX-B */
            zgemm_(trans, "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldb, &z_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((imatrix == 'O') && (scal != NULL))
            {
                zscal_(&n, scal, X, &i_one);
            }

            if((imatrix == 'U') && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            residual = ((norm / norm_x) / (n * eps));
            break;
        }
        default:
            residual = err_thresh;
            break;
    }

    if(imatrix == 'U')
    {
        free_vector(Y);
    }

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}
