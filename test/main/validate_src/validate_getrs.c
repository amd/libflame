/*
    Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_getrs.c
 *  @brief Defines validate function of GETRS() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_getrs(char *tst_api, char *trans, integer n, integer nrhs, void *A, integer lda,
                    void *A_fact, void *A_fact_save, integer *IPIV, integer *IPIV_save, void *B,
                    integer ldb, void *X, integer datatype, double err_thresh, char imatrix,
                    void *scal, void *params)
{
    void *work = NULL, *Y = NULL;
    char NORM = '1';
    integer i;
    double residual, resid1 = 0., resid2 = 0., resid3 = 0.;

    /* Early return conditions */
    if(n == 0 || nrhs == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    if(same_char(imatrix, 'U'))
    {
        create_vector(datatype, &Y, i_one);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_x, norm;

            /* Test 1 */
            /* Compute AX-B */
            sgemm_(trans, "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldb, &s_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((same_char(imatrix, 'O')) && (scal != NULL))
            {
                sscal_(&n, scal, X, &i_one);
            }

            if((same_char(imatrix, 'U')) && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        case DOUBLE:
        {
            double norm_x, norm;

            /* Test 1 */
            /* Compute AX-B */
            dgemm_(trans, "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldb, &d_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((same_char(imatrix, 'O')) && (scal != NULL))
            {
                dscal_(&n, scal, X, &i_one);
            }

            if((same_char(imatrix, 'U')) && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        case COMPLEX:
        {
            float norm_x, norm;

            /* Test 1 */
            /* Compute AX-B */
            cgemm_(trans, "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldb, &c_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((same_char(imatrix, 'O')) && (scal != NULL))
            {
                cscal_(&n, scal, X, &i_one);
            }

            if((same_char(imatrix, 'U')) && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_x, norm;

            /* Test 1 */
            /* Compute AX-B */
            zgemm_(trans, "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldb, &z_n_one, B, &ldb);

            /* Scaling the X for overflow/underflow cases */
            if((same_char(imatrix, 'O')) && (scal != NULL))
            {
                zscal_(&n, scal, X, &i_one);
            }

            if((same_char(imatrix, 'U')) && (scal != NULL))
            {
                get_reciprocal_real_vector(get_realtype(datatype), scal, i_one, Y, i_one);
                scal_matrix(datatype, scal, X, n, nrhs, ldb, i_one);
            }

            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);

            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        default:
            resid1 = err_thresh;
            break;
    }

    if(same_char(imatrix, 'U'))
    {
        free_vector(Y);
    }

    /* Test 2: Ensure factored input matrix A was not modified by GETRS */
    resid2 = compare_matrix(datatype, "full", n, n, A_fact_save, lda, A_fact, lda);

    /* Test 3: Ensure IPIV was not modified by GETRS */
    if(IPIV != NULL && IPIV_save != NULL)
    {
        for(i = 0; i < n; i++)
        {
            if(IPIV[i] != IPIV_save[i])
            {
                resid3 = DBL_MAX;
                break;
            }
        }
    }

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(residual, resid3);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
}
