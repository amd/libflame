/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_orgqr.c
 *  @brief Defines validate function of ORGQR() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_orgqr(char *tst_api, integer m, integer n, void *A, integer lda, void *Q, void *R,
                    integer datatype, double err_thresh, char imatrix)
{
    integer k;
    double residual, resid1 = 0., resid2 = 0.;
    void *work = NULL;

    /* Early return conditions */
    if(m == 0 || n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    k = m;
    char NORM = '1';

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            sgemm_("T", "N", &n, &n, &k, &s_n_one, Q, &lda, A, &lda, &s_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, n, lda);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            dgemm_("T", "N", &n, &n, &k, &d_n_one, Q, &lda, A, &lda, &d_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, n, lda);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            cgemm_("C", "N", &n, &n, &k, &c_n_one, Q, &lda, A, &lda, &c_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, n, lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            zgemm_("C", "N", &n, &n, &k, &z_n_one, Q, &lda, A, &lda, &z_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, n, lda);
            break;
        }
    }

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
