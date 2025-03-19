/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_sytrd.c
 *  @brief Defines validate function of SYTRD() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

void validate_sytrd(char *tst_api, integer datatype, char uplo, integer n, void *A_test, void *A_in,
                    integer lda, void *D, void *E, void *tau, double err_thresh)
{
    double residual = 0., resid1 = 0., resid2 = 0.;
    void *A_save = NULL, *work = NULL, *Q_temp = NULL;
    void *Q = NULL, *T = NULL;
    integer lwork = -1;
    integer info = 0;

    /* Early return conditions */
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, n);
    reset_matrix(datatype, n, n, A_save, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_temp, n);
    reset_matrix(datatype, n, n, Q_temp, n);
    create_vector(datatype, &work, 1);

    /* Create Q and T matrices. */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &T, n);

    reset_matrix(datatype, n, n, Q, n);
    reset_matrix(datatype, n, n, T, n);

    /* Copying the matrices needed */
    copy_matrix(datatype, &uplo, n, n, A_test, lda, Q, n);
    copy_sym_tridiag_matrix(datatype, D, E, n, n, T, n);
    copy_matrix(datatype, "Full", n, n, T, n, A_save, n);

    lwork = -1;
    invoke_orgtr(datatype, &uplo, &n, NULL, &n, tau, work, &lwork, &info);

    /* Get work size for ORGTR */
    lwork = get_work_value(datatype, work);
    free_vector(work);
    create_vector(datatype, &work, lwork);

    /* Generate orthogonal matrix in Q */
    invoke_orgtr(datatype, &uplo, &n, Q, &n, tau, work, &lwork, &info);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps;
            eps = fla_lapack_slamch("P");

            /* Test 1: Generating symmetric matrix A using  T - Q**T * A * Q | / ( |T| n ulp ) */
            norm_A = fla_lapack_slange("1", &n, &n, A_save, &n, NULL);
            sgemm_("T", "N", &n, &n, &n, &s_one, Q, &n, A_in, &lda, &s_zero, Q_temp, &n);
            sgemm_("N", "N", &n, &n, &n, &s_one, Q_temp, &n, Q, &n, &s_n_one, A_save, &n);
            norm = fla_lapack_slange("1", &n, &n, A_save, &n, NULL);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2: Compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n, n, n);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps;
            eps = fla_lapack_dlamch("P");

            /* Test 1: Generating symmetric matrix A using  T - Q**T * A * Q | / ( |T| n ulp ) */
            norm_A = fla_lapack_dlange("1", &n, &n, A_save, &n, NULL);
            dgemm_("T", "N", &n, &n, &n, &d_one, Q, &n, A_in, &lda, &d_zero, Q_temp, &n);
            dgemm_("N", "N", &n, &n, &n, &d_one, Q_temp, &n, Q, &n, &d_n_one, A_save, &n);
            norm = fla_lapack_dlange("1", &n, &n, A_save, &n, NULL);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2: Compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n, n, n);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps;
            eps = fla_lapack_slamch("P");

            /* Test 1: Generating symmetric matrix A using  T - Q**T * A * Q | / ( |T| n ulp ) */
            norm_A = fla_lapack_clange("1", &n, &n, A_save, &n, work);
            cgemm_("C", "N", &n, &n, &n, &c_one, Q, &n, A_in, &lda, &c_zero, Q_temp, &n);
            cgemm_("N", "N", &n, &n, &n, &c_one, Q_temp, &n, Q, &n, &c_n_one, A_save, &n);
            norm = fla_lapack_clange("1", &n, &n, A_save, &n, work);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2: Compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n, n, n);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps;
            eps = fla_lapack_dlamch("P");

            /* Test 1: Generating symmetric matrix A using  T - Q**T * A * Q | / ( |T| n ulp ) */
            norm_A = fla_lapack_zlange("1", &n, &n, A_save, &n, work);
            zgemm_("C", "N", &n, &n, &n, &z_one, Q, &n, A_in, &lda, &z_zero, Q_temp, &n);
            zgemm_("N", "N", &n, &n, &n, &z_one, Q_temp, &n, Q, &n, &z_n_one, A_save, &n);
            norm = fla_lapack_zlange("1", &n, &n, A_save, &n, work);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2: Compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n, n, n);
            break;
        }
    }

    free_matrix(T);
    free_matrix(A_save);
    free_matrix(Q);
    free_matrix(Q_temp);
    free_vector(work);

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
