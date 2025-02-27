/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_geqp3.c
 *  @brief Defines validate function of GEQP3() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_geqp3(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    integer *jpvt, void *T_test, integer datatype, double err_thresh, char imatrix)
{
    void *Q = NULL, *R = NULL, *work = NULL;
    char NORM = '1';
    integer min_A;
    integer lwork = -1, FLA_TRUE = 1;
    integer info = 0;
    double residual, resid1 = 0., resid2 = 0.;

    /* Early return conditions */
    if(m_A == 0 || n_A == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m_A, n_A, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m_A, n_A, err_thresh);

    min_A = fla_min(m_A, n_A);

    // Create Q and R matrices.
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, m_A, &Q, m_A);
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &R, m_A);

    reset_matrix(datatype, m_A, m_A, Q, m_A);
    reset_matrix(datatype, m_A, n_A, R, m_A);

    // Extract R matrix and elementary reflectors from the input/output matrix parameter A_test.
    copy_matrix(datatype, "full", m_A, min_A, A_test, lda, Q, m_A);
    copy_matrix(datatype, "Upper", min_A, n_A, A_test, lda, R, m_A);

    switch(datatype)
    {
        case FLOAT:
        {
            float twork;
            float norm, norm_A, eps;
            eps = fla_lapack_slamch("P");

            /* permute A using the permuted vector jpvt to get (A * P) */
            fla_lapack_slapmt(&FLA_TRUE, &m_A, &n_A, A, &lda, jpvt);
            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm_A, imatrix, work);

            /* sorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_sorgqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;
            lwork = twork;
            create_vector(datatype, &work, lwork);
            fla_lapack_sorgqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm(((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            sgemm_("N", "N", &m_A, &n_A, &m_A, &s_one, Q, &m_A, R, &m_A, &s_n_one, A, &lda);
            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm, imatrix, work);

            resid1 = norm / (eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (double)check_orthogonality(datatype, Q, m_A, m_A, m_A);
            break;
        }
        case DOUBLE:
        {
            double twork;
            double norm, norm_A, eps;

            eps = fla_lapack_dlamch("P");
            /* permute A using the permuted vector jpvt to get (A * P) */
            fla_lapack_dlapmt(&FLA_TRUE, &m_A, &n_A, A, &lda, jpvt);
            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm_A, imatrix, work);

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_dorgqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;
            lwork = twork;
            create_vector(datatype, &work, lwork);

            fla_lapack_dorgqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            dgemm_("N", "N", &m_A, &n_A, &m_A, &d_one, Q, &m_A, R, &m_A, &d_n_one, A, &lda);

            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm, imatrix, work);
            resid1 = norm / (eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m_A, m_A, m_A);
            break;
        }
        case COMPLEX:
        {
            scomplex twork;
            float norm, norm_A, eps;

            eps = fla_lapack_slamch("P");
            /* permute A using the permuted vector jpvt to get (A * P) */
            fla_lapack_clapmt(&FLA_TRUE, &m_A, &n_A, A, &lda, jpvt);
            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm_A, imatrix, work);

            /* corgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_cungqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_cungqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            cgemm_("N", "N", &m_A, &n_A, &m_A, &c_one, Q, &m_A, R, &m_A, &c_n_one, A, &lda);

            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm, imatrix, work);
            resid1 = norm / (eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m_A, m_A, m_A);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex twork;
            double norm, norm_A, eps;

            eps = fla_lapack_dlamch("P");
            /* permute A using the permuted vector jpvt to get (A * P) */
            fla_lapack_zlapmt(&FLA_TRUE, &m_A, &n_A, A, &lda, jpvt);
            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm_A, imatrix, work);

            /* zorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_zungqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_zungqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            zgemm_("N", "N", &m_A, &n_A, &m_A, &z_n_one, Q, &m_A, R, &m_A, &z_one, A, &lda);

            compute_matrix_norm(datatype, NORM, m_A, n_A, A, lda, &norm, imatrix, work);
            resid1 = norm / (eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m_A, m_A, m_A);
            break;
        }
    }
    // Free up buffers
    free_matrix(R);
    free_matrix(Q);
    free_vector(work);

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(m_A, n_A, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
