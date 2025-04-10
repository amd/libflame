/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gerqf.c
 *  @brief Defines validate function of GERQF() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_gerqf(char *tst_api, integer m_A, integer n_A, void *A, void *A_test, integer lda,
                    void *T_test, integer datatype, double err_thresh)
{
    void *R, *Q, *work = NULL;
    integer min_A, diff_A;
    integer lwork = -1;
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

    if(m_A <= n_A)
        diff_A = n_A - m_A;
    else
        diff_A = m_A - n_A;

    // Create R and Q matrices.
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &R, m_A);

    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &Q, n_A);
    reset_matrix(datatype, m_A, n_A, R, m_A);
    reset_matrix(datatype, n_A, n_A, Q, n_A);

    // Extract R matrix and elementary reflectors from the input/output matrix parameter A_test.
    if(m_A <= n_A)
    {
        copy_matrix(datatype, "Upper", m_A, m_A, get_m_ptr(datatype, A_test, 0, diff_A, lda), lda,
                    get_m_ptr(datatype, R, 0, diff_A, m_A), m_A);
        copy_matrix(datatype, "full", m_A, n_A, A_test, lda, get_m_ptr(datatype, Q, diff_A, 0, n_A),
                    n_A);
    }
    else
    {
        copy_matrix(datatype, "full", diff_A, n_A, A_test, lda, R, m_A);
        copy_matrix(datatype, "Upper", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, lda), lda,
                    get_m_ptr(datatype, R, diff_A, 0, m_A), m_A);
        copy_matrix(datatype, "full", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, lda), lda, Q,
                    n_A);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float twork;
            float norm, norm_A, eps;

            /* sorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_sorgrq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;

            lwork = twork;
            create_vector(datatype, &work, lwork);

            fla_lapack_sorgrq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            sgemm_("N", "T", &m_A, &n_A, &n_A, &s_n_one, A, &lda, Q, &n_A, &s_one, R, &m_A);

            norm_A = fla_lapack_slange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_slange("1", &m_A, &n_A, R, &m_A, work);

            // Get machine precision
            eps = fla_lapack_slamch("P");

            resid1 = norm / (eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n_A, n_A, n_A);
            break;
        }
        case DOUBLE:
        {
            double twork;
            double norm, norm_A, eps;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_dorgrq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;

            lwork = twork;
            create_vector(datatype, &work, lwork);

            fla_lapack_dorgrq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            dgemm_("N", "T", &m_A, &n_A, &n_A, &d_n_one, A, &lda, Q, &n_A, &d_one, R, &m_A);

            norm_A = fla_lapack_dlange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_dlange("1", &m_A, &n_A, R, &m_A, work);

            eps = fla_lapack_dlamch("P");

            resid1 = norm / (eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n_A, n_A, n_A);
            break;
        }
        case COMPLEX:
        {
            scomplex twork;
            float norm, norm_A, eps;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_cungrq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_cungrq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            cgemm_("N", "C", &m_A, &n_A, &n_A, &c_n_one, A, &lda, Q, &n_A, &c_one, R, &m_A);
            norm_A = fla_lapack_clange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_clange("1", &m_A, &n_A, R, &m_A, work);

            eps = fla_lapack_slamch("P");

            resid1 = norm / (eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n_A, n_A, n_A);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex twork;
            double norm, norm_A, eps;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_zungrq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &info);
            if(info < 0)
                break;

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_zungrq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &info);
            if(info < 0)
                break;

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            zgemm_("N", "C", &m_A, &n_A, &n_A, &z_n_one, A, &lda, Q, &n_A, &z_one, R, &m_A);
            norm_A = fla_lapack_zlange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_zlange("1", &m_A, &n_A, R, &m_A, work);

            eps = fla_lapack_dlamch("P");

            resid1 = norm / (eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n_A, n_A, n_A);
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
