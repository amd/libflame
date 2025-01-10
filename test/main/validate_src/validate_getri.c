/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_getri.c
 *  @brief Defines validate function of GETRI() to use in test suite.
 *  */

#include "test_common.h"

extern double perf;
extern double time_min;

void validate_getri(char *tst_api, integer m_A, integer n_A, void *A, void *A_inv, integer lda,
                    integer *IPIV, integer datatype, double err_thresh, char imatrix)
{
    void *a_temp, *work;
    char NORM = '1';
    double residual;

    /* Early return conditions */
    if(m_A == 0 || n_A == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m_A, n_A, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m_A, n_A, err_thresh);

    /* Create Identity matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &a_temp, n_A);
    create_vector(datatype, &work, 2 * m_A);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_I, eps;

            eps = fla_lapack_slamch("Epsilon");
            /* compute I - A' * A */
            fla_lapack_slaset("full", &m_A, &m_A, &s_zero, &s_one, a_temp, &m_A);
            norm_I = sqrt(m_A);
            /* compute I - A' * A */
            sgemm_("N", "N", &m_A, &n_A, &m_A, &s_n_one, A_inv, &lda, A, &lda, &s_one, a_temp,
                   &m_A);

            compute_matrix_norm(datatype, NORM, m_A, m_A, a_temp, m_A, &norm, imatrix, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            residual = (double)(norm / (norm_I * eps * n_A));
            break;
        }

        case DOUBLE:
        {
            double norm_I, norm, eps;

            eps = fla_lapack_dlamch("Epsilon");
            /* compute I - A' * A */
            fla_lapack_dlaset("full", &m_A, &m_A, &d_zero, &d_one, a_temp, &m_A);
            norm_I = sqrt(m_A);
            /* compute I - A' * A */
            dgemm_("N", "N", &m_A, &n_A, &m_A, &d_n_one, A_inv, &lda, A, &lda, &d_one, a_temp,
                   &m_A);

            compute_matrix_norm(datatype, NORM, m_A, m_A, a_temp, m_A, &norm, imatrix, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            residual = (double)(norm / (norm_I * eps * n_A));
            break;
        }
        case COMPLEX:
        {
            float norm, norm_I, eps;

            eps = fla_lapack_slamch("Epsilon");
            /* compute I - A' * A */
            fla_lapack_claset("full", &m_A, &m_A, &c_zero, &c_one, a_temp, &m_A);
            norm_I = sqrt(m_A);
            /* compute I - A' * A */
            cgemm_("N", "N", &m_A, &n_A, &m_A, &c_n_one, A_inv, &lda, A, &lda, &c_one, a_temp,
                   &m_A);
            compute_matrix_norm(datatype, NORM, m_A, m_A, a_temp, m_A, &norm, imatrix, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            residual = (double)(norm / (norm_I * eps * n_A));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_I, eps;

            eps = fla_lapack_dlamch("Epsilon");
            /* compute I - A' * A */
            fla_lapack_zlaset("full", &m_A, &m_A, &z_zero, &z_one, a_temp, &m_A);
            norm_I = sqrt(m_A);
            /* compute I - A' * A */
            zgemm_("N", "N", &m_A, &n_A, &m_A, &z_n_one, A_inv, &lda, A, &lda, &z_one, a_temp,
                   &m_A);

            compute_matrix_norm(datatype, NORM, m_A, m_A, a_temp, m_A, &norm, imatrix, work);
            /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
            residual = (double)(norm / (norm_I * eps * n_A));
            break;
        }
        default:
            residual = err_thresh;
            break;
    }

    // Free up buffers
    free_vector(work);
    free_vector(a_temp);

    FLA_PRINT_TEST_STATUS(m_A, n_A, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}
