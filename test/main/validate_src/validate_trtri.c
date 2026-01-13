/*
    Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_trtri.c
 *  @brief Defines validate function of TRTRI() to use in test suite.
 */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_trtri(char *tst_api, char uplo, char diag, integer n, void *A, void *A_inv,
                    integer lda, integer datatype, double err_thresh, char imatrix, void *params)
{
    void *I_mat = NULL;
    char NORM = '1';
    double residual = 0;
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    /* Create Identity matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &I_mat, lda);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_I;
            /* compute I - A' * A */
            fla_lapack_slaset("full", &n, &n, &s_zero, &s_one, I_mat, &lda);
            norm_I = sqrt(n);
            sgemm_("N", "N", &n, &n, &n, &s_n_one, A_inv, &lda, A, &lda, &s_one, I_mat, &lda);
            /* Mask irrelevant triangle */
            if(uplo == 'U' || uplo == 'u')
            {
                /* Mask below the diagonal (lower triangle) */
                fla_lapack_slaset("L", &n, &n, &s_zero, &s_zero, I_mat, &lda);
            }
            else
            {
                /* Mask above the diagonal (upper triangle) */
                fla_lapack_slaset("U", &n, &n, &s_zero, &s_zero, I_mat, &lda);
            }
            compute_matrix_norm(datatype, NORM, n, n, I_mat, lda, &norm, imatrix, NULL);
            residual = fla_compute_residual(datatype, 'E', norm, norm_I, n, params);
            break;
        }
        case DOUBLE:
        {
            double norm_I, norm;
            /* compute I - A' * A */
            fla_lapack_dlaset("full", &n, &n, &d_zero, &d_one, I_mat, &lda);
            norm_I = sqrt(n);
            dgemm_("N", "N", &n, &n, &n, &d_n_one, A_inv, &lda, A, &lda, &d_one, I_mat, &lda);
            /* Mask irrelevant triangle */
            if(uplo == 'U' || uplo == 'u')
            {
                /* Mask below the diagonal (lower triangle) */
                fla_lapack_dlaset("L", &n, &n, &d_zero, &d_zero, I_mat, &lda);
            }
            else
            {
                /* Mask above the diagonal (upper triangle) */
                fla_lapack_dlaset("U", &n, &n, &d_zero, &d_zero, I_mat, &lda);
            }
            compute_matrix_norm(datatype, NORM, n, n, I_mat, lda, &norm, imatrix, NULL);
            residual = fla_compute_residual(datatype, 'E', norm, norm_I, n, params);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_I;
            /* compute I - A' * A */
            fla_lapack_claset("full", &n, &n, &c_zero, &c_one, I_mat, &lda);
            norm_I = sqrt(n);
            cgemm_("N", "N", &n, &n, &n, &c_n_one, A_inv, &lda, A, &lda, &c_one, I_mat, &lda);
            /* Mask irrelevant triangle */
            if(uplo == 'U' || uplo == 'u')
            {
                /* Mask below the diagonal (lower triangle) */
                fla_lapack_claset("L", &n, &n, &c_zero, &c_zero, I_mat, &lda);
            }
            else
            {
                /* Mask above the diagonal (upper triangle) */
                fla_lapack_claset("U", &n, &n, &c_zero, &c_zero, I_mat, &lda);
            }
            compute_matrix_norm(datatype, NORM, n, n, I_mat, lda, &norm, imatrix, NULL);
            residual = fla_compute_residual(datatype, 'E', norm, norm_I, n, params);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_I;
            /* compute I - A' * A */
            fla_lapack_zlaset("full", &n, &n, &z_zero, &z_one, I_mat, &lda);
            norm_I = sqrt(n);
            zgemm_("N", "N", &n, &n, &n, &z_n_one, A_inv, &lda, A, &lda, &z_one, I_mat, &lda);
            /* Mask irrelevant triangle */
            if(uplo == 'U' || uplo == 'u')
            {
                /* Mask below the diagonal (lower triangle) */
                fla_lapack_zlaset("L", &n, &n, &z_zero, &z_zero, I_mat, &lda);
            }
            else
            {
                /* Mask above the diagonal (upper triangle) */
                fla_lapack_zlaset("U", &n, &n, &z_zero, &z_zero, I_mat, &lda);
            }
            compute_matrix_norm(datatype, NORM, n, n, I_mat, lda, &norm, imatrix, NULL);
            residual = fla_compute_residual(datatype, 'E', norm, norm_I, n, params);
            break;
        }
    }

    free_vector(I_mat);

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}