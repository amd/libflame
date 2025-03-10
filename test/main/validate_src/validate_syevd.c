/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_syevd.c
 *  @brief Defines validate function of SYEVD() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_syevd(char *tst_api, char *jobz, integer n, void *A, void *A_test, integer lda,
                    void *w, integer datatype, double err_thresh)
{
    double residual, resid1 = 0., resid2 = 0.;

    /* Early return conditions */
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    if(!same_char(*jobz, 'N'))
    {
        void *lambda = NULL, *zlambda = NULL, *Z = NULL;
        void *work = NULL;

        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &lambda, n);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &zlambda, n);
        create_matrix(datatype, LAPACK_COL_MAJOR, lda, n, &Z, n);

        reset_matrix(datatype, n, n, zlambda, n);
        reset_matrix(datatype, n, n, Z, lda);

        copy_matrix(datatype, "full", n, n, A_test, lda, Z, lda);

        diagonalize_realtype_vector(datatype, w, lambda, n, n, n);

        switch(datatype)
        {
            case FLOAT:
            {
                float norm, norm_A, eps;
                eps = fla_lapack_slamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Z, &lda, lambda, &n, &s_zero, zlambda, &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, Z, &lda, &s_n_one, A, &lda);
                norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = (float)check_orthogonality(datatype, Z, n, n, lda);
                break;
            }

            case DOUBLE:
            {
                double norm, norm_A, eps;
                eps = fla_lapack_dlamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Z, &lda, lambda, &n, &d_zero, zlambda, &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, Z, &lda, &d_n_one, A, &lda);
                norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = check_orthogonality(datatype, Z, n, n, lda);
                break;
            }

            case COMPLEX:
            {
                float norm, norm_A, eps;
                eps = fla_lapack_slamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Z, &lda, lambda, &n, &c_zero, zlambda, &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, Z, &lda, &c_n_one, A, &lda);
                norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = (float)check_orthogonality(datatype, Z, n, n, lda);
                break;
            }

            case DOUBLE_COMPLEX:
            {
                double norm, norm_A, eps;
                eps = fla_lapack_dlamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Z, &lda, lambda, &n, &z_zero, zlambda, &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, Z, &lda, &z_n_one, A, &lda);
                norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = check_orthogonality(datatype, Z, n, n, lda);
                break;
            }
        }
        free_matrix(lambda);
        free_matrix(zlambda);
        free_matrix(Z);
    }

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
