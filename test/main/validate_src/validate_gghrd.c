/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gghrd.c
 *  @brief Defines validate function of GGHRD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gghrd(char *compq, char *compz, integer n, void *A, void *A_test, integer lda,
                    void *B, void *B_test, integer ldb, void *Q, void *Q_test, integer ldq, void *Z,
                    void *Z_test, integer ldz, integer datatype, double *residual)
{
    void *work = NULL, *lambda = NULL, *alambda = NULL, *Q_tmp = NULL, *Z_tmp = NULL;

    if(n == 0)
        return;

    /* If compq=N or/and compz=N, just compare A/A_test and B/B_test matrices and return */
    if(*compq == 'N' || *compz == 'N')
    {
        switch(datatype)
        {
            case FLOAT:
            {
                /* Compare A and A_test matrices */
                float eps, norm_H, norm, resid1, norm_T, resid2;
                eps = fla_lapack_slamch("P");
                create_vector(datatype, &work, n);
                norm_H = fla_lapack_slange("1", &n, &n, A, &lda, work);
                matrix_difference(datatype, n, n, A, lda, A_test, lda);
                norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
                resid1 = norm / (float)n / norm_H / eps;

                /* Compare B and B_test matrices */
                norm_T = fla_lapack_slange("1", &n, &n, B, &ldb, work);
                matrix_difference(datatype, n, n, B, lda, B_test, lda);
                norm = fla_lapack_slange("1", &n, &n, B, &ldb, work);
                resid2 = norm / (float)n / norm_T / eps;

                *residual = (double)fla_max(resid1, resid2);
                break;
            }
            case DOUBLE:
            {
                /* Compare A and A_test matrices */
                double eps, norm_H, norm, resid1, norm_T, resid2;
                eps = fla_lapack_dlamch("P");
                create_vector(datatype, &work, n);
                norm_H = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                matrix_difference(datatype, n, n, A, lda, A_test, lda);
                norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                resid1 = norm / (double)n / norm_H / eps;

                /* Compare B and B_test matrices */
                norm_T = fla_lapack_dlange("1", &n, &n, B, &ldb, work);
                matrix_difference(datatype, n, n, B, lda, B_test, lda);
                norm = fla_lapack_dlange("1", &n, &n, B, &ldb, work);
                resid2 = norm / (double)n / norm_T / eps;

                *residual = (double)fla_max(resid1, resid2);
                break;
            }
            case COMPLEX:
            {
                /* Compare A and A_test matrices */
                float eps, norm_H, norm, resid1, norm_T, resid2;
                eps = fla_lapack_slamch("P");
                create_vector(datatype, &work, n);
                norm_H = fla_lapack_clange("1", &n, &n, A, &lda, work);
                matrix_difference(datatype, n, n, A, lda, A_test, lda);
                norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
                resid1 = norm / (float)n / norm_H / eps;

                /* Compare B and B_test matrices */
                norm_T = fla_lapack_clange("1", &n, &n, B, &ldb, work);
                matrix_difference(datatype, n, n, B, lda, B_test, lda);
                norm = fla_lapack_clange("1", &n, &n, B, &ldb, work);
                resid2 = norm / (float)n / norm_T / eps;

                *residual = (double)fla_max(resid1, resid2);
                break;
            }
            case DOUBLE_COMPLEX:
            {
                /* Compare A and A_test matrices */
                double eps, norm_H, norm, resid1, norm_T, resid2;
                eps = fla_lapack_dlamch("P");
                create_vector(datatype, &work, n);
                norm_H = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                matrix_difference(datatype, n, n, A, lda, A_test, lda);
                norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                resid1 = norm / (double)n / norm_H / eps;

                /* Compare B and B_test matrices */
                norm_T = fla_lapack_zlange("1", &n, &n, B, &ldb, work);
                matrix_difference(datatype, n, n, B, lda, B_test, lda);
                norm = fla_lapack_zlange("1", &n, &n, B, &ldb, work);
                resid2 = norm / (double)n / norm_T / eps;

                *residual = (double)fla_max(resid1, resid2);
                break;
            }
        }
        free_vector(work);
        return;
    }

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &lambda, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &alambda, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_tmp, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_tmp, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, norm_T, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_slamch("P");

            /* Test 1
                | A - Q**T * Q_test * H  * Z_test**T * Z  | / ( |A| n ulp ) */
            norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);

            if(*compq == 'V')
            {
                /* Compute Q_tmp = Q**T * Q_test
                   Compute lambda = Q**T * Q_test * H */
                sgemm_("T", "N", &n, &n, &n, &s_one, Q, &ldq, Q_test, &ldq, &s_zero, Q_tmp, &n);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_tmp, &n, A_test, &lda, &s_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * H */
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, A_test, &lda, &s_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute Z_tmp = Z_test**T * Z
                   Compute A = A - lambda * Z_test**T * Z */
                sgemm_("T", "N", &n, &n, &n, &s_one, Z_test, &ldz, Z, &ldz, &s_zero, Z_tmp, &n);
                sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z_tmp, &n, &s_n_one, A, &lda);
            }
            else
            {
                /* Compute A = A - lambda * Z_test**T */
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_n_one, A, &lda);
            }
            norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2
                | B - Q**T * Q_test *  T * Z_test**T * Z  | / ( |B| n ulp ) */
            norm_T = fla_lapack_slange("1", &n, &n, B, &ldb, work);
            if(*compq == 'V')
            {
                /* Compute lambda = Q_tmp * T */
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_tmp, &n, B_test, &ldb, &s_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * T */
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, B_test, &ldb, &s_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute B = B - lambda * Z_tmp */
                sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z_tmp, &n, &s_n_one, B, &ldb);
            }
            else
            {
                /* Compute B = B - lambda * Z_test**T */
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_n_one, B, &ldb);
            }
            norm = fla_lapack_slange("1", &n, &n, B, &ldb, work);
            resid2 = norm / (eps * norm_T * (float)n);

            /* Test 3
                compute norm(I - Z_test'*Z_test) / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 4
                compute norm(I - Q_test'*Q_test) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, norm_B, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                | A - Q**T * Q_test * H  * Z_test**T * Z  | / ( |A| n ulp ) */
            norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);

            if(*compq == 'V')
            {
                /* Compute Q_tmp = Q**T * Q_test
                   Compute lambda = Q**T * Q_test * H */
                dgemm_("T", "N", &n, &n, &n, &d_one, Q, &ldq, Q_test, &ldq, &d_zero, Q_tmp, &n);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_tmp, &n, A_test, &lda, &d_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * H */
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, A_test, &lda, &d_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute Z_tmp = Z_test**T * Z
                   Compute A = A - lambda * Z_test**T * Z */
                dgemm_("T", "N", &n, &n, &n, &d_one, Z_test, &ldz, Z, &ldz, &d_zero, Z_tmp, &n);
                dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z_tmp, &n, &d_n_one, A, &lda);
            }
            else
            {
                /* Compute A = A - lambda * Z_test**T */
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_n_one, A, &lda);
            }
            norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2
                | B - Q**H * Q_test *  T * Z_test**H * Z  | / ( |B| n ulp ) */
            norm_B = fla_lapack_dlange("1", &n, &n, B, &ldb, work);
            if(*compq == 'V')
            {
                /* Compute lambda = Q_tmp * T */
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_tmp, &n, B_test, &ldb, &d_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * T */
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, B_test, &ldb, &d_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute B = B - lambda * Z_tmp */
                dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z_tmp, &n, &d_n_one, B, &ldb);
            }
            else
            {
                /* Compute B = B - lambda * Z_test**T */
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_n_one, B, &ldb);
            }
            norm = fla_lapack_dlange("1", &n, &n, B, &ldb, work);
            resid2 = norm / (eps * norm_B * (float)n);

            /* Test 3
                compute norm(I - Z_test'*Z_test) / (N * EPS)*/
            resid3 = check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 4
                compute norm(I - Q_test'*Q_test) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, norm_B, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_slamch("P");

            /* Test 1
                | A - Q**H * Q_test * H  * Z_test**H * Z  | / ( |A| n ulp ) */
            norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
            if(*compq == 'V')
            {
                /* Compute Q_tmp = Q**H * Q_test
                   Compute lambda = Q**H * Q_test * H */
                cgemm_("C", "N", &n, &n, &n, &c_one, Q, &ldq, Q_test, &ldq, &c_zero, Q_tmp, &n);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_tmp, &n, A_test, &lda, &c_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * H */
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, A_test, &lda, &c_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute Z_tmp = Z_test**H * Z
                   Compute A = A - lambda * Z_test**H * Z */
                cgemm_("C", "N", &n, &n, &n, &c_one, Z_test, &ldz, Z, &ldz, &c_zero, Z_tmp, &n);
                cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z_tmp, &n, &c_n_one, A, &lda);
            }
            else
            {
                /* Compute A = A - lambda * Z_test**H */
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_n_one, A, &lda);
            }
            norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2
                | B - Q**H * Q_test *  T * Z_test**H * Z  | / ( |B| n ulp ) */
            norm_B = fla_lapack_clange("1", &n, &n, B, &ldb, work);
            if(*compq == 'V')
            {
                /* Compute lambda = Q_tmp * T */
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_tmp, &n, B_test, &ldb, &c_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * T */
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, B_test, &ldb, &c_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute B = B - lambda * Z_tmp */
                cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z_tmp, &n, &c_n_one, B, &ldb);
            }
            else
            {
                /* Compute B = B - lambda * Z_test**H */
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_n_one, B, &ldb);
            }
            norm = fla_lapack_clange("1", &n, &n, B, &ldb, work);
            resid2 = norm / (eps * norm_B * (float)n);

            /* Test 3
                compute norm(I - Z_test'*Z_test) / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 4
                compute norm(I - Q_test'*Q_test) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, norm_B, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                | A - Q**H * Q_test * H  * Z_test**H * Z  | / ( |A| n ulp ) */
            norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
            if(*compq == 'V')
            {
                /* Compute Q_tmp = Q**H * Q_test
                   Compute lambda = Q**H * Q_test * H */
                zgemm_("C", "N", &n, &n, &n, &z_one, Q, &ldq, Q_test, &ldq, &z_zero, Q_tmp, &n);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_tmp, &n, A_test, &lda, &z_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * H */
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, A_test, &lda, &z_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute Z_tmp = Z_test**H * Z
                   Compute A = A - lambda * Z_test**H * Z */
                zgemm_("C", "N", &n, &n, &n, &z_one, Z_test, &ldz, Z, &ldz, &z_zero, Z_tmp, &n);
                zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z_tmp, &n, &z_n_one, A, &lda);
            }
            else
            {
                /* Compute A = A - lambda * Z_test**H */
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_n_one, A, &lda);
            }
            norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
            resid1 = norm / (eps * norm_A * (float)n);

            /* Test 2
                | B - Q**H * Q_test *  T * Z_test**H * Z  | / ( |B| n ulp ) */
            norm_B = fla_lapack_zlange("1", &n, &n, B, &ldb, work);
            if(*compq == 'V')
            {
                /* Compute lambda = Q_tmp * T */
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_tmp, &n, B_test, &ldb, &z_zero, lambda, &n);
            }
            else
            {
                /* Compute lambda = Q_test * T */
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, B_test, &ldb, &z_zero, lambda,
                       &n);
            }
            if(*compz == 'V')
            {
                /* Compute B = B - lambda * Z_tmp */
                zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z_tmp, &n, &z_n_one, B, &ldb);
            }
            else
            {
                /* Compute B = B - lambda * Z_test**H */
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_n_one, B, &ldb);
            }
            norm = fla_lapack_zlange("1", &n, &n, B, &ldb, work);
            resid2 = norm / (eps * norm_B * (float)n);

            /* Test 3
                compute norm(I - Z_test'*Z_test) / (N * EPS)*/
            resid3 = check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 4
                compute norm(I - Q_test'*Q_test) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(alambda);
    free_matrix(Q_tmp);
    free_matrix(Z_tmp);
}
