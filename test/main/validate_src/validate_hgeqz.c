/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_hgeqz.c
 *  @brief Defines validate function of HGEQZ() to use in test suite.
 *  */

#include "test_common.h"

void validate_hgeqz(char *job, char *compq, char *compz, integer n, void *H, void *H_test, void *A,
                    integer ldh, void *T, void *T_test, void *B, integer ldt, void *Q, void *Q_test,
                    void *Q_A, integer ldq, void *Z, void *Z_test, void *Z_A, integer ldz,
                    integer datatype, double *residual, char imatrix, integer *info)
{
    if((*job == 'E') || (*compz == 'N') || (*compq == 'N') || (n <= 0))
        return;

    void *work = NULL, *lambda = NULL, *alambda = NULL, *h_input = NULL, *t_input = NULL,
         *Ilambda = NULL;
    *info = 0;
    char NORM = '1';

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &lambda, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &alambda, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &h_input, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &t_input, ldt);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Ilambda, n);
    set_identity_matrix(datatype, n, n, Ilambda, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_H, norm_T, eps, resid1 = 0.f, resid2 = 0.f, resid3, resid4;
            double res_max, res_max1, *res_gghrd = residual;
            eps = fla_lapack_slamch("P");

            if(*compq == 'V' && *compz == 'V')
            {
                /* Test 1
                    | A - Q S Z**T  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, H_test, &ldh, &s_zero, lambda,
                       &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, alambda,
                       &n);
                sgemm_("T", "N", &n, &n, &n, &s_one, Q_A, &ldq, alambda, &n, &s_zero, lambda, &n);
                sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z_A, &ldz, &s_n_one, A, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);

                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q P Z**T  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, T_test, &ldt, &s_zero, lambda,
                       &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, alambda,
                       &n);
                sgemm_("T", "N", &n, &n, &n, &s_one, Q_A, &ldq, alambda, &n, &s_zero, lambda, &n);
                sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z_A, &ldz, &s_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }
            if(*compq == 'I' && *compz == 'I')
            {
                /* Test 1
                    | H - Q S Z**T  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, H_test, &ldh, &s_zero, lambda,
                       &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, h_input,
                       &ldh);
                sgemm_("N", "N", &n, &n, &n, &s_one, h_input, &ldh, Ilambda, &n, &s_n_one, H, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q P Z**T  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, T_test, &ldt, &s_zero, lambda,
                       &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, t_input,
                       &ldt);
                sgemm_("N", "N", &n, &n, &n, &s_one, t_input, &ldt, Ilambda, &n, &s_n_one, T, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd, info);
                if(*info < 0)
                    break;
            }

            /* Test 4
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 5
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            if(*compq == 'I' && *compz == 'I')
                *residual = (double)fla_max(*residual, *res_gghrd);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_H, norm_T, eps, resid1 = 0., resid2 = 0., resid3, resid4;
            double res_max, res_max1, *res_gghrd = residual;
            eps = fla_lapack_dlamch("P");

            if(*compq == 'V' && *compz == 'V')
            {
                /* Test 1
                    | A - Q S Z**T  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, H_test, &ldh, &d_zero, lambda,
                       &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, alambda,
                       &n);
                dgemm_("T", "N", &n, &n, &n, &d_one, Q_A, &ldq, alambda, &n, &d_zero, lambda, &n);
                dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z_A, &ldz, &d_n_one, A, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q P Z**T  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, T_test, &ldt, &d_zero, lambda,
                       &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, alambda,
                       &n);
                dgemm_("T", "N", &n, &n, &n, &d_one, Q_A, &ldq, alambda, &n, &d_zero, lambda, &n);
                dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z_A, &ldz, &d_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }
            if(*compq == 'I' && *compz == 'I')
            {
                /* Test 1
                    | H - Q S Z**T  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, H_test, &ldh, &d_zero, lambda,
                       &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, h_input,
                       &ldh);
                dgemm_("N", "N", &n, &n, &n, &d_one, h_input, &ldh, Ilambda, &n, &d_n_one, H, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q P Z**T  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, T_test, &ldt, &d_zero, lambda,
                       &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, t_input,
                       &ldt);
                dgemm_("N", "N", &n, &n, &n, &d_one, t_input, &ldt, Ilambda, &n, &d_n_one, T, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd, info);
                if(*info < 0)
                    break;
            }

            /* Test 4
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 5
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            if(*compq == 'I' && *compz == 'I')
                *residual = (double)fla_max(*residual, *res_gghrd);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_H, norm_T, eps, resid1 = 0.f, resid2 = 0.f, resid3, resid4;
            double res_max, res_max1, *res_gghrd = residual;
            eps = fla_lapack_slamch("P");

            if(*compq == 'V' && *compz == 'V')
            {
                /* Test 1
                    | A - Q H Z**T  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, H_test, &ldh, &c_zero, lambda,
                       &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, alambda,
                       &n);
                cgemm_("C", "N", &n, &n, &n, &c_one, Q_A, &ldq, alambda, &n, &c_zero, lambda, &n);
                cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z_A, &ldz, &c_n_one, A, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q T Z**T  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, T_test, &ldt, &c_zero, lambda,
                       &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, alambda,
                       &n);
                cgemm_("C", "N", &n, &n, &n, &c_one, Q_A, &ldq, alambda, &n, &c_zero, lambda, &n);
                cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z_A, &ldz, &c_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }

            if(*compq == 'I' && *compz == 'I')
            {
                /* Test 1
                    | H - Q S Z**T  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, H_test, &ldh, &c_zero, lambda,
                       &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, h_input,
                       &ldh);
                cgemm_("N", "N", &n, &n, &n, &c_one, h_input, &ldh, Ilambda, &n, &c_n_one, H, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q P Z**T  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, T_test, &ldt, &c_zero, lambda,
                       &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, t_input,
                       &ldt);
                cgemm_("N", "N", &n, &n, &n, &c_one, t_input, &ldt, Ilambda, &n, &c_n_one, T, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd, info);
                if(*info < 0)
                    break;
            }

            /* Test 4
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 5
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            if(*compq == 'I' && *compz == 'I')
                *residual = (double)fla_max(*residual, *res_gghrd);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_H, norm_T, eps, resid1 = 0., resid2 = 0., resid3, resid4;
            double res_max, res_max1, *res_gghrd = residual;
            eps = fla_lapack_dlamch("P");

            if(*compq == 'V' && *compz == 'V')
            {
                /* Test 1
                    | A - Q S Z**T  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, H_test, &ldh, &z_zero, lambda,
                       &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, alambda,
                       &n);
                zgemm_("C", "N", &n, &n, &n, &z_one, Q_A, &ldq, alambda, &n, &z_zero, lambda, &n);
                zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z_A, &ldz, &z_n_one, A, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q P Z**T  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, T_test, &ldt, &z_zero, lambda,
                       &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, alambda,
                       &n);
                zgemm_("C", "N", &n, &n, &n, &z_one, Q_A, &ldq, alambda, &n, &z_zero, lambda, &n);
                zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z_A, &ldz, &z_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }

            if(*compq == 'I' && *compz == 'V')
            {
                /* Test 1
                    | H - Q S Z**T  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, H_test, &ldh, &z_zero, lambda,
                       &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, h_input,
                       &ldh);
                zgemm_("N", "N", &n, &n, &n, &z_one, h_input, &ldh, Ilambda, &n, &z_n_one, H, &ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q P Z**T  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, T_test, &ldt, &z_zero, lambda,
                       &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, t_input,
                       &ldt);
                zgemm_("N", "N", &n, &n, &n, &z_one, t_input, &ldt, Ilambda, &n, &z_n_one, T, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd, info);
                if(*info < 0)
                    break;
            }
            /* Test 4
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = check_orthogonality(datatype, Z_test, n, n, ldz);

            /* Test 5
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, Q_test, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            if(*compq == 'I' && *compz == 'I')
                *residual = (double)fla_max(*residual, *res_gghrd);
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(alambda);
    free_matrix(Ilambda);
    free_matrix(h_input);
    free_matrix(t_input);
}
