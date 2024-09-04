/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_hgeqz.c
 *  @brief Defines validate function of HGEQZ() to use in test suite.
 *  */

#include "test_common.h"

void validate_hgeqz_comp_n(integer datatype, integer n, void *H_test, void *H_ntest, integer ldh,
                           void *T_test, void *T_ntest, integer ldt, void *alpha, void *alphar,
                           void *alphai, void *beta, void *alphan, void *alphanr, void *alphani,
                           void *betan, double *residual)
{
    if(n <= 0)
        return;
    /* Validate the eigen values when compq=N or compz=N */
    validate_hgeqz_eigen_values(datatype, n, alpha, alphar, alphai, beta, alphan, alphanr, alphani,
                                betan, residual);

    void *work = NULL;
    switch(datatype)
    {
        case FLOAT:
        {
            /* Compare the H_test and H_ntest matrices */
            float eps, norm_original, norm, resid;
            eps = fla_lapack_slamch("P");
            norm_original = fla_lapack_slange("1", &n, &n, H_test, &ldh, work);
            matrix_difference(datatype, n, n, H_test, ldh, H_ntest, ldh);
            norm = fla_lapack_slange("1", &n, &n, H_test, &ldh, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, (double)resid);

            /* Compare the T_test and T_ntest matrices */
            norm_original = fla_lapack_slange("1", &n, &n, T_test, &ldt, work);
            matrix_difference(datatype, n, n, T_test, ldt, T_ntest, ldt);
            norm = fla_lapack_slange("1", &n, &n, T_test, &ldt, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, (double)resid);
            break;
        }
        case DOUBLE:
        {
            /* Compare the H_test and H_ntest matrices */
            double eps, norm_original, norm, resid;
            eps = fla_lapack_dlamch("P");
            norm_original = fla_lapack_dlange("1", &n, &n, H_test, &ldh, work);
            matrix_difference(datatype, n, n, H_test, ldh, H_ntest, ldh);
            norm = fla_lapack_dlange("1", &n, &n, H_test, &ldh, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, resid);

            /* Compare the T_test and T_ntest matrices */
            norm_original = fla_lapack_dlange("1", &n, &n, T_test, &ldt, work);
            matrix_difference(datatype, n, n, T_test, ldt, T_ntest, ldt);
            norm = fla_lapack_dlange("1", &n, &n, T_test, &ldt, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, resid);
            break;
        }
        case COMPLEX:
        {
            /* Compare the H_test and H_ntest matrices */
            float eps, norm_original, norm, resid;
            eps = fla_lapack_slamch("P");
            norm_original = fla_lapack_clange("1", &n, &n, H_test, &ldh, work);
            matrix_difference(datatype, n, n, H_test, ldh, H_ntest, ldh);
            norm = fla_lapack_clange("1", &n, &n, H_test, &ldh, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, (double)resid);

            /* Compare the T_test and T_ntest matrices */
            norm_original = fla_lapack_clange("1", &n, &n, T_test, &ldt, work);
            matrix_difference(datatype, n, n, T_test, ldt, T_ntest, ldt);
            norm = fla_lapack_clange("1", &n, &n, T_test, &ldt, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, (double)resid);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Compare the H_test and H_ntest matrices */
            double eps, norm_original, norm, resid;
            eps = fla_lapack_dlamch("P");
            norm_original = fla_lapack_zlange("1", &n, &n, H_test, &ldh, work);
            matrix_difference(datatype, n, n, H_test, ldh, H_ntest, ldh);
            norm = fla_lapack_zlange("1", &n, &n, H_test, &ldh, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, resid);

            /* Compare the T_test and T_ntest matrices */
            norm_original = fla_lapack_zlange("1", &n, &n, T_test, &ldt, work);
            matrix_difference(datatype, n, n, T_test, ldt, T_ntest, ldt);
            norm = fla_lapack_zlange("1", &n, &n, T_test, &ldt, work);

            resid = norm / (norm_original * n * eps);
            *residual = fla_max(*residual, resid);
            break;
        }
    }
}

void validate_hgeqz_eigen_values(integer datatype, integer n, void *alpha, void *alphar,
                                 void *alphai, void *beta, void *alphae, void *alphaer,
                                 void *alphaei, void *betae, double *residual)
{
    if(n <= 0)
        return;
    void *work = NULL;
    *residual = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm_original, norm, resid;

            /* Compare alphar and alphaer vectors */
            eps = fla_lapack_slamch("P");
            norm_original = fla_lapack_slange("1", &n, &i_one, alphar, &i_one, work);
            saxpy_(&n, &s_n_one, alphaer, &i_one, alphar, &i_one);
            norm = fla_lapack_slange("1", &n, &i_one, alphar, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, (double)resid);

            /* Compare alphai and alphaei vectors */
            norm_original = fla_lapack_slange("1", &n, &i_one, alphai, &i_one, work);
            saxpy_(&n, &s_n_one, alphaei, &i_one, alphai, &i_one);
            norm = fla_lapack_slange("1", &n, &i_one, alphai, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, (double)resid);

            /* Compare beta and betae vectors */
            norm_original = fla_lapack_slange("1", &n, &i_one, beta, &i_one, work);
            saxpy_(&n, &s_n_one, betae, &i_one, beta, &i_one);
            norm = fla_lapack_slange("1", &n, &i_one, beta, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, (double)resid);

            break;
        }

        case DOUBLE:
        {
            double eps, norm_original, norm, resid;

            /* Compare alphar and alphaer vectors */
            eps = fla_lapack_dlamch("P");
            norm_original = fla_lapack_dlange("1", &n, &i_one, alphar, &i_one, work);
            daxpy_(&n, &d_n_one, alphaer, &i_one, alphar, &i_one);
            norm = fla_lapack_dlange("1", &n, &i_one, alphar, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, resid);

            /* Compare alphai and alphaei vectors */
            norm_original = fla_lapack_dlange("1", &n, &i_one, alphai, &i_one, work);
            daxpy_(&n, &d_n_one, alphaei, &i_one, alphai, &i_one);
            norm = fla_lapack_dlange("1", &n, &i_one, alphai, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, resid);

            /* Compare beta and betae vectors */
            norm_original = fla_lapack_dlange("1", &n, &i_one, beta, &i_one, work);
            daxpy_(&n, &d_n_one, betae, &i_one, beta, &i_one);
            norm = fla_lapack_dlange("1", &n, &i_one, beta, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, resid);

            break;
        }
        case COMPLEX:
        {
            float eps, norm_original, norm, resid;

            /* Compare alpha and alphae vectors */
            eps = fla_lapack_slamch("P");
            norm_original = fla_lapack_clange("1", &n, &i_one, alpha, &i_one, work);
            caxpy_(&n, &c_n_one, alphae, &i_one, alpha, &i_one);
            norm = fla_lapack_clange("1", &n, &i_one, alpha, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, (double)resid);

            /* Compare beta and betae vectors */
            norm_original = fla_lapack_clange("1", &n, &i_one, beta, &i_one, work);
            caxpy_(&n, &c_n_one, betae, &i_one, beta, &i_one);
            norm = fla_lapack_clange("1", &n, &i_one, beta, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, (double)resid);

            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm_original, norm, resid;

            /* Compare alpha and alphae vectors */
            eps = fla_lapack_dlamch("P");
            norm_original = fla_lapack_zlange("1", &n, &i_one, alpha, &i_one, work);
            zaxpy_(&n, &z_n_one, alphae, &i_one, alpha, &i_one);
            norm = fla_lapack_zlange("1", &n, &i_one, alpha, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, resid);

            /* Compare beta and betae vectors */
            norm_original = fla_lapack_zlange("1", &n, &i_one, beta, &i_one, work);
            zaxpy_(&n, &z_n_one, betae, &i_one, beta, &i_one);
            norm = fla_lapack_zlange("1", &n, &i_one, beta, &i_one, work);

            resid = norm / (eps * norm_original * n);
            *residual = fla_max(*residual, resid);

            break;
        }
    }
}

void validate_hgeqz(char *job, char *compq, char *compz, integer n, void *H, void *H_test, void *A,
                    integer ldh, void *T, void *T_test, void *B, integer ldt, void *Q, void *Q_test,
                    void *Q_A, integer ldq, void *Z, void *Z_test, void *Z_A, integer ldz,
                    integer datatype, double *residual, char imatrix, integer *info)
{
    if(n <= 0)
        return;

    void *work = NULL, *lambda = NULL, *h_input = NULL, *t_input = NULL, *Q_tmp = NULL,
         *Z_tmp = NULL;
    *info = 0;
    char NORM = '1';

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &lambda, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &h_input, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &t_input, ldt);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_tmp, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_tmp, n);

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
                    | A - Q_A**T * Q_test * S * Z_test**T * Z_A**T  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                sgemm_("T", "N", &n, &n, &n, &s_one, Q_A, &ldq, Q_test, &ldq, &s_zero, Q_tmp, &n);
                sgemm_("T", "N", &n, &n, &n, &s_one, Z_test, &ldz, Z_A, &ldz, &s_zero, Z_tmp, &n);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_tmp, &n, H_test, &ldh, &s_zero, lambda, &n);
                sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z_tmp, &n, &s_n_one, A, &ldh);

                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q_A**T * Q_test * P * Z_test**T * Z_A**T  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_tmp, &n, T_test, &ldt, &s_zero, lambda, &n);
                sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z_tmp, &n, &s_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }
            if(*compq == 'I' && *compz == 'I')
            {
                /* Test 1
                    | H - Q_test S Z_test**T  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, H_test, &ldh, &s_zero, lambda,
                       &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, h_input,
                       &ldh);
                matrix_difference(datatype, n, n, H, ldh, h_input, ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q_test P Z_test**T  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, T_test, &ldt, &s_zero, lambda,
                       &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, t_input,
                       &ldt);
                matrix_difference(datatype, n, n, T, ldt, t_input, ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd);
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
            double norm, norm_H, norm_T, eps, resid1 = 0.f, resid2 = 0.f, resid3, resid4;
            double res_max, res_max1, *res_gghrd = residual;
            eps = fla_lapack_dlamch("P");

            if(*compq == 'V' && *compz == 'V')
            {
                /* Test 1
                    | A - Q_A**T * Q_test * S * Z_test**T * Z_A**T  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                dgemm_("T", "N", &n, &n, &n, &d_one, Q_A, &ldq, Q_test, &ldq, &d_zero, Q_tmp, &n);
                dgemm_("T", "N", &n, &n, &n, &d_one, Z_test, &ldz, Z_A, &ldz, &d_zero, Z_tmp, &n);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_tmp, &n, H_test, &ldh, &d_zero, lambda, &n);
                dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z_tmp, &n, &d_n_one, A, &ldh);

                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q_A**T * Q_test * P * Z_test**T * Z_A**T  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_tmp, &n, T_test, &ldt, &d_zero, lambda, &n);
                dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z_tmp, &n, &d_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }
            if(*compq == 'I' && *compz == 'I')
            {
                /* Test 1
                    | H - Q_test S Z_test**T  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, H_test, &ldh, &d_zero, lambda,
                       &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, h_input,
                       &ldh);
                matrix_difference(datatype, n, n, H, ldh, h_input, ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q_test P Z_test**T  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, T_test, &ldt, &d_zero, lambda,
                       &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, t_input,
                       &ldt);
                matrix_difference(datatype, n, n, T, ldt, t_input, ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd);
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
                    | A - Q_A**H * Q_test * S * Z_test**H * Z_A**H  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                cgemm_("C", "N", &n, &n, &n, &c_one, Q_A, &ldq, Q_test, &ldq, &c_zero, Q_tmp, &n);
                cgemm_("C", "N", &n, &n, &n, &c_one, Z_test, &ldz, Z_A, &ldz, &c_zero, Z_tmp, &n);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_tmp, &n, H_test, &ldh, &c_zero, lambda, &n);
                cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z_tmp, &n, &c_n_one, A, &ldh);

                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q_A**H * Q_test * P * Z_test**H * Z_A**H  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_tmp, &n, T_test, &ldt, &c_zero, lambda, &n);
                cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z_tmp, &n, &c_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }
            if(*compq == 'I' && *compz == 'I')
            {
                /* Test 1
                    | H - Q_test S Z_test**H  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, H_test, &ldh, &c_zero, lambda,
                       &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, h_input,
                       &ldh);
                matrix_difference(datatype, n, n, H, ldh, h_input, ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q_test P Z_test**H  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, T_test, &ldt, &c_zero, lambda,
                       &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, t_input,
                       &ldt);
                matrix_difference(datatype, n, n, T, ldt, t_input, ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd);
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
            double norm, norm_H, norm_T, eps, resid1 = 0.f, resid2 = 0.f, resid3, resid4;
            double res_max, res_max1, *res_gghrd = residual;
            eps = fla_lapack_dlamch("P");

            if(*compq == 'V' && *compz == 'V')
            {
                /* Test 1
                    | A - Q_A**H * Q_test * S * Z_test**H * Z_A**H  | / ( |A| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm_H, imatrix, work);
                zgemm_("C", "N", &n, &n, &n, &z_one, Q_A, &ldq, Q_test, &ldq, &z_zero, Q_tmp, &n);
                zgemm_("C", "N", &n, &n, &n, &z_one, Z_test, &ldz, Z_A, &ldz, &z_zero, Z_tmp, &n);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_tmp, &n, H_test, &ldh, &z_zero, lambda, &n);
                zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z_tmp, &n, &z_n_one, A, &ldh);

                compute_matrix_norm(datatype, NORM, n, n, A, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | B - Q_A**H * Q_test * P * Z_test**H * Z_A**H  | / ( |B| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm_T, imatrix, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_tmp, &n, T_test, &ldt, &z_zero, lambda, &n);
                zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z_tmp, &n, &z_n_one, B, &ldt);
                compute_matrix_norm(datatype, NORM, n, n, B, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);
            }
            if(*compq == 'I' && *compz == 'I')
            {
                /* Test 1
                    | H - Q_test S Z_test**H  | / ( |H| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, H_test, &ldh, &z_zero, lambda,
                       &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, h_input,
                       &ldh);
                matrix_difference(datatype, n, n, H, ldh, h_input, ldh);
                compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
                resid1 = norm / (eps * norm_H * (float)n);

                /* Test 2
                    | T - Q_test P Z_test**H  | / ( |T| n ulp ) */
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm_T, imatrix, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, T_test, &ldt, &z_zero, lambda,
                       &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, t_input,
                       &ldt);
                matrix_difference(datatype, n, n, T, ldt, t_input, ldt);
                compute_matrix_norm(datatype, NORM, n, n, T, ldt, &norm, imatrix, work);
                resid2 = norm / (eps * norm_T * (float)n);

                /* Test 3 */
                validate_gghrd(compq, compz, n, A, h_input, ldh, B, t_input, ldt, Q_A, Q, ldq, Z_A,
                               Z, ldz, datatype, res_gghrd);
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
    free_matrix(h_input);
    free_matrix(t_input);
    free_matrix(Q_tmp);
    free_matrix(Z_tmp);
}
