/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_hseqr.c
 *  @brief Defines validate function of HSEQR() to use in test suite.
 *  */

#include "test_common.h"

extern double perf;
extern double time_min;

void validate_hseqr(char *tst_api, char *job, char *compz, integer n, void *H, void *H_test,
                    integer ldh, void *Z, void *Z_test, integer ldz, void *wr, void *wr_in,
                    void *wi, void *wi_in, void *w, integer datatype, double err_thresh,
                    integer *ilo, integer *ihi, char imatrix, void *scal_H)
{
    char NORM = '1';
    void *Y = NULL;
    double residual, resid1 = 0., resid2 = 0.;
    double resid3 = 0., resid4 = 0.;

    /* Early return conditions */
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    if((imatrix == 'O' || imatrix == 'U') && (scal_H != NULL))
    {
        create_realtype_vector(datatype, &Y, 1);
        get_reciprocal_real_vector(get_realtype(datatype), scal_H, 1, Y, 1);

        if(datatype == FLOAT || datatype == DOUBLE)
        {
            scal_matrix(datatype, Y, wr, n, i_one, i_one, i_one);
            scal_matrix(datatype, Y, wi, n, i_one, i_one, i_one);
        }
        else
        {
            scal_matrix(datatype, Y, w, n, i_one, i_one, i_one);
        }
        free_vector(Y);
    }

    if(datatype == FLOAT || datatype == DOUBLE)
    {
        /* Find negative value of every 2nd element (starting from ilo-1 till ihi-2) and store in
         * next location. Used to store imaginary parts of complex conjuate pair of eigen
         * values */
        add_negative_values_ilo_ihi(datatype, wi_in, *ilo, *ihi);

        sort_vector(datatype, "A", n, wi_in, 1);
        sort_vector(datatype, "A", n, wr_in, 1);
        sort_vector(datatype, "A", n, wr, 1);
        sort_vector(datatype, "A", n, wi, 1);
    }
    else
    {
        sort_vector(datatype, "A", n, w, 1);
        sort_vector(datatype, "A", n, wr_in, 1);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            /* Validate the eigen values returned by the api */
            void *work = NULL;
            float norm, norm1, norm2, eps;
            eps = fla_lapack_slamch("P");
            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm1, imatrix, work);
            saxpy_(&n, &s_n_one, wr, &i_one, wr_in, &i_one);
            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm, imatrix, work);
            resid1 = norm / (eps * norm1 * n);

            if(*ilo != *ihi)
            {
                compute_matrix_norm(datatype, NORM, n, i_one, wi_in, i_one, &norm2, imatrix, work);
                saxpy_(&n, &s_n_one, wi, &i_one, wi_in, &i_one);
                compute_matrix_norm(datatype, NORM, n, i_one, wi_in, i_one, &norm, imatrix, work);
                resid2 = norm / (eps * norm2 * n);
            }
            break;
        }

        case DOUBLE:
        {
            /* Validate the eigen values returned by the api */
            void *work = NULL;
            double norm, norm1, norm2, eps;
            eps = fla_lapack_dlamch("P");
            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm1, imatrix, work);
            daxpy_(&n, &d_n_one, wr, &i_one, wr_in, &i_one);
            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm, imatrix, work);
            resid1 = norm / (eps * norm1 * n);

            if(*ilo != *ihi)
            {
                compute_matrix_norm(datatype, NORM, n, i_one, wi_in, i_one, &norm2, imatrix, work);
                daxpy_(&n, &d_n_one, wi, &i_one, wi_in, &i_one);
                compute_matrix_norm(datatype, NORM, n, i_one, wi_in, i_one, &norm, imatrix, work);
                resid2 = norm / (eps * norm2 * n);
            }

            break;
        }
        case COMPLEX:
        {
            /* Validate the eigen values returned by the api */
            void *work = NULL;
            float norm, norm1, eps;
            eps = fla_lapack_slamch("P");

            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm1, imatrix, work);
            caxpy_(&n, &c_n_one, w, &i_one, wr_in, &i_one);
            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm, imatrix, work);
            norm = fla_lapack_clange("1", &n, &i_one, wr_in, &i_one, work);
            resid1 = norm / (eps * norm1 * n);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Validate the eigen values returned by the api */
            void *work = NULL;
            double norm, norm1, eps;
            eps = fla_lapack_dlamch("P");

            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm1, imatrix, work);
            zaxpy_(&n, &z_n_one, w, &i_one, wr_in, &i_one);
            compute_matrix_norm(datatype, NORM, n, i_one, wr_in, i_one, &norm, imatrix, work);
            resid1 = norm / (eps * norm1 * n);
            break;
        }
    }
    residual = fla_test_max(resid1, resid2);

    if(*job == 'E' || *compz == 'N')
    {
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
        FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
        FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
        return;
    }

    void *zlambda = NULL, *work = NULL, *lambda = NULL;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &zlambda, n);
    reset_matrix(datatype, n, n, zlambda, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &lambda, n);
    reset_matrix(datatype, n, n, lambda, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_H, eps;
            eps = fla_lapack_slamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
            sgemm_("T", "N", &n, &n, &n, &s_one, Z, &ldz, Z_test, &ldz, &s_zero, lambda, &n);
            sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, H_test, &ldh, &s_zero, zlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, lambda, &n, &s_n_one, H, &ldh);
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
            resid3 = norm / (eps * norm_H * (float)n);
            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, lambda, n, n, n);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_H, eps;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
            dgemm_("T", "N", &n, &n, &n, &d_one, Z, &ldz, Z_test, &ldz, &d_zero, lambda, &n);
            dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, H_test, &ldh, &d_zero, zlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, lambda, &n, &d_n_one, H, &ldh);
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
            resid3 = norm / (norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, lambda, n, n, n);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_H, eps;
            eps = fla_lapack_slamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
            cgemm_("C", "N", &n, &n, &n, &c_one, Z, &ldz, Z_test, &ldz, &c_zero, lambda, &n);
            cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, H_test, &ldh, &c_zero, zlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, lambda, &n, &c_n_one, H, &ldh);
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
            resid3 = norm / (norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, lambda, n, n, n);

            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_H, eps;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm_H, imatrix, work);
            zgemm_("C", "N", &n, &n, &n, &z_one, Z, &ldz, Z_test, &ldz, &z_zero, lambda, &n);
            zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, H_test, &ldh, &z_zero, zlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, lambda, &n, &z_n_one, H, &ldh);
            compute_matrix_norm(datatype, NORM, n, n, H, ldh, &norm, imatrix, work);
            resid3 = norm / (norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, lambda, n, n, n);
            break;
        }
    }
    free_matrix(zlambda);
    free_matrix(lambda);

    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
}
