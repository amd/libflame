/******************************************************************************
 * Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_hseqr.c
 *  @brief Defines validate function of HSEQR() to use in test suite.
 *  */

#include "test_common.h"

void validate_hseqr(char *job, char *compz, integer n, void *H, void *H_test, integer ldh, void *Z,
                    void *Z_test, integer ldz, void *wr, void *wr_in, void *wi, void *wi_in,
                    void *w, integer datatype, double *residual, integer *info, integer *ilo,
                    integer *ihi)
{
    *residual = 0;
    if(n == 0)
        return;

    if(datatype == FLOAT || datatype == DOUBLE)
    {
        /* Find negative value of every 2nd element (starting from ilo-1 till ihi-2) and s tore in
         * next location. Used to store imaginary parts of complex conjuate pair o           f eigen
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
            float norm, norm1, norm2, eps, resid;
            eps = fla_lapack_slamch("P");
            norm1 = fla_lapack_slange("1", &n, &i_one, wr_in, &i_one, work);
            saxpy_(&n, &s_n_one, wr, &i_one, wr_in, &i_one);
            norm = fla_lapack_slange("1", &n, &i_one, wr_in, &i_one, work);
            resid = norm / (eps * norm1 * n);

            *residual = fla_max(*residual, (double)resid);

            if(*ilo != *ihi)
            {
                norm2 = fla_lapack_slange("1", &n, &i_one, wi_in, &i_one, work);
                saxpy_(&n, &s_n_one, wi, &i_one, wi_in, &i_one);
                norm = fla_lapack_slange("1", &n, &i_one, wi_in, &i_one, work);
                resid = norm / (eps * norm2 * n);

                *residual = fla_max(*residual, (double)resid);
            }

            break;
        }

        case DOUBLE:
        {
            /* Validate the eigen values returned by the api */
            void *work = NULL;
            double norm, norm1, norm2, eps, resid;
            eps = fla_lapack_dlamch("P");
            norm1 = fla_lapack_dlange("1", &n, &i_one, wr_in, &i_one, work);
            daxpy_(&n, &d_n_one, wr, &i_one, wr_in, &i_one);
            norm = fla_lapack_dlange("1", &n, &i_one, wr_in, &i_one, work);
            resid = norm / (eps * norm1 * n);

            *residual = fla_max(*residual, resid);

            if(*ilo != *ihi)
            {
                norm2 = fla_lapack_dlange("1", &n, &i_one, wi_in, &i_one, work);
                daxpy_(&n, &d_n_one, wi, &i_one, wi_in, &i_one);
                norm = fla_lapack_dlange("1", &n, &i_one, wi_in, &i_one, work);
                resid = norm / (eps * norm2 * n);
                *residual = fla_max(*residual, resid);
            }

            break;
        }
        case COMPLEX:
        {
            /* Validate the eigen values returned by the api */
            void *work = NULL;
            float norm, norm1, eps, resid;
            eps = fla_lapack_slamch("P");
            norm1 = fla_lapack_clange("1", &n, &i_one, wr_in, &i_one, work);
            caxpy_(&n, &c_n_one, w, &i_one, wr_in, &i_one);
            norm = fla_lapack_clange("1", &n, &i_one, wr_in, &i_one, work);
            resid = norm / (eps * norm1 * n);

            *residual = fla_max(*residual, (double)resid);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* Validate the eigen values returned by the api */
            void *work = NULL;
            double norm, norm1, eps, resid;
            eps = fla_lapack_dlamch("P");
            norm1 = fla_lapack_zlange("1", &n, &i_one, wr_in, &i_one, work);
            zaxpy_(&n, &z_n_one, w, &i_one, wr_in, &i_one);
            norm = fla_lapack_zlange("1", &n, &i_one, wr_in, &i_one, work);
            resid = norm / (eps * norm1 * n);

            *residual = fla_max(*residual, resid);
            break;
        }
    }

    if(*job == 'E' || *compz == 'N')
        return;

    void *zlambda = NULL, *work = NULL, *lambda = NULL;
    *info = 0;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &zlambda, n);
    reset_matrix(datatype, n, n, zlambda, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &lambda, n);
    reset_matrix(datatype, n, n, lambda, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_slange("1", &n, &n, H, &ldh, work);
            sgemm_("T", "N", &n, &n, &n, &s_one, Z, &ldz, Z_test, &ldz, &s_zero, lambda, &n);
            sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, H_test, &ldh, &s_zero, zlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, lambda, &n, &s_n_one, H, &ldh);
            norm = fla_lapack_slange("1", &n, &n, H, &ldh, work);
            resid1 = norm / (eps * norm_H * (float)n);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, lambda, n, n, n);
            *residual = fla_max(*residual, fla_max((double)resid1, (double)resid2));
            break;
        }
        case DOUBLE:
        {
            double norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_dlange("1", &n, &n, H, &ldh, work);
            dgemm_("T", "N", &n, &n, &n, &d_one, Z, &ldz, Z_test, &ldz, &d_zero, lambda, &n);
            dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, H_test, &ldh, &d_zero, zlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, lambda, &n, &d_n_one, H, &ldh);
            norm = fla_lapack_dlange("1", &n, &n, H, &ldh, work);
            resid1 = norm / (norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = check_orthogonality(datatype, lambda, n, n, n);
            *residual = fla_max(*residual, fla_max(resid1, resid2));
            break;
        }
        case COMPLEX:
        {
            float norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_clange("1", &n, &n, H, &ldh, work);
            cgemm_("C", "N", &n, &n, &n, &c_one, Z, &ldz, Z_test, &ldz, &c_zero, lambda, &n);
            cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, H_test, &ldh, &c_zero, zlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, lambda, &n, &c_n_one, H, &ldh);
            norm = fla_lapack_clange("1", &n, &n, H, &ldh, work);
            resid1 = norm / (norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, lambda, n, n, n);
            *residual = fla_max(*residual, fla_max((double)resid1, (double)resid2));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_zlange("1", &n, &n, H, &ldh, work);
            zgemm_("C", "N", &n, &n, &n, &z_one, Z, &ldz, Z_test, &ldz, &z_zero, lambda, &n);
            zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, H_test, &ldh, &z_zero, zlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, lambda, &n, &z_n_one, H, &ldh);
            norm = fla_lapack_zlange("1", &n, &n, H, &ldh, work);
            resid1 = norm / (norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = check_orthogonality(datatype, lambda, n, n, n);
            *residual = fla_max(*residual, fla_max(resid1, resid2));
            break;
        }
    }
    free_matrix(zlambda);
    free_matrix(lambda);
}
