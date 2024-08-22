/******************************************************************************
 * Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_orgqr.c
 *  @brief Defines validate function of ORGQR() to use in test suite.
 *  */

#include "test_common.h"

void validate_orgqr(integer m, integer n, void *A, integer lda, void *Q, void *R, void *work,
                    integer datatype, double *residual, integer *info, char imatrix)
{
    if(m == 0 || n == 0)
        return;
    integer k;
    *info = 0;

    k = m;
    char NORM = '1';

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            sgemm_("T", "N", &n, &n, &k, &s_n_one, Q, &lda, A, &lda, &s_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, n, lda);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            dgemm_("T", "N", &n, &n, &k, &d_n_one, Q, &lda, A, &lda, &d_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, n, lda);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            cgemm_("C", "N", &n, &n, &k, &c_n_one, Q, &lda, A, &lda, &c_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, n, lda);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm_A, imatrix, work);
            zgemm_("C", "N", &n, &n, &k, &z_n_one, Q, &lda, A, &lda, &z_one, R, &n);

            compute_matrix_norm(datatype, NORM, n, n, R, n, &norm, imatrix, work);
            resid1 = (norm / norm_A) / (eps * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, n, lda);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
}
