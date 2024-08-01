/******************************************************************************
 * Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gesv.c
 *  @brief Defines validate function of GESV() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesv(integer n, integer nrhs, void *A, integer lda, void *B, integer ldb, void *X,
                   integer datatype, double *residual, integer *info, char imatrix, void *scal)
{
    if(n == 0 || nrhs == 0)
        return;
    void *work = NULL;
    char NORM = '1';
    integer ldx;
    *info = 0;
    ldx = ldb;

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_x, norm, eps, resid;

            /* Test 1 */
            /* Compute AX-B */
            sgemm_("N", "N", &n, &nrhs, &n, &s_one, A, &lda, X, &ldx, &s_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    sscal_(&n, scal, X, &i_one);
                }
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");

            resid = norm / (norm_x * n * eps);

            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
            double norm_x, norm, eps, resid;

            /* Test 1 */
            /* Compute AX-B */
            dgemm_("N", "N", &n, &nrhs, &n, &d_one, A, &lda, X, &ldx, &d_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    dscal_(&n, scal, X, &i_one);
                }
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid = norm / (norm_x * n * eps);

            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
            float norm_x, norm, eps, resid;

            /* Test 1 */
            /* Compute AX-B */
            cgemm_("N", "N", &n, &nrhs, &n, &c_one, A, &lda, X, &ldx, &c_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    cscal_(&n, scal, X, &i_one);
                }
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_slamch("E");
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid = norm / (norm_x * n * eps);

            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_x, norm, eps, resid;

            /* Test 1 */
            /* Compute AX-B */
            zgemm_("N", "N", &n, &nrhs, &n, &z_one, A, &lda, X, &ldx, &z_n_one, B, &ldb);
            if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    zscal_(&n, scal, X, &i_one);
                }
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
            eps = fla_lapack_dlamch("E");
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);

            resid = norm / (norm_x * n * eps);

            *residual = (double)resid;
            break;
        }
    }
}