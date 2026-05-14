/*
    Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gbtrs.c
 *  @brief Defines validate function of GBTRS() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_gbtrs(char *tst_api, char *trans, integer n, integer kl, integer ku, integer nrhs,
                    void *AB_orig, void *AB, void *AB_save, integer ldab, integer *IPIV,
                    integer *IPIV_save, void *B, integer ldb, void *X, integer datatype,
                    double err_thresh, char imatrix, void *params)
{
    void *A = NULL, *work = NULL;
    char NORM = '1';
    double residual, resid1 = 0., resid2 = 0., resid3 = 0.;
    integer i;

    /* Early return conditions */
    if(n == 0 || nrhs == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    /* Test 1: Compute residual ||AX - B|| / (||X|| * n * eps)
       Convert original band storage to full matrix for residual computation.
       Note: Unlike GETRS, GBTRS does not perform overflow/underflow scaling
       on X, so no scaling correction is applied before computing norm_x.
       The denominator uses only ||X|| (consistent with GETRS), not ||A||. */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, n);
    reset_matrix(datatype, n, n, A, n);
    get_band_matrix_from_band_storage(datatype, n, n, kl, ku, AB_orig, ldab, A, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_x, norm;

            sgemm_(trans, "N", &n, &nrhs, &n, &s_one, A, &n, X, &ldb, &s_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        case DOUBLE:
        {
            double norm_x, norm;

            dgemm_(trans, "N", &n, &nrhs, &n, &d_one, A, &n, X, &ldb, &d_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        case COMPLEX:
        {
            float norm_x, norm;

            cgemm_(trans, "N", &n, &nrhs, &n, &c_one, A, &n, X, &ldb, &c_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_x, norm;

            zgemm_(trans, "N", &n, &nrhs, &n, &z_one, A, &n, X, &ldb, &z_n_one, B, &ldb);
            compute_matrix_norm(datatype, NORM, n, nrhs, X, ldb, &norm_x, imatrix, work);
            compute_matrix_norm(datatype, NORM, n, nrhs, B, ldb, &norm, imatrix, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_x, n, params);
            break;
        }
        default:
            resid1 = err_thresh;
            break;
    }

    free_matrix(A);

    /* Test 2: Ensure factored band matrix AB was not modified by GBTRS */
    resid2 = compare_matrix(datatype, "full", ldab, n, AB_save, ldab, AB, ldab);

    /* Test 3: Ensure IPIV was not modified by GBTRS */
    if(IPIV_save != NULL)
    {
        for(i = 0; i < n; i++)
        {
            if(IPIV[i] != IPIV_save[i])
            {
                resid3 = DBL_MAX;
                break;
            }
        }
    }

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
}
