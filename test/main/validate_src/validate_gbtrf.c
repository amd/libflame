/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gbtrf.c
 *  @brief Defines validate function of GBTRF() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_gbtrf(char *tst_api, integer m_A, integer n_A, integer kl, integer ku, void *AB,
                    void *AB_test, integer ldab, integer *IPIV, integer datatype, double err_thresh)
{
    void *work = NULL;
    double residual;

    /* Early return conditions */
    if(m_A == 0 || n_A == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m_A, n_A, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m_A, n_A, err_thresh);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_B, eps;

            eps = fla_lapack_slamch("Epsilon");

            /* Test 1 - Check for input AB. */
            norm_B = fla_lapack_slange("1", &ldab, &n_A, AB, &ldab, work);
            // Get reconstructed matrix and perform AB_test-AB
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, ldab, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_slange("1", &ldab, &n_A, AB_test, &ldab, work);
            residual = norm / (float)m_A / norm_B / eps;
            break;
        }
        case DOUBLE:
        {
            double norm, norm_B, eps;

            eps = fla_lapack_dlamch("Epsilon");

            /* Test 1 - Check for input AB. */
            norm_B = fla_lapack_dlange("1", &ldab, &n_A, AB, &ldab, work);
            // Get reconstructed matrix and perform AB_test-AB
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, ldab, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_dlange("1", &ldab, &n_A, AB_test, &ldab, work);
            residual = norm / (double)m_A / norm_B / eps;
            break;
        }
        case COMPLEX:
        {
            float norm, norm_B, eps;

            eps = fla_lapack_slamch("Epsilon");

            /* Test 1 - Check for input AB. */
            norm_B = fla_lapack_clange("1", &ldab, &n_A, AB, &ldab, work);
            // Get reconstructed matrix and perform AB_test-AB
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, ldab, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_clange("1", &ldab, &n_A, AB_test, &ldab, work);
            residual = norm / (float)m_A / norm_B / eps;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_B, eps;

            eps = fla_lapack_dlamch("Epsilon");

            /* Test 1 - Check for input AB. */
            norm_B = fla_lapack_zlange("1", &ldab, &n_A, AB, &ldab, work);
            // Get reconstructed matrix and perform AB_test-AB
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, ldab, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_zlange("1", &ldab, &n_A, AB_test, &ldab, work);
            residual = norm / (double)m_A / norm_B / eps;
            break;
        }
        default:
            residual = err_thresh;
            break;
    }

    FLA_PRINT_TEST_STATUS(m_A, n_A, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}
