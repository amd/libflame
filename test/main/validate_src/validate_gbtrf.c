/*
    Copyright (C) 2024 Advanced Micro Devices, Inc.  All rights reserved. Portions of this file
   consist of AI-generated content.
*/

/*! @file validate_gbtrf.c
 *  @brief Defines validate function of GBTRF() to use in test suite.
 *  */

#include "test_common.h"

void validate_gbtrf(integer m_A, integer n_A, integer kl, integer ku, void *AB, void *AB_test,
                    integer ldab, integer *IPIV, integer datatype, double *residual, integer *info)
{
    void *work = NULL;

    *info = 0;

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
            *residual = norm / (float)m_A / norm_B / eps;
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
            *residual = norm / (double)m_A / norm_B / eps;
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
            *residual = norm / (float)m_A / norm_B / eps;
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
            *residual = norm / (double)m_A / norm_B / eps;
            break;
        }
    }
}
