/*
    Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gbtrf.c
 *  @brief Defines validate function of GBTRF() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_gbtrf(char *tst_api, integer m_A, integer n_A, integer kl, integer ku, void *AB,
                    void *AB_test, integer ldab, integer *IPIV, integer datatype, double err_thresh,
                    void *params)
{
    void *work = NULL;
    double residual, resid1 = 0., resid2 = 0.;
    integer m_band = 2 * kl + ku + 1;

    /* Early return conditions */
    if(m_A == 0 || n_A == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m_A, n_A, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m_A, n_A, err_thresh);

    /* Test 1 - Check for input AB.
       Use m_band (= 2*kl+ku+1) as the row count instead of ldab so that the
       norm/difference operate only on actual band data and never touch the
       padding sentinels in rows m_band..ldab-1 (which Test 2 verifies below). */
    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_B;

            norm_B = fla_lapack_slange("1", &m_band, &n_A, AB, &ldab, work);
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, m_band, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_slange("1", &m_band, &n_A, AB_test, &ldab, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_B, m_A, params);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_B;

            norm_B = fla_lapack_dlange("1", &m_band, &n_A, AB, &ldab, work);
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, m_band, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_dlange("1", &m_band, &n_A, AB_test, &ldab, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_B, m_A, params);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_B;

            norm_B = fla_lapack_clange("1", &m_band, &n_A, AB, &ldab, work);
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, m_band, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_clange("1", &m_band, &n_A, AB_test, &ldab, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_B, m_A, params);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_B;

            norm_B = fla_lapack_zlange("1", &m_band, &n_A, AB, &ldab, work);
            reconstruct_band_storage_matrix(datatype, m_A, n_A, kl, ku, AB_test, ldab, IPIV);
            matrix_difference(datatype, m_band, n_A, AB_test, ldab, AB, ldab);
            norm = fla_lapack_zlange("1", &m_band, &n_A, AB_test, &ldab, work);
            resid1 = fla_compute_residual(datatype, 'E', norm, norm_B, m_A, params);
            break;
        }
        default:
            resid1 = err_thresh;
            break;
    }

    /* Test 2: Check padding rows not modified.
       Band storage uses m_band = 2*kl+ku+1 rows for actual data; anything in
       rows m_band..ldab-1 must still hold the sentinel pattern. */
    resid2 = check_padding(datatype, m_band, n_A, AB_test, ldab);

    residual = fla_test_max(resid1, resid2);

    FLA_PRINT_TEST_STATUS(m_A, n_A, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
