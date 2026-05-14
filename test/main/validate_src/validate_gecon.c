/******************************************************************************
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gels.c
 *  @brief Defines validate function of GECON() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

void validate_gecon(char *tst_api, integer datatype, char norm, integer n, void *A, void *A_save,
                    integer lda, double err_thresh, char imatrix_char, void *params)
{
    double residual, resid1 = 0., resid2 = 0.;

    /* Early return conditions */
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    /* Test 1: Check if matrix A was not modified by the call to GECON */
    resid1 = compare_matrix(datatype, "Full", n, n, A, lda, A_save, lda);

    /* Test 2: Check padding rows not modified */
    resid2 = check_padding(datatype, n, n, A, lda);

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
