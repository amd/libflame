/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gels.c
 *  @brief Defines validate function of GECON() to use in test suite.
 *  */

#include "test_common.h"

void validate_gecon(char *tst_api, integer datatype, char norm, integer n, void *A, void *A_save,
                    integer lda, double err_thresh, char imatrix_char)
{
    double residual;

    /* Early return conditions */
    if(n == 0 || n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    /* Check if matrix A was not modified by the call to GECON
     */
    if(compare_matrix(datatype, "Full", n, n, A, lda, A_save, lda) == 0)
    {
        residual = DBL_MAX;
    }
    else
    {
        residual = 0;
    }

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}
