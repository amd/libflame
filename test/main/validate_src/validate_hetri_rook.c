/******************************************************************************
 * Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_hetri_rook.c
 *  @brief Defines validate function of HETRI_ROOK() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

void validate_hetri_rook(char *tst_api, char uplo, integer n, void *A, void *A_inv, integer lda,
                         integer *ipiv, integer datatype, double err_thresh, char imatrix,
                         void *params)
{
    double resid_tri = 0.;

    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
        return;
    }

    /* Test 3: Ensure unused triangle was not modified.
       Only run for normal validation (positive threshold); negative/error
       test cases use special err_thresh values handled by validate_getri. */
    if(err_thresh > 0)
    {
        resid_tri = compare_matrix(datatype, same_char(uplo, 'U') ? "L" : "U", n, n, A, lda, A_inv, lda);
    }

    /*If UPLO = 'U', the upper triangular part of the
      inverse is formed in upper triangular part of A_inv;
      if UPLO = 'L' the lower triangular part of the
      inverse is formed in lower triangular part of A_inv.
     */
    form_symmetric_matrix(datatype, n, A_inv, lda, "C", uplo);

    /* validate_getri prints the overall status and its own subtests (01, 02).
       Padding on A_inv is already checked by validate_getri as its Test 2. */
    validate_getri(tst_api, n, n, A, A_inv, lda, ipiv, datatype, err_thresh, imatrix, params);
    if(err_thresh > 0)
        FLA_PRINT_SUBTEST_STATUS(resid_tri, err_thresh, "03");
}
