/******************************************************************************
 * Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_hetri_rook.c
 *  @brief Defines validate function of HETRI_ROOK() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

void validate_hetri_rook(char *tst_api, char uplo, integer n, void *A, void *A_inv, integer lda, integer *ipiv,
                         integer datatype, double err_thresh, char imatrix)
{ 
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
        return;
    }
    /*If UPLO = 'U', the upper triangular part of the
      inverse is formed in upper triangular part of A_inv;
      if UPLO = 'L' the lower triangular part of the 
      inverse is formed in lower triangular part of A_inv.
     */
    form_symmetric_matrix(datatype, n, A_inv, lda, "C", uplo);
    validate_getri(tst_api, n, n, A, A_inv, lda, ipiv, datatype, err_thresh, imatrix);
}
