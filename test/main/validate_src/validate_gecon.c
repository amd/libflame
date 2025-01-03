/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gels.c
 *  @brief Defines validate function of GECON() to use in test suite.
 *  */

#include "test_common.h"

void validate_gecon(integer datatype, char norm, integer n, void *A, void *A_save, integer lda,
                    double *residual, char imatrix_char)
{
    /* Check if matrix A was not modified by the call to GECON
     */
    if(compare_matrix(datatype, "Full", n, n, A, lda, A_save, lda) == 0)
    {
        *residual = DBL_MAX;
    }
    else
    {
        *residual = 0;
    }
}
