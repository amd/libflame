/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file test_overflow_underflow.h
 *  @brief Defines function declarations for overflow and underflow test to use in APIs of test
 * suite.
 *  */

#ifndef TEST_OVERFLOW_UNDERFLOW_H
#define TEST_OVERFLOW_UNDERFLOW_H

#include "test_common.h"

#include <float.h>

/* Initializing matrix with values around overflow underflow */
void init_matrix_overflow_underflow_svd(integer datatype, integer m, integer n, void *A,
                                        integer lda, char imatrix, void *scal);
/* Calculating the scaling value with respect to max and min for SVD */
void calculate_svd_scale_value(integer datatype, integer m, integer n, void *A, integer lda,
                               char imatrix, void *scal);
/* Finding ratio between two real datatype values */
void compute_ratio(integer datatype, void *scal, float flt_quotient, double dbl_quotient,
                   void *denominator);
/* Scaling matrix with values around overflow underflow for gerq2 */
void scale_matrix_overflow_underflow_gerq2(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
#endif // TEST_OVERFLOW_UNDERFLOW_H