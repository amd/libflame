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
/* Initializing asymmetrix matrix with values around overflow underflow (GEEV)*/
void init_matrix_overflow_underflow_asym(integer datatype, integer m, integer n, void *A,
                                         integer lda, char imatrix, void *scal);
/* Calculating the scaling value with respect to max and min for SVD */
void calculate_svd_scale_value(integer datatype, integer m, integer n, void *A, integer lda,
                               char imatrix, void *scal);
/* Calculating the scaling value with respect to max and min for ASYM */
void calculate_asym_scale_value(integer datatype, integer m, integer n, void *A, integer lda,
                                char imatrix, void *scal);
/* Finding ratio between two real datatype values */
void compute_ratio(integer datatype, void *scal, float flt_quotient, double dbl_quotient,
                   void *denominator);
/* Scaling matrix with values around overflow underflow for gerq2 */
void scale_matrix_overflow_underflow_gerq2(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scale matrix with values around overflow underflow for potrf */
void scale_matrix_overflow_underflow_potrf(integer datatype, integer m, void *A, integer lda,
                                           char imatrix);
/* Scaling matrix with values around overflow underflow for gelqf */
void scale_matrix_underflow_overflow_gelqf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for geqrf */
void scale_matrix_underflow_overflow_geqrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for larfg */
void scale_matrix_underflow_overflow_larfg(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
void scale_matrix_underflow_overflow_getrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
#endif // TEST_OVERFLOW_UNDERFLOW_H