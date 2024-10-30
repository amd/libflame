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
/* Scaling matrix with values around overflow underflow for larf */
void scale_matrix_underflow_overflow_larf(integer datatype, integer m, integer n, void *A,
                                          integer lda, char imatrix_char);
void scale_matrix_underflow_overflow_getrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scale matrix with values around overflow underflow for hetrf_rook */
void scale_matrix_overflow_underflow_hetrf_rook(integer datatype, integer m, void *A, integer lda,
                                                char imatrix);
/* Scale matrix with values around overflow underflow for sytrf */
void scale_matrix_overflow_underflow_sytrf(integer datatype, integer m, void *A, integer lda,
                                           char imatrix);
/* Scale matrix with values around overflow underflow for org2r */
void scale_matrix_underflow_overflow_org2r(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for geqrf */
void scale_matrix_underflow_overflow_geqrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for larfg */
void scale_matrix_underflow_overflow_larfg(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
void scale_matrix_underflow_overflow_getrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scale matrix with values around overflow underflow for sytrf */
void scale_matrix_overflow_underflow_sytrf(integer datatype, integer m, void *A, integer lda,
                                           char imatrix);
/* Scale matrix with values around overflow underflow for org2r */
void scale_matrix_underflow_overflow_org2r(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for gtsv */
void init_matrix_overflow_underflow_gtsv(integer datatype, integer m, integer n, void *A,
                                         integer lda, char imatrix, void *scal);
/* Scaling matrix with values around overflow underflow for gels */
void scale_matrix_underflow_overflow_gels(integer datatype, char *trans, integer m, integer n,
                                          void *A, integer lda, char imatrix_char, integer sysmat);
/* Scaling matrix with values around overflow underflow for gelsd */
void scale_matrix_overflow_underflow_gelsd(integer datatype, integer m, integer n, integer nrhs,
                                           void *A, integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for gbtrf */
void scale_matrix_underflow_overflow_gbtrf(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for gbtrs */
void scale_matrix_underflow_overflow_gbtrs(integer datatype, integer m, integer nrhs, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow underflow for gelss */
void scale_matrix_overflow_underflow_gelss(integer datatype, integer m, integer n, integer nrhs,
                                           void *A, integer lda, char imatrix_char);
/* Scaling matrix with values around overflow, underflow for STEDC */
void scale_matrix_underflow_overflow_stedc(integer datatype, integer n, void *A, integer lda,
                                           char *imatrix_char, char *scal);
/* Scaling matrix with values around overflow, underflow for STEVD */
void scale_matrix_underflow_overflow_stevd(integer datatype, integer n, void *A, integer lda,
                                           char *imatrix_char, char *scal);
/* Scaling matrix with values around overflow, underflow for SYEV/HEEV */
void scale_matrix_underflow_overflow_syev(integer datatype, integer n, void *A, integer lda,
                                          char *imatrix_char, char *scal);
/* Scaling matrix with values around overflow, underflow for STEQR */
void scale_matrix_underflow_overflow_steqr(integer datatype, integer n, void *A, integer lda,
                                           char *imatrix_char, char *scal);
/* Scale matrix with values around overflow underflow for ggev */
void scale_matrix_overflow_underflow_ggev(integer datatype, integer m, void *A, integer lda,
                                          char imatrix);
/* Scaling matrix with values around overflow underflow for syevx */
void scale_matrix_underflow_overflow_syevx(integer datatype, integer n, void *A, integer lda,
                                           char imatrix_char, void *scal);
/* Scale matrix with values around overflow underflow for gesv */
void scale_matrix_underflow_overflow_gesv(integer datatype, integer n, void *A, integer lda,
                                          char imatrix_char, void *scal);
/* Scaling matrix with values around overflow, underflow for GETRS */
void scale_matrix_underflow_overflow_getrs(integer datatype, char *trans, integer m, integer n,
                                           void *A, integer lda, char imatrix_char, void *scal);
/* Scaling matrix with values around overflow underflow for getri */
void scale_matrix_underflow_overflow_getri(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow, underflow for GEQP3 */
void scale_matrix_underflow_overflow_geqp3(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow, underflow for SYEVD/HEEVD */
void scale_matrix_underflow_overflow_syevd(integer datatype, integer n, void *A, integer lda,
                                           char *imatrix_char, void *scal);
/* Scale matrix with values around overflow underflow for gesv */
void scale_matrix_underflow_overflow_gesdd(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char, void *scal);
/* Scaling matrix with values around overflow, underflow for ORGQR/UNGQR */
void scale_matrix_underflow_overflow_orgqr(integer datatype, integer m, integer n, void *A,
                                           integer lda, char imatrix_char);
/* Scaling matrix with values around overflow, underflow for POTRS */
void scale_matrix_underflow_overflow_potrs(integer datatype, integer n, void *A, integer lda,
                                           char imatrix_char, void *scal);
/* Scale matrix with values around overflow underflow for hgeqz A matrix */
void scale_matrix_overflow_underflow_hgeqz_A(integer datatype, integer n, void *A, integer lda,
                                             char imatrix_char, void *scal);
/* Scale matrix with values around overflow underflow for hgeqz B matrix */
void scale_matrix_overflow_underflow_hgeqz_B(integer datatype, integer n, void *A, integer lda,
                                             char imatrix_char, void *scal);
/* Scale matrix with values around overflow underflow for hseqr */
void scale_matrix_overflow_underflow_hseqr(integer datatype, integer n, void *A, integer lda,
                                           char imatrix_char, void *scal);
/* Scaling matrix with values around overflow underflow for GEHRD */
void scale_matrix_underflow_overflow_gehrd(integer datatype, integer n, void *A, integer lda,
                                           char imatrix_char);
/* Scaling matrix with values around overflow underflow for GGHRD */
void scale_matrix_underflow_overflow_gghrd(integer datatype, integer n, void *A, integer lda,
                                           char imatrix_char);
/* Scaling matrix with values around overflow, underflow for SYGVD/HEGVD */
void scale_matrix_underflow_overflow_sygvd(integer datatype, integer n, void *A, integer lda,
                                           void *B, integer ldb, integer itype, char imatrix_char,
                                           void *scal);
#endif // TEST_OVERFLOW_UNDERFLOW_H
