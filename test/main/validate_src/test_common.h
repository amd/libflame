/******************************************************************************
 * Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file test_common.h
 *  @brief Defines function declarations to use in APIs of test suite.
 *  */

#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include "blis.h"
#include "test_prototype.h"
#include "validate_common.h"
#include "test_overflow_underflow.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

// global variables
extern integer i_zero, i_one, i_n_one;
extern float s_zero, s_one, s_n_one;
extern double d_zero, d_one, d_n_one;
extern scomplex c_zero, c_one, c_n_one;
extern dcomplex z_zero, z_one, z_n_one;

#define DRAND() ((double)rand() / ((double)RAND_MAX / 2.0F)) - 1.0F
#define SRAND() (float)((double)rand() / ((double)RAND_MAX / 2.0F)) - 1.0F
#define SRAND_IN_RANGE(lower, upper) \
    (float)(lower + (upper - lower) * ((double)rand() / (double)RAND_MAX))
#define DRAND_IN_RANGE(lower, upper) \
    (double)(lower + (upper - lower) * ((double)rand() / (double)RAND_MAX))
#define FLA_FABS(x) ((x) >= 0) ? (x) : (-1 * x)

// Datatype
#define CONSTANT 101
#define INTEGER 102
#define FLOAT 103
#define DOUBLE 104
#define COMPLEX 105
#define DOUBLE_COMPLEX 106
#define INVALID_TYPE -106

#define MAX_FLT_DIFF 0.00001 // Maximum allowed difference for float comparision
#define MAX_DBL_DIFF 0.0000000001 // Maximum allowed difference for double comparision

#if defined(FLA_ENABLE_ILP64)
#ifdef _WIN32
#define FT_IS "lld"
#else
#define FT_IS "ld"
#endif
#else
#define FT_IS "d"
#endif

/* Integer absolute function */
integer fla_i_abs(integer *x);
/* vector functions*/
void create_vector(integer datatype, void **A, integer M);
void create_realtype_vector(integer datatype, void **A, integer M);
void free_vector(void *A);
void reset_vector(integer datatype, void *A, integer M, integer incA);
/* Initialize vector,
 * if range = V then initialize real type vector with random values between given range (VL, VU)
 * if range = U then initialize real type vector within specific range with uniform initialization
 * if range = R initialize vector with random values
 */
void rand_vector(integer datatype, integer M, void *A, integer LDA, double VL, double VU,
                 char range);
void copy_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB);
void copy_subvector(integer datatype, integer m, void *A, integer lda, void *B, integer ldb,
                    integer srow, integer scol, integer drow, integer dcol);
void swap_row_col(integer datatype, integer *m, void *A, integer lda, integer *incx, integer *incy,
                  integer srow, integer scol, integer drow, integer dcol);
void copy_realtype_vector(integer datatype, integer M, void *A, integer LDA, void *B, integer LDB);

/* matrix functions*/
void create_matrix(integer datatype, void **A, integer M, integer N);
void create_realtype_matrix(integer datatype, void **A, integer M, integer N);
integer get_datatype(char stype);
integer get_realtype(integer datatype);
double get_realtype_value(integer datatype, void *value);
void create_block_diagonal_matrix(integer datatype, void *wr, void *wi, void *lambda, integer m,
                                  integer n, integer lda);
void *get_m_ptr(integer datatype, void *A, integer M, integer N, integer LDA);
void free_matrix(void *A);
void rand_matrix(integer datatype, void *A, integer M, integer N, integer LDA);
void rand_sym_matrix(integer datatype, void *A, integer M, integer N, integer LDA);
void rand_spd_matrix(integer datatype, char *uplo, void *A, integer m, integer lda);
void rand_hermitian_matrix(integer datatype, integer n, void **A, integer lda);
void copy_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA, void *B,
                 integer LDB);
void copy_realtype_matrix(integer datatype, char *uplo, integer M, integer N, void *A, integer LDA,
                          void *B, integer LDB);
void reset_matrix(integer datatype, integer M, integer N, void *A, integer LDA);
void set_identity_matrix(integer datatype, integer M, integer N, void *A, integer LDA);
void copy_submatrix(integer datatype, integer m, integer n, void *A, integer lda, void *B,
                    integer ldb, integer srow, integer scol, integer drow, integer dcol);
void copy_realtype_subvector(integer datatype, integer m, void *A, void *B, integer index);
void assign_value(integer datatype, void *x, double data_real, double data_imag);

void matrix_difference(integer datatype, integer m, integer n, void *A, integer lda, void *B,
                       integer ldb);

/* GEMM implementation for  C := alpha*op( A )*op( B ) + beta*C
 * Where alpha = 1, beta = 0
 */
void fla_invoke_gemm(integer datatype, char *transA, char *transB, integer *m, integer *n,
                     integer *k, void *A, integer *lda, void *B, integer *ldb, void *C,
                     integer *ldc);
/* orthgonality property of matrix */
double check_orthogonal_matrix(char trn, integer datatype, void *A, integer m, integer n, integer k,
                               integer lda);
double check_orthogonality(integer datatype, void *A, integer m, integer n, integer lda);
void get_diagonal(integer datatype, void *A, integer m, integer n, integer lda, void *Diag);
void get_subdiagonal(integer datatype, void *A, integer m, integer n, integer lda, void *Subdiag);
/*Tridiagonal matrix functions*/
void rand_sym_tridiag_matrix(integer datatype, void *A, integer M, integer N, integer LDA);
void copy_sym_tridiag_matrix(integer datatype, void *D, void *E, integer M, integer N, void *A,
                             integer LDA);
void copy_tridiag_matrix(integer datatype, void *dl, void *d, void *du, integer M, integer N,
                         void *A, integer LDA);
void copy_tridiag_vector(integer datatype, void *dl, void *d, void *du, integer M, integer N,
                         void *A, integer LDA);
void tridiag_matrix_multiply(integer datatype, integer n, integer nrhs, void *dl, void *d, void *du,
                             void *B, integer ldb, void *C, integer ldc);
void copy_sym_tridiag_matrix(integer datatype, void *D, void *E, integer M, integer N, void *B,
                             integer LDA);

/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp);
void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp);

/* work value calculation */
integer get_work_value(integer datatype, void *work);

/* Diagonal Scaling*/
void diagmv(integer datatype, integer m, integer n, void *x, integer incx, void *a, integer a_rs,
            integer a_cs);
void scalv(integer datatype, integer n, void *x, integer incx, void *y, integer incy);

/* set Transpose based on uplo */
void set_transpose(integer datatype, char *uplo, char *trans_A, char *trans_B);

/* Create diagonal matrix by copying elements from vector to matrix */
void diagonalize_realtype_vector(integer datatype, void *s, void *sigma, integer m, integer n,
                                 integer LDA);
/* To calculate matrix multiplication with real and complex datatypes */
void scgemv(char TRANS, integer real_alpha, integer m, integer n, scomplex *alpha, float *a,
            integer lda, scomplex *v, integer incv, float beta, scomplex *c, integer inc);

/* To find the maximum from the array */
void get_max_from_array(integer datatype, void *arr, void *max_val, integer n);
/* To find the minimum from the array */
void get_min_from_array(integer datatype, void *arr, void *min_val, integer n);
/* Reading matrix input data from a file */
void init_matrix_from_file(integer datatype, void *A, integer m, integer n, integer lda,
                           FILE *fptr);
/* Reading vector input data from a file */
void init_vector_from_file(integer datatype, void *A, integer m, integer inc, FILE *fptr);
/* Intialize vector with special values */
void init_vector_spec_in(integer datatype, void *A, integer M, integer incx, char type);
/* Allocate dynamic memory. If FLA_MEM_UNALIGNED is set, unaligned memory is allocated */
char *fla_mem_alloc(size_t size);
/* Generate Hessenberg matrix */
void get_hessenberg_matrix(integer datatype, integer n, void *A, integer lda, void *Z, integer ldz,
                           integer *ilo, integer *ihi, integer *info, bool AInitialized);
/* Generate Hessenberg matrix from eigen values */
void get_hessenberg_matrix_from_EVs(integer datatype, integer n, void *A, integer lda, void *Z,
                                    integer ldz, integer *ilo, integer *ihi, integer *info,
                                    bool AInitialized, void *wr_in, void *wi_in);
/* Convert matrix to upper hessenberg form */
void convert_upper_hessenberg(integer datatype, integer n, void *A, integer lda);
/* Pack a symmetric matrix in column first order */
void pack_matrix_lt(integer datatype, void *A, void *B, integer N, integer lda);
/* Convert matrix to upper hessenberg form */
void extract_upper_hessenberg_matrix(integer datatype, integer n, void *A, integer lda);
/* Convert matrix according to ILO and IHI values */
void get_generic_triangular_matrix(integer datatype, integer N, void *A, integer LDA, integer ilo,
                                   integer ihi, bool AInitialized);
/* Decompose matrix A in to QR and store orthogonal matrix in Q and R in A*/
void get_orthogonal_matrix_from_QR(integer datatype, integer n, void *A, integer lda, void *Q,
                                   integer ldq, integer *info);
/* Print matrix contents for visual inspection
   if order == 'C' matrix will be printed in columns first order
   else if order == 'R' matrix will be printed in rows first order*/
void print_matrix(char *desc, char *order, integer datatype, integer M, integer N, void *A,
                  integer lda);
/* Get upper triangular matrix or lower triangular matrix based on UPLO */
void get_triangular_matrix(char *uplo, integer datatype, integer m, integer n, void *A,
                           integer lda);
/*To Check order of Singular values of SVD (positive and non-decreasing)*/
double svd_check_order(integer datatype, void *s, integer m, integer n, double residual);
/*Generate Matrix for SVD*/
void create_svd_matrix(integer datatype, char range, integer m, integer n, void *A_input,
                       integer lda, void *S, double vl, double vu, integer il, integer iu,
                       integer info);
void get_abs_vector_value(integer datatype, void *S, integer M, integer inc);
/* Intialize matrix with special values*/
void init_matrix_spec_in(integer datatype, void *A, integer M, integer N, integer LDA, char type);
/*Intialize matrix according to given input*/
void init_matrix(integer datatype, void *A, integer M, integer N, integer LDA, FILE *g_ext_fptr,
                 char imatrix_char);
/* Intialize matrix with special values in random locations */
void init_matrix_spec_rand_in(integer datatype, void *A, integer M, integer N, integer LDA,
                              char type);
/*Test to check the extreme values propagation in output matrix */
bool check_extreme_value(integer datatype, integer M, integer N, void *A, integer LDA, char type);
/*Intialize vector according to given input*/
void init_vector(integer datatype, void *A, integer M, integer incx, FILE *g_ext_fptr,
                 char ivector_char);
/*Initialize vector with special values*/
void init_vector_spec_in(integer datatype, void *A, integer M, integer incx, char type);
/* Checks whether the value is zero or not */
double is_value_zero(integer datatype, void *value, double residual);
/* Multiply general m * n matrix with diagonal vector (of an n * n diagonal matrix) of size n */
void multiply_matrix_diag_vector(integer datatype, integer m, integer n, void *A, integer lda,
                                 void *X, integer incx);
/* Generate square matrix of size n x n using Eigen decomposition(ED) */
void generate_matrix_from_ED(integer datatype, integer n, void *A, integer lda, void *Q,
                             void *lambda);
/* Sort the real type vector in given order */
void sort_realtype_vector(integer datatype, char *order, integer vect_len, void *w, integer incw);
/* Compare two vectors starting from offset_A in A vector with B vector (starting from offset 0 in
 * B) */
integer compare_realtype_vector(integer datatype, integer vect_len, void *A, integer inca,
                                integer offset_A, void *B, integer incb);
/* Create input matrix A(symmetric) by randomly generating eigen values(EVs) in given range (vl,vu)
 */
void generate_matrix_from_EVs(integer datatype, char range, integer n, void *A, integer lda,
                              void *L, double vl, double vu);
/* Initialize band matrix with random values.
Note: Input buffer A has to be allocated by caller.*/
void rand_band_matrix(integer datatype, integer M, integer N, integer kl, integer ku, void *A,
                      integer LDA);
/* Initialize band storage for given band matrix.
Note: Input buffer A has to be allocated, initialized with band matrix by caller.*/
void get_band_storage_matrix(integer datatype, integer M, integer N, integer kl, integer ku,
                             void *A, integer LDA, void *AB, integer LDAB);
/* Initialize band storage with random band matrix.
Note: Input buffer AB has to be allocated by caller.*/
void rand_band_storage_matrix(integer datatype, integer M, integer N, integer kl, integer ku,
                              void *AB, integer LDAB);
/* Get band matrix from band storage matrix.
Note: Input buffer A has to be allocated, initialized with band storage matrix by caller.*/
void get_band_matrix_from_band_storage(integer datatype, integer M, integer N, integer kl,
                                       integer ku, void *AB, integer LDAB, void *A, integer LDA);
/* On input, AB is the output of GBTRF().
   On output, AB is the reconstructed band storage matrix same as input of GBTRF().*/
void reconstruct_band_storage_matrix(integer datatype, integer m, integer n, integer kl, integer ku,
                                     void *AB, integer ldab, integer *ipiv);
/* Test for checking whether solution x of Ax = B from least square api belongs to row space of A*/
void check_vector_in_rowspace(integer datatype, char *trans, integer m, integer n, integer nrhs,
                              void *A, integer lda, void *x, integer ldb, void *resid);
/* Compute norm for matrix/vectors
 * Choose 2-norm or 1-norm based on type of test:
 *     Overflow/Underflow vs Normal.
 * ntype: Choose norm type for normal tests (refer lange API).
 * imatrix if 'U', 2-norm is computed.
 */
void compute_matrix_norm(integer datatype, char ntype, integer m, integer n, void *A, integer lda,
                         void *nrm2, char imatrix, void *work);
/* To calculate the resudial sum of squares of solution for solution x of Ax = b and m < n*/
void residual_sum_of_squares(int datatype, integer m, integer n, integer nrhs, void *x, integer ldx,
                             double *resid);
/* Generate a symmetric or hermitian matrix from existing matrix A
 * If type = "C" hermitian matrix formed.
 * If type = "S" symmetric matrix is formed.
 */
void form_symmetric_matrix(integer datatype, integer n, void *A, integer lda, char *type);
/* Scaling the matrix by x scalar */
void scal_matrix(integer datatype, void *x, void *A, integer m, integer n, integer lda,
                 integer inc);
/* Get the maximum value from the matrix */
void get_max_from_matrix(integer datatype, void *A, void *max_val, integer m, integer n,
                         integer lda);
/* Get the minimum value from the matrix */
void get_min_from_matrix(integer datatype, void *A, void *min_val, integer m, integer n,
                         integer lda);
/* Sort the given vector in specified order */
void sort_vector(integer datatype, char *order, integer vect_len, void *w, integer incw);
/* Generate a block diagonal matrix with complex conjugate eigen value pairs as
   2 * 2 blocks along the diagonal. This is used for generating asymmetric matrix */
void create_realtype_block_diagonal_matrix(integer datatype, void *A, integer n, integer lda);
/* Create input matrix A(Asymmetric) by randomly generating eigen values(EVs) */
void generate_asym_matrix_from_EVs(integer datatype, integer n, void *A, integer lda, void *L,
                                   char imatrix, void *scal);
/* Generate asymmetric square matrix of size n x n using Eigen decomposition(ED) */
void generate_asym_matrix_from_ED(integer datatype, integer n, void *A, integer lda, void *Q,
                                  void *lambda);
/* Compare two vectors starting from offset_A in A vector with B vector */
integer compare_vector(integer datatype, integer vect_len, void *A, integer inca, integer offset_A,
                       void *B, integer incb);
/* Create diagonal matrix by copying elements from a vector to matrix */
void diagonalize_vector(integer datatype, void *s, void *sigma, integer m, integer n, integer LDA);
/* Find negative value of each element and store in next location
   Used to store imaginary parts of complex conjuate pair of eigen values
   in asymmetric matrix eigen decomposition APIs
   Ex: input vector {a, 0, -b, 0 ...}
       output vector {a, -a, -b, b, ...} */
void add_negative_values(integer datatype, void *vect, integer n);
void add_negative_values_ilo_ihi(integer datatype, void *vect, integer ilo, integer ihi);
/* Convert the given matrix from column major layout to row major layout and vice versa */
void convert_matrix_layout(integer matrix_layout, integer datatype, integer m, integer n, void *a,
                           integer lda, integer *lda_t);
#endif // TEST_COMMON_H
