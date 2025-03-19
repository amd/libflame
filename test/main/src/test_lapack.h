/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#ifndef TEST_LAPACK_H
#define TEST_LAPACK_H

#include <ctype.h>
#include <float.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#ifdef _WIN32
#include <Windows.h>
#endif
#include <inttypes.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#ifndef _WIN32
#include <sys/time.h>
#include <unistd.h>
#endif
#include "test_common.h"

#define LAPACK_OPERATIONS_FILENAME "input.global.operations"
#define LAPACKE_OPERATIONS_FILENAME "input.global.operations.lapacke"

#define COMMENT_CHAR '#'
#define MAX_BINARY_NAME_LENGTH 256
#define INPUT_BUFFER_SIZE 256
#define MAX_DT_STRING_LENGTH 32
#define MAX_STOR_STRING_LENGTH 32
#define MAX_FUNC_STRING_LENGTH 20
#define MAX_NUM_STORAGE 4
#define MAX_NUM_DATATYPES 4
#define FLOPS_PER_UNIT_PERF 1e9

#define MAX_NUM_NORMTYPES 4

#define DISABLE_ALL 0
#define SPECIFY 1
#define DISABLE 0
#define ENABLE 1

#define RECT_INPUT 0
#define SQUARE_INPUT 1

#define CLI_NORM_THRESH 30.0
#define CLI_DECIMAL_BASE 10

#define MAX_PASS_STRING_LENGTH 32

#define NUM_STORAGE_CHARS 3
#define STORAGE_SCHEME_CHARS "crg"

#define NUM_SUB_TESTS (4)

// API categories
#define LIN (1)
#define EIG_SYM (2)
#define EIG_NSYM (3)
#define SVD (4)
#define AUX (5)

// pass 1 to test  standard AOCL_FLA_PROGRESS fucntion,
// pass 2 to test  register callback function
#define AOCL_FLA_SET_PROGRESS_ENABLE 0

#if AOCL_FLA_SET_PROGRESS_ENABLE == 2
int test_progress(const char *const api, const integer lenapi, const integer *const progress,
                  const integer *const current_thread, const integer *const total_threads);
#endif

/* Flag to indicate lwork/liwork/lrwork availability status
 * <= 0 - To be calculated
 * > 0  - Use the value
 * */
extern integer g_lwork;
extern integer g_liwork;
extern integer g_lrwork;
/* Variable to indicate the source of inputs
 * = 0 - Inputs are from command line
 * = 1 - Inputs are from config file
 * */
extern integer config_data;
/* File pointer for external file which is used
 * to pass the input matrix values
 * */
extern FILE *g_ext_fptr;

#define FLA_TEST_PARSE_LAST_ARG(argv)                                                             \
    integer i;                                                                                    \
    char *info;                                                                                   \
    char info_value[2][MAX_PASS_STRING_LENGTH];                                                   \
                                                                                                  \
    i = 0;                                                                                        \
    if(strstr(argv, "--einfo") != NULL)                                                           \
    {                                                                                             \
        info = strtok(argv, "=");                                                                 \
        while(info != NULL && i < 2)                                                              \
        {                                                                                         \
            strcpy(info_value[i], info);                                                          \
            i++;                                                                                  \
            info = strtok(NULL, "=");                                                             \
        }                                                                                         \
        einfo = atoi(info_value[1]);                                                              \
    }                                                                                             \
    else if(strstr(argv, "--imatrix") != NULL)                                                    \
    {                                                                                             \
        info = strtok(argv, "=");                                                                 \
        while(info != NULL && i < 2)                                                              \
        {                                                                                         \
            strcpy(info_value[i], info);                                                          \
            i++;                                                                                  \
            info = strtok(NULL, "=");                                                             \
        }                                                                                         \
        if(strlen(info_value[1]) != 1                                                             \
           || (info_value[1][0] != 'I' && info_value[1][0] != 'N' && info_value[1][0] != 'A'      \
               && info_value[1][0] != 'F' && info_value[1][0] != 'U' && info_value[1][0] != 'O')) \
        {                                                                                         \
            printf("\n Invalid input for imatrix \n");                                            \
            return;                                                                               \
        }                                                                                         \
        params->imatrix_char = info_value[1][0];                                                  \
    }                                                                                             \
    else                                                                                          \
    {                                                                                             \
        g_ext_fptr = fopen(argv, "r");                                                            \
        if(g_ext_fptr == NULL)                                                                    \
        {                                                                                         \
            printf("\n Invalid input file argument \n");                                          \
            return;                                                                               \
        }                                                                                         \
    }

#define FLA_TEST_CHECK_EINFO(residual, info, einfo) \
    {                                               \
        if(einfo == 0 && info < 0)                  \
            *residual = DBL_MIN;                    \
        else if(info != einfo)                      \
            *residual = DBL_MAX;                    \
    }

#define FLA_EXTREME_CASE_TEST                                                                  \
    (params->imatrix_char == 'A' || params->imatrix_char == 'F' || params->imatrix_char == 'N' \
     || params->imatrix_char == 'I')

#define FLA_OVERFLOW_UNDERFLOW_TEST (params->imatrix_char == 'O' || params->imatrix_char == 'U')

/* Macro to check if a LAPACK API have different names for its
   (precision)variants and modify API display string
   according to datatype/precision

   Ex: Eigen API SYEVX - FLOAT          - SSYEVX
                         DOUBLE         - DSYEVX
                         COMPLEX        - CHEEVX
                         DOUBLE_COMPLEX - ZHEEVX                 */
#define FLA_MAP_API_NAME(datatype, func_str)                  \
    if((datatype == COMPLEX) || (datatype == DOUBLE_COMPLEX)) \
    {                                                         \
        if(strcmp(func_str, "SYEVX") == 0)                    \
            func_str = "HEEVX";                               \
        else if(strcmp(func_str, "SYEV") == 0)                \
            func_str = "HEEV";                                \
        else if(strcmp(func_str, "SYEVD") == 0)               \
            func_str = "HEEVD";                               \
        else if(strcmp(func_str, "ORGQR") == 0)               \
            func_str = "UNGQR";                               \
        else if(strcmp(func_str, "ORG2R") == 0)               \
            func_str = "UNG2R";                               \
        else if(strcmp(func_str, "SYGVD") == 0)               \
            func_str = "HEGVD";                               \
    }

/* Macro to skip complex and double complex tests of not supported APIs */
#define FLA_SKIP_TEST(datatype_char, func_str)                                           \
    ((((datatype_char == 'c' || datatype_char == 'z') && strcmp(func_str, "STEVD") == 0) \
      || ((datatype_char == 's' || datatype_char == 'd')                                 \
          && (strcmp(func_str, "HETRF") == 0 || strcmp(func_str, "HETRF_ROOK") == 0    \
        || strcmp(func_str, "HETRI_ROOK") == 0)))   \
         ? TRUE                                                                          \
         : FALSE)

/* Assign leading dimension value based on matrix layout and cmdline/config inputs */
#define SELECT_LDA(g_ext_fptr, config_data, layout, n, rm_lda, lda_t) \
    if((g_ext_fptr == NULL) && (layout == LAPACK_ROW_MAJOR))          \
    {                                                                 \
        if(config_data)                                               \
        {                                                             \
            lda_t = fla_max(1, n);                                    \
        }                                                             \
        else                                                          \
        {                                                             \
            lda_t = rm_lda;                                           \
        }                                                             \
    }

typedef struct Lin_solver_paramlist_t
{
    // below params are used only by Lin solver driver APIs.
    doublereal rcond; // used to determine the effective rank of matrix
    float solver_threshold; // threshold to verify PASS/FAIL criteria
    integer n_err_bnds_porfsx;
    integer nparams_porfsx;
    integer kl_gbcon; // number of subdiagonals
    integer ku_gbcon; // number of superdiagonals
    integer ldab_gbcon; // leading dimension of the array ab

    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer num_tests;
    integer num_repeats;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    integer matrix_layout;
    integer nrhs; // number of rhight hand sides
    integer ncolm; // number of columns to factor
    integer lda; //  leading dimension of the array a
    integer ldb; //  leading dimension of the array b
    integer ldq; //  leading dimension of the array q
    integer ldz; //  leading dimension of the array z
    integer ldab; //  leading dimension of the array ab. For GBTRF, GBTRS, LDAB >= 2*KL+KU+1
    integer kl; // number of subdiagonals
    integer ku; // number of superdiagonals
    integer kd; // number of super or sub diagonals
    integer ilo;
    integer ihi;
    char data_types_char[MAX_NUM_DATATYPES];
    char diag; // flag to indicate unit diagonal
    char Uplo;
    char transr; // Must be 'N' or 'T' or 'C'.
    char compq_gghrd;
    char compz_gghrd;
    // below params are used only by Lin solver driver APIs.
    char fact; // Must be 'F', 'N', or 'E'.
    char equed; // Must be 'N', 'R'. 'C', 'B'
    char symm; // if symmetric 'S' or Hermitian 'H'
    char equed_porfsx; // Must be 'N', 'Y'.
    char norm_gbcon; // norm param for gbcon API
} Lin_solver_paramlist;

/* struct to hold eigen parameters */
typedef struct EIG_paramlist_t
{
    real VL;
    real VU;
    real abstol;
    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer num_tests;
    integer num_repeats;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    integer matrix_layout;
    integer threshold_value; // threshold value for EIG
    integer nrhs; // number of rhight hand sides
    integer lda; //  leading dimension of the array a
    integer ldb; //  leading dimension of the array b
    integer ldz;
    integer ldq;
    integer nb; //  leading dimension of the array ab
    integer ldt; // number of subdiagonals
    integer k;
    integer isgn;
    integer kb;
    integer itype;
    integer tsize;
    integer ilo;
    integer ihi;
    integer IL;
    integer IU;
    char data_types_char[MAX_NUM_DATATYPES];
    char trans; // Must be 'N' or 'T' or 'C'.
    char uplo; // Must be 'U' or 'L'
    char job; // Must be 'N', 'P', 'S' or 'B'
    char jobz; // Must be 'N' or 'V'
    char vect; // Vector must be 'Q' or  'P'
    char compq_hgeqz;
    char compz_hgeqz;
    char compz;
    char compz_hseqr;
    char vect_rd;
    char side;
    char job_seqr; // Must be 'E', 'S'
    char eigsrc; // Must be 'Q' or  'N'.
    char initv; // Must be 'N' or 'U'.
    char norm;
    char diag;
    char storev;
    char range_x; // range must be 'A', 'V' or 'I'
} EIG_paramlist;

/* struct to hold eigen parameters */
typedef struct EIG_Non_symmetric_paramlist_t
{
    /* Thresholds for the APIs  */
    float gghrd_threshold; // threshold for the gghrd API
    float ggbal_threshold; // threshold for the ggbal API
    float GenNonSymEigProblem_threshold; // threshold for the ggbal API

    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer lda; // The leading dimension of A. LDA >= fla_max(1,N).
    integer ldb; // The leading dimension of B.  LDB >= fla_max(1,N).
    integer ldvl; // The leading dimension of the matrix VL. LDVL >= 1, and
                  // if JOBVL = 'V', LDVL >= N.
    integer ldvr; // The leading dimension of the matrix VR. LDVR >= 1, and
                  // if JOBVR = 'V', LDVR >= N.
    integer num_repeats;
    integer num_tests;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    integer matrix_layout;
    integer wantz; // Must be 1 or 0
    integer wantq; // Must be 1 or 0
    integer tgsen_ijob; // Must be between 0 to 5
    integer ilo;
    integer ihi;

    /* used params for trsyl API  */
    integer isgn; //  +1 or -1

    char data_types_char[MAX_NUM_DATATYPES];
    char howmny; // Must be 'A' or 'B' or 'S'.
    char initv; // Must be 'N' or 'U'.
    char job_seqr; // Must be 'E', 'S'
    char eigsrc; // Must be 'Q' or  'N'.
    char side; // Must be 'R' or 'L' or 'B'.
    char job; //  Must be 'E' or 'V' or 'B'.
    char howmny_trsna; // Must be 'A' or 'S'.

    /* used params for trsen API  */
    char job_trsen; // must bie 'N' or 'E' or 'V' or 'B'.
    char compq; // Must be 'V' or 'N' .

    /* used params for trsyl API  */
    char trana_real; //  Must be 'N' or 'T' or 'C'.
    char trana_complex; //  Must be 'N' or 'T' or 'C'.
    char tranb_real; //  Must be 'N' or 'T' or 'C'.
    char tranb_complex; //  Must be 'N' or 'T' or 'C'.

    char compq_hgeqz; // Must be 'I' or 'V' or 'N'
    char compz_hgeqz; // Must be 'I' or 'V' or 'N'

    char side_tgevc; // Must be 'R', 'L', or 'B'.
    char jobvsl; // Must be 'N', or 'V'.
    char jobvsr; // Must be 'N', or 'V'.
    char sort; // Must be 'N', or 'S'.
    char sense_ggesx; // must be 'N' or 'E' or 'V' or 'B'.
    char balance_ggevx; // must be 'N' or 'P' or 'S' or 'B'.
    char sense_ggevx; // must be 'N' or 'E' or 'V' or 'B'.
    char sort_gees; // Must be 'N', or 'S'.
    char unmhr_trans; // Must be N or C
} EIG_Non_symmetric_paramlist;

/* struct to hold SVD parameters */
typedef struct SVD_paramlist_t
{
    /* Parameters for 'gesvj' API  */
    float ctol_gesvj; // convergence of threshold

    /* Parameters for 'gesvdx' API  */
    float vl, vu; //  the lower and upper bounds of the interval.

    /* Thresholds for the APIs  */
    float svd_threshold; // threshold for the gghrd API

    float tola;
    float tolb;
    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer lda; // Leading dimension of Array A. LDA >= fla_max(1, n)
    integer ldu; // The leading dimension of the array U.  LDU >= 1; if JOBZ = 'S' or 'A' or JOBZ =
                 // 'O' and M < N, LDU >= M.
    integer ldvt; // The leading dimension of the array VT.  LDVT >= 1; if JOBZ = 'A' or JOBZ = 'O'
                  // and M >= N, LDVT >= N;if JOBZ = 'S', LDVT >= min(M,N).
    integer num_repeats;
    integer num_tests;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    integer matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
    integer m; // The number of rows of the matrix A
    integer p; // The number of rows of the matrix B
    integer n; // The number of columns of the matrices A and B
    integer m_gejsv; // The number of rows of the matrix A
    integer n_gejsv; // The number of rows of the matrix B

    /* Parameters for 'gesvj' API  */
    integer m_gesvj; // The number of rows of the matrix A
    integer n_gesvj; // The number of rows of the matrix B
    integer mv_gesvj;

    /* Parameters for 'gesvdx' API  */
    integer il, iu; // the indices of the smallest and largest singular values.

    char data_types_char[MAX_NUM_DATATYPES];
    char jobu; // Must be 'U' or 'N'.
    char jobv; // Must be 'V' or 'N'.
    char jobq; // Must be 'Q' or 'N'.
    char jobu_gesvd; // Must be 'A', 'S', 'O', or 'N'.
    char jobvt_gesvd; // Must be 'A', 'S', 'O', or 'N'.

    char joba_gejsv; //  Must be 'C', 'E', 'F', 'G', 'A', or 'R'.
    char jobu_gejsv; // Must be 'U', 'F', 'W', or 'N'.
    char jobv_gejsv; // Must be 'V', 'J', 'W', or 'N'.
    char jobr_gejsv; // Must be 'N' or 'R'.
    char jobt_gejsv; // Must be 'T' or 'N'.
    char jobp_gejsv; //  Must be 'P' or 'N'.

    /* Parameters for 'gesvj' API  */
    char joba_gesvj; //  Must be 'L', 'U' or 'G'.
    char jobu_gesvj; // Must be 'U', 'C' or 'N'.
    char jobv_gesvj; // Must be 'V', 'A' or 'N'.

    /* Parameters for 'gesvdx' API  */
    char jobu_gesvdx; //  Must be 'V', or 'N'.
    char jobvt_gesvdx; // Must be 'V', or 'N'.
    char range_gesvdx; // Must be 'A', 'V', 'I'.

    /* Parameters for 'gesvdq' API  */
    char joba_gesvdq; //  Must be 'A', 'H', 'M' , 'E'
    char jobu_gesvdq; // Must be 'A', 'S', 'R' , 'N'
    char jobv_gesvdq; // Must be 'A', 'V', 'R' , 'N'.
} SVD_paramlist;

/* struct to hold AUX parameters */
typedef struct AUX_paramlist_t
{
    /* Parameter for 'larfg' API */
    double alpha_real;
    double alpha_imag; // The alpha values for larfg

    /* Thresholds for the APIs  */
    float aux_threshold; // threshold for the aux API

    /* Parameter for 'larfg' API */
    integer incx_larfg; // The increment between successive values of X in larfg(incx > 0)
    integer incv; // The increment between elements of V for larf
    integer ldc; // The leading dimension of the array C for larf
    integer num_repeats;
    integer num_tests;
    integer num_data_types;
    integer data_types[MAX_NUM_DATATYPES];
    integer matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR

    integer num_ranges; // number of ranges to run
    integer m_range_start;
    integer m_range_end;
    integer m_range_step_size;
    integer n_range_start;
    integer n_range_end;
    integer n_range_step_size;
    integer lda; // Leading dimension of Array A. LDA >= fla_max(1, n)
    integer incx; // The increment between successive values of CX
    integer incy; // The increment between successive values of CY
    /* Parameter for 'larfg' API */
    char side; // The side (either L or R) for larf. L means left and R means right
    char data_types_char[MAX_NUM_DATATYPES];
    /* Norm types to be tested for lange APIs
        M: max absolute value
        1/O: 1-norm
        I: infinity norm
        F/E: Frobenius norm  */
    char norm_types_str[MAX_NUM_NORMTYPES];
    /* The display name of the API being tested */
    char *front_str;
} AUX_paramlist;

typedef enum
{
    LAPACK_TEST,
    LAPACKE_ROW_TEST,
    LAPACKE_COLUMN_TEST,
    LAPACK_CPP_TEST
} test_interface;

typedef struct
{
    char *datatype_char;
    integer *datatype;
    integer n_repeats;
    integer n_datatypes;
    integer p_first;
    integer p_max;
    integer p_inc;
    integer p_nfact;
    integer test_lapacke_interface;
    test_interface interfacetype;
    int matrix_major;
    char imatrix_char;

    struct SVD_paramlist_t svd_paramslist[NUM_SUB_TESTS];
    struct EIG_Non_symmetric_paramlist_t eig_non_sym_paramslist[NUM_SUB_TESTS];
    struct EIG_paramlist_t eig_sym_paramslist[NUM_SUB_TESTS];
    struct Lin_solver_paramlist_t lin_solver_paramslist[NUM_SUB_TESTS];
    struct AUX_paramlist_t aux_paramslist[NUM_SUB_TESTS];

} test_params_t;

typedef struct
{
    double failwarn_s;
    double warnpass_s;
    double failwarn_d;
    double warnpass_d;
    double failwarn_c;
    double warnpass_c;
    double failwarn_z;
    double warnpass_z;
} test_thresh_t;

typedef struct
{
    integer type; /* 0 for LIN, 1 for EIG, 2 for SVD, 3 for AUX */
    char *ops;
    void (*fp)(integer argc, char **argv, test_params_t *);
} OPERATIONS;

// Prototypes.
char *fla_test_get_string_for_result(double residual, integer datatype, double thresh);
void fla_test_init_strings(void);
void fla_test_execute_cli_api(integer argc, char **argv, test_params_t *params);
void fla_test_output_op_struct(char *op_str, integer op);
void fla_test_output_info(char *message, ...);
void fla_test_output_error(char *message, ...);
void fla_test_parse_message(FILE *output_stream, char *message, va_list args);
void fla_test_read_next_line(char *buffer, FILE *input_stream);
integer fla_test_check_run_only(FILE *input_stream, integer *op, char *buffer);
integer fla_test_read_tests_for_op(FILE *input_stream, integer *op, char *buffer);

/*Read Linear API parameters from config file */
void fla_test_read_linear_param(const char *input_filename, test_params_t *params);

/*Function to read EIG parametes from config file */
void fla_test_read_sym_eig_params(const char *input_filename, test_params_t *params);

/*Function to read EIG non-symmetric parametes from config file */
void fla_test_read_non_sym_eig_params(const char *input_filename, test_params_t *params);

/*Function to read SVD parametes from config file */
void fla_test_read_svd_params(const char *input_filename, test_params_t *params);

/*Function to read AUX parametes from config file */
void fla_test_read_aux_params(const char *input_filename, test_params_t *params);

void fla_test_lapack_suite(char *input_filename, test_params_t *params);

void fla_test_op_driver(char *func_str, integer square_inp, test_params_t *params, integer api_type,
                        void (*f_exp)(test_params_t *, // params
                                      integer, // datatype
                                      integer, // p_cur
                                      integer, // q_cur
                                      integer, // pci (param combo counter)
                                      integer, // n_repeats
                                      integer, // einfo
                                      double *, // perf
                                      double *, // time
                                      double *)); // residual
void fla_test_print_status(char *func_str, char datatype_char, integer sqr_inp, integer p_cur,
                           integer q_cur, double residual, double thresh, double time, double perf);
void fla_test_print_summary();
void fla_test_build_function_string(char *func_base_str, char *impl_var_str, char *func_str);
void fill_string_with_n_spaces(char *str, integer n_spaces);
double fla_test_clock(void);
void fla_test_get_time_unit(char *scale, double *time);
integer fla_test_get_group_id(char *buffer);
void fla_check_lapacke_interface(integer *arg_count, char **argv, test_params_t *params);
bool fla_check_interface(integer *arg_count, char **argv, test_params_t *params);
#endif // TEST_LAPACK_H
