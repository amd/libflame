/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_overflow_underflow.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_trtrs_lda;
integer row_major_trtrs_ldb;

/* Local prototypes */
void fla_test_trtrs_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_trtrs_run(char *uplo, char *trans, char *diag, integer m, integer nrhs, void *A,
                       integer lda, integer datatype, void *b, integer ldb, integer *info,
                       integer interfacetype, int matrix_layout, test_params_t *params);
void invoke_trtrs(char *uplo, char *trans, char *diag, integer datatype, integer *m, void *A,
                  integer *lda, integer *nrhs, void *b, integer *ldb, integer *info);
double prepare_lapacke_trtrs_run(integer datatype, int layout, char *uplo, char *trans, char *diag,
                                 integer m, integer nrhs, void *A, integer lda, void *b,
                                 integer ldb, integer *info);

void fla_test_trtrs(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Triangular solve";
    char *front_str = "TRTRS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_trtrs_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        params->lin_solver_paramslist[0].transr = argv[4][0];
        params->lin_solver_paramslist[0].diag = argv[5][0];
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_trtrs_lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_trtrs_ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
            params->lin_solver_paramslist[0].ldb = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        }

        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        if(n_repeats > 0)
        {
            params->lin_solver_paramslist[0].solver_threshold = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalid datatype */
                if(datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }

                /* Check for duplicate datatype presence */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_trtrs_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error message if help is required */
    if(tests_not_run == 1)
    {
        fla_test_output_info("\n");
        if(invalid_dtype)
        {
            fla_test_output_info("Invalid data type. Supported datatypes are s, d, c, z\n");
        }
        fla_test_output_info("Usage: trtrs <precisions - sdcz> <uplo> <trans> <diag> <n> <nrhs> "
                             "<lda> <ldb> <repeats>\n");
        fla_test_output_info("\n");
    }
}

void fla_test_trtrs_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, nrhs, lda, ldb, info = 0, matrix_layout;
    char uplo, trans, diag;
    double residual, err_thresh;
    void *A = NULL, *B = NULL, *B_test = NULL;
    void *scal = NULL;

    /* Get input matrix dimensions */
    n = p_cur;
    nrhs = params->lin_solver_paramslist[pci].nrhs;
    uplo = params->lin_solver_paramslist[pci].Uplo;
    trans = params->lin_solver_paramslist[pci].transr;
    diag = params->lin_solver_paramslist[pci].diag;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;
    matrix_layout
        = (params->interfacetype == LAPACKE_ROW_TEST) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
        if(ldb == -1)
        {
            ldb = fla_max(1, n);
        }
    }
    /* Allocate memory for matrices */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &B, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &B_test, ldb);

    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        /* allocate scaling factor for overflow/underflow tests */
        create_vector(datatype, &scal, 1);
    }

    /* Generate input matrix based on the input matrix type */
    /* Generate random matrix A */
    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        get_non_singular_triangular_matrix(&uplo, datatype, n, n, A, lda,
                                           diag == 'U' ? UNIT_DIAG : NON_UNIT_DIAG);
    }

    /* Generate random RHS matrix B */
    init_matrix(datatype, B, n, nrhs, ldb, g_ext_fptr, params->imatrix_char);
    /* Oveflow or underflow test initialization */
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_overflow_underflow_trtrs(datatype, n, nrhs, A, B, lda, ldb, trans,
                                              params->imatrix_char, diag, uplo);
    }

    /* Make a copy of the input matrices for testing */
    copy_matrix(datatype, "full", n, nrhs, B, ldb, B_test, ldb);

    /* Prepare and run trtrs call */
    prepare_trtrs_run(&uplo, &trans, &diag, n, nrhs, A, lda, datatype, B_test, ldb, &info,
                      params->interfacetype, matrix_layout, params);

    /* Compute performance */
    perf = (double)(n * n * nrhs) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;
    /* Validate trtrs call by computing residual */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        /* Simple validation - compute basic residual norm */
        /* For a complete implementation, would need proper validation function */
        validate_trtrs(tst_api, datatype, &uplo, &trans, &diag, n, nrhs, A, lda, B_test, B, ldb,
                       residual, params->imatrix_char, params);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, n, A, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(B_test);
    free_matrix(B);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal);
    }
}

void prepare_trtrs_run(char *uplo, char *trans, char *diag, integer n, integer nrhs, void *A,
                       integer lda, integer datatype, void *B, integer ldb, integer *info,
                       integer interfacetype, int layout, test_params_t *params)
{
    void *A_save = NULL, *B_test = NULL;
    double exe_time;

    /* Make a copy of the input matrix A and B. Same input values will be passed in
       each iteration.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &B_test, ldb);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrix A and B values for each iteration */
        copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
        copy_matrix(datatype, "full", n, nrhs, B, ldb, B_test, ldb);
        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_trtrs_run(datatype, layout, uplo, trans, diag, n, nrhs,
                                                 A_save, lda, B_test, ldb, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP trtrs API */
            invoke_cpp_trtrs(uplo, trans, diag, datatype, &n, A_save, &lda, &nrhs, B_test, &ldb,
                             info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK trtrs API */
            invoke_trtrs(uplo, trans, diag, datatype, &n, A_save, &lda, &nrhs, B_test, &ldb, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /* Copy the output matrices */
    copy_matrix(datatype, "full", n, nrhs, B_test, ldb, B, ldb);
    /* Free up the saved matrices */
    free_matrix(A_save);
    free_matrix(B_test);
}

double prepare_lapacke_trtrs_run(integer datatype, int layout, char *uplo, char *trans, char *diag,
                                 integer n, integer nrhs, void *A, integer lda, void *b,
                                 integer ldb, integer *info)
{
    double exe_time = 0;
    integer lda_t = lda;
    integer ldb_t = ldb;
    void *A_t = NULL;
    void *b_t = NULL;
    A_t = A;
    b_t = b;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_trtrs_lda, lda_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_trtrs_ldb, ldb_t);

    /* In case of row_major matrix layout,
    convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        create_matrix(datatype, layout, n, n, &A_t, lda_t);
        create_matrix(datatype, layout, n, nrhs, &b_t, ldb_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, nrhs, b, ldb, b_t,
                              fla_max(ldb_t, nrhs));
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE sytrf API */
    *info = invoke_lapacke_trtrs(datatype, layout, uplo, trans, diag, n, nrhs, A_t, lda_t, b_t,
                                 ldb_t);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
    to column_major layout */
    if((layout == LAPACK_ROW_MAJOR))
    {
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        convert_matrix_layout(layout, datatype, n, nrhs, b_t, ldb_t, b, ldb);
        free_matrix(A_t);
        free_matrix(b_t);
    }

    return exe_time;
}

void invoke_trtrs(char *uplo, char *trans, char *diag, integer datatype, integer *n, void *A,
                  integer *lda, integer *nrhs, void *b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_strtrs(uplo, trans, diag, n, nrhs, A, lda, b, ldb, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dtrtrs(uplo, trans, diag, n, nrhs, A, lda, b, ldb, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_ctrtrs(uplo, trans, diag, n, nrhs, A, lda, b, ldb, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_ztrtrs(uplo, trans, diag, n, nrhs, A, lda, b, ldb, info);
            break;
        }
    }
}