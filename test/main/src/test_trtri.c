/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_trtri_lda;

/* Local prototypes */
void fla_test_trtri_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_trtri_run(integer n, void *A, char uplo, char diag, integer lda, integer datatype,
                       integer n_repeats, integer *info, integer interfacetype, int layout,
                       test_params_t *params);
void invoke_trtri(integer datatype, char *uplo, char *diag, integer *n, void *a, integer *lda,
                  integer *info);
double prepare_lapacke_trtri_run(integer datatype, int layout, char uplo, char diag, integer n,
                                 void *A, integer lda, integer *info);

void fla_test_trtri(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Triangular matrix inversion";
    char *front_str = "TRTRI";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_trtri_experiment);
        tests_not_run = 0;
    }
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }
    if(argc >= 8 && argc <= 9)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        params->lin_solver_paramslist[0].diag = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_trtri_lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        }

        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        if(n_repeats > 0)
        {
            params->lin_solver_paramslist[0].solver_threshold = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalid dataype */
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
                fla_test_trtri_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for trtri\n");
        printf("./<EXE> trtri <precisions - sdcz> <UPLO> <DIAG> <N> <LDA> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
}

void fla_test_trtri_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, lda;
    integer info = 0;
    void *A = NULL, *A_test = NULL;
    double residual, err_thresh;
    void *filename = NULL;
    char uplo = params->lin_solver_paramslist[pci].Uplo;
    char diag = params->lin_solver_paramslist[pci].diag;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);

    /* This code path is run to generate the matrix to be passed to the API. This is the default
     * input generation logic accessed both when BRT is run in Ground truth mode and for non BRT
     * Test cases. For verification runs the input is loaded from the input generated during Ground
     * truth run */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
        {
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            get_non_singular_triangular_matrix(&uplo, datatype, n, n, A, lda,
                                               (diag == 'U') ? UNIT_DIAG : NON_UNIT_DIAG);
        }

        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_trtri(datatype, diag, n, A, lda, params->imatrix_char);
        }
    }

    /* This macro is used in the BRT test cases for the following purposes:
     *    - In the Ground truth runs (BRT_char => G, F), the input is stored in a file for future
     * reference
     *    - In the verification runs (BRT_char => V, M), the input is loaded from the file and
     * passed as input to the API
     * */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, n, n, A, lda, "ccdd", uplo, diag, n, lda)

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_trtri_run(n, A_test, uplo, diag, lda, datatype, n_repeats, &info, interfacetype, layout,
                      params);

    /* performance computation */
    perf = (double)((1.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    IF_FLA_BRT_VALIDATION(
        n, n, store_outputs_base(filename, params, 1, 0, datatype, n, n, A_test, lda),
        validate_trtri(tst_api, uplo, diag, n, A, A_test, lda, datatype, residual,
                       params->imatrix_char, params),
        check_reproducibility_base(filename, params, 1, 0, datatype, n, n, A_test, lda))
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_trtri(tst_api, uplo, diag, n, A, A_test, lda, datatype, residual,
                       params->imatrix_char, params);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if(!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char))
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
    free_matrix(A_test);
free_buffers:
    FLA_FREE_FILENAME(filename);
    free_matrix(A);
}

void prepare_trtri_run(integer n, void *A, char uplo, char diag, integer lda, integer datatype,
                       integer n_repeats, integer *info, integer interfacetype, int layout,
                       test_params_t *params)
{
    void *A_save = NULL;
    double exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);

    *info = 0;

    FLA_EXEC_LOOP_BEGIN
    {
        /* Copy original input */
        copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time
                = prepare_lapacke_trtri_run(datatype, layout, uplo, diag, n, A_save, lda, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST) /* Call CPP TRTRI API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_trtri(datatype, &uplo, &diag, &n, A_save, &lda, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* call to API */
            invoke_trtri(datatype, &uplo, &diag, &n, A_save, &lda, info);
            exe_time = fla_test_clock() - exe_time;
        }

        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /* Save the output to vector A */
    copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
    free_matrix(A_save);
}

double prepare_lapacke_trtri_run(integer datatype, int layout, char uplo, char diag, integer n,
                                 void *A, integer lda, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_trtri_lda, lda_t);

    A_t = A;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* Call LAPACKE trtri API */
    *info = invoke_lapacke_trtri(datatype, layout, uplo, diag, n, A_t, lda_t);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        /* free temporary buffers */
        free_matrix(A_t);
    }

    return exe_time;
}

void invoke_trtri(integer datatype, char *uplo, char *diag, integer *n, void *a, integer *lda,
                  integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_strtri(uplo, diag, n, (float *)a, lda, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dtrtri(uplo, diag, n, (double *)a, lda, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_ctrtri(uplo, diag, n, (scomplex *)a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_ztrtri(uplo, diag, n, (dcomplex *)a, lda, info);
            break;
        }
    }
}