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
integer row_major_potri_lda;

/* Local prototypes.*/
void fla_test_potri_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_potri_run(char *uplo, integer n, void *A, integer lda, integer datatype, integer *info,
                       integer interfacetype, int matrix_layout, test_params_t *params);
void invoke_potri(char *uplo, integer datatype, integer *n, void *a, integer *lda, integer *info);
double prepare_lapacke_potri_run(integer datatype, int matrix_layout, char *uplo, integer n,
                                 void *A, integer lda, integer *info);

#define VALIDATE_POTRI                                                                            \
    form_symmetric_matrix(datatype, n, A, lda, "C", uplo);                                        \
    form_symmetric_matrix(datatype, n, A_test, lda, "C", uplo);                                   \
    validate_getri(tst_api, n, n, A, A_test, lda, NULL, datatype, residual, params->imatrix_char, \
                   params);

void fla_test_potri(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Inverse of Cholesky factorization";
    char *front_str = "POTRI";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potri_experiment);
        tests_not_run = 0;
    }
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if(argc >= 7 && argc <= 8)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_potri_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_potri_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for potri\n");
        printf("./<EXE> potri <precisions - sdcz> <Uplo> <N> <LDA> <repeats>\n");
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
    return;
}

void fla_test_potri_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, lda;
    integer info = 0;
    void *A = NULL, *A_test = NULL;
    char uplo = params->lin_solver_paramslist[pci].Uplo;
    double residual, err_thresh;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* Get input matrix dimensions */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;

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

    /* Skip input generation during BRT verification runs since FLA_BRT_PROCESS_ macros load inputs
     * from files */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST && !FLA_OVERFLOW_UNDERFLOW_TEST))
        {
            /* Initialize input matrix with custom data */
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
            if((params->imatrix_char != '\0'))
            {
                form_symmetric_matrix(datatype, n, A, lda, "C", uplo);
            }
        }
        else
        {
            rand_spd_matrix(datatype, &uplo, A, n, lda);
            /* Oveflow or underflow test initialization */
            if(FLA_OVERFLOW_UNDERFLOW_TEST)
            {
                scale_matrix_overflow_underflow_potri(datatype, n, A, lda, params->imatrix_char);
            }
        }
    }

    /* BRT input processing: store input matrix during ground truth runs, load during verification
     * runs */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, n, n, A, lda, "cdd", uplo, n, lda)

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    invoke_potrf(&uplo, datatype, &n, A_test, &lda, &info);

    prepare_potri_run(&uplo, n, A_test, lda, datatype, &info, interfacetype, layout, params);

    /* Compute the performance of the best experiment repeat */
    /* (1/3)n^3 for real and (4/3)n^3 for scomplex*/
    perf = (double)(1.0 / 3.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    IF_FLA_BRT_VALIDATION(
        n, n, store_outputs_base(filename, params, 1, 0, datatype, n, n, A_test, lda),
        VALIDATE_POTRI,
        check_reproducibility_base(filename, params, 1, 0, datatype, n, n, A_test, lda))
    else if(FLA_RANDOM_INIT_MODE)
    {
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }
    else if(!FLA_EXTREME_CASE_TEST)
    {
        /* Form full matrices before calling validate code of GETRI */
        form_symmetric_matrix(datatype, n, A, lda, "C", uplo);
        form_symmetric_matrix(datatype, n, A_test, lda, "C", uplo);
        validate_getri(tst_api, n, n, A, A_test, lda, NULL, datatype, residual,
                       params->imatrix_char, params);
    }
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(A);
    free_matrix(A_test);
}

void prepare_potri_run(char *uplo, integer n, void *A, integer lda, integer datatype, integer *info,
                       integer interfacetype, int layout, test_params_t *params)
{
    void *A_save = NULL;
    double exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrix A value and allocate memory to output buffers
        for each iteration */
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_potri_run(datatype, layout, uplo, n, A, lda, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP potri API */
            invoke_cpp_potri(uplo, datatype, &n, A, &lda, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK potri API */
            invoke_potri(uplo, datatype, &n, A, &lda, info);

            exe_time = fla_test_clock() - exe_time;
        }
        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    free_matrix(A_save);
}

double prepare_lapacke_potri_run(integer datatype, int layout, char *uplo, integer n, void *A,
                                 integer lda, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_potri_lda, lda_t);

    A_t = A;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, fla_max(n, lda_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    *info = invoke_lapacke_potri(datatype, layout, *uplo, n, A_t, lda_t);

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

void invoke_potri(char *uplo, integer datatype, integer *n, void *a, integer *lda, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_spotri(uplo, n, a, lda, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dpotri(uplo, n, a, lda, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cpotri(uplo, n, a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zpotri(uplo, n, a, lda, info);
            break;
        }
    }
}
