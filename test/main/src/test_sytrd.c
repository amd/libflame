/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

integer row_major_sytrd_lda;

void fla_test_sytrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_sytrd_run(integer datatype, char uplo, integer n, void *A, integer lda, void *D,
                       void *E, void *tau, integer *info, integer interfacetype, integer layout,
                       test_params_t *params);
void invoke_sytrd(integer datatype, char *uplo, integer *n, void *A, integer *lda, void *D, void *E,
                  void *tau, void *work, integer *lwork, integer *info);
double prepare_lapacke_sytrd_run(integer datatype, integer layout, char uplo, integer n, void *A,
                                 integer lda, void *D, void *E, void *tau, integer *info);

void fla_test_sytrd(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Reduction to tridiagonal matrix - Eigen routine";
    char *front_str = "SYTRD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_sytrd_experiment);
        tests_not_run = 0;
    }
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }
    if((argc == 8) || (argc == 9))
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].uplo = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_sytrd_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].lda = N;
        }
        else
        {
            params->eig_sym_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }

        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        if(n_repeats > 0)
        {
            params->eig_sym_paramslist[0].threshold_value = CLI_NORM_THRESH;

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
                fla_test_sytrd_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("Invalid arguments for sytrd\n");
        printf("Usage: ./<EXE> sytrd <precisions - sd> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
        printf("Usage: ./<EXE> hetrd <precisions - cz> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
    }
    else if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
}

void fla_test_sytrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, info = 0, realtype, lda;
    void *A = NULL, *A_test = NULL, *D = NULL, *E = NULL, *tau = NULL;
    char uplo;
    double residual, err_thresh;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;

    /* Get input matrix dimensions. */
    n = q_cur;

    /* Initialize parameter needed for sytrd() call. */
    uplo = params->eig_sym_paramslist[pci].uplo;
    lda = params->eig_sym_paramslist[pci].lda;
    err_thresh = params->eig_sym_paramslist[pci].threshold_value;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
    }

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);

    reset_matrix(datatype, lda, n, A, lda);
    reset_matrix(datatype, lda, n, A_test, lda);

    realtype = get_realtype(datatype);
    create_vector(realtype, &D, n);
    create_vector(realtype, &E, n - 1);

    /* This code path is run to generate the matrix to be passed to the API. This is the default
     * input generation logic accessed both when BRT is run in Ground truth mode and for non BRT
     * Test cases. For verification runs the input is loaded from the input generated during Ground
     * truth run */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr || FLA_EXTREME_CASE_TEST)
        {
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
            if(FLA_EXTREME_CASE_TEST)
            {
                /* Get the symmetric/hermitian matrix.*/
                form_symmetric_matrix(datatype, n, A, lda, "C", uplo);
            }
        }
        else
        {
            /* TODO: changes to create A from known inputs Q, T */
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
            form_symmetric_matrix(datatype, n, A, lda, "C", uplo);
            /* Scaling matrix with values around overflow, underflow for SYTRD/HETRD */
            if(FLA_OVERFLOW_UNDERFLOW_TEST)
            {
                scale_matrix_underflow_overflow_sytrd(datatype, n, A, lda, params->imatrix_char);
            }
        }
    }

    /* This macro is used in the BRT test cases for the following purposes:
     *    - In the Ground truth runs (BRT_char => G, F), the input is stored in a file for future
     * reference
     *    - In the verification runs (BRT_char => V, M), the input is loaded from the file and
     * passed as input to the API
     * */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, n, n, A, lda, "cddd", uplo, n, lda, g_lwork)

    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    create_vector(datatype, &tau, n - 1);

    prepare_sytrd_run(datatype, uplo, n, A_test, lda, D, E, tau, &info, interfacetype, layout,
                      params);

    /* Performance computation (4/3)n^3 flops. */
    perf = (double)((4.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;

    /* Output validation. */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    IF_FLA_BRT_VALIDATION(
        n, n,
        store_outputs_base(filename, params, 1, 3, datatype, n, n, A_test, lda,
                           get_realtype(datatype), n, D, get_realtype(datatype), n - 1, E, datatype,
                           n - 1, tau),
        validate_sytrd(tst_api, datatype, uplo, n, A_test, A, lda, D, E, tau, residual, params),
        check_reproducibility_base(filename, params, 1, 3, datatype, n, n, A_test, lda,
                                   get_realtype(datatype), n, D, get_realtype(datatype), n - 1, E,
                                   datatype, n - 1, tau))
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_sytrd(tst_api, datatype, uplo, n, A_test, A, lda, D, E, tau, residual, params);
    }
    /* Check for output matrix & vectors when inputs are extreme values */
    else
    {
        residual = err_thresh;
        if(!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char))
        {
            residual = DBL_MAX;
        }
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

    /* Free up buffers. */
    free_vector(tau);
free_buffers:
    FLA_FREE_FILENAME(filename);
    free_matrix(A);
    free_matrix(A_test);
    free_vector(D);
    free_vector(E);
}

void prepare_sytrd_run(integer datatype, char uplo, integer n, void *A, integer lda, void *D,
                       void *E, void *tau, integer *info, integer interfacetype, integer layout,
                       test_params_t *params)
{
    integer lwork;
    void *A_save = NULL;
    void *work = NULL;
    double exe_time;

    /* Make a copy of the input matrices. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    /* Call to sytrd() API to get work buffers size.
       NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork == -1))
    {
        /* Make a workspace query the first time. This will provide us with
        and ideal workspace size based on internal block size.*/
        create_vector(datatype, &work, 1);
        lwork = g_lwork;
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_sytrd(datatype, &uplo, &n, A, &lda, D, E, tau, work, &lwork, info);
        }
        else
#endif
        {
            invoke_sytrd(datatype, &uplo, &n, A, &lda, D, E, tau, work, &lwork, info);
        }

        /* Get work buffers size. */
        if(*info == 0)
        {
            if(g_lwork == -1)
            {
                lwork = get_work_value(datatype, work);
            }
        }

        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }
    create_vector(datatype, &work, lwork);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrices and allocate memory to output buffers
           for each iteration. */
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        reset_vector(datatype, work, lwork, 1);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time
                = prepare_lapacke_sytrd_run(datatype, layout, uplo, n, A, lda, D, E, tau, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP sytrd API */
            invoke_cpp_sytrd(datatype, &uplo, &n, A, &lda, D, E, tau, work, &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();

            /* Call LAPACK sytrd API */
            invoke_sytrd(datatype, &uplo, &n, A, &lda, D, E, tau, work, &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /* Free up buffers. */
    free_vector(A_save);
    free_vector(work);
}

double prepare_lapacke_sytrd_run(integer datatype, integer layout, char uplo, integer n, void *A,
                                 integer lda, void *D, void *E, void *tau, integer *info)
{
    double exe_time = 0;
    void *A_t = A;
    integer lda_t = lda;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_sytrd_lda, lda_t);

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        create_matrix(datatype, layout, n, n, &A_t, fla_max(n, lda_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE sytrd API */
    *info = invoke_lapacke_sytrd(datatype, layout, uplo, n, A_t, lda_t, D, E, tau);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */
    if((layout == LAPACK_ROW_MAJOR))
    {
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        free_matrix(A_t);
    }

    return exe_time;
}

/* LAPACK SYTRD API invoke function to form Symmetric tridiagonal matrix. */
void invoke_sytrd(integer datatype, char *uplo, integer *n, void *A, integer *lda, void *D, void *E,
                  void *tau, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssytrd(uplo, n, A, lda, D, E, tau, work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsytrd(uplo, n, A, lda, D, E, tau, work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_chetrd(uplo, n, A, lda, D, E, tau, work, lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhetrd(uplo, n, A, lda, D, E, tau, work, lwork, info);
            break;
        }
    }
}
