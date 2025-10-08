/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

extern double perf;
extern double time_min;

/* Local prototypes.*/
void fla_test_potf2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_potf2_run(char *uplo, integer m, void *A, integer lda, integer datatype, integer *info,
                       test_params_t *params);
void invoke_potf2(char *uplo, integer datatype, integer *m, void *a, integer *lda, integer *info);

void fla_test_potf2(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Cholesky factorization (unblocked algorithm)";
    char *front_str = "POTF2";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potf2_experiment);
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
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_potf2_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for potf2\n");
        printf("./<EXE> potf2 <precisions - sdcz> <Uplo> <N> <LDA> <repeats>\n");
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

void fla_test_potf2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, lda;
    integer info = 0;
    void *A = NULL, *A_test = NULL;
    char uplo = params->lin_solver_paramslist[pci].Uplo;
    double residual, err_thresh;
    void *filename = NULL;

    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* Get input matrix dimensions */
    m = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A, lda);

    /* Skip input generation for BRT verification runs */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr != NULL)
        {
            /* Initialize input matrix with custom data */
            init_matrix(datatype, A, m, m, lda, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            rand_spd_matrix(datatype, &uplo, A, m, lda);
        }
    }

    /* BRT input processing */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, m, m, A, lda, "cdd", uplo, m, lda)

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_test, lda);
    copy_matrix(datatype, "full", m, m, A, lda, A_test, lda);

    prepare_potf2_run(&uplo, m, A_test, lda, datatype, &info, params);

    /* Compute the performance of the best experiment repeat */
    /* (1/3)m^3 for real and (4/3)m^3 for scomplex*/
    perf = (double)(1.0 / 3.0 * m * m * m) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* BRT validation */
    IF_FLA_BRT_VALIDATION(
        m, m, store_outputs_base(filename, params, 1, 0, datatype, m, m, A_test, lda),
        FLA_PRINT_TEST_STATUS(m, m, residual, err_thresh),
        check_reproducibility_base(filename, params, 1, 0, datatype, m, m, A_test, lda))
    else
    {
        FLA_PRINT_TEST_STATUS(m, m, residual, err_thresh);
    }

    /* Free up the buffers */
free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(A);
    free_matrix(A_test);
}

void prepare_potf2_run(char *uplo, integer m, void *A, integer lda, integer datatype, integer *info,
                       test_params_t *params)
{
    void *A_save = NULL;
    double exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each iteration.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_save, lda);
    copy_matrix(datatype, "full", m, m, A, lda, A_save, lda);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrix A value for each iteration */
        copy_matrix(datatype, "full", m, m, A_save, lda, A, lda);

        exe_time = fla_test_clock();
        /* Call LAPACK potf2 API */
        invoke_potf2(uplo, datatype, &m, A, &lda, info);
        exe_time = fla_test_clock() - exe_time;

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    free_matrix(A_save);
}

void invoke_potf2(char *uplo, integer datatype, integer *m, void *a, integer *lda, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_spotf2(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dpotf2(uplo, m, a, lda, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cpotf2(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zpotf2(uplo, m, a, lda, info);
            break;
        }
    }
}