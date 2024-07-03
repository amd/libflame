/******************************************************************************
 * Copyright (C) 2023-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"

#define SYTRF_VU 10.0 // Maximum eigen value for condition number.
#define SYTRF_VL 0.01 // Minimum eigen value for condition number.

void invoke_sytrf(integer datatype, char *uplo, integer *n, void *a, integer *lda, integer *ipiv,
                  void *work, integer *lwork, integer *info);
void fla_test_sytrf_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_sytrf_run(integer datatype, integer n, void *A, char uplo, integer lda, integer *ipiv,
                       void *work, integer lwork, integer n_repeats, double *time_min_,
                       integer *info);

void fla_test_sytrf(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LDL**T / U**TDU factorization";
    char *front_str = "SYTRF";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        config_data = 1;
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_sytrf_experiment);
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
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

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

                /* Check for duplicate datatype */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_sytrf_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);

                /* Print the result */
                fla_test_print_status(front_str, stype, SQUARE_INPUT, N, N, residual,
                                      params->lin_solver_paramslist[0].solver_threshold, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }
    /* Print error message */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for sytrf\n");
        printf("./<EXE> sytrf <precisions - sdcz> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
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

void fla_test_sytrf_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer n, lda, lwork = -1, info = 0;
    void *A = NULL, *A_test = NULL, *ipiv = NULL, *work = NULL, *L = NULL;
    char uplo;

    /* Determine the dimensions */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    uplo = params->lin_solver_paramslist[pci].Uplo;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
    }

    /* Create the matrices for the current operation */
    create_matrix(datatype, &A, lda, n);
    create_vector(INTEGER, &ipiv, n);
    create_matrix(datatype, &A_test, lda, n);

    /* Initialize the test matrices */
    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST && !FLA_OVERFLOW_UNDERFLOW_TEST))
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        if(params->imatrix_char != NULL)
        {
            form_symmetric_matrix(datatype, n, A, lda, "S");
        }
    }
    else
    {
        /* Generating input matrix with condition number < 1000 */
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, 'V', n, A, lda, L, SYTRF_VL, SYTRF_VU);
        form_symmetric_matrix(datatype, n, A, lda, "S");
        free_vector(L);
        /* Oveflow or underflow test initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_overflow_underflow_sytrf(datatype, n, A, lda, params->imatrix_char);
        }
    }

    /* Save the original matrix */
    copy_matrix(datatype, "full", lda, n, A, lda, A_test, lda);

    /* call to API */
    prepare_sytrf_run(datatype, n, A_test, uplo, lda, ipiv, work, lwork, n_repeats, t, &info);

    /* Performance computation */
    *perf = (double)(n * n * n) * (1.0 / 3.0) / *t / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Output validataion */
    if((!FLA_EXTREME_CASE_TEST) && info >= 0)
    {
        validate_sytrf(&uplo, n, lda, A_test, datatype, ipiv, residual, &info, A);
        info = 0;
    }
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up buffers */
    free_vector(ipiv);
    free_matrix(A);
    free_matrix(A_test);
}

void prepare_sytrf_run(integer datatype, integer n, void *A, char uplo, integer lda, integer *ipiv,
                       void *work, integer lwork, integer n_repeats, double *time_min_,
                       integer *info)
{
    integer i;
    void *A_save = NULL;
    double time_min = 1e9, exe_time;

    create_matrix(datatype, &A_save, lda, n);

    if(g_lwork == -1)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* Getting lwork from api by passing lwork = -1 */
        invoke_sytrf(datatype, &uplo, &n, NULL, &lda, ipiv, work, &lwork, info);
        if(*info == 0)
        {
            lwork = get_work_value(datatype, work);
        }
        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }

    *info = 0;

    for(i = 0; i < n_repeats && *info == 0; i++)
    {
        /* Copy original input */
        copy_matrix(datatype, "full", lda, n, A, lda, A_save, lda);

        /* Create work buffer */
        create_vector(datatype, &work, lwork);

        exe_time = fla_test_clock();

        /*  call to API */
        invoke_sytrf(datatype, &uplo, &n, A_save, &lda, ipiv, work, &lwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        free_vector(work);
    }

    *time_min_ = time_min;

    /* Save the output to vector A */
    copy_matrix(datatype, "full", lda, n, A_save, lda, A, lda);
    free_matrix(A_save);
}

/*
SYTRF_API calls LAPACK interface for factorization
of a real symmetric matrix A using the Bunch-Kaufman
diagonal pivoting method(A = LDL*T or A = U*TDU)
*/
void invoke_sytrf(integer datatype, char *uplo, integer *n, void *a, integer *lda, integer *ipiv,
                  void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssytrf(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dsytrf(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_csytrf(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zsytrf(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }
    }
}