/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

integer row_major_potrs_lda;
integer row_major_potrs_ldb;

/* Local prototypes */
void fla_test_potrs_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *time_min, double *residual);
void prepare_potrs_run(char *uplo, integer m, integer nrhs, void *A, integer lda, integer datatype,
                       void *b, integer ldb, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int matrix_layout);
void invoke_potrs(char *uplo, integer datatype, integer *m, void *A, integer *lda, integer *nrhs,
                  void *b, integer *ldb, integer *info);
double prepare_lapacke_potrs_run(integer datatype, int matrix_layout, char *uplo, integer m,
                                 integer nrhs, void *A, integer lda, void *b, integer ldb,
                                 integer *info);
integer invoke_lapacke_potrs(integer datatype, int matrix_layout, char uplo, integer n,
                             integer nrhs, const void *a, integer lda, void *b, integer ldb);

void fla_test_potrs(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Cholesky factorization";
    char *front_str = "POTRS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potrs_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }
    if(argc >= 9 && argc <= 10)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_potrs_lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            row_major_potrs_ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
            params->lin_solver_paramslist[0].ldb = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        }

        n_repeats = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->lin_solver_paramslist[0].solver_threshold = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalide dataype */
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
                fla_test_potrs_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, SQUARE_INPUT, N, N, residual,
                                      params->lin_solver_paramslist[0].solver_threshold, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for potrs\n");
        printf("./<EXE> potrs <precisions - sdcz> <UPLO> <N> <NRHS> <LDA> <LDB> <repeats>\n");
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

void fla_test_potrs_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer n, info = 0, nrhs, lda, ldb, vinfo = 0;
    void *A = NULL, *A_test = NULL, *scal = NULL;
    void *B = NULL, *X = NULL;
    void *B_test = NULL;
    double time_min = 1e9;
    char uplo = params->lin_solver_paramslist[pci].Uplo;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    nrhs = params->lin_solver_paramslist[pci].nrhs;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    /* Get input matrix dimensions. */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
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

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &B, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &X, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &B_test, ldb);

    /* Initialize input symmetric positive definite matrix A */
    reset_matrix(datatype, n, n, A, lda);
    if((!FLA_EXTREME_CASE_TEST) && (g_ext_fptr == NULL))
    {
        rand_spd_matrix(datatype, &uplo, A, n, lda);
        /* Initialize input matrix with random numbers */
        rand_matrix(datatype, B, n, nrhs, ldb);

        /* Scaling matrix with values around overflow, underflow for POTRS */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_realtype_vector(get_datatype(datatype), &scal, 1);
            scale_matrix_underflow_overflow_potrs(datatype, n, A, lda, params->imatrix_char, scal);
        }
    }
    else
    {
        /* Initialize input matrix with custom data */
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        init_matrix(datatype, B, n, nrhs, ldb, g_ext_fptr, params->imatrix_char);
    }
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    /* cholesky factorisation of A as input to potrs */
    invoke_potrf(&uplo, datatype, &n, A, &lda, &info);

    copy_matrix(datatype, "full", n, nrhs, B, ldb, B_test, ldb);

    /* Invoke potrs API to find x using Ax-b */
    prepare_potrs_run(&uplo, n, nrhs, A, lda, datatype, B_test, ldb, n_repeats, &time_min, &info,
                      interfacetype, layout);
    copy_matrix(datatype, "full", n, nrhs, B_test, ldb, X, n);
    /* execution time */
    *t = time_min;
    /* Compute the performance of the best experiment repeat. */
    /* (2.0)m^2 flops for Ax=b computation. */
    *perf
        = (double)((2.0 * n * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Validate potrs call by computing Ax-b */
    if((info == 0) && (!FLA_EXTREME_CASE_TEST))
        validate_potrs(n, nrhs, A_test, lda, X, B, ldb, datatype, residual, &vinfo,
                       params->imatrix_char);
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, nrhs, B_test, ldb, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B_test);
    free_matrix(B);
    free_matrix(X);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal);
    }
}

void prepare_potrs_run(char *uplo, integer n, integer nrhs, void *A, integer lda, integer datatype,
                       void *B, integer ldb, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int layout)
{
    void *A_save = NULL, *B_test = NULL;
    double time_min = 1e9, exe_time;
    integer i;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nrhs, &B_test, ldb);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
        for each iteration */
        copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
        copy_matrix(datatype, "full", n, nrhs, B, ldb, B_test, ldb);
        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_potrs_run(datatype, layout, uplo, n, nrhs, A_save, lda,
                                                 B_test, ldb, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP potrs API */
            invoke_cpp_potrs(uplo, datatype, &n, A_save, &lda, &nrhs, B_test, &ldb, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK potrs API */
            invoke_potrs(uplo, datatype, &n, A_save, &lda, &nrhs, B_test, &ldb, info);

            exe_time = fla_test_clock() - exe_time;
        }
        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }
    copy_matrix(datatype, "full", n, nrhs, B_test, ldb, B, ldb);
    *time_min_ = time_min;
    free_matrix(A_save);
    free_vector(B_test);
}

double prepare_lapacke_potrs_run(integer datatype, int layout, char *uplo, integer n, integer nrhs,
                                 void *A_save, integer lda, void *B_test, integer ldb,
                                 integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldb_t = ldb;
    void *A_t = NULL, *B_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_potrs_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, nrhs, row_major_potrs_ldb, ldb_t);

    A_t = A_save;
    B_t = B_test;
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, fla_max(n, lda_t));
        create_matrix(datatype, layout, n, nrhs, &B_t, fla_max(nrhs, ldb_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A_save, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, nrhs, B_test, ldb, B_t, ldb_t);
    }
    exe_time = fla_test_clock();

    *info = invoke_lapacke_potrs(datatype, layout, *uplo, n, nrhs, A_t, lda_t, B_t, ldb_t);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A_save, lda);
        convert_matrix_layout(layout, datatype, n, nrhs, B_t, ldb_t, B_test, ldb);
        /* free temporary buffers */
        free_matrix(A_t);
        free_matrix(B_t);
    }
    return exe_time;
}

void invoke_potrs(char *uplo, integer datatype, integer *n, void *A, integer *lda, integer *nrhs,
                  void *B, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_spotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dpotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cpotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zpotrs(uplo, n, nrhs, A, lda, B, ldb, info);
            break;
        }
    }
}

integer invoke_lapacke_potrs(integer datatype, int layout, char uplo, integer n, integer nrhs,
                             const void *A, integer lda, void *B, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_spotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dpotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cpotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zpotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
    }
    return info;
}
