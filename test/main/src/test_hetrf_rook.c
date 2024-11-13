/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

#define HETRF_ROOK_VU 10.0 // Maximum eigen value for condition number.
#define HETRF_ROOK_VL 0.01 // Minimum eigen value for condition number.

integer row_major_hetrf_rook_lda;

void invoke_hetrf_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                       integer *ipiv, void *work, integer *lwork, integer *info);
void fla_test_hetrf_rook_experiment(test_params_t *params, integer datatype, integer p_cur,
                                    integer q_cur, integer pci, integer n_repeats, integer einfo,
                                    double *perf, double *t, double *residual);
void prepare_hetrf_rook_run(integer datatype, integer n, void *A, char uplo, integer lda,
                            integer *ipiv, void *work, integer lwork, integer n_repeats,
                            double *time_min_, integer *info, integer interfacetype,
                            integer mlayout);
double prepare_lapacke_hetrf_rook_run(integer datatype, integer layout, char uplo, integer n,
                                      void *A, integer lda, void *ipiv, integer *info);
integer invoke_lapacke_hetrf_rook(integer datatype, integer layout, char uplo, integer n, void *a,
                                  integer lda, integer *ipiv);

void fla_test_hetrf_rook(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LDL**H / UDU**H factorization";
    char *front_str = "HETRF_ROOK";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        config_data = 1;
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_hetrf_rook_experiment);
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
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_hetrf_rook_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
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
                if(datatype == FLOAT || datatype == DOUBLE ||  datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }

                /* Check for duplicate datatype */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_hetrf_rook_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
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
        printf("\nIllegal arguments for hetrf_rook\n");
        printf("./<EXE> hetrf_rook <precisions - cz> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'cz'\n\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
}

void fla_test_hetrf_rook_experiment(test_params_t *params, integer datatype, integer p_cur,
                                    integer q_cur, integer pci, integer n_repeats, integer einfo,
                                    double *perf, double *t, double *residual)
{
    integer n, lda, lwork = -1, info = 0;
    void *A = NULL, *A_test = NULL, *ipiv = NULL, *work = NULL, *L = NULL;
    char uplo;
    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;

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
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_vector(INTEGER, &ipiv, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);

    /* Initialize the test matrices */
    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST && !FLA_OVERFLOW_UNDERFLOW_TEST))
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        if(params->imatrix_char != NULL)
        {
            form_symmetric_matrix(datatype, n, A, lda, "C");
        }
    }
    else
    {
        /* Create input matrix (hermitian) */
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, 'V', n, A, lda, L, HETRF_ROOK_VL, HETRF_ROOK_VU,
                                 USE_SIGNED_EIGEN_VALUES);
        form_symmetric_matrix(datatype, n, A, lda, "C");
        free_vector(L);
        /* Oveflow or underflow test initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_overflow_underflow_hetrf_rook(datatype, n, A, lda, params->imatrix_char);
        }
    }

    /* Save the original matrix */
    copy_matrix(datatype, "full", lda, n, A, lda, A_test, lda);

    /* call to API */
    prepare_hetrf_rook_run(datatype, n, A_test, uplo, lda, ipiv, work, lwork, n_repeats, t, &info,
                           interfacetype, layout);

    /* Performance computation */
    *perf = (double)(n * n * n) * (1.0 / 3.0) / *t / FLOPS_PER_UNIT_PERF;
    *perf *= 4.0;

    /* Output validataion */
    if((!FLA_EXTREME_CASE_TEST) && info >= 0)
    {
        validate_hetrf_rook(&uplo, n, lda, A_test, datatype, ipiv, residual, &info, A);
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

void prepare_hetrf_rook_run(integer datatype, integer n, void *A, char uplo, integer lda,
                            integer *ipiv, void *work, integer lwork, integer n_repeats,
                            double *time_min_, integer *info, integer interfacetype,
                            integer layout)
{
    integer i;
    void *A_save = NULL;
    double time_min = 1e9, exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.
     NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST) && (g_lwork == -1))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* Getting lwork from api by passing lwork = -1 */
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            invoke_cpp_hetrf_rook(datatype, &uplo, &n, NULL, &lda, ipiv, work, &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
        else
#endif
        {
            invoke_hetrf_rook(datatype, &uplo, &n, NULL, &lda, ipiv, work, &lwork, info);
        }
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
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_hetrf_rook_run(datatype, layout, uplo, n, A_save, lda, ipiv,
                                                      info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)   /* Call CPP hetrf_rook API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_hetrf_rook(datatype, &uplo, &n, A_save, &lda, ipiv, work, &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();

            /*  call to API */
            invoke_hetrf_rook(datatype, &uplo, &n, A_save, &lda, ipiv, work, &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }
        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        free_vector(work);
    }

    *time_min_ = time_min;

    /* Save the output to vector A */
    copy_matrix(datatype, "full", lda, n, A_save, lda, A, lda);
    free_matrix(A_save);
}

double prepare_lapacke_hetrf_rook_run(integer datatype, integer layout, char uplo, integer n,
                                      void *A, integer lda, void *ipiv, integer *info)
{
    double exe_time = 0;
    integer lda_t = lda;
    void *A_t = NULL;
    A_t = A;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_hetrf_rook_lda, lda_t);

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        create_matrix(datatype, layout, n, n, &A_t, fla_max(lda_t, n));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, fla_max(lda_t, n));
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE hetrf_rook API */
    *info = invoke_lapacke_hetrf_rook(datatype, layout, uplo, n, A_t, lda_t, ipiv);

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

/*
HETRF_ROOK_API calls LAPACK interface for factorization
of a complex hermitian matrix A using the bounded
Bunch-Kaufman("rook") diagonal pivoting method
(A = L*D*L**H or A = U*D*U**H)
*/
void invoke_hetrf_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                       integer *ipiv, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case COMPLEX:
        {
            fla_lapack_chetrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhetrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }
    }
}

/*
LAPACKE HETRF_ROOK API invoke function
*/
integer invoke_lapacke_hetrf_rook(integer datatype, integer layout, char uplo, integer n, void *a,
                                  integer lda, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case COMPLEX:
        {
            info = LAPACKE_chetrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zhetrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }
    }
    return info;
}
