/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include <invoke_lapacke.h>

#define SYTRF_ROOK_VU 10.0 // Maximum eigen value for condition number.
#define SYTRF_ROOK_VL 0.01 // Minimum eigen value for condition number.

extern double perf;
extern double time_min;
integer row_major_sytrf_rook_lda;

void invoke_sytrf_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                       integer *ipiv, void *work, integer *lwork, integer *info);
void fla_test_sytrf_rook_experiment(char *tst_api, test_params_t *params, integer datatype,
                                    integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                    integer einfo);
void prepare_sytrf_rook_run(integer datatype, integer n, void *A, char uplo, integer lda,
                            integer *ipiv, void *work, integer lwork, integer n_repeats,
                            double *time_min_, integer *info, integer test_lapacke_interface,
                            integer mlayout);
double prepare_lapacke_sytrf_rook_run(integer datatype, integer layout, char uplo, integer n,
                                      void *A, integer lda, void *ipiv, integer *info);

void fla_test_sytrf_rook(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LDL**T / U**TDU factorization";
    char *front_str = "SYTRF_ROOK";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        config_data = 1;
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_sytrf_rook_experiment);
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
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        if((g_ext_fptr == NULL) && params->test_lapacke_interface
           && (params->matrix_major == LAPACK_ROW_MAJOR))
        {
            row_major_sytrf_rook_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_sytrf_rook_experiment(front_str, params, datatype, N, N, 0, n_repeats,
                                               einfo);
                tests_not_run = 0;
            }
        }
    }
    /* Print error message */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for sytrf_rook\n");
        printf("./<EXE> sytrf_rook <precisions - sdcz> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
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

void fla_test_sytrf_rook_experiment(char *tst_api, test_params_t *params, integer datatype,
                                    integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                    integer einfo)
{
    integer n, lda, lwork = -1, info = 0;
    void *A = NULL, *A_test = NULL, *ipiv = NULL, *work = NULL, *L = NULL;
    char uplo;
    double residual, err_thresh;

    integer test_lapacke_interface = params->test_lapacke_interface;
    integer layout = params->matrix_major;

    /* Determine the dimensions */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
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
        if(params->imatrix_char != '\0')
        {
            form_symmetric_matrix(datatype, n, A, lda, "S", 'U');
        }
    }
    else
    {
        /* Generating input matrix with condition number < 1000 */
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, 'V', n, A, lda, L, SYTRF_ROOK_VL, SYTRF_ROOK_VU,
                                 USE_ABS_EIGEN_VALUES);
        form_symmetric_matrix(datatype, n, A, lda, "S", 'U');
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
    prepare_sytrf_rook_run(datatype, n, A_test, uplo, lda, ipiv, work, lwork, n_repeats, &time_min,
                           &info, test_lapacke_interface, layout);

    /* Performance computation */
    perf = (double)(n * n * n) * (1.0 / 3.0) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* Output validataion */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_sytrf(tst_api, &uplo, n, lda, A_test, datatype, ipiv, residual, A);
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

    /* Free up buffers */
    free_vector(ipiv);
    free_matrix(A);
    free_matrix(A_test);
}

void prepare_sytrf_rook_run(integer datatype, integer n, void *A, char uplo, integer lda,
                            integer *ipiv, void *work, integer lwork, integer n_repeats,
                            double *time_min_, integer *info, integer test_lapacke_interface,
                            integer layout)
{
    integer i;
    void *A_save = NULL;
    double t_min = 1e9, exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.
     NOTE: LAPACKE interface handles workspace query internally */
    if((test_lapacke_interface == 0) && (g_lwork == -1))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* Getting lwork from api by passing lwork = -1 */
        invoke_sytrf_rook(datatype, &uplo, &n, NULL, &lda, ipiv, work, &lwork, info);
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
        if(test_lapacke_interface == 1)
        {
            exe_time = prepare_lapacke_sytrf_rook_run(datatype, layout, uplo, n, A_save, lda, ipiv,
                                                      info);
        }
        else
        {
            exe_time = fla_test_clock();

            /*  call to API */
            invoke_sytrf_rook(datatype, &uplo, &n, A_save, &lda, ipiv, work, &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }
        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        free_vector(work);
    }

    *time_min_ = t_min;

    /* Save the output to vector A */
    copy_matrix(datatype, "full", lda, n, A_save, lda, A, lda);
    free_matrix(A_save);
}

double prepare_lapacke_sytrf_rook_run(integer datatype, integer layout, char uplo, integer n,
                                      void *A, integer lda, void *ipiv, integer *info)
{
    double exe_time = 0;
    integer lda_t = lda;
    void *A_t = NULL;
    A_t = A;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_sytrf_rook_lda, lda_t);

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        create_matrix(datatype, layout, n, n, &A_t, fla_max(lda_t, n));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, fla_max(lda_t, n));
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE sytrf_rook API */
    *info = invoke_lapacke_sytrf_rook(datatype, layout, uplo, n, A_t, lda_t, ipiv);

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
SYTRF_ROOK_API calls LAPACK interface for factorization
of a real symmetric matrix A using the bounded Bunch-Kaufman ("rook")
diagonal pivoting method(A = LDL**T or A = U**TDU)
*/
void invoke_sytrf_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                       integer *ipiv, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dsytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_csytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zsytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
            break;
        }
    }
}