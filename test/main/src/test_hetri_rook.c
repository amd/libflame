/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

#define HETRI_ROOK_VU 10.0 // Maximum eigen value for condition number.
#define HETRI_ROOK_VL 0.1 // Minimum eigen value for condition number.

extern double perf;
extern double time_min;
integer row_major_hetri_rook_lda;

void invoke_hetri_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                       integer *ipiv, void *work, integer *info);
void fla_test_hetri_rook_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                                    integer q_cur, integer pci, integer n_repeats, integer einfo);
void prepare_hetri_rook_run(integer datatype, integer n, void *A, char uplo, integer lda,
                            integer *ipiv, void *work, integer n_repeats, double *time_min_,
                            integer *info, integer interfacetype, integer mlayout);

void fla_test_hetri_rook(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LDL**H / UDU**H Inverse";
    char *front_str = "HETRI_ROOK";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_hetri_rook_experiment);
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
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        if((g_ext_fptr == NULL) && params->interfacetype
           && (params->matrix_major == LAPACK_ROW_MAJOR))
        {
            row_major_hetri_rook_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->lin_solver_paramslist[0].solver_threshold = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalid datatype */
                if(datatype == FLOAT || datatype == DOUBLE || datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }

                /* Check for duplicate datatype */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_hetri_rook_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }
    /* Print error message */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for hetri_rook\n");
        printf("./<EXE> hetri_rook <precisions - cz> <UPLO> <N> <LDA> <repeats>\n");
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

void fla_test_hetri_rook_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                                    integer q_cur, integer pci, integer n_repeats, integer einfo)
{
    integer n, lda, info = 0;
    void *A = NULL, *A_test = NULL, *A_original = NULL, *ipiv = NULL, *work = NULL, *L = NULL;
    char uplo;
    double residual, err_thresh;
    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major, lwork;

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
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_original, lda);

    /* Initialize the test matrices */
    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST && !FLA_OVERFLOW_UNDERFLOW_TEST))
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        for(integer i = 0; i < n; i++)
        {
            ((integer*)ipiv)[i] = i + 1;
        }
        if(params->imatrix_char != '\0')
        {
            form_symmetric_matrix(datatype, n, A, lda, "C", 'U');
        }
    }
    else
    {
        /* Create input matrix (hermitian) */
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, 'V', n, A, lda, L, HETRI_ROOK_VL, HETRI_ROOK_VU,
                                 USE_SIGNED_EIGEN_VALUES);
        /* Oveflow or underflow test initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_overflow_underflow_hetri_rook(datatype, n, A, lda, params->imatrix_char);
        }
        form_symmetric_matrix(datatype, n, A, lda, "C", 'U');
        copy_matrix(datatype, "full", lda, n, A, lda, A_original, lda);
        free_vector(L);
        create_vector(datatype, &work, 1);
        lwork = -1;
        invoke_hetrf_rook(datatype, &uplo, &n, NULL, &lda, ipiv, work, &lwork, &info);
        if(info == 0)
        {
            lwork = get_work_value(datatype, work);
        }
        free_vector(work);
        info = 0;
        create_vector(datatype, &work, lwork);
        invoke_hetrf_rook(datatype, &uplo, &n, A, &lda, ipiv, work, &lwork, &info);
        free_vector(work);
    }
    /* Save the original matrix */
    copy_matrix(datatype, "full", lda, n, A, lda, A_test, lda);

    /* call to API */
    prepare_hetri_rook_run(datatype, n, A_test, uplo, lda, ipiv, work, n_repeats, &time_min, &info,
                           interfacetype, layout);

    /* Performance computation */
    perf = (double)(n * n * n) * (1.0 / 3.0) / time_min / FLOPS_PER_UNIT_PERF;
    perf *= 4.0;

    /* Output validataion */
        FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_hetri_rook(tst_api, uplo, n, A_original, A_test, lda, ipiv, datatype, residual,
                            params->imatrix_char);
    }
    else
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
    free_matrix(A_original);
}

void prepare_hetri_rook_run(integer datatype, integer n, void *A, char uplo, integer lda,
                            integer *ipiv, void *work, integer n_repeats, double *time_min_,
                            integer *info, integer interfacetype, integer layout)
{
    integer i, lwork;
    void *A_save = NULL;
    double t_min = 1e9, exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);

    /* Work buffer allocation */
    lwork = fla_max(1, n);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; i++)
    {
        /* Copy original input */
        copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

        /* Create work buffer */
        create_vector(datatype, &work, lwork);

#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP gecon API */
            invoke_cpp_hetri_rook(datatype, &uplo, &n, A_save, &lda, ipiv, work, info);

            exe_time = fla_test_clock() - exe_time;
        }
        else
        {
#endif
            exe_time = fla_test_clock();

            /*  call to API */
            invoke_hetri_rook(datatype, &uplo, &n, A_save, &lda, ipiv, work, info);

            exe_time = fla_test_clock() - exe_time;

#if ENABLE_CPP_TEST
        }
#endif
        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);
        free_vector(work);
    }

    *time_min_ = t_min;

    /* Save the output to vector A */
    copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
    free_matrix(A_save);
}

/*
HETRI_ROOK_API calls LAPACK interface for factorization
of a complex hermitian matrix A using the bounded
Bunch-Kaufman("rook") diagonal pivoting method
(A = L*D*L**H or A = U*D*U**H)
*/
void invoke_hetri_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                       integer *ipiv, void *work, integer *info)
{
    switch(datatype)
    {
        case COMPLEX:
        {
            fla_lapack_chetri_rook(uplo, n, a, lda, ipiv, work, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhetri_rook(uplo, n, a, lda, ipiv, work, info);
            break;
        }
    }
}
