/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

#define GETRS_VL 0.1
#define GETRS_VU 10

integer row_major_getrs_lda;
integer row_major_getrs_ldb;

/* Local prototypes */
void fla_test_getrs_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_getrs_run(char *trans, integer m_A, integer n_A, void *A, integer lda, void *B,
                       integer ldb, integer *ipiv, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer interfacetype,
                       int matrix_layout);
void invoke_getrs(integer datatype, char *trans, integer *nrhs, integer *n, void *a, integer *lda,
                  integer *ipiv, void *b, integer *ldb, integer *info);
double prepare_lapacke_getrs_run(integer datatype, int matrix_layout, char *trans, integer n,
                                 integer nrhs, void *A, integer lda, void *B, integer ldb,
                                 integer *ipiv, integer *info);
integer invoke_lapacke_getrs(integer datatype, int matrix_layout, char trans, integer n,
                             integer nrhs, const void *a, integer lda, const integer *ipiv, void *b,
                             integer ldb);

void fla_test_getrs(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LU factorization";
    char *front_str = "GETRS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_getrs_experiment);
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
        params->lin_solver_paramslist[0].transr = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_getrs_lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            row_major_getrs_ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_getrs_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
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
        printf("\nIllegal arguments for getrs\n");
        printf("./<EXE> getrs <precisions - sdcz> <TRANS> <N> <NRHS> <LDA> <LDB> <repeats>\n");
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

void fla_test_getrs_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer n, lda, ldb, NRHS;
    integer info = 0, vinfo = 0;
    void *IPIV, *s_test = NULL;
    char range = 'U';
    void *A, *A_test, *B, *B_save, *X, *scal = NULL;
    double time_min = 1e9;
    char TRANS = params->lin_solver_paramslist[pci].transr;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    NRHS = params->lin_solver_paramslist[pci].nrhs;

    /* Determine the dimensions*/
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

    /* Create the matrices for the current operation*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_vector(INTEGER, &IPIV, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &B, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &B_save, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &X, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    create_realtype_vector(datatype, &s_test, n);

    /* Initialize the test matrices*/
    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, n, n, A, lda, s_test, GETRS_VL, GETRS_VU, i_zero, i_zero,
                          info);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_realtype_vector(get_datatype(datatype), &scal, n);
            scale_matrix_underflow_overflow_getrs(datatype, &TRANS, n, n, A, lda,
                                                  params->imatrix_char, scal);
        }
    }
    init_matrix(datatype, B, n, NRHS, ldb, g_ext_fptr, params->imatrix_char);

    /* Save the original matrix*/

    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    copy_matrix(datatype, "full", n, NRHS, B, ldb, B_save, ldb);

    /*  call to API getrf to get AFACT */
    invoke_getrf(datatype, &n, &n, A_test, &lda, IPIV, &info);

    /* call to API */
    prepare_getrs_run(&TRANS, n, NRHS, A_test, lda, B, ldb, IPIV, datatype, n_repeats, &time_min,
                      &info, interfacetype, layout);
    copy_matrix(datatype, "full", n, NRHS, B, ldb, X, ldb);
    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2*n^2 * nrhs flops */
    *perf = (double)(2.0 * n * n * NRHS) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if(info == 0 && !FLA_EXTREME_CASE_TEST)
        validate_getrs(&TRANS, n, NRHS, A, lda, B_save, ldb, X, datatype, residual, &vinfo,
                       params->imatrix_char, scal);
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if(!check_extreme_value(datatype, n, NRHS, B, ldb, params->imatrix_char))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
    free_matrix(B);
    free_matrix(X);
    free_matrix(B_save);
    free_vector(s_test);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal);
    }
}

void prepare_getrs_run(char *TRANS, integer n_A, integer nrhs, void *A, integer lda, void *B,
                       integer ldb, integer *IPIV, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer interfacetype,
                       int layout)
{
    integer i;
    void *A_save, *B_test;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, nrhs, &B_test, ldb);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
        copy_matrix(datatype, "full", n_A, nrhs, B, ldb, B_test, ldb);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_getrs_run(datatype, layout, TRANS, n_A, nrhs, A_save, lda,
                                                 B_test, ldb, IPIV, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)   /* Call CPP getrs API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_getrs(datatype, TRANS, &n_A, &nrhs, A_save, &lda, IPIV, B_test, &ldb, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK getrs API */
            invoke_getrs(datatype, TRANS, &n_A, &nrhs, A_save, &lda, IPIV, B_test, &ldb, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the final result to B matrix*/
    copy_matrix(datatype, "full", n_A, nrhs, B_test, ldb, B, ldb);

    free_matrix(A_save);
    free_vector(B_test);
}

double prepare_lapacke_getrs_run(integer datatype, int layout, char *trans, integer n_A,
                                 integer nrhs, void *A, integer lda, void *B, integer ldb,
                                 integer *ipiv, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldb_t = ldb;
    void *A_t = NULL, *B_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_getrs_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, nrhs, row_major_getrs_ldb, ldb_t);

    A_t = A;
    B_t = B;
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n_A, n_A, &A_t, fla_max(n_A, lda_t));
        create_matrix(datatype, layout, n_A, nrhs, &B_t, fla_max(nrhs, ldb_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n_A, n_A, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n_A, nrhs, B, ldb, B_t, ldb_t);
    }

    exe_time = fla_test_clock();

    /* Call LAPACKE getrs API */
    *info = invoke_lapacke_getrs(datatype, layout, *trans, n_A, nrhs, A_t, lda_t, ipiv, B_t, ldb_t);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n_A, n_A, A_t, lda_t, A, lda);
        convert_matrix_layout(layout, datatype, n_A, nrhs, B_t, ldb_t, B, ldb);
        /* free temporary buffers */
        free_matrix(A_t);
        free_matrix(B_t);
    }
    return exe_time;
}

/*
 *  GETRS_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
void invoke_getrs(integer datatype, char *trans, integer *n, integer *nrhs, void *a, integer *lda,
                  integer *ipiv, void *b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }
    }
}

integer invoke_lapacke_getrs(integer datatype, int layout, char trans, integer n, integer nrhs,
                             const void *a, integer lda, const integer *ipiv, void *b, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgetrs(layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }
    }
    return info;
}