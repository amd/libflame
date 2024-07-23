/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#define GESV_VL 0.1
#define GESV_VU 10

/* Local prototypes */
void fla_test_gesv_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual);
void prepare_gesv_run(integer n_A, integer nrhs, void *A, integer lda, void *B, integer ldb,
                      integer *ipiv, integer datatype, integer n_repeats, double *time_min_,
                      integer *info, integer test_lapacke_interface, integer matrix_layout);

void invoke_gesv(integer datatype, integer *nrhs, integer *n, void *a, integer *lda, integer *ipiv,
                 void *b, integer *ldb, integer *info);
double prepare_lapacke_gesv_run(integer datatype, integer matrix_layout, integer n_A, integer nrhs,
                                void *A, integer *lda, void *B, integer *ldb, integer *ipiv,
                                integer *info);
integer invoke_lapacke_gesv(integer datatype, integer matrix_layout, integer n, integer nrhs,
                            void *a, integer lda, integer *ipiv, void *b, integer ldb);

void fla_test_gesv(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Linear Solve using LU", *front_str = "GESV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gesv_experiment);
        tests_not_run = 0;
    }
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }
    if(argc >= 8 && argc <= 9)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_gesv_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
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
        printf("\nIllegal arguments for gesv\n");
        printf("./<EXE> gesv <precisions - sdcz>  <N> <NRHS> <LDA> <LDB> <repeats>\n");
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

void fla_test_gesv_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual)
{
    integer n, lda, ldb, NRHS, info = 0, vinfo = 0;
    void *IPIV = NULL, *A = NULL, *A_save = NULL, *B = NULL, *B_save = NULL, *s_test = NULL;
    double time_min = 1e9;
    char range = 'U';
    integer test_lapacke_interface = params->test_lapacke_interface;
    integer matrix_layout = params->matrix_major;
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
    create_matrix(datatype, matrix_layout, n, n, &A, lda);
    create_matrix(datatype, matrix_layout, n, n, &A_save, lda);
    create_vector(INTEGER, &IPIV, n);
    create_matrix(datatype, matrix_layout, n, NRHS, &B, ldb);
    create_matrix(datatype, matrix_layout, n, NRHS, &B_save, ldb);
    create_realtype_vector(datatype, &s_test, n);

    /* Initialize the test matrices*/
    if(params->imatrix_char == NULL && g_ext_fptr == NULL)
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, n, n, A, lda, s_test, GESV_VL, GESV_VU, i_zero, i_zero,
                          info);
    }
    else
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    init_matrix(datatype, B, n, NRHS, ldb, g_ext_fptr, params->imatrix_char);

    /* Save the original matrix*/
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
    copy_matrix(datatype, "full", n, NRHS, B, ldb, B_save, ldb);
    /* call to API */
    prepare_gesv_run(n, NRHS, A_save, lda, B_save, ldb, IPIV, datatype, n_repeats, &time_min, &info,
                     test_lapacke_interface, matrix_layout);

    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    *perf
        = (double)((2.0 * n * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if(!params->imatrix_char && info == 0)
        validate_gesv(n, NRHS, A, lda, B, ldb, B_save, datatype, residual, &vinfo);
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, n, A_save, lda, params->imatrix_char))
           && (!check_extreme_value(datatype, n, NRHS, B_save, ldb, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_save);
    free_vector(IPIV);
    free_matrix(B);
    free_matrix(B_save);
    free_vector(s_test);
}

void prepare_gesv_run(integer n_A, integer nrhs, void *A, integer lda, void *B, integer ldb,
                      integer *IPIV, integer datatype, integer n_repeats, double *time_min_,
                      integer *info, integer test_lapacke_interface, integer matrix_layout)
{
    integer i;
    void *A_test, *B_test;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, matrix_layout, n_A, n_A, &A_test, lda);
    create_matrix(datatype, matrix_layout, n_A, nrhs, &B_test, ldb);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {

        /* Copy original input data */
        copy_matrix(datatype, "full", n_A, n_A, A, lda, A_test, lda);
        copy_matrix(datatype, "full", n_A, nrhs, B, ldb, B_test, ldb);
        /* Check if LAPACKE interface is enabled */
        if(test_lapacke_interface == 1)
        {
            exe_time = prepare_lapacke_gesv_run(datatype, matrix_layout, n_A, nrhs, A_test, &lda,
                                                B_test, &ldb, IPIV, info);
        }
        else
        {
            exe_time = fla_test_clock();
            /*  call  gesv API  */
            invoke_gesv(datatype, &n_A, &nrhs, A_test, &lda, IPIV, B_test, &ldb, info);
            exe_time = fla_test_clock() - exe_time;
        }
        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the final result to B matrix*/
    copy_matrix(datatype, "full", n_A, nrhs, B_test, ldb, B, ldb);

    free_matrix(A_test);
    free_matrix(B_test);
}

double prepare_lapacke_gesv_run(integer datatype, integer matrix_layout, integer n_A, integer nrhs,
                                void *A, integer *lda, void *B, integer *ldb, integer *ipiv,
                                integer *info)
{
    double exe_time = 0;
    integer lda_t = *lda;
    integer ldb_t = *ldb;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(matrix_layout == LAPACK_ROW_MAJOR)
    {
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, *lda, n_A, A, *lda, &lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, *ldb, nrhs, B, *ldb, &ldb_t);
    }

    exe_time = fla_test_clock();

    /*  call LAPACKE gesv API  */
    *info = invoke_lapacke_gesv(datatype, matrix_layout, n_A, nrhs, A, lda_t, ipiv, B, ldb_t);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */
    if(matrix_layout == LAPACK_ROW_MAJOR)
    {
        convert_matrix_layout(matrix_layout, datatype, *lda, n_A, A, lda_t, lda);
        convert_matrix_layout(matrix_layout, datatype, *ldb, nrhs, B, ldb_t, ldb);
    }
    return exe_time;
}

/*
 *  gesv_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
void invoke_gesv(integer datatype, integer *n, integer *nrhs, void *a, integer *lda, integer *ipiv,
                 void *b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgesv(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgesv(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgesv(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgesv(n, nrhs, a, lda, ipiv, b, ldb, info);
            break;
        }
    }
}

integer invoke_lapacke_gesv(integer datatype, integer matrix_layout, integer n, integer nrhs,
                            void *a, integer lda, integer *ipiv, void *b, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgesv(matrix_layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgesv(matrix_layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgesv(matrix_layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgesv(matrix_layout, n, nrhs, a, lda, ipiv, b, ldb);
            break;
        }
    }
    return info;
}