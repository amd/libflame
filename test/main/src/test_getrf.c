/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#define GETRF_VL 0.1
#define GETRF_VU 10

/* Local prototypes */
void fla_test_getrf_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_getrf_run(integer m_A, integer n_A, void *A, integer lda, integer *ipiv,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer test_lapacke_interface, int matrix_layout);
void invoke_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv,
                  integer *info);
double prepare_lapacke_getrf_run(integer datatype, int matrix_layout, integer m_A, integer n_A,
                                 void *A, integer lda, integer *ipiv, integer *info);
integer invoke_lapacke_getrf(integer datatype, int matrix_layout, integer m, integer n, void *a,
                             integer lda, integer *ipiv);

void fla_test_getrf(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LU factorization";
    char *front_str = "GETRF";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_getrf_experiment);
        tests_not_run = 0;
    }
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if(argc >= 7 && argc <= 8)
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_getrf_experiment(params, datatype, M, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, RECT_INPUT, M, N, residual,
                                      params->lin_solver_paramslist[0].solver_threshold, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for getrf\n");
        printf("./<EXE> getrf <precisions - sdcz> <M> <N> <LDA> <repeats>\n");
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

void fla_test_getrf_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer m, n, lda, info = 0, vinfo = 0;
    void *IPIV = NULL, *A = NULL, *A_test = NULL, *s_test = NULL;
    char range = 'U';
    double time_min = 1e9;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Determine the dimensions*/
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_vector(INTEGER, &IPIV, fla_min(m, n));
    create_realtype_vector(datatype, &s_test, fla_min(m, n));

    /* Initialize the test matrices*/
    if(g_ext_fptr != NULL || FLA_EXTREME_CASE_TEST)
    {
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, m, n, A, lda, s_test, GETRF_VL, GETRF_VU, i_zero, i_zero,
                          info);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_getrf(datatype, m, n, A, lda, params->imatrix_char);
        }
    }
    /* Save the original matrix*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    /* call to API */
    prepare_getrf_run(m, n, A_test, lda, IPIV, datatype, n_repeats, &time_min, &info,
                      test_lapacke_interface, layout);

    /* execution time */
    *t = time_min;

    /* performance computation */
    if(m == n)
    {
        *perf = (2.0 / 3.0) * n * n * n / time_min / FLOPS_PER_UNIT_PERF;
    }
    else if(m > n)
    {
        *perf = (1.0 / 3.0) * n * n * (3 * m - n) / time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        *perf = (1.0 / 3.0) * m * m * (3 * n - m) / time_min / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if((!FLA_EXTREME_CASE_TEST) && info == 0)
    {
        validate_getrf(m, n, A, A_test, lda, IPIV, datatype, residual, &vinfo,
                       params->imatrix_char);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char)))
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
    free_vector(s_test);
}

void prepare_getrf_run(integer m_A, integer n_A, void *A, integer lda, integer *IPIV,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer test_lapacke_interface, int layout)
{
    integer i;
    void *A_save;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {

        /* Copy original input data */
        copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

        /* Check if LAPACKE interface is enabled */
        if(test_lapacke_interface == 1)
        {
            exe_time
                = prepare_lapacke_getrf_run(datatype, layout, m_A, n_A, A_save, lda, IPIV, info);
        }
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK getrf API */
            invoke_getrf(datatype, &m_A, &n_A, A_save, &lda, IPIV, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the AFACT to matrix A */
    copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);
    free_matrix(A_save);
}

double prepare_lapacke_getrf_run(integer datatype, int layout, integer m_A, integer n_A,
                                 void *A, integer lda, integer *ipiv, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;
    A_t = A;

    if(layout == LAPACK_ROW_MAJOR)
    {
        lda_t = fla_max(1, n_A);
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m_A, n_A, &A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m_A, n_A, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /*  call LAPACKE getrf API */
    *info = invoke_lapacke_getrf(datatype, layout, m_A, n_A, A_t, lda_t, ipiv);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, m_A, n_A, A_t, lda_t, A, lda);
        /* free temporary buffers */
        free_matrix(A_t);
    }
    return exe_time;
}

/*
 *  GETRF_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
void invoke_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv,
                  integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgetrf(m, n, a, lda, ipiv, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgetrf(m, n, a, lda, ipiv, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgetrf(m, n, a, lda, ipiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgetrf(m, n, a, lda, ipiv, info);
            break;
        }
    }
}

integer invoke_lapacke_getrf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgetrf(layout, m, n, a, lda, ipiv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgetrf(layout, m, n, a, lda, ipiv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgetrf(layout, m, n, a, lda, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgetrf(layout, m, n, a, lda, ipiv);
            break;
        }
    }
    return info;
}