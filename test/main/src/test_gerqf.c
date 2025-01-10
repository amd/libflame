/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;
integer row_major_gerqf_lda;

// Local prototypes.
void fla_test_gerqf_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_gerqf_run(integer m_A, integer n_A, void *A, integer lda, void *T, integer datatype,
                       integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                       int matrix_layout);
void invoke_gerqf(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                  void *work, integer *lwork, integer *info);
double prepare_lapacke_gerqf_run(integer datatype, int matrix_layout, integer m_A, integer n_A,
                                 void *A, integer lda, void *T, integer *info);
integer invoke_lapacke_gerqf(integer datatype, int matrix_layout, integer m, integer n, void *a,
                             integer lda, void *tau);

void fla_test_gerqf(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "RQ factorization";
    char *front_str = "GERQF";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gerqf_experiment);
        tests_not_run = 0;
    }
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }
    if(argc >= 8 && argc <= 9)
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gerqf_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_gerqf_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gerqf\n");
        printf("./<EXE> gerqf <precisions - sdcz> <M> <N> <LDA> <LWORK> <repeats>\n");
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

void fla_test_gerqf_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, lda;
    integer info = 0;
    void *A, *A_test, *T;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    // Get input matrix dimensions.
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
    }

    // Create input matrix parameters
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_vector(datatype, &T, fla_min(m, n));

    init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_underflow_overflow_geqrf(datatype, m, n, A, lda, params->imatrix_char);
    }

    // Make a copy of input matrix A. This is required to validate the API functionality.
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_gerqf_run(m, n, A_test, lda, T, datatype, n_repeats, &time_min, &info, interfacetype,
                      layout);

    /* performance computation */
    if(m >= n)
        perf = (double)((2.0 * m * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min
               / FLOPS_PER_UNIT_PERF;
    else
        perf = (double)((2.0 * n * m * m) - ((2.0 / 3.0) * m * m * m)) / time_min
               / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_gerqf(tst_api, m, n, A, A_test, lda, T, datatype, residual);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if(!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    }

    // Free up the buffers
    free_matrix(A);
    free_matrix(A_test);
    free_vector(T);
}

void prepare_gerqf_run(integer m_A, integer n_A, void *A, integer lda, void *T, integer datatype,
                       integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                       int layout)
{
    integer min_A, i;
    void *A_save, *T_test, *work;
    integer lwork = -1;
    double t_min = 1e9, exe_time;

    min_A = fla_min(m_A, n_A);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    /* Make a workspace query the first time. This will provide us with
       and ideal workspace size based on internal block size.
       NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call to gerqf API */
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_gerqf(datatype, &m_A, &n_A, NULL, &lda, NULL, work, &lwork, info);
        }
        else
#endif
        {
            invoke_gerqf(datatype, &m_A, &n_A, NULL, &lda, NULL, work, &lwork, info);
        }
        if(*info == 0)
        {
            /* Get work size */
            lwork = get_work_value(datatype, work);
        }

        /* Output buffers will be freshly allocated for each iterations, free up
       the current output buffers.*/
        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);

        // T_test vector will hold the scalar factors of the elementary reflectors.
        create_vector(datatype, &T_test, min_A);

        // create work buffer
        create_vector(datatype, &work, lwork);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gerqf_run(datatype, layout, m_A, n_A, A, lda, T_test, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP gerqf API */
            invoke_cpp_gerqf(datatype, &m_A, &n_A, A, &lda, T_test, work, &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK gerqf API */
            invoke_gerqf(datatype, &m_A, &n_A, A, &lda, T_test, work, &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        // Get the best execution time
        t_min = fla_min(t_min, exe_time);

        // Make a copy of the output buffers. This is required to validate the API functionality.
        copy_vector(datatype, min_A, T_test, 1, T, 1);

        // Free up the output buffers
        free_vector(work);
        free_vector(T_test);
    }

    *time_min_ = t_min;

    free_matrix(A_save);
}

double prepare_lapacke_gerqf_run(integer datatype, int layout, integer m_A, integer n_A, void *A,
                                 integer lda, void *T, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_gerqf_lda, lda_t);

    A_t = A;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m_A, n_A, &A_t, fla_max(n_A, lda_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m_A, n_A, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* Call LAPACKE gerqf API */
    *info = invoke_lapacke_gerqf(datatype, layout, m_A, n_A, A_t, lda_t, T);

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

void invoke_gerqf(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                  void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgerqf(m, n, a, lda, tau, work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dgerqf(m, n, a, lda, tau, work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cgerqf(m, n, a, lda, tau, work, lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgerqf(m, n, a, lda, tau, work, lwork, info);
            break;
        }
    }
}

integer invoke_lapacke_gerqf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgerqf(layout, m, n, a, lda, tau);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dgerqf(layout, m, n, a, lda, tau);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cgerqf(layout, m, n, a, lda, tau);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgerqf(layout, m, n, a, lda, tau);
            break;
        }
    }
    return info;
}
