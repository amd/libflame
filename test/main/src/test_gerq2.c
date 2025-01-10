/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;

// Local prototypes.
void fla_test_gerq2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_gerq2_run(integer m_A, integer n_A, void *A, integer lda, void *T, integer datatype,
                       integer n_repeats, double *time_min_, integer *info, integer interfacetype);
void invoke_gerq2(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                  void *work, integer *info);

void fla_test_gerq2(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "RQ factorization with unblocked algorithm";
    char *front_str = "GERQ2";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gerq2_experiment);
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
                fla_test_gerq2_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gerq2\n");
        printf("./<EXE> gerq2 <precisions - sdcz> <M> <N> <LDA> <repeats>\n");
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

void fla_test_gerq2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, lda;
    integer info = 0;
    void *A = NULL, *A_test = NULL, *T = NULL;
    double residual, err_thresh;
    integer interfacetype = params->interfacetype;

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
        scale_matrix_overflow_underflow_gerq2(datatype, m, n, A, lda, params->imatrix_char);
    }

    // Make a copy of input matrix A. This is required to validate the API functionality.
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_gerq2_run(m, n, A_test, lda, T, datatype, n_repeats, &time_min, &info, interfacetype);

    /* performance computation */
    if(m >= n)
        perf = (double)((2.0 * m * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min
               / FLOPS_PER_UNIT_PERF;
    else
        perf = (double)((2.0 * n * m * m) - ((2.0 / 3.0) * m * m * m)) / time_min
               / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    // output validation
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_gerq2(tst_api, m, n, A, A_test, lda, T, datatype, residual);
    }
    else
    {
        if((!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    }

    // Free up buffers
    free_matrix(A);
    free_matrix(A_test);
    free_vector(T);
}

void prepare_gerq2_run(integer m_A, integer n_A, void *A, integer lda, void *T, integer datatype,
                       integer n_repeats, double *time_min_, integer *info, integer interfacetype)
{
    integer min_A, i;
    void *A_save = NULL, *T_test = NULL, *work = NULL;
    double t_min = 1e9, exe_time;

    min_A = fla_min(m_A, n_A);

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);
        create_vector(datatype, &T_test, min_A);
        create_vector(datatype, &work, m_A);

#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST) /* Call CPP gerq2 API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_gerq2(datatype, &m_A, &n_A, A, &lda, T_test, work, info);
            exe_time = fla_test_clock() - exe_time;
        }
        else
#endif
        {
            exe_time = fla_test_clock();

            // call to gerq2 API
            invoke_gerq2(datatype, &m_A, &n_A, A, &lda, T_test, work, info);

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

void invoke_gerq2(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                  void *work, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgerq2(m, n, a, lda, tau, work, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgerq2(m, n, a, lda, tau, work, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgerq2(m, n, a, lda, tau, work, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgerq2(m, n, a, lda, tau, work, info);
            break;
        }
    }
}
