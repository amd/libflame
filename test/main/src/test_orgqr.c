/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_orgqr_lda;

/* Local prototypes.*/
void fla_test_orgqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_orgqr_run(integer m, integer n, void *A, integer lda, void *T, integer datatype,
                       integer *info, integer interfacetype, int matrix_layout,
                       test_params_t *params);
void invoke_orgqr(integer datatype, integer *m, integer *n, integer *min_A, void *a, integer *lda,
                  void *tau, void *work, integer *lwork, integer *info);
double prepare_lapacke_orgqr_run(integer datatype, int matrix_layout, integer m, integer n, void *A,
                                 integer lda, void *T, integer *info);

void fla_test_orgqr(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "QR factorization";
    char *front_str = "ORGQR";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_orgqr_experiment);
        tests_not_run = 0;
    }
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }
    if(argc >= 8 && argc <= 9)
    {
        integer i, num_types, N, M;
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
            row_major_orgqr_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

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
                fla_test_orgqr_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for orgqr/ungqr\n");
        printf("./<EXE> orgqr <precisions - sd> <M> <N> <lda> <LWORK> <repeats>\n");
        printf("./<EXE> ungqr <precisions - cz> <M> <N> <lda> <LWORK> <repeats>\n");
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

void fla_test_orgqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, lda;
    void *A = NULL, *A_test = NULL, *T_test = NULL;
    void *work = NULL;
    void *Q = NULL, *R = NULL;
    integer lwork = -1, info = 0;
    double residual, err_thresh;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    time_min = 0.;
    perf = 0.;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* When inputs are from config file,
       1) if m < n(invalid case), interchange m, n
       2) if leading dimensions = -1, set them to default value */
    if(g_config_data)
    {
        if(p_cur < q_cur)
        {
            m = q_cur;
            n = p_cur;
        }
        if(lda == -1 || lda < m)
        {
            lda = fla_max(1, m);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);

    /* create Q matrix to check orthogonality */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &Q, lda);
    reset_matrix(datatype, m, n, Q, lda);

    /* create tau vector */
    create_vector(datatype, &T_test, fla_min(m, n));

    if(!FLA_BRT_VERIFICATION_RUN)
    {
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);

        /* Scaling matrix with values around overflow, underflow for ORGQR/UNGQR */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_orgqr(datatype, m, n, A, lda, params->imatrix_char);
        }

        /* Make a copy of input matrix A. This is required to validate the API functionality.*/
        copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

        /* Make a workspace query the first time. This will provide us with
        and ideal workspace size based on internal block size for geqrf.*/
        create_vector(datatype, &work, 1);
        /* call to  geqrf API */
        invoke_geqrf(datatype, &m, &n, NULL, &lda, NULL, work, &lwork, &info);
        if(info == 0)
        {
            /* Get work size */
            lwork = get_work_value(datatype, work);
        }
        else
        {
            lwork = fla_max(1, n);
        }
        /* Free the query workspace */
        free_vector(work);
        work = NULL;

        /* create work buffer */
        create_vector(datatype, &work, lwork);

        /* QR Factorisation on matrix A to generate Q and R */
        invoke_geqrf(datatype, &m, &n, A_test, &lda, T_test, work, &lwork, &info);

        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &R, n);
        reset_matrix(datatype, n, n, R, n);
        copy_matrix(datatype, "Upper", n, n, A_test, lda, R, n);

        copy_matrix(datatype, "full", m, n, A_test, lda, Q, lda);
    }
    FLA_BRT_PROCESS_TWO_INPUT(datatype, m, n, Q, lda, datatype, 1, fla_min(m, n), T_test, 1, "dddd",
                              m, n, lda, g_lwork)
    /*invoke orgqr API */
    prepare_orgqr_run(m, n, Q, lda, T_test, datatype, &info, interfacetype, layout, params);

    /* performance computation
    (2/3)*n2*(3m - n) */
    perf = (double)((2.0 * m * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    IF_FLA_BRT_VALIDATION(
        m, n, store_outputs_base(filename, params, 1, 0, datatype, m, n, Q, lda),
        validate_orgqr(tst_api, m, n, A, lda, Q, R, datatype, residual, params->imatrix_char,
                       params),
        check_reproducibility_base(filename, params, 1, 0, datatype, m, n, Q, lda))
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_orgqr(tst_api, m, n, A, lda, Q, R, datatype, residual, params->imatrix_char,
                       params);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if((!check_extreme_value(datatype, m, n, T_test, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    }

    /* Free up the buffers */
free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(work);
    free_vector(T_test);
    free_matrix(Q);
    free_matrix(R);
}

void prepare_orgqr_run(integer m, integer n, void *A, integer lda, void *T, integer datatype,
                       integer *info, integer interfacetype, int layout, test_params_t *params)
{
    integer lwork;
    void *A_save = NULL, *work = NULL;
    double exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_save, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_save, lda);
    /* Make a workspace query the first time. This will provide us with
       and ideal workspace size based on internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* call to  geqrf API */
        invoke_orgqr(datatype, &m, &n, &n, NULL, &lda, NULL, work, &lwork, info);

        if(*info == 0)
        {
            /* Get work size */
            lwork = get_work_value(datatype, work);
        }
        /* Output buffers will be freshly allocated for each iterations, free up
           the current output buffers.*/
        free_vector(work);
        work = NULL;
    }
    else
    {
        lwork = g_lwork;
    }

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m, n, A_save, lda, A, lda);
        create_vector(datatype, &work, lwork);
        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_orgqr_run(datatype, layout, m, n, A, lda, T, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP orgqr API */
            invoke_cpp_orgqr(datatype, &m, &n, &n, A, &lda, T, work, &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK orgqr API */
            invoke_orgqr(datatype, &m, &n, &n, A, &lda, T, work, &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
        free_vector(work);
    }

    free_matrix(A_save);
}

double prepare_lapacke_orgqr_run(integer datatype, int layout, integer m, integer n, void *A,
                                 integer lda, void *T, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_orgqr_lda, lda_t);

    A_t = A;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m, n, &A_t, fla_max(n, lda_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m, n, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* Call to LAPACKE orgqr API */
    *info = invoke_lapacke_orgqr(datatype, layout, m, n, n, A_t, lda_t, T);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, m, n, A_t, lda_t, A, lda);
        /* free temporary buffers */
        free_matrix(A_t);
    }

    return exe_time;
}

void invoke_orgqr(integer datatype, integer *m, integer *n, integer *min_A, void *a, integer *lda,
                  void *tau, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sorgqr(m, n, n, a, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dorgqr(m, n, n, a, lda, tau, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cungqr(m, n, n, a, lda, tau, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zungqr(m, n, n, a, lda, tau, work, lwork, info);
            break;
        }
    }
}
