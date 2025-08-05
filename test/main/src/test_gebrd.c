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
integer row_major_gebrd_lda;

/* Local prototypes */
void fla_test_gebrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_gebrd_run(integer m_A, integer n_A, void *A, integer lda, void *d, void *e, void *tauq,
                       void *taup, void *work, integer lwork, integer datatype, integer n_repeats,
                       double *time_min_, integer interfacetype, integer *info,
                       test_params_t *params);
void invoke_gebrd(integer datatype, integer *m, integer *n, void *a, integer *lda, void *d, void *e,
                  void *tauq, void *taup, void *work, integer *lwork, integer *info);
double prepare_lapacke_gebrd_run(integer datatype, int layout, integer m_A, integer n_A, void *A,
                                 integer lda, void *d, void *e, void *tauq, void *taup,
                                 integer *info);
/*
 * Test driver for the GEBRD (bidiagonal reduction) routine.
 * This function parses command-line/config arguments, sets up the test,
 * and calls the experiment runner for each datatype.
 */
void fla_test_gebrd(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Bidiagonal reduction";
    char *front_str = "GEBRD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    /* If no arguments, run all tests from config file */
    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gebrd_experiment);
        tests_not_run = 0;
    }
    /* Parse extra argument for error info or matrix input */
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }
    /* Parse command-line arguments for direct invocation */
    if(argc >= 8 && argc <= 9)
    {
        integer i, num_types, N, M;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gebrd_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = M;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        /* Run tests for each requested datatype */
        if(n_repeats > 0)
        {
            params->lin_solver_paramslist[0].solver_threshold = CLI_NORM_THRESH;
            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);
                if(datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;
                fla_test_gebrd_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }
    /* Print error/help messages if needed */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gebrd\n");
        printf("./<EXE> gebrd <precisions - sdcz> <M> <N> <LDA> <lwork> <repeats>\n");
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

/*
 * Run a single GEBRD experiment for the given datatype and matrix size.
 * Allocates and initializes matrices, performs workspace query, times the routine,
 * and (optionally) validates the result.
 */
void fla_test_gebrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, lda, lwork = -1, info = 0;
    void *A = NULL, *A_test = NULL, *d = NULL, *e = NULL, *tauq = NULL, *taup = NULL, *work = NULL;
    double residual, err_thresh;
    integer interfacetype = params->interfacetype;
    double time_min_local = 0.0;
    integer realtype = get_realtype(datatype);

    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    if(g_config_data && lda == -1)
        lda = fla_max(1, m);

    /* Create input/output arrays for the test */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_realtype_vector(datatype, &d, fla_min(m, n));
    create_realtype_vector(datatype, &e, fla_min(m, n) - 1);
    create_vector(datatype, &tauq, fla_min(m, n));
    create_vector(datatype, &taup, fla_min(m, n));

    /* Initialize input matrix with random or file-based values */
    init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
        scale_matrix_underflow_overflow_labrd(datatype, m, n, A, lda, params->imatrix_char);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    /* Workspace query: call GEBRD with lwork=-1 to get optimal size if g_lwork <= 0 */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork <= 0))
    {
        integer lwork_query = -1;
        integer info_query = 0;
        create_vector(datatype, &work, 1); // minimal for query
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_gebrd(datatype, &m, &n, NULL, &lda, NULL, NULL, NULL, NULL, work,
                             &lwork_query, &info_query);
        }
        else
#endif
        {
            invoke_gebrd(datatype, &m, &n, NULL, &lda, NULL, NULL, NULL, NULL, work, &lwork_query,
                         &info_query);
        }
        if(info_query == 0)
        {
            /* Get work size */
            lwork = get_work_value(datatype, work);
        }
    }
    else
    {
        lwork = g_lwork;
    }
    free_vector(work);
    create_vector(datatype, &work, lwork);

    /* Run the actual test and measure minimum execution time */
    prepare_gebrd_run(m, n, A_test, lda, d, e, tauq, taup, work, lwork, datatype, n_repeats,
                      &time_min_local, interfacetype, &info, params);
    time_min = time_min_local;

    /* Performance computation (see LABRD for formula)
     * The number of floating point operations in GEBRD is 4n^2(3m - n)/3 if m>=n else 4m^2(3n-m)/3
     * See: https://support.nag.com/numeric/nl/nagdoc_latest/clhtml/f08/f08kec.html
     */
    perf = (double)(((4.0 * n * n) * ((3.0 * m) - n)) / 3.0) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        /* Validate the GEBRD result */
        validate_gebrd(datatype, tst_api, m, n, A, lda, A_test, lda, d, e, tauq, taup, residual,
                       params);
    }
    else
    {
        /* Check for output matrix and vectors when inputs are extreme values */
        int extreme_fail = 0;
        if(!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char))
            extreme_fail = 1;
        if(!check_extreme_value(realtype, fla_min(m, n), 1, d, 1, params->imatrix_char))
            extreme_fail = 1;
        if(!check_extreme_value(realtype, fla_min(m, n) - 1, 1, e, 1, params->imatrix_char))
            extreme_fail = 1;
        if(!check_extreme_value(datatype, fla_min(m, n), 1, tauq, 1, params->imatrix_char))
            extreme_fail = 1;
        if(!check_extreme_value(datatype, fla_min(m, n), 1, taup, 1, params->imatrix_char))
            extreme_fail = 1;

        if(extreme_fail)
            residual = DBL_MAX;
        else
            residual = err_thresh;

        /* Print the test status */
        FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    }

    /* Free all allocated memory */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(d);
    free_vector(e);
    free_vector(tauq);
    free_vector(taup);
    free_vector(work);
}

/*
 * Helper to run GEBRD multiple times and record the minimum execution time.
 * Allocates fresh output buffers for each repeat, and copies results to output.
 */
void prepare_gebrd_run(integer m_A, integer n_A, void *A, integer lda, void *d, void *e, void *tauq,
                       void *taup, void *work, integer lwork, integer datatype, integer n_repeats,
                       double *time_min_, integer interfacetype, integer *info,
                       test_params_t *params)
{
    void *A_save = NULL, *d_test = NULL, *e_test = NULL, *tauq_test = NULL, *taup_test = NULL,
         *work_test = NULL;
    double exe_time;
    integer min_mn = fla_min(m_A, n_A);
    int layout = params->matrix_major;

    /* Make a copy of the input matrix A. Same input values will be passed in each iteration. */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrix and allocate output buffers for each repeat */
        copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);
        create_realtype_vector(datatype, &d_test, min_mn);
        create_realtype_vector(datatype, &e_test, min_mn - 1);
        create_vector(datatype, &tauq_test, min_mn);
        create_vector(datatype, &taup_test, min_mn);
        create_vector(datatype, &work_test, lwork);

        exe_time = fla_test_clock();
        /* Check if LAPACKE interface is enabled
           and call gebrd API based on interface type */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gebrd_run(datatype, layout, m_A, n_A, A, lda, d_test, e_test,
                                                 tauq_test, taup_test, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            invoke_cpp_gebrd(datatype, &m_A, &n_A, A, &lda, d_test, e_test, tauq_test, taup_test,
                             work_test, &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
        else
#endif
        {
            exe_time = fla_test_clock();
            /* Call the GEBRD routine */
            invoke_gebrd(datatype, &m_A, &n_A, A, &lda, d_test, e_test, tauq_test, taup_test,
                         work_test, &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Copy results to output buffers for validation/performance reporting */
        copy_realtype_vector(datatype, min_mn, d_test, 1, d, 1);
        copy_realtype_vector(datatype, min_mn - 1, e_test, 1, e, 1);
        copy_vector(datatype, min_mn, tauq_test, 1, tauq, 1);
        copy_vector(datatype, min_mn, taup_test, 1, taup, 1);

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO

        free_vector(d_test);
        free_vector(e_test);
        free_vector(tauq_test);
        free_vector(taup_test);
        free_vector(work_test);
    }

    *time_min_ = time_min;

    /* Free the saved input matrix */
    free_matrix(A_save);
}

double prepare_lapacke_gebrd_run(integer datatype, int layout, integer m_A, integer n_A, void *A,
                                 integer lda, void *d, void *e, void *tauq, void *taup,
                                 integer *info)
{
    double exe_time;
    void *A_t = A;
    integer lda_t = lda;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n_A, row_major_gebrd_lda, lda_t);

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m_A, n_A, &A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m_A, n_A, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* call LAPACKE gesdd API */
    *info = invoke_lapacke_gebrd(datatype, layout, m_A, n_A, A_t, lda_t, d, e, tauq, taup);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Convert output matrix back to column_major layout if needed */
        convert_matrix_layout(LAPACK_ROW_MAJOR, datatype, m_A, n_A, A_t, lda_t, A, lda);
        free_matrix(A_t);
    }

    return exe_time;
}

/*
 * Helper to call the correct GEBRD routine for the given datatype.
 */
void invoke_gebrd(integer datatype, integer *m, integer *n, void *a, integer *lda, void *d, void *e,
                  void *tauq, void *taup, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
            break;
        }
    }
}

