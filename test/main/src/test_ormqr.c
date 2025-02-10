/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;
integer row_major_ormqr_lda;
integer row_major_ormqr_ldc;

/* Local prototypes.*/
void fla_test_ormqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_ormqr_run(char side, char trans, integer m, integer n, integer k, integer m_A,
                       integer n_A, void *A, integer lda, void *tau, void *c, integer ldc,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int matrix_layout);
void invoke_ormqr(integer datatype, char *side, char *trans, integer *m, integer *n, integer *k,
                  void *a, integer *lda, void *tau, void *c, integer *ldc, void *work,
                  integer *lwork, integer *info);
double prepare_lapacke_ormqr_run(integer datatype, int matrix_layout, char side, char trans,
                                 integer m, integer n, integer k, integer m_A, integer n_A, void *A,
                                 integer lda, void *tau, void *c, integer ldc, integer *info);
integer invoke_lapacke_ormqr(integer datatype, int matrix_layout, char side, char trans, integer m,
                             integer n, integer k, void *a, integer lda, const void *tau, void *c,
                             integer ldc);

void fla_test_ormqr(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Orthogonal matrix times matrix multiplied with a QR factorization";
    char *front_str = "ORMQR";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_ormqr_experiment);
        tests_not_run = 0;
    }
    if(argc == 13)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[12]);
    }
    if(argc >= 12 && argc <= 13)
    {
        integer i, num_types, N, M;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].side = argv[3][0];
        params->lin_solver_paramslist[0].transr = argv[4][0];
        M = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].kl = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && params->test_lapacke_interface
           && (params->matrix_major == LAPACK_ROW_MAJOR))
        {
            row_major_ormqr_lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_ormqr_ldc = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = params->lin_solver_paramslist[0].kl;
            params->lin_solver_paramslist[0].ldb = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);

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

                /* Check for duplicate datatype presence */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_ormqr_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for ormqr/unmqr\n");
        printf("./<EXE> ormqr <precisions - sd> <SIDE> <TRANS> <M> <N> <K> <LDA> <LDC> <LWORK> "
               "<repeats>\n");
        printf("./<EXE> unmqr <precisions - cz> <SIDE> <TRANS> <M> <N> <K> <LDA> <LDC> <LWORK> "
               "<repeats>\n");
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

void fla_test_ormqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, k, lda, ldc, m_A, n_A;
    void *A = NULL, *A_test = NULL, *T_test = NULL, *C = NULL, *C_test = NULL;
    void *work = NULL, *qwork = NULL, *tau = NULL;
    integer lwork = -1, info = 0;
    char side, trans;
    double residual, err_thresh;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    m = p_cur;
    n = q_cur;
    k = params->lin_solver_paramslist[0].kl;
    side = params->lin_solver_paramslist[pci].side;
    trans = params->lin_solver_paramslist[pci].transr;
    if((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && trans == 'T')
        trans = 'C';
    lda = params->lin_solver_paramslist[pci].lda;
    ldc = params->lin_solver_paramslist[pci].ldb;

    time_min = 0.;
    perf = 0.;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    if(side == 'L')
    {
        m_A = m;
        n_A = n;
    }
    else
    {
        m_A = n;
        n_A = m;
    }
    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(k < fla_min(m, n))
        {
            k = fla_min(m, n);
        }
        if(lda == -1 || lda < m_A)
        {
            lda = fla_max(1, m_A);
        }
        if(ldc == -1 || ldc < m)
        {
            ldc = fla_max(1, m);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, k, &A, lda);
    init_matrix(datatype, A, m_A, k, lda, g_ext_fptr, params->imatrix_char);

    /* Make a copy of input matrix A to validate the API functionality.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, k, &A_test, lda);
    copy_matrix(datatype, "full", m_A, k, A, lda, A_test, lda);

    /* Create and initialize matrix C */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &C, ldc);
    init_matrix(datatype, C, m, n, ldc, g_ext_fptr, params->imatrix_char);
    /* Scaling matrix with values around overflow, underflow for ORGQR/UNGQR */
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_underflow_overflow_ormqr(datatype, m_A, k, A, lda, params->imatrix_char);
    }
    /* Make a copy of matrix C to validate the API functionality.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &C_test, ldc);
    copy_matrix(datatype, "full", m, n, C, ldc, C_test, ldc);

    /* Create tau vector */
    create_vector(datatype, &T_test, fla_min(m_A, n_A));
    lwork = -1;

    create_vector(datatype, &qwork, 1);
    invoke_geqrf(datatype, &m_A, &k, NULL, &lda, NULL, qwork, &lwork, &info);
    if(info == 0)
    {
        lwork = get_work_value(datatype, qwork);
    }
    else
    {
        lwork = fla_max(1, n_A);
    }

    /* create work buffer */
    create_vector(datatype, &work, lwork);

    /* QR Factorisation on matrix A to generate Q and R */
    invoke_geqrf(datatype, &m_A, &k, A_test, &lda, T_test, work, &lwork, &info);
    create_vector(datatype, &tau, k);
    copy_vector(datatype, fla_min(n_A, k), T_test, 1, tau, 1);
    /*invoke ormqr API */
    prepare_ormqr_run(side, trans, m, n, k, m_A, n_A, A_test, lda, tau, C, ldc, datatype,
                      n_repeats, &time_min, &info, test_lapacke_interface, layout);

    /* performance computation
       perf = 2nk(2m-k) if side = L
            = 2mk(2n-k) if side = R */
    perf = (double)((2.0 * n_A * k) * (2 * m_A - k)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_ormqr(tst_api, side, trans, m, n, k, A_test, lda, C, tau, ldc, C_test, datatype,
                       residual, params->imatrix_char);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if((!check_extreme_value(datatype, m, n, C, ldc, params->imatrix_char))
           && (!check_extreme_value(datatype, m_A, k, A_test, lda, params->imatrix_char)))
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
    free_matrix(A);
    free_matrix(A_test);
    free_vector(work);
    free_vector(T_test);
    free_vector(tau);
    free_matrix(C);
    free_matrix(C_test);
}

void prepare_ormqr_run(char side, char trans, integer m, integer n, integer k, integer m_A,
                       integer n_A, void *A, integer lda, void *tau, void *c, integer ldc,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int layout)
{
    integer i, lwork;
    void *C_save = NULL, *work = NULL;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix C. Same input values will be passed in
       each iteration.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &C_save, ldc);
    copy_matrix(datatype, "full", m, n, c, ldc, C_save, ldc);
    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.
    NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* Getting lwork from api by passing lwork = -1 */
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_ormqr(datatype, &side, &trans, &m, &n, &k, NULL, &lda, tau, NULL, &ldc, work,
                             &lwork, info);
        }
        else
#endif
        {
            invoke_ormqr(datatype, &side, &trans, &m, &n, &k, NULL, &lda, tau, NULL, &ldc, work,
                         &lwork, info);
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

    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix C value for each iteration*/
        copy_matrix(datatype, "full", m, n, C_save, ldc, c, ldc);
        create_vector(datatype, &work, lwork);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_ormqr_run(datatype, layout, side, trans, m, n, k, m_A, n_A,
                                                 A, lda, tau, c, ldc, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP ormqr API */
            invoke_cpp_ormqr(datatype, &side, &trans, &m, &n, &k, A, &lda, tau, c, &ldc, work,
                             &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK ormqr API */
            invoke_ormqr(datatype, &side, &trans, &m, &n, &k, A, &lda, tau, c, &ldc, work, &lwork,
                         info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
        free_vector(work);
    }

    *time_min_ = time_min;
}

double prepare_lapacke_ormqr_run(integer datatype, int layout, char side, char trans, integer m,
                                 integer n, integer k, integer m_A, integer n_A, void *A,
                                 integer lda, void *tau, void *C, integer ldc, integer *info)
{
    double exe_time;
    integer lda_t = lda, ldc_t = ldc;
    void *A_t = NULL, *C_t = NULL;

    if (lda >= m_A)
    {
        /* Configure leading dimensions as per the input matrix layout */
        SELECT_LDA(g_ext_fptr, config_data, layout, fla_max(m, n), row_major_ormqr_lda, lda_t);
        SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_ormqr_ldc, ldc_t);

        A_t = A;
        C_t = C;

        /* In case of row_major matrix layout,
           convert input matrix to row_major */
        if(layout == LAPACK_ROW_MAJOR)
        {
            /* Create temporary buffers for converting matrix layout */
            create_matrix(datatype, layout, m_A, k, &A_t, lda_t);
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m_A, k, A, lda, A_t, lda_t);

            create_matrix(datatype, layout, m, n, &C_t, fla_max(n, ldc_t));
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m, n, C, ldc, C_t, ldc_t);
        }
    }
    exe_time = fla_test_clock();

    /* Call to LAPACKE ormqr API */
    *info
        = invoke_lapacke_ormqr(datatype, layout, side, trans, m, n, k, A_t, lda_t, tau, C_t, ldc_t);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR && *info == 0)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, m, n, C_t, ldc_t, C, ldc);
        /* free temporary buffers */
        free_matrix(A_t);
        free_matrix(C_t);
    }

    return exe_time;
}

void invoke_ormqr(integer datatype, char *side, char *trans, integer *m, integer *n, integer *k,
                  void *a, integer *lda, void *tau, void *c, integer *ldc, void *work,
                  integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
            break;
        }
    }
}

integer invoke_lapacke_ormqr(integer datatype, int layout, char side, char trans, integer m,
                             integer n, integer k, void *a, integer lda, const void *tau, void *c,
                             integer ldc)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sormqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dormqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cunmqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zunmqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc);
            break;
        }
    }
    return info;
}