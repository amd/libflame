#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

integer row_major_ormlq_lda, row_major_ormlq_ldc;
extern double perf;
extern double time_min;

void fla_test_ormlq_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_ormlq_run(integer datatype, char side, char trans, integer m, integer n, integer k,
                       void *A, integer lda, integer m_A, integer n_A, void *tau, void *C,
                       integer ldc, void *work, integer *info, integer interfacetype,
                       integer layout, test_params_t *params);
void invoke_ormlq(integer datatype, char *side, char *trans, integer *m, integer *n, integer *k,
                  void *A, integer *lda, void *tau, void *C, integer *ldc, void *work,
                  integer *lwork, integer *info);
double prepare_lapacke_ormlq_run(integer datatype, integer layout, char side, char trans, integer m,
                                 integer n, integer k, void *A, integer lda, integer m_A,
                                 integer n_A, void *tau, void *C, integer ldc, integer *info);
void fla_test_ormlq(integer argc, char **argv, test_params_t *params)
{
    srand(1.0);
    char *op_str
        = "ORMLQ overwrites the general real M-by-N matrix C with real orthogonal matrix Q";
    char *front_str = "ORMLQ";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_ormlq_experiment);
        tests_not_run = 0;
    }
    if(argc == 13)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[12]);
    }
    if(argc >= 12 && argc <= 13)
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].kl = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].transr = argv[4][0];
        params->lin_solver_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldc = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].side = argv[3][0];

        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_ormlq_lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_ormlq_ldc = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = fla_max(M, N);
            params->lin_solver_paramslist[0].ldc = M;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldc = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

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
                fla_test_ormlq_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }
    /* Print error message */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gecon\n");
        printf("./<EXE> ormlq <precisions - sd> <SIDE> <TRANS> <M> <N> <K> <LDA> <LDC> <LWORK> "
               "<repeats>\n");
        printf("./<EXE> unmlq <precisions - cz> <SIDE> <TRANS> <M> <N> <K> <LDA> <LDC> <LWORK> "
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

void fla_test_ormlq_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, k, lda, ldc, info = 0, m_A, n_A, gelqf_info = 0, gelqf_lwork = -1;
    void *A = NULL, *C = NULL, *work = NULL, *C_save = NULL, *tau = NULL;
    char side, trans;
    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;

    double residual, err_thresh;

    /* Determine the dimensions */
    m = p_cur;
    n = q_cur;
    k = params->lin_solver_paramslist[0].kl;
    lda = params->lin_solver_paramslist[pci].lda;
    ldc = params->lin_solver_paramslist[pci].ldc;
    side = params->lin_solver_paramslist[pci].side;
    trans = params->lin_solver_paramslist[pci].transr;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    time_min = 0.;
    perf = 0.;
    /* If leading dimensions = -1, set them to default value
    when inputs are from config files */
    if(side == 'L' || side == 'l')
    {
        m_A = m;
        n_A = n;
    }
    else
    {
        m_A = n;
        n_A = m;
    }
    if(g_config_data)
    {
        if(lda == -1 || lda < m_A)
        {
            lda = fla_max(m_A, n_A);
        }
        if(ldc == -1 || ldc < m)
        {
            ldc = fla_max(1, m);
        }
    }
    if((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && (trans == 'T'))
    {
        trans = 'C';
    }
    /* Create the matrices for the current operation */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, fla_max(m_A, n_A), &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &C, ldc);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &C_save, ldc);
    create_vector(datatype, &tau, k);

    init_matrix(datatype, A, k, m_A, lda, g_ext_fptr, params->imatrix_char);
    init_matrix(datatype, C, m, n, ldc, g_ext_fptr, params->imatrix_char);

    /* Scaling matrix with values around overflow, underflow for ORMLQ/UNMLQ */
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_underflow_overflow_ormlq(datatype, m_A, k, A, lda, params->imatrix_char);
    }
    /*Generate input matrix A for ormlq, where the i-th row must contain the vector which defines
    the elementary reflector H(i), for i = 1,2,...,k, as returned by GELQF
    */
    if(g_ext_fptr == NULL && !(FLA_EXTREME_CASE_TEST))
    {
        create_vector(datatype, &work, i_one);
        invoke_gelqf(datatype, &k, &m_A, NULL, &lda, NULL, work, &gelqf_lwork, &gelqf_info);
        if(gelqf_info == 0)
        {
            gelqf_lwork = get_work_value(datatype, work);
        }
        else
        {
            gelqf_lwork = fla_max(m_A, 1);
        }
        free_vector(work);
        create_vector(datatype, &work, gelqf_lwork);
        gelqf_info = 0;
        invoke_gelqf(datatype, &k, &m_A, A, &lda, tau, work, &gelqf_lwork, &gelqf_info);
        free_vector(work);
    }
    /* Save the original matrix */
    copy_matrix(datatype, "Full", m, n, C, ldc, C_save, ldc);

    /* call to API */
    prepare_ormlq_run(datatype, side, trans, m, n, k, A, lda, m_A, n_A, tau, C, ldc, work, &info,
                      interfacetype, layout, params);

    /* performance computation
       perf = 2nk(2m-k) if side = L
            = 2mk(2n-k) if side = R */
    perf = (double)(2.0 * (m_A * k) * ((2.0 * n_A) - k)) / time_min / FLOPS_PER_UNIT_PERF;

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    /* Output validataion */
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_ormqr(tst_api, side, trans, m, n, k, A, lda, C, tau, ldc, C_save, datatype,
                       residual, params->imatrix_char, params);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if(!check_extreme_value(datatype, m, n, C, ldc, params->imatrix_char))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    }

    /* Free up buffers */
    free_matrix(A);
    free_matrix(C);
    free_matrix(C_save);
    free_vector(tau);
}

void prepare_ormlq_run(integer datatype, char side, char trans, integer m, integer n, integer k,
                       void *A, integer lda, integer m_A, integer n_A, void *tau, void *C,
                       integer ldc, void *work, integer *info, integer interfacetype,
                       integer layout, test_params_t *params)
{
    integer lwork;
    void *C_save = NULL;
    double exe_time;
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &C_save, ldc);

    /* Workspace size calculations for complex datatype */
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
            invoke_cpp_ormlq(datatype, &side, &trans, &m, &n, &k, NULL, &lda, tau, NULL, &ldc, work,
                             &lwork, info);
        }
        else
#endif
        {
            invoke_ormlq(datatype, &side, &trans, &m, &n, &k, NULL, &lda, tau, NULL, &ldc, work,
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
    FLA_EXEC_LOOP_BEGIN
    {
        /* Copy original input */
        copy_matrix(datatype, "full", m, n, C, ldc, C_save, ldc);

        /* Create work buffer */
        create_vector(datatype, &work, lwork);
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            /*call to lapacke gecon*/
            exe_time = prepare_lapacke_ormlq_run(datatype, layout, side, trans, m, n, k, A, lda,
                                                 m_A, n_A, tau, C_save, ldc, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP gecon API */
            invoke_cpp_ormlq(datatype, &side, &trans, &m, &n, &k, A, &lda, tau, C_save, &ldc, work,
                             &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();

            /*  call to ormlq API */
            invoke_ormlq(datatype, &side, &trans, &m, &n, &k, A, &lda, tau, C_save, &ldc, work,
                         &lwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO

        free_vector(work);
    }

    /* Save the output to vector C */
    copy_matrix(datatype, "full", m, n, C_save, ldc, C, ldc);

    /* Free up buffers */
    free_matrix(C_save);
}

double prepare_lapacke_ormlq_run(integer datatype, integer layout, char side, char trans, integer m,
                                 integer n, integer k, void *A, integer lda, integer m_A,
                                 integer n_A, void *tau, void *C, integer ldc, integer *info)
{
    double exe_time = 0;
    void *A_t = NULL, *C_t = NULL;
    integer lda_t = lda, ldc_t = ldc;

    A_t = A;
    C_t = C;
    /* Configure leading dimensions as per the input matrix layout */

    SELECT_LDA(g_ext_fptr, g_config_data, layout, fla_max(m, n), row_major_ormlq_lda, lda_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, fla_max(m, n), row_major_ormlq_ldc, ldc_t);
    /* In case of row_major matrix layout,
    convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m_A, k, &A_t, lda_t);
        create_matrix(datatype, layout, fla_max(m, n), n, &C_t, fla_max(n, ldc_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, k, m_A, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m, n, C, ldc, C_t, ldc_t);
    }
    exe_time = fla_test_clock();

    /* call to LAPACKE ormlq API */
    *info
        = invoke_lapacke_ormlq(datatype, layout, side, trans, m, n, k, A_t, lda_t, tau, C_t, ldc_t);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
    to column_major layout */

    if((layout == LAPACK_ROW_MAJOR))
    {
        convert_matrix_layout(layout, datatype, m, n, C_t, ldc_t, C, ldc);

        free_matrix(A_t);
        free_matrix(C_t);
    }
    return exe_time;
}

void invoke_ormlq(integer datatype, char *side, char *trans, integer *m, integer *n, integer *k,
                  void *A, integer *lda, void *tau, void *C, integer *ldc, void *work,
                  integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sormlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dormlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cunmlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zunmlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info);
            break;
        }
    }
}
