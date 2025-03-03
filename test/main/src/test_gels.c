/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

#define GELS_VL 0.1
#define GELS_VU 10

extern double perf;
extern double time_min;
integer row_major_gels_lda;
integer row_major_gels_ldb;

void invoke_gels(integer datatype, char *trans, integer *m, integer *n, integer *nrhs, void *A,
                 integer *lda, void *B, integer *ldb, void *work, integer *lwork, integer *info);
void fla_test_gels_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo);
void prepare_gels_run(integer datatype, char trans, integer m, integer n, integer m_b, integer nrhs,
                      void *A, integer lda, void *B, integer ldb, void *work, integer lwork,
                      integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                      integer layout);
double prepare_lapacke_gels_run(integer datatype, integer layout, char trans, integer m, integer n,
                                integer nrhs, integer m_b, void *A, integer lda, void *B,
                                integer ldb, integer *info);

void fla_test_gels(integer argc, char **argv, test_params_t *params)
{
    srand(14);
    char *op_str = "Solves overdetermined or underdetermined systems for GE matrices";
    char *front_str = "GELS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gels_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].transr = argv[3][0];
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gels_lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            row_major_gels_ldb = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = M;
            params->lin_solver_paramslist[0].ldb = fla_max(M, N);
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldb = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_gels_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }
    /* Print error message */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gels\n");
        printf("./<EXE> gels <precisions - sdcz> <TRANS> <M> <N> <NRHS> <LDA> <LDB> <LWORK> "
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

void fla_test_gels_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo)
{
    integer m, n, m_b, nrhs, lda, ldb, lwork = -1, info = 0;
    void *A = NULL, *A_test = NULL, *B = NULL, *B_test = NULL, *work = NULL, *s_test = NULL;
    char trans, range = 'U';
    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;
    double residual, err_thresh;

    /* Determine the dimensions */
    m = p_cur;
    n = q_cur;
    nrhs = params->lin_solver_paramslist[pci].nrhs;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    trans = params->lin_solver_paramslist[pci].transr;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
        if(ldb == -1)
        {
            ldb = fla_max(fla_max(1, m), n);
        }
    }

    /* Based on the value of TRANS the dimension of B changes.
     * Dimension of B is (m, nrhs) if TRANS = "N"
     * Dimension of B is (n, nrhs) if TRANS = "T" */

    m_b = n;
    if(same_char(trans, 'N'))
    {
        trans = 'N';
        m_b = m;
    }

    /* trans for complex number should be equal to 'C' (or 'c') while passing to the GEL api
     */
    if((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && (same_char(trans, 'T')))
    {
        trans = 'C';
    }

    /* Create the matrices for the current operation */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m_b, nrhs, &B, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, m_b, nrhs, &B_test, ldb);
    reset_matrix(datatype, ldb, nrhs, B, ldb);
    reset_matrix(datatype, ldb, nrhs, B_test, ldb);
    create_realtype_vector(datatype, &s_test, fla_min(m, n));

    /* Initialize the test matrices */
    init_matrix(datatype, B, m_b, nrhs, ldb, g_ext_fptr, params->imatrix_char);

    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, m, n, A, lda, s_test, GELS_VL, GELS_VU, i_zero, i_zero,
                          info);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_gels(datatype, &trans, m, n, A, lda,
                                                 params->imatrix_char, 1);
        }
    }

    /* Save the original matrix */
    copy_matrix(datatype, "full", lda, n, A, lda, A_test, lda);
    copy_matrix(datatype, "full", ldb, nrhs, B, ldb, B_test, ldb);

    /* call to API */
    prepare_gels_run(datatype, trans, m, n, m_b, nrhs, A_test, lda, B_test, ldb, work, lwork,
                     n_repeats, &time_min, &info, interfacetype, layout);

    /* Performance computation */
    if(m >= n)
    {
        perf = (double)((n * n) * (2.0 / 3.0) * ((3 * m) - n)) / time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        perf = (double)((m * m) * (2.0 / 3.0) * ((3 * n) - m)) / time_min / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* Output validataion */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_gels(tst_api, &trans, m, n, nrhs, A, lda, B, ldb, B_test, datatype, residual,
                      params->imatrix_char);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if((!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char))
           && (!check_extreme_value(datatype, m, n, B_test, ldb, params->imatrix_char)))
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
    free_matrix(A_test);
    free_matrix(B);
    free_matrix(B_test);
    free_vector(s_test);
}

void prepare_gels_run(integer datatype, char trans, integer m, integer n, integer m_b, integer nrhs,
                      void *A, integer lda, void *B, integer ldb, void *work, integer lwork,
                      integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                      integer layout)
{
    integer i;
    void *A_save = NULL, *B_save = NULL;
    double t_min = 1e9, exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_save, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m_b, nrhs, &B_save, ldb);

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
            invoke_cpp_gels(datatype, &trans, &m, &n, &nrhs, NULL, &lda, NULL, &ldb, work, &lwork,
                            info);
        }
        else
#endif
        {
            invoke_gels(datatype, &trans, &m, &n, &nrhs, NULL, &lda, NULL, &ldb, work, &lwork,
                        info);
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
    for(i = 0; i < n_repeats && *info == 0; i++)
    {
        /* Copy original input */
        copy_matrix(datatype, "full", lda, n, A, lda, A_save, lda);
        copy_matrix(datatype, "full", ldb, nrhs, B, ldb, B_save, ldb);

        /* Create work buffer */
        create_vector(datatype, &work, lwork);
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gels_run(datatype, layout, trans, m, n, nrhs, m_b, A_save,
                                                lda, B_save, ldb, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP gels API */
            invoke_cpp_gels(datatype, &trans, &m, &n, &nrhs, A_save, &lda, B_save, &ldb, work,
                            &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /*  call to LAPACK gels API */
            invoke_gels(datatype, &trans, &m, &n, &nrhs, A_save, &lda, B_save, &ldb, work, &lwork,
                        info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        free_vector(work);
    }
    *time_min_ = t_min;

    /* Save the output to vector A */
    copy_matrix(datatype, "full", lda, n, A_save, lda, A, lda);
    copy_matrix(datatype, "full", ldb, nrhs, B_save, ldb, B, ldb);

    /* Free up buffers */
    free_matrix(A_save);
    free_matrix(B_save);
}

double prepare_lapacke_gels_run(integer datatype, integer layout, char trans, integer m, integer n,
                                integer nrhs, integer m_b, void *A, integer lda, void *B,
                                integer ldb, integer *info)
{
    double exe_time = 0;
    void *A_t = NULL, *B_t = NULL;
    integer lda_t = lda, ldb_t = ldb, m_x = m;
    A_t = A;
    B_t = B;

    if(same_char(trans, 'N'))
    {
        m_x = n;
    }

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_gels_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, nrhs, row_major_gels_ldb, ldb_t);

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m, n, &A_t, fla_max(n, lda_t));
        create_matrix(datatype, layout, fla_max(m, n), nrhs, &B_t, fla_max(nrhs, ldb_t));
        reset_matrix(datatype, fla_max(m, n), nrhs, B_t, fla_max(m, n));

        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m, n, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m_b, nrhs, B, ldb, B_t, ldb_t);
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE gels API */
    *info = invoke_lapacke_gels(datatype, layout, trans, m, n, nrhs, A_t, lda_t, B_t, ldb_t);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */

    if((layout == LAPACK_ROW_MAJOR))
    {
        convert_matrix_layout(layout, datatype, m, n, A_t, lda_t, A, lda);
        convert_matrix_layout(layout, datatype, m_x, nrhs, B_t, ldb_t, B, ldb);

        free_matrix(A_t);
        free_matrix(B_t);
    }
    return exe_time;
}

/*
LAPACK GELS API invoke function
*/
void invoke_gels(integer datatype, char *trans, integer *m, integer *n, integer *nrhs, void *A,
                 integer *lda, void *B, integer *ldb, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
    }
}