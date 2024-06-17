/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#define GELSD_VL 0.1
#define GELSD_VU 10

/* Local prototypes */
void fla_test_gelsd_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_gelsd_run(integer m_A, integer n_A, integer nrhs, void *A, integer lda, void *B,
                       integer ldb, void *s, void *rcond, integer *rank, integer datatype,
                       integer n_repeats, double *time_min_, integer *info);
void invoke_gelsd(integer datatype, integer *m, integer *n, integer *nrhs, void *a, integer *lda,
                  void *b, integer *ldb, void *s, void *rcond, integer *rank, void *work,
                  integer *lwork, void *rwork, integer *iwork, integer *info);

void fla_test_gelsd(integer argc, char **argv, test_params_t *params)
{
    srand(1);
    char *op_str = "Linear least squares";
    char *front_str = "GELSD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gelsd_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        /* Test with parameters from commandline */
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

        params->lin_solver_paramslist[0].rcond = atof(argv[8]);
        g_lwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_gelsd_experiment(params, datatype, M, N, 0, n_repeats, einfo, &perf,
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
        printf("\nIllegal arguments for gelsd\n");
        printf("./<EXE> gelsd <precisions - sdcz>  <M> <N> <NRHS> <LDA>"
               " <LDB> <RCOND> <LWORK> <repeats>\n");
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

void fla_test_gelsd_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer m, n, lda, ldb, NRHS;
    integer info = 0, rank;
    void *A = NULL, *A_save = NULL, *B = NULL, *B_save = NULL;
    void *S = NULL, *rcond = NULL, *s_test = NULL;
    double time_min = 1e9;
    char range = 'U';
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    NRHS = params->lin_solver_paramslist[pci].nrhs;
    m = p_cur;
    n = q_cur;

    rank = fla_min(m, n);
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;

    create_realtype_vector(datatype, &rcond, 1);

    if(datatype == FLOAT || datatype == COMPLEX)
    {
        *(real *)rcond = params->lin_solver_paramslist[pci].rcond;
    }
    else
    {
        *(doublereal *)rcond = params->lin_solver_paramslist[pci].rcond;
    }

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
            ldb = fla_max(1, fla_max(m, n));
        }
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &A_save, lda, n);
    create_matrix(datatype, &B, ldb, NRHS);
    create_matrix(datatype, &B_save, ldb, NRHS);
    create_realtype_vector(datatype, &s_test, fla_min(m, n));

    /* Initialize the test matrices */
    init_matrix(datatype, B, m, NRHS, ldb, g_ext_fptr, params->imatrix_char);
    if(params->imatrix_char == NULL && g_ext_fptr == NULL)
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, m, n, A, lda, s_test, GELSD_VL, GELSD_VU, i_zero, i_zero,
                          '\0', NULL, info);
    }
    else
    {
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    }

    /* Save the original matrix*/
    copy_matrix(datatype, "full", m, n, A, lda, A_save, lda);
    copy_matrix(datatype, "full", m, NRHS, B, ldb, B_save, ldb);
    create_realtype_vector(datatype, &S, fla_min(m, n));

    /* call to API */
    prepare_gelsd_run(m, n, NRHS, A_save, lda, B_save, ldb, S, rcond, &rank, datatype, n_repeats,
                      &time_min, &info);
    /* execution time */
    *t = time_min;

    /* performance computation */
    if(m >= n)
    {
        *perf = (double)(4.0 * m * n * (n + NRHS) + NRHS) / *t / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        *perf = (double)((2.0 * m * NRHS * (m + n) + m * (4.0 * n * n + 1.0) + NRHS)) / *t
                / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Output validation, accuracy test */
    if(!params->imatrix_char && info == 0)
        validate_gelsd(m, n, NRHS, A, lda, B, ldb, S, B_save, rcond, &rank, datatype, residual);
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, m, n, A_save, lda, params->imatrix_char))
           && (!check_extreme_value(datatype, m, n, B_save, ldb, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_save);
    free_matrix(B);
    free_matrix(B_save);
    free_vector(S);
    free_vector(rcond);
    free_vector(s_test);
}

void prepare_gelsd_run(integer m_A, integer n_A, integer nrhs, void *A, integer lda, void *B,
                       integer ldb, void *s, void *rcond, integer *rank, integer datatype,
                       integer n_repeats, double *time_min_, integer *info)
{
    integer i, lwork, liwork = 1, lrwork = 1, realtype;
    void *A_test, *B_test;
    void *work = NULL, *rwork = NULL, *iwork = NULL;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, &A_test, lda, n_A);
    create_matrix(datatype, &B_test, ldb, nrhs);
    realtype = get_realtype(datatype);

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.*/
    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        create_vector(realtype, &rwork, 1);
        create_vector(INTEGER, &iwork, 1);

        /* call to  gelsd API */
        invoke_gelsd(datatype, &m_A, &n_A, &nrhs, NULL, &lda, NULL, &ldb, NULL, rcond, rank, work,
                     &lwork, rwork, iwork, info);
        /* Get work size */
        if(*info == 0)
        {
            lwork = get_work_value(datatype, work);
            liwork = get_work_value(INTEGER, iwork);
            lrwork = get_work_value(realtype, rwork);
        }
        free_vector(work);
        free_vector(iwork);
        free_vector(rwork);
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
        /* Copy original input data */
        copy_matrix(datatype, "full", m_A, n_A, A, lda, A_test, lda);
        copy_matrix(datatype, "full", m_A, nrhs, B, ldb, B_test, ldb);

        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, fla_max(1, liwork));

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, fla_max(1, lrwork));
        else
            rwork = NULL;

        exe_time = fla_test_clock();

        /* call to API */
        invoke_gelsd(datatype, &m_A, &n_A, &nrhs, A_test, &lda, B_test, &ldb, s, rcond, rank, work,
                     &lwork, rwork, iwork, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);
    }

    *time_min_ = time_min;

    /*  Save the final result to B matrix*/
    copy_matrix(datatype, "full", n_A, nrhs, B_test, ldb, B, ldb);
    free_matrix(A_test);
    free_matrix(B_test);
}

/*
 *  Call to LAPACK interface of Linear least square - gelsd
 *  */
void invoke_gelsd(integer datatype, integer *m, integer *n, integer *nrhs, void *a, integer *lda,
                  void *b, integer *ldb, void *s, void *rcond, integer *rank, void *work,
                  integer *lwork, void *rwork, integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork,
                              info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork,
                              info);
            break;
        }
    }
}
