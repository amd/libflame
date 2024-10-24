/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"

integer row_major_gesdd_lda;
integer row_major_gesdd_ldu;
integer row_major_gesdd_ldvt;

/* Local prototypes.*/
void fla_test_gesdd_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_gesdd_run(char *jobz, integer m_A, integer n_A, void *A, integer lda, void *s, void *U,
                       integer ldu, void *V, integer ldvt, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer test_lapacke_interface,
                       int matrix_layout);
void invoke_gesdd(integer datatype, char *jobz, integer *m, integer *n, void *a, integer *lda,
                  void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                  integer *lwork, void *rwork, integer *iwork, integer *info);
double prepare_lapacke_gesdd_run(integer datatype, int matrix_layout, char *jobz, integer m_A,
                                 integer n_A, void *A, integer lda, void *s, void *U, integer ldu,
                                 void *V, integer ldvt, integer *info);
integer invoke_lapacke_gesdd(integer datatype, int matrix_layout, char jobz, integer m, integer n,
                             void *a, integer lda, void *s, void *u, integer ldu, void *vt,
                             integer ldvt);

void fla_test_gesdd(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Singular value decomposition";
    char *front_str = "GESDD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, SVD, fla_test_gesdd_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        integer i, num_types, N, M;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->svd_paramslist[0].jobu_gesvd = argv[3][0];
        M = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && params->test_lapacke_interface
           && (params->matrix_major == LAPACK_ROW_MAJOR))
        {
            row_major_gesdd_lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            row_major_gesdd_ldu = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            row_major_gesdd_ldvt = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].lda = N;
            params->svd_paramslist[0].ldu = N;
            params->svd_paramslist[0].ldvt = N;
        }
        else
        {
            params->svd_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldu = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldvt = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        }

        g_lwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->svd_paramslist[0].svd_threshold = CLI_NORM_THRESH;

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
                fla_test_gesdd_experiment(params, datatype, M, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, RECT_INPUT, M, N, residual,
                                      params->svd_paramslist[0].svd_threshold, time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gesdd\n");
        printf("./<EXE> gesdd <precisions - sdcz> <JOBZ> <M> <N> <LDA> <LDU> <LDVT> <LWORK> "
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

void fla_test_gesdd_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *time_min, double *residual)
{
    integer m, n, lda, ldu, ldvt, n_U, m_V, ns;
    integer info = 0, vinfo = 0;
    char jobz;
    void *A = NULL, *U = NULL, *V = NULL, *s = NULL, *A_test = NULL, *s_in = NULL, *scal = NULL;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    jobz = params->svd_paramslist[pci].jobu_gesvd;
    *residual = params->svd_paramslist[pci].svd_threshold;

    m = p_cur;
    n = q_cur;
    ns = fla_min(m, n);
    lda = params->svd_paramslist[pci].lda;
    ldu = params->svd_paramslist[pci].ldu;
    ldvt = params->svd_paramslist[pci].ldvt;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
        /* LDU >= 1;
           if JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M. */
        if(ldu == -1)
        {
            if(((jobz == 'S') || (jobz == 'A')) || ((jobz == 'O') && (m < n)))
            {
                ldu = m;
            }
            else
            {
                ldu = 1;
            }
        }
        /* LDVT >= 1;
           if JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
           if JOBZ = 'S', LDVT >= min(M,N). */
        if(ldvt == -1)
        {
            if((jobz == 'A') || ((jobz == 'O') && (m >= n)))
            {
                ldvt = n;
            }
            else if(jobz == 'S')
            {
                ldvt = ns;
            }
            else
            {
                ldvt = 1;
            }
        }
    }

    n_U = (jobz == 'S') ? ns : m;
    m_V = (jobz == 'S') ? ns : n;

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    if (jobz != 'N')
    {
        if (jobz == 'A' || jobz == 'S')
        {
            create_matrix(datatype, LAPACK_COL_MAJOR, m, n_U, &U, ldu);
            create_matrix(datatype, LAPACK_COL_MAJOR, m_V, n, &V, ldvt);
        }
        else if (jobz == 'O' && m >= n)
        {
            create_matrix(datatype, LAPACK_COL_MAJOR, m_V, n, &V, ldvt);
        }
        else
        {
            create_matrix(datatype, LAPACK_COL_MAJOR, m, n_U, &U, ldu);
        }
    }

    create_realtype_vector(datatype, &s, fla_min(m, n));

    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        create_realtype_vector(datatype, &s_in, fla_min(m, n));
        /* Generate matrix A from known singular values */
        create_svd_matrix(datatype, 'A', m, n, A, lda, s_in, s_one, s_one, i_one, i_one, info);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_vector(get_realtype(datatype), &scal, 1);
            scale_matrix_underflow_overflow_gesdd(datatype, m, n, A, lda, params->imatrix_char,
                                                  scal);
        }
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_gesdd_run(&jobz, m, n, A_test, lda, s, U, ldu, V, ldvt, datatype, n_repeats, time_min,
                      &info, test_lapacke_interface, layout);

    /* performance computation
       6mn^2 + 8n^3 flops */
    if(m >= n)
        *perf = (double)((6.0 * m * n * n) + (8.0 * n * n * n)) / *time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)((6.0 * n * m * m) + ((8.0) * m * m * m)) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if(info == 0 && (FLA_OVERFLOW_UNDERFLOW_TEST || (!FLA_EXTREME_CASE_TEST)))
        validate_gesdd(&jobz, m, n, A, A_test, lda, s, s_in, U, ldu, V, ldvt, datatype, residual,
                       &vinfo, params->imatrix_char, scal);
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if(!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    if(jobz == 'A' || jobz == 'S' || (jobz == 'O' && m < n))
    {
        free_matrix(U);
    }
    if(jobz == 'A' || jobz == 'S' || (jobz == 'O' && m >= n))
    {
        free_matrix(V);
    }
    free_vector(s);
    if((g_ext_fptr == NULL) && !(params->imatrix_char))
    {
        free_vector(s_in);
    }
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal);
    }
}

void prepare_gesdd_run(char *jobz, integer m_A, integer n_A, void *A, integer lda, void *s, void *U,
                       integer ldu, void *V, integer ldvt, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer test_lapacke_interface, int layout)
{
    integer min_m_n, max_m_n, n_U, m_V;
    void *A_save, *s_test, *work, *iwork, *rwork;
    void *U_test, *V_test;
    integer lwork, liwork, lrwork;
    integer i;
    double time_min = 1e9, exe_time;

    min_m_n = fla_min(m_A, n_A);
    max_m_n = fla_max(m_A, n_A);
    n_U = (*jobz == 'S') ? min_m_n : m_A;
    m_V = (*jobz == 'S') ? min_m_n : n_A;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = fla_max((5 * min_m_n * min_m_n + 5 * min_m_n),
                     (2 * max_m_n * min_m_n + 2 * min_m_n * min_m_n + min_m_n));
    liwork = 8 * min_m_n;

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.
     NOTE: LAPACKE interface handles workspace query internally */
    if((test_lapacke_interface == 0) && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call to  gesdd API */
        invoke_gesdd(datatype, jobz, &m_A, &n_A, NULL, &lda, NULL, NULL, &ldu, NULL, &ldvt, work,
                     &lwork, NULL, NULL, info);
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
        if (*jobz != 'N')
        {
            if (*jobz == 'A' || *jobz == 'S')
            {
                create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_U, &U_test, ldu);
                create_matrix(datatype, LAPACK_COL_MAJOR, m_V, n_A, &V_test, ldvt);
            }
            else if (*jobz == 'O' && m_A >= n_A)
            {
                create_matrix(datatype, LAPACK_COL_MAJOR, m_V, n_A, &V_test, ldvt);
            }
            else
            {
                create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_U, &U_test, ldu);
            }
        }
        create_realtype_vector(datatype, &s_test, min_m_n);
        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, liwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, lrwork);
        else
            rwork = NULL;

        /* Check if LAPACKE interface is enabled */
        if(test_lapacke_interface == 1)
        {
            exe_time = prepare_lapacke_gesdd_run(datatype, layout, jobz, m_A, n_A, A, lda, s_test,
                                                 U_test, ldu, V_test, ldvt, info);
        }
        else
        {
            exe_time = fla_test_clock();
            /* call to API */
            invoke_gesdd(datatype, jobz, &m_A, &n_A, A, &lda, s_test, U_test, &ldu, V_test, &ldvt,
                         work, &lwork, rwork, iwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Make a copy of the output buffers. This is required to validate the API functionality.*/
        if (*jobz != 'N')
        {
            if (*jobz == 'A' || *jobz == 'S')
            {
                copy_matrix(datatype, "full", m_A, n_U, U_test, ldu, U, ldu);
                copy_matrix(datatype, "full", m_V, n_A, V_test, ldvt, V, ldvt);
            }
            else if (*jobz == 'O' && m_A >= n_A)
            {
                copy_matrix(datatype, "full", m_V, n_A, V_test, ldvt, V, ldvt);
            }
            else
            {
                copy_matrix(datatype, "full", m_A, n_U, U_test, ldu, U, ldu);
            }
        }
        copy_realtype_vector(datatype, min_m_n, s_test, 1, s, 1);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
        if(*jobz == 'A' || *jobz == 'S' || (*jobz == 'O' && m_A < n_A))
        {
            free_matrix(U_test);
        }
        if(*jobz == 'A' || *jobz == 'S' || (*jobz == 'O' && m_A >= n_A))
        {
            free_matrix(V_test);
        }
        free_vector(s_test);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
}

double prepare_lapacke_gesdd_run(integer datatype, int layout, char *jobz, integer m_A, integer n_A,
                                 void *A, integer lda, void *s, void *U, integer ldu, void *V,
                                 integer ldvt, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldu_t = ldu;
    integer ldvt_t = ldvt;
    void *A_t = NULL, *U_t = NULL, *V_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_gesdd_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, m_A, row_major_gesdd_ldu, ldu_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_gesdd_ldvt, ldvt_t);

    A_t = A;
    U_t = U;
    V_t = V;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m_A, n_A, &A_t, fla_max(n_A, lda_t));
        switch(*jobz)
        {
            case 'A':
            {
                create_matrix(datatype, layout, m_A, m_A, &U_t, fla_max(m_A, ldu_t));
                create_matrix(datatype, layout, n_A, n_A, &V_t, fla_max(n_A, ldvt_t));
                break;
            }
            case 'S':
            {
                create_matrix(datatype, layout, m_A, fla_min(m_A, n_A), &U_t, fla_max(m_A, ldu_t));
                create_matrix(datatype, layout, fla_min(m_A, n_A), n_A, &V_t, fla_max(n_A, ldvt_t));
                break;
            }
            case 'O':
            {
                if(m_A < n_A)
                {
                    create_matrix(datatype, layout, m_A, m_A, &U_t, fla_max(m_A, ldu_t));
                }
                else
                {
                    create_matrix(datatype, layout, n_A, n_A, &V_t, fla_max(n_A, ldvt_t));
                }
                break;
            }
            default:
            {
                break;
            }
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m_A, n_A, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* call LAPACKE gesdd API */
    *info = invoke_lapacke_gesdd(datatype, layout, *jobz, m_A, n_A, A_t, lda_t, s, U_t, ldu_t, V_t,
                                 ldvt_t);

    exe_time = fla_test_clock() - exe_time;
    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */
    if(layout == LAPACK_ROW_MAJOR)
    {
        switch(*jobz)
        {
            case 'A':
            {
                convert_matrix_layout(layout, datatype, m_A, m_A, U_t, ldu_t, U, ldu);
                convert_matrix_layout(layout, datatype, n_A, n_A, V_t, ldvt_t, V, ldvt);
                free_matrix(U_t);
                free_matrix(V_t);
                break;
            }
            case 'S':
            {
                convert_matrix_layout(layout, datatype, m_A, fla_min(m_A, n_A), U_t, ldu_t, U, ldu);
                convert_matrix_layout(layout, datatype, fla_min(m_A, n_A), n_A, V_t, ldvt_t, V,
                                      ldvt);
                free_matrix(U_t);
                free_matrix(V_t);
                break;
            }
            case 'O':
            {
                convert_matrix_layout(layout, datatype, m_A, n_A, A_t, lda_t, A, lda);
                if(m_A < n_A)
                {
                    convert_matrix_layout(layout, datatype, m_A, m_A, U_t, fla_max(m_A, ldu_t), U,
                                          ldu);
                    free_matrix(U_t);
                }
                else
                {
                    convert_matrix_layout(layout, datatype, n_A, n_A, V_t, fla_max(n_A, ldvt_t), V,
                                          ldvt);
                    free_matrix(V_t);
                }
                break;
            }
            default:
                break;
        }
        free_matrix(A_t);
    }
    return exe_time;
}

void invoke_gesdd(integer datatype, char *jobz, integer *m, integer *n, void *a, integer *lda,
                  void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                  integer *lwork, void *rwork, integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork,
                              info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork,
                              info);
            break;
        }
    }
}

integer invoke_lapacke_gesdd(integer datatype, int layout, char jobz, integer m, integer n, void *a,
                             integer lda, void *s, void *u, integer ldu, void *vt, integer ldvt)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgesdd(layout, jobz, m, n, a, lda, s, u, ldu, vt, ldvt);
            break;
        }
    }
    return info;
}