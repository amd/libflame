/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

integer row_major_ggev_lda;
integer row_major_ggev_ldb;
integer row_major_ggev_ldvl;
integer row_major_ggev_ldvr;

/* Local prototypes */
void fla_test_ggev_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual);
void prepare_ggev_run(char *jobvl, char *jobvr, integer n, void *a, integer lda, void *b,
                      integer ldb, void *alpha, void *alphar, void *alphai, void *beta, void *vl,
                      integer ldvl, void *vr, integer ldvr, integer datatype, integer n_repeats,
                      double *time_min_, integer *info, integer interfacetype,
                      int matrix_layout);
void invoke_ggev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                 void *b, integer *ldb, void *alpha, void *alphar, void *alphai, void *beta,
                 void *vl, integer *ldvl, void *vr, integer *ldvr, void *work, integer *lwork,
                 void *rwork, integer *info);
double prepare_lapacke_ggev_run(integer datatype, int matrix_layout, char *jobvl, char *jobvr,
                                integer n, void *a, integer lda, void *b, integer ldb, void *alpha,
                                void *alphar, void *alphai, void *beta, void *vl, integer ldvl,
                                void *vr, integer ldvr, integer *info);
integer invoke_lapacke_ggev(integer datatype, int matrix_layout, char jobvl, char jobvr,
                            integer n, void *a, integer lda, void *b, integer ldb, void *alpha,
                            void *alphar, void *alphai, void *beta, void *vl, integer ldvl,
                            void *vr, integer ldvr);

void fla_test_ggev(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Computing Eigen value and Eigen vectors";
    char *front_str = "GGEV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_ggev_experiment);
        tests_not_run = 0;
    }
    if(argc == 13)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[12]);
    }
    if(argc >= 12 && argc <= 13)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_non_sym_paramslist[0].jobvsl = argv[3][0];
        params->eig_non_sym_paramslist[0].jobvsr = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_ggev_lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            row_major_ggev_ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            row_major_ggev_ldvl = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_ggev_ldvr = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].lda = N;
            params->eig_non_sym_paramslist[0].ldb = N;
            params->eig_non_sym_paramslist[0].ldvl = N;
            params->eig_non_sym_paramslist[0].ldvr = N;
        }
        else
        {
            params->eig_non_sym_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldb = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->eig_non_sym_paramslist[0].GenNonSymEigProblem_threshold = CLI_NORM_THRESH;

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
                fla_test_ggev_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
                                         &time_min, &residual);
                /* Print the results */
                fla_test_print_status(
                    front_str, stype, SQUARE_INPUT, N, N, residual,
                    params->eig_non_sym_paramslist[0].GenNonSymEigProblem_threshold, time_min,
                    perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for ggev\n");
        printf("./<EXE> ggev <precisions - sdcz> <jobvl> <jobvr> <N> <LDA> <LDB> <LDVL> <LDVR> "
               "<LWORK> <repeats>\n");
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

void fla_test_ggev_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual)
{
    integer m, lda, ldvl, ldvr, ldb;
    integer info = 0, vinfo = 0;
    void *A = NULL, *B = NULL, *VL = NULL, *VR = NULL;
    void *alpha = NULL, *alphar = NULL, *alphai = NULL, *beta, *A_test, *B_test;
    double time_min = 1e9;
    *residual = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    char JOBVL = params->eig_non_sym_paramslist[pci].jobvsl;
    char JOBVR = params->eig_non_sym_paramslist[pci].jobvsr;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions */
    m = p_cur;

    lda = params->eig_non_sym_paramslist[pci].lda;
    ldb = params->eig_non_sym_paramslist[pci].ldb;
    ldvl = params->eig_non_sym_paramslist[pci].ldvl;
    ldvr = params->eig_non_sym_paramslist[pci].ldvr;

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
            ldb = fla_max(1, m);
        }
        /* LDVL >= 1, and
           if JOBVL = 'V', LDVL >= M */
        if(ldvl == -1)
        {
            if(JOBVL == 'V')
            {
                ldvl = m;
            }
            else
            {
                ldvl = 1;
            }
        }
        /* LDVR >= 1, and
           if JOBVR = 'V', LDVR >= M */
        if(ldvr == -1)
        {
            if(JOBVR == 'V')
            {
                ldvr = m;
            }
            else
            {
                ldvr = 1;
            }
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &B, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &VL, ldvl);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &VR, ldvr);
    if(datatype == FLOAT || datatype == DOUBLE)
    {
        create_vector(datatype, &alphar, m);
        create_vector(datatype, &alphai, m);
    }
    else
    {
        create_vector(datatype, &alpha, m);
    }
    create_vector(datatype, &beta, m);

    init_matrix(datatype, A, m, m, lda, g_ext_fptr, params->imatrix_char);
    init_matrix(datatype, B, m, m, lda, g_ext_fptr, params->imatrix_char);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_overflow_underflow_ggev(datatype, m, A, lda, params->imatrix_char);
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_test, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &B_test, ldb);
    copy_matrix(datatype, "full", m, m, A, lda, A_test, lda);
    copy_matrix(datatype, "full", m, m, B, ldb, B_test, ldb);

    prepare_ggev_run(&JOBVL, &JOBVR, m, A_test, lda, B_test, ldb, alpha, alphar, alphai, beta, VL,
                     ldvl, VR, ldvr, datatype, n_repeats, &time_min, &info, interfacetype,
                     layout);

    /* execution time */
    *t = time_min;

    /* performance computation */
    /* 2m^3 - (2/3)m^3 flops */
    *perf
        = (double)((2.0 * m * m * m) - ((2.0 / 3.0) * m * m * m)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if((!FLA_EXTREME_CASE_TEST) && info == 0)
    {
        if(JOBVL == 'V' || JOBVR == 'V')
        {
            validate_ggev(&JOBVL, &JOBVR, m, A, lda, B, ldb, alpha, alphar, alphai, beta, VL, ldvl,
                          VR, ldvr, datatype, residual, &vinfo);
        }
        else
        { /* For JOBVL = JOBVR = N, eigen values are validated by comparing them with
             the eigen values generated from JOBVL = V, JOBVR = N case */
            if(JOBVL == 'N' && JOBVR == 'N')
            {
                void *A_copy = NULL, *B_copy = NULL, *beta_copy = NULL, *VL_copy = NULL;
                void *alphar_copy = NULL, *alphai_copy = NULL, *alpha_copy = NULL;
                double time_min_copy = 1e9;

                g_lwork = -1;
                create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &VL_copy, m);
                create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_copy, lda);
                create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &B_copy, ldb);
                if(datatype == FLOAT || datatype == DOUBLE)
                {
                    create_vector(datatype, &alphar_copy, m);
                    create_vector(datatype, &alphai_copy, m);
                }
                else
                {
                    create_vector(datatype, &alpha_copy, m);
                }
                create_vector(datatype, &beta_copy, m);
                copy_matrix(datatype, "full", m, m, A, lda, A_copy, lda);
                copy_matrix(datatype, "full", m, m, B, ldb, B_copy, ldb);

                prepare_ggev_run("V", &JOBVR, m, A_copy, lda, B_copy, ldb, alpha_copy, alphar_copy,
                                 alphai_copy, beta_copy, VL_copy, m, NULL, ldvr, datatype,
                                 n_repeats, &time_min_copy, &info, interfacetype, layout);
                /* Valdiate eigen values from both the runs
                  (JOBVL = JOBVR = N with that of JOBVL = V and JOBVR = N)*/
                if(info == 0)
                {
                    validate_ggev_EVs(m, alpha, alphar, alphai, beta, alpha_copy, alphar_copy,
                                      alphai_copy, beta_copy, datatype, residual);
                }
                free_matrix(VL_copy);
                free_matrix(A_copy);
                free_matrix(B_copy);
                if(datatype == FLOAT || datatype == DOUBLE)
                {
                    free_vector(alphar_copy);
                    free_vector(alphai_copy);
                }
                else
                {
                    free_vector(alpha_copy);
                }
                free_vector(beta_copy);
            }
        }
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, m, m, A_test, lda, params->imatrix_char))
           && (!check_extreme_value(datatype, m, m, B_test, ldb, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(VL);
    free_matrix(VR);
    free_matrix(B);
    free_matrix(B_test);
    if(datatype == FLOAT || datatype == DOUBLE)
    {
        free_vector(alphar);
        free_vector(alphai);
    }
    else
    {
        free_vector(alpha);
    }

    free_vector(beta);
}

void prepare_ggev_run(char *jobvl, char *jobvr, integer n_A, void *A, integer lda, void *B,
                      integer ldb, void *alpha, void *alphar, void *alphai, void *beta, void *VL,
                      integer ldvl, void *VR, integer ldvr, integer datatype, integer n_repeats,
                      double *time_min_, integer *info, integer interfacetype,
                      int layout)
{
    void *A_save = NULL, *B_save = NULL, *work = NULL, *rwork = NULL;
    integer i;
    integer lwork;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &B_save, ldb);
    copy_matrix(datatype, "full", n_A, n_A, B, ldb, B_save, ldb);

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.
       NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST) && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 8 * n_A);

        /* call to ggev API to get work query */
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_ggev(datatype, jobvl, jobvr, &n_A, NULL, &lda, NULL, &ldb, NULL, NULL, NULL, NULL,
                            NULL, &ldvl, NULL, &ldvr, work, &lwork, rwork, info);
        }
        else
#endif
        {
            invoke_ggev(datatype, jobvl, jobvr, &n_A, NULL, &lda, NULL, &ldb, NULL, NULL, NULL, NULL,
                        NULL, &ldvl, NULL, &ldvr, work, &lwork, rwork, info);
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
        /* Restore input matrix A value and allocate memory to output buffers for each iteration */
        copy_matrix(datatype, "full", n_A, n_A, A_save, lda, A, lda);
        copy_matrix(datatype, "full", n_A, n_A, B_save, ldb, B, ldb);
        create_vector(datatype, &work, lwork);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            create_realtype_vector(datatype, &rwork, 8 * n_A);
        }
        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time
                = prepare_lapacke_ggev_run(datatype, layout, jobvl, jobvr, n_A, A, lda, B, ldb,
                                           alpha, alphar, alphai, beta, VL, ldvl, VR, ldvr, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)   /* Call CPP ggev API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_ggev(datatype, jobvl, jobvr, &n_A, A, &lda, B, &ldb, alpha, alphar, alphai, beta,
                            VL, &ldvl, VR, &ldvr, work, &lwork, rwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK ggev API */
            invoke_ggev(datatype, jobvl, jobvr, &n_A, A, &lda, B, &ldb, alpha, alphar, alphai, beta,
                        VL, &ldvl, VR, &ldvr, work, &lwork, rwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
    }

    *time_min_ = time_min;
    copy_matrix(datatype, "full", n_A, n_A, A_save, lda, A, lda);
    copy_matrix(datatype, "full", n_A, n_A, B_save, ldb, B, ldb);

    free_matrix(A_save);
    free_matrix(B_save);
}

double prepare_lapacke_ggev_run(integer datatype, int layout, char *jobvl, char *jobvr,
                                integer n_A, void *A, integer lda, void *B, integer ldb,
                                void *alpha, void *alphar, void *alphai, void *beta, void *vl,
                                integer ldvl, void *vr, integer ldvr, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldb_t = ldb;
    integer ldvl_t = ldvl;
    integer ldvr_t = ldvr;
    void *A_t = NULL, *B_t = NULL, *vl_t = NULL, *vr_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggev_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggev_ldb, ldb_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggev_ldvl, ldvl_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggev_ldvr, ldvr_t);

    A_t = A;
    B_t = B;
    vl_t = vl;
    vr_t = vr;
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n_A, n_A, &A_t, fla_max(n_A, lda_t));
        create_matrix(datatype, layout, n_A, n_A, &B_t, fla_max(n_A, ldb_t));
        if(*jobvl == 'V')
        {
            create_matrix(datatype, layout, n_A, n_A, &vl_t, fla_max(n_A, ldvl_t));
        }
        if(*jobvr == 'V')
        {
            create_matrix(datatype, layout, n_A, n_A, &vr_t, fla_max(n_A, ldvr_t));
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n_A, n_A, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n_A, n_A, B, ldb, B_t, ldb_t);
    }

    exe_time = fla_test_clock();

    /* call LAPACKE ggev API */
    *info = invoke_lapacke_ggev(datatype, layout, *jobvl, *jobvr, n_A, A_t, lda_t, B_t, ldb_t,
                                alpha, alphar, alphai, beta, vl_t, ldvl_t, vr_t, ldvr_t);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n_A, n_A, A_t, lda_t, A, lda);
        convert_matrix_layout(layout, datatype, n_A, n_A, B_t, ldb_t, B, ldb);

        if(*jobvl == 'V')
        {
            convert_matrix_layout(layout, datatype, n_A, n_A, vl_t, ldvl_t, vl, ldvl);
            free_matrix(vl_t);
        }
        if(*jobvr == 'V')
        {
            convert_matrix_layout(layout, datatype, n_A, n_A, vr_t, ldvr_t, vr, ldvr);
            free_matrix(vr_t);
        }
        /* free temporary buffers */
        free_matrix(A_t);
        free_matrix(B_t);
    }
    return exe_time;
}

void invoke_ggev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                 void *b, integer *ldb, void *alpha, void *alphar, void *alphai, void *beta,
                 void *vl, integer *ldvl, void *vr, integer *ldvr, void *work, integer *lwork,
                 void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr,
                             ldvr, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr,
                             ldvr, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work,
                             lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work,
                             lwork, rwork, info);
            break;
        }
    }
}

integer invoke_lapacke_ggev(integer datatype, int layout, char jobvl, char jobvr, integer n,
                            void *a, integer lda, void *b, integer ldb, void *alpha, void *alphar,
                            void *alphai, void *beta, void *vl, integer ldvl, void *vr,
                            integer ldvr)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                                 ldvl, vr, ldvr);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                                 ldvl, vr, ldvr);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr,
                                 ldvr);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zggev(layout, jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr,
                                 ldvr);
            break;
        }
    }
    return info;
}