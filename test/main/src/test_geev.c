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
integer row_major_geev_lda;
integer row_major_geev_ldvl;
integer row_major_geev_ldvr;

/* Local prototypes.*/
void fla_test_geev_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo);
void prepare_geev_run(char *jobvl, char *jobvr, integer n, void *a, integer lda, void *wr, void *wi,
                      void *w, void *vl, integer ldvl, void *vr, integer ldvr, integer datatype,
                      integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                      int matrix_layout);
void invoke_geev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                 void *wr, void *wi, void *w, void *vl, integer *ldvl, void *vr, integer *ldvr,
                 void *work, integer *lwork, void *rwork, integer *info);
double prepare_lapacke_geev_run(integer datatype, int matrix_layout, char *jobvl, char *jobvr,
                                integer n, void *a, integer lda, void *wr, void *wi, void *w,
                                void *vl, integer ldvl, void *vr, integer ldvr, integer *info);

void fla_test_geev(integer argc, char **argv, test_params_t *params)
{
    srand(5); /* Setting the seed for random input genetation values */
    char *op_str = "Eigen Decomposition of non symmetric matrix";
    char *front_str = "GEEV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_geev_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
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
            row_major_geev_lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            row_major_geev_ldvl = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            row_major_geev_ldvr = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].lda = N;
            params->eig_non_sym_paramslist[0].ldvl = N;
            params->eig_non_sym_paramslist[0].ldvr = N;
        }
        else
        {
            params->eig_non_sym_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        }

        g_lwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_geev_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for geev\n");
        printf("./<EXE> geev <precisions - sdcz> <jobvl> <jobvr> <N> <LDA> <LDVL> <LDVR> <LWORK> "
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

void fla_test_geev_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo)
{
    integer n, lda, ldvl, ldvr;
    integer info = 0;
    void *A = NULL, *wr = NULL, *wi = NULL, *w = NULL, *VL = NULL, *VR = NULL;
    void *A_test = NULL, *L = NULL, *wr_in = NULL, *wi_in = NULL, *scal = NULL;
    char jobvl, jobvr;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    n = p_cur;
    lda = params->eig_non_sym_paramslist[pci].lda;
    ldvl = params->eig_non_sym_paramslist[pci].ldvl;
    ldvr = params->eig_non_sym_paramslist[pci].ldvr;

    err_thresh = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    jobvl = params->eig_non_sym_paramslist[pci].jobvsl;
    jobvr = params->eig_non_sym_paramslist[pci].jobvsr;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
        /* LDVL >= 1; if JOBVL = 'V', LDVL >= N */
        if(ldvl == -1)
        {
            if(jobvl == 'V')
            {
                ldvl = n;
            }
            else
            {
                ldvl = 1;
            }
        }
        /* LDVR >= 1; if JOBVR = 'V', LDVR >= N */
        if(ldvr == -1)
        {
            if(jobvr == 'V')
            {
                ldvr = n;
            }
            else
            {
                ldvr = 1;
            }
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &VL, ldvl);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &VR, ldvr);
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_vector(datatype, &w, n);
    }
    else
    {
        create_vector(datatype, &wr, n);
        create_vector(datatype, &wi, n);
    }
    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        /*  Creating input matrix A by generating random eigen values */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &L, n);
        generate_asym_matrix_from_EVs(datatype, n, A, lda, L);

        /* Overflow/underflow initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_vector(get_realtype(datatype), &scal, 1);
            init_matrix_overflow_underflow_asym(datatype, n, n, A, lda, params->imatrix_char, scal);
        }
        /* Diagonal and sub-diagonals(upper and lower sub-diagonal together
           contains imaginary parts) contain real and imaginary parts
           of eigen values respectively. Storing them for valiation purpose */
        create_vector(datatype, &wr_in, n);
        get_diagonal(datatype, L, n, n, n, wr_in);

        if(datatype == FLOAT || datatype == DOUBLE)
        {
            create_vector(datatype, &wi_in, n);
            reset_vector(datatype, wi_in, n, 1);
            get_subdiagonal(datatype, L, n, n, n, wi_in);
        }
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_geev_run(&jobvl, &jobvr, n, A_test, lda, wr, wi, w, VL, ldvl, VR, ldvr, datatype,
                     n_repeats, &time_min, &info, interfacetype, layout);

    /* performance computation
       4/3 n^3 flops if job = 'N'
       8/3 n^3 flops if job = 'V' */

    if(jobvl == 'N' && jobvr == 'N')
        perf = (double)((4.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    else
        perf = (double)((8.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_geev(tst_api, &jobvl, &jobvr, n, A, A_test, lda, VL, ldvl, VR, ldvr, w, wr, wi,
                      datatype, params->imatrix_char, scal, residual, wr_in, wi_in);
    }
    else
    {
        printf("Extreme Value tests not supported for xGEEV APIs\n");
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(VL);
    free_matrix(VR);
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        free_vector(w);
    }
    else
    {
        free_vector(wr);
        free_vector(wi);
    }
    if((g_ext_fptr == NULL) && !(FLA_EXTREME_CASE_TEST))
    {
        free_matrix(L);
        free_vector(wr_in);
        if(datatype == FLOAT || datatype == DOUBLE)
        {
            free_vector(wi_in);
        }
    }
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
        free_vector(scal);
}

void prepare_geev_run(char *jobvl, char *jobvr, integer n_A, void *A, integer lda, void *wr,
                      void *wi, void *w, void *VL, integer ldvl, void *VR, integer ldvr,
                      integer datatype, integer n_repeats, double *time_min_, integer *info,
                      integer interfacetype, int layout)
{
    void *A_save = NULL, *rwork = NULL, *work = NULL;
    integer lwork, lrwork;
    integer i;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = 2 * n_A;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_realtype_vector(datatype, &rwork, lrwork);
    }

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.
     NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

#if ENABLE_CPP_TEST
        /* call to  geev API */
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_geev(datatype, jobvl, jobvr, &n_A, NULL, &lda, NULL, NULL, NULL, NULL, &ldvl,
                            NULL, &ldvr, work, &lwork, rwork, info);
        }
        else
#endif
        {
            invoke_geev(datatype, jobvl, jobvr, &n_A, NULL, &lda, NULL, NULL, NULL, NULL, &ldvl,
                        NULL, &ldvr, work, &lwork, rwork, info);
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
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        free_vector(rwork);
    }

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n_A, n_A, A_save, lda, A, lda);

        create_vector(datatype, &work, lwork);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            create_realtype_vector(datatype, &rwork, lrwork);
        }
        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_geev_run(datatype, layout, jobvl, jobvr, n_A, A, lda, wr, wi,
                                                w, VL, ldvl, VR, ldvr, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST) /* Call CPP geev API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_geev(datatype, jobvl, jobvr, &n_A, A, &lda, wr, wi, w, VL, &ldvl, VR, &ldvr,
                            work, &lwork, rwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else /* Call LAPACK geev API */
        {
            exe_time = fla_test_clock();
            invoke_geev(datatype, jobvl, jobvr, &n_A, A, &lda, wr, wi, w, VL, &ldvl, VR, &ldvr,
                        work, &lwork, rwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
    }

    *time_min_ = t_min;

    free_matrix(A_save);
}

double prepare_lapacke_geev_run(integer datatype, int layout, char *jobvl, char *jobvr, integer n,
                                void *a, integer lda, void *wr, void *wi, void *w, void *vl,
                                integer ldvl, void *vr, integer ldvr, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldvl_t = ldvl;
    integer ldvr_t = ldvr;
    void *a_t = NULL, *vl_t = NULL, *vr_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_geev_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_geev_ldvl, ldvl_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_geev_ldvr, ldvr_t);

    a_t = a;
    vl_t = vl;
    vr_t = vr;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &a_t, fla_max(n, lda_t));
        if(*jobvl == 'V')
        {
            create_matrix(datatype, layout, n, n, &vl_t, fla_max(n, ldvl_t));
        }
        if(*jobvr == 'V')
        {
            create_matrix(datatype, layout, n, n, &vr_t, fla_max(n, ldvr_t));
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, a, lda, a_t, lda_t);
    }
    exe_time = fla_test_clock();

    /* call to LAPACKE geev API */
    *info = invoke_lapacke_geev(datatype, layout, *jobvl, *jobvr, n, a_t, lda_t, wr, wi, w, vl_t,
                                ldvl_t, vr_t, ldvr_t);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */
    if(layout == LAPACK_ROW_MAJOR)
    {
        convert_matrix_layout(layout, datatype, n, n, a_t, lda_t, a, lda);
        if(*jobvl == 'V')
        {
            convert_matrix_layout(layout, datatype, n, n, vl_t, ldvl_t, vl, ldvl);
            free_matrix(vl_t);
        }
        if(*jobvr == 'V')
        {
            convert_matrix_layout(layout, datatype, n, n, vr_t, ldvr_t, vr, ldvr);
            free_matrix(vr_t);
        }
        free_matrix(a_t);
    }
    return exe_time;
}

void invoke_geev(integer datatype, char *jobvl, char *jobvr, integer *n, void *a, integer *lda,
                 void *wr, void *wi, void *w, void *vl, integer *ldvl, void *vr, integer *ldvr,
                 void *work, integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork,
                             info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork,
                             info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork,
                             info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork,
                             info);
            break;
        }
    }
}
