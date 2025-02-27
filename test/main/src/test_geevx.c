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
integer row_major_geevx_lda;
integer row_major_geevx_ldvl;
integer row_major_geevx_ldvr;

/* Local prototypes.*/
void fla_test_geevx_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_geevx_run(char *balanc, char *jobvl, char *jobvr, char *sense, integer n, void *a,
                       integer lda, void *wr, void *wi, void *w, void *vl, integer ldvl, void *vr,
                       integer ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                       void *rconde, void *rcondv, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer interfacetype, int matrix_layout);
void invoke_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                  void *a, integer *lda, void *wr, void *wi, void *w, void *vl, integer *ldvl,
                  void *vr, integer *ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                  void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                  integer *iwork, integer *info);
double prepare_lapacke_geevx_run(integer datatype, int matrix_layout, char *balanc, char *jobvl,
                                 char *jobvr, char *sense, integer n, void *a, integer lda,
                                 void *wr, void *wi, void *w, void *vl, integer ldvl, void *vr,
                                 integer ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                                 void *rconde, void *rcondv, integer *info);


void fla_test_geevx(integer argc, char **argv, test_params_t *params)
{
    srand(1); /* Setting the seed for random input genetation values */
    char *op_str = "Eigen Decomposition of non symmetric matrix";
    char *front_str = "GEEVX";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_geevx_experiment);
        tests_not_run = 0;
    }
    if(argc == 14)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[13]);
    }
    if(argc >= 13 && argc <= 14)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_non_sym_paramslist[0].balance_ggevx = argv[3][0];
        params->eig_non_sym_paramslist[0].jobvsl = argv[4][0];
        params->eig_non_sym_paramslist[0].jobvsr = argv[5][0];
        params->eig_non_sym_paramslist[0].sense_ggevx = argv[6][0];
        N = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_geevx_lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_geevx_ldvl = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            row_major_geevx_ldvr = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].lda = N;
            params->eig_non_sym_paramslist[0].ldvl = N;
            params->eig_non_sym_paramslist[0].ldvr = N;
        }
        else
        {
            params->eig_non_sym_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_geevx_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for geevx\n");
        printf("./<EXE> geevx <precisions - sdcz> <balanc> <jobvl> <jobvr> <sense> <N> <LDA> "
               "<LDVL> <LDVR> <LWORK> <repeats>\n");
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

void fla_test_geevx_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, lda, ldvl, ldvr;
    integer info = 0;
    integer ilo, ihi;
    void *A = NULL, *wr = NULL, *wi = NULL, *w = NULL, *VL = NULL, *VR = NULL;
    void *scale = NULL, *abnrm = NULL, *rconde = NULL, *rcondv = NULL;
    void *A_test = NULL, *L = NULL, *wr_in = NULL, *wi_in = NULL, *scal = NULL;
    char balanc, jobvl, jobvr, sense;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    m = p_cur;
    lda = params->eig_non_sym_paramslist[pci].lda;
    ldvl = params->eig_non_sym_paramslist[pci].ldvl;
    ldvr = params->eig_non_sym_paramslist[pci].ldvr;

    err_thresh = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    balanc = params->eig_non_sym_paramslist[pci].balance_ggevx;
    jobvl = params->eig_non_sym_paramslist[pci].jobvsl;
    jobvr = params->eig_non_sym_paramslist[pci].jobvsr;
    sense = params->eig_non_sym_paramslist[pci].sense_ggevx;
    if(sense == 'B' || sense == 'E')
    {
        jobvl = 'V';
        jobvr = 'V';
    }
    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
        /* LDVL >= 1; if JOBVL = 'V', LDVL >= M */
        if(ldvl == -1)
        {
            if(jobvl == 'V')
            {
                ldvl = m;
            }
            else
            {
                ldvl = 1;
            }
        }
        /* LDVR >= 1; if JOBVR = 'V', LDVR >= M */
        if(ldvr == -1)
        {
            if(jobvr == 'V')
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

    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &VL, ldvl);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &VR, ldvr);
    create_realtype_vector(datatype, &scale, m);
    create_realtype_vector(datatype, &abnrm, 1);
    create_realtype_vector(datatype, &rconde, m);
    create_realtype_vector(datatype, &rcondv, m);

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_vector(datatype, &w, m);
    }
    else
    {
        create_vector(datatype, &wr, m);
        create_vector(datatype, &wi, m);
    }

    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        init_matrix(datatype, A, m, m, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        /*  Creating input matrix A by generating random eigen values */
        create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &L, m);
        generate_asym_matrix_from_EVs(datatype, m, A, lda, L);

        /* Overflow/underflow initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_vector(get_realtype(datatype), &scal, 1);
            init_matrix_overflow_underflow_asym(datatype, m, m, A, lda, params->imatrix_char, scal);
        }
        /* Diagonal and sub-diagonals(upper and lower sub-diagonal together
           contains imaginary parts) contain real and imaginary parts
           of eigen values respectively. Storing them for valiation purpose */
        create_vector(datatype, &wr_in, m);
        get_diagonal(datatype, L, m, m, m, wr_in);

        if(datatype == FLOAT || datatype == DOUBLE)
        {
            create_vector(datatype, &wi_in, m);
            reset_vector(datatype, wi_in, m, 1);
            get_subdiagonal(datatype, L, m, m, m, wi_in);
        }
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_test, lda);
    copy_matrix(datatype, "full", m, m, A, lda, A_test, lda);

    prepare_geevx_run(&balanc, &jobvl, &jobvr, &sense, m, A_test, lda, wr, wi, w, VL, ldvl, VR,
                      ldvr, &ilo, &ihi, scale, abnrm, rconde, rcondv, datatype, n_repeats,
                      &time_min, &info, interfacetype, layout);

    /* performance computation
       4/3 m^3 flops if job = 'N'
       8/3 m^3 + m^2 flops if job = 'V' */
    if(jobvl == 'N' && jobvr == 'N')
        perf = (double)((4.0 / 3.0) * m * m * m) / time_min / FLOPS_PER_UNIT_PERF;
    else
        perf = (double)(((8.0 / 3.0) * m * m * m) + (m * m)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_geevx(tst_api, &jobvl, &jobvr, &sense, &balanc, m, A, A_test, lda, VL, ldvl, VR,
                       ldvr, w, wr, wi, scale, abnrm, rconde, rcondv, datatype,
                       params->imatrix_char, scal, residual, wr_in, wi_in);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if(!check_extreme_value(datatype, m, m, A_test, lda, params->imatrix_char))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, m, residual, err_thresh);
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(VL);
    free_matrix(VR);
    free_vector(scale);
    free_vector(abnrm);
    free_vector(rconde);
    free_vector(rcondv);
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

void prepare_geevx_run(char *balanc, char *jobvl, char *jobvr, char *sense, integer m_A, void *A,
                       integer lda, void *wr, void *wi, void *w, void *VL, integer ldvl, void *VR,
                       integer ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                       void *rconde, void *rcondv, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer interfacetype, int layout)
{
    void *A_save = NULL, *rwork = NULL, *iwork = NULL, *work = NULL;
    integer lwork, liwork, lrwork;
    integer i;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, m_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, m_A, A, lda, A_save, lda);

    /* Get rwork and iwork array size since it is not depedent on internal blocks*/
    lrwork = 2 * m_A;
    liwork = 2 * m_A - 2;

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
        /* call to  geevx API */
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, NULL, &lda, NULL, NULL,
                             NULL, NULL, &ldvl, NULL, &ldvr, ilo, ihi, NULL, NULL, NULL, NULL, work,
                             &lwork, rwork, NULL, info);
        }
        else
#endif
        {
            invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, NULL, &lda, NULL, NULL, NULL,
                         NULL, &ldvl, NULL, &ldvr, ilo, ihi, NULL, NULL, NULL, NULL, work, &lwork,
                         rwork, NULL, info);
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
        copy_matrix(datatype, "full", m_A, m_A, A_save, lda, A, lda);

        create_vector(datatype, &work, lwork);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            create_realtype_vector(datatype, &rwork, lrwork);
        }
        else
        {
            create_vector(INTEGER, &iwork, liwork);
        }
        /* Check if LAPACKE interface enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_geevx_run(datatype, layout, balanc, jobvl, jobvr, sense, m_A,
                                                 A, lda, wr, wi, w, VL, ldvl, VR, ldvr, ilo, ihi,
                                                 scale, abnrm, rconde, rcondv, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST) /* Call CPP geevx API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, A, &lda, wr, wi, w, VL,
                             &ldvl, VR, &ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, &lwork,
                             rwork, iwork, info);

            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else /* Call LAPACK geevx API */
        {
            exe_time = fla_test_clock();
            invoke_geevx(datatype, balanc, jobvl, jobvr, sense, &m_A, A, &lda, wr, wi, w, VL, &ldvl,
                         VR, &ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work, &lwork, rwork,
                         iwork, info);

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
        else
        {
            free_vector(iwork);
        }
    }
    *time_min_ = t_min;

    free_matrix(A_save);
}

double prepare_lapacke_geevx_run(integer datatype, int layout, char *balanc, char *jobvl,
                                 char *jobvr, char *sense, integer n, void *a, integer lda,
                                 void *wr, void *wi, void *w, void *vl, integer ldvl, void *vr,
                                 integer ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                                 void *rconde, void *rcondv, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldvl_t = ldvl;
    integer ldvr_t = ldvr;
    void *a_t = NULL, *vl_t = NULL, *vr_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_geevx_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_geevx_ldvl, ldvl_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_geevx_ldvr, ldvr_t);

    a_t = a;
    vl_t = vl;
    vr_t = vr;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &a_t, lda_t);
        if(*jobvl == 'V')
        {
            create_matrix(datatype, layout, n, n, &vl_t, ldvl_t);
        }
        if(*jobvr == 'V')
        {
            create_matrix(datatype, layout, n, n, &vr_t, ldvr_t);
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, a, lda, a_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* Call LAPACKE geevx API */
    *info = invoke_lapacke_geevx(datatype, layout, *balanc, *jobvl, *jobvr, *sense, n, a_t, lda_t,
                                 wr, wi, w, vl_t, ldvl_t, vr_t, ldvr_t, ilo, ihi, scale, abnrm,
                                 rconde, rcondv);

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

void invoke_geevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                  void *a, integer *lda, void *wr, void *wi, void *w, void *vl, integer *ldvl,
                  void *vr, integer *ldvr, integer *ilo, integer *ihi, void *scale, void *abnrm,
                  void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                  integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgeevx(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr,
                              ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgeevx(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr,
                              ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo,
                              ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo,
                              ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info);
            break;
        }
    }
}
