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
integer row_major_ggevx_lda;
integer row_major_ggevx_ldb;
integer row_major_ggevx_ldvl;
integer row_major_ggevx_ldvr;

/* Local prototypes */
void fla_test_ggevx_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_ggevx_run(char *balanc, char *jobvl, char *jobvr, char *sense, integer n_A, void *A,
                       integer lda, void *B, integer ldb, void *alpha, void *alphar, void *alphai,
                       void *beta, void *VL, integer ldvl, void *VR, integer ldvr, integer *ilo,
                       integer *ihi, void *lscale, void *rscale, void *abnrm, void *bbnrm,
                       void *rconde, void *rcondv, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer interfacetype, int matrix_layout);
void invoke_ggevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                  void *a, integer *lda, void *b, integer *ldb, void *alpha, void *alphar,
                  void *alphai, void *beta, void *vl, integer *ldvl, void *vr, integer *ldvr,
                  integer *ilo, integer *ihi, void *lscale, void *rscale, void *abnrm, void *bbnrm,
                  void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                  integer *iwork, integer *bwork, integer *info);
double prepare_lapacke_ggevx_run(integer datatype, int matrix_layout, char *balanc, char *jobvl,
                                 char *jobvr, char *sense, integer n_A, void *A, integer lda,
                                 void *B, integer ldb, void *alpha, void *alphar, void *alphai,
                                 void *beta, void *VL, integer ldvl, void *VR, integer ldvr,
                                 integer *ilo, integer *ihi, void *lscale, void *rscale,
                                 void *abnrm, void *bbnrm, void *rconde, void *rcondv,
                                 integer *info);

void fla_test_ggevx(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Computing Eigen value and Eigen vectors with condition numbers";
    char *front_str = "GGEVX";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        /* Test with parameters from config */
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_ggevx_experiment);
        tests_not_run = 0;
    }
    if(argc == 15)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[14]);
    }
    if(argc >= 14 && argc <= 15)
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
            row_major_ggevx_lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_ggevx_ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            row_major_ggevx_ldvl = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            row_major_ggevx_ldvr = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].lda = N;
            params->eig_non_sym_paramslist[0].ldb = N;
            params->eig_non_sym_paramslist[0].ldvl = N;
            params->eig_non_sym_paramslist[0].ldvr = N;
        }
        else
        {
            params->eig_non_sym_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvl = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].ldvr = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_ggevx_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for GGEVX\n");
        printf("./<EXE> ggevx <precisions - sdcz> <balanc> <jobvl> <jobvr> <sense> <N> <LDA> <LDB> "
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

void fla_test_ggevx_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, lda, ldb, ldvl, ldvr;
    integer info = 0;
    integer ilo, ihi;
    void *A = NULL, *B = NULL, *VL = NULL, *VR = NULL;
    void *rscale = NULL, *lscale = NULL, *alpha = NULL, *alphar = NULL, *alphai = NULL,
         *beta = NULL, *A_test = NULL, *B_test = NULL;
    void *abnrm = NULL, *bbnrm = NULL, *rconde = NULL, *rcondv = NULL;
    char JOBVL = params->eig_non_sym_paramslist[pci].jobvsl;
    char JOBVR = params->eig_non_sym_paramslist[pci].jobvsr;
    char BALANC = params->eig_non_sym_paramslist[pci].balance_ggevx;
    char SENSE = params->eig_non_sym_paramslist[pci].sense_ggevx;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    err_thresh = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    /* Get input matrix dimensions */
    n = p_cur;
    lda = params->eig_non_sym_paramslist[pci].lda;
    ldvl = params->eig_non_sym_paramslist[pci].ldvl;
    ldvr = params->eig_non_sym_paramslist[pci].ldvr;
    ldb = params->eig_non_sym_paramslist[pci].ldb;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
        if(ldb == -1)
        {
            ldb = fla_max(1, n);
        }
        /* LDVL >= 1, and
           if JOBVL = 'V', LDVL >= N */
        if(ldvl == -1)
        {
            if(JOBVL == 'V')
            {
                ldvl = n;
            }
            else
            {
                ldvl = 1;
            }
        }
        /* LDVR >= 1, and
           if JOBVR = 'V', LDVR >= N */
        if(ldvr == -1)
        {
            if(JOBVR == 'V')
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
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &VL, ldvl);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &VR, ldvr);
    create_vector(datatype, &beta, n);

    if(datatype == FLOAT || datatype == DOUBLE)
    {
        create_vector(datatype, &alphar, n);
        create_vector(datatype, &alphai, n);
    }
    else
    {
        create_vector(datatype, &alpha, n);
    }

    create_realtype_vector(datatype, &lscale, n);
    create_realtype_vector(datatype, &rscale, n);
    create_realtype_vector(datatype, &rconde, n);
    create_realtype_vector(datatype, &rcondv, n);
    create_realtype_vector(datatype, &abnrm, 1);
    create_realtype_vector(datatype, &bbnrm, 1);

    init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    init_matrix(datatype, B, n, n, ldb, g_ext_fptr, params->imatrix_char);

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_test, ldb);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    copy_matrix(datatype, "full", n, n, B, ldb, B_test, ldb);

    prepare_ggevx_run(&BALANC, &JOBVL, &JOBVR, &SENSE, n, A, lda, B, ldb, alpha, alphar, alphai,
                      beta, VL, ldvl, VR, ldvr, &ilo, &ihi, lscale, rscale, abnrm, bbnrm, rconde,
                      rcondv, datatype, n_repeats, &time_min, &info, interfacetype, layout);

    /* performance computation */
    /* 2mn^2 - (2/3)n^3 flops */
    perf = (double)((2.0 * n * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        if(JOBVL == 'V' || JOBVR == 'V')
        {
            validate_ggevx(tst_api, &BALANC, &JOBVL, &JOBVR, &SENSE, n, A_test, lda, B_test, ldb,
                           alpha, alphar, alphai, beta, VL, ldvl, VR, ldvr, datatype, residual);
        }
    }
    else
    {
        printf("Extreme Value tests not supported for xGGEVX APIs\n");
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B);
    free_matrix(B_test);
    free_matrix(VL);
    free_matrix(VR);
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
    free_vector(lscale);
    free_vector(rscale);
    free_vector(rconde);
    free_vector(rcondv);
    free_vector(abnrm);
    free_vector(bbnrm);
}

void prepare_ggevx_run(char *balanc, char *jobvl, char *jobvr, char *sense, integer n_A, void *A,
                       integer lda, void *B, integer ldb, void *alpha, void *alphar, void *alphai,
                       void *beta, void *VL, integer ldvl, void *VR, integer ldvr, integer *ilo,
                       integer *ihi, void *lscale, void *rscale, void *abnrm, void *bbnrm,
                       void *rconde, void *rcondv, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer interfacetype, int layout)
{
    void *A_save = NULL, *B_save = NULL, *work = NULL, *rwork = NULL, *iwork = NULL, *bwork = NULL;
    ;
    integer i;
    integer lwork;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &A_save, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &B_save, ldb);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);
    copy_matrix(datatype, "full", n_A, n_A, B, ldb, B_save, ldb);

    if(*sense != 'E')
    {
        create_vector(INTEGER, &iwork, 6 + n_A);
    }
    if(*sense != 'N')
    {
        create_vector(INTEGER, &bwork, n_A);
    }

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.
     NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* call to ggevx API */
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_ggevx(datatype, balanc, jobvl, jobvr, sense, &n_A, A, &lda, B, &ldb, alpha,
                             alphar, alphai, beta, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale,
                             abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, info);
        }
        else
#endif
        {
            invoke_ggevx(datatype, balanc, jobvl, jobvr, sense, &n_A, A, &lda, B, &ldb, alpha,
                         alphar, alphai, beta, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale,
                         abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, info);
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

    create_realtype_vector(datatype, &rwork, 8 * n_A);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers for each iteration */
        copy_matrix(datatype, "full", n_A, n_A, A_save, lda, A, lda);
        copy_matrix(datatype, "full", n_A, n_A, B_save, ldb, B, ldb);
        create_vector(datatype, &work, lwork);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_ggevx_run(datatype, layout, balanc, jobvl, jobvr, sense, n_A,
                                                 A, lda, B, ldb, alpha, alphar, alphai, beta, VL,
                                                 ldvl, VR, ldvr, ilo, ihi, lscale, rscale, abnrm,
                                                 bbnrm, rconde, rcondv, info);
        }
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP ggevx API */
            invoke_cpp_ggevx(datatype, balanc, jobvl, jobvr, sense, &n_A, A, &lda, B, &ldb, alpha,
                             alphar, alphai, beta, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale,
                             abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK ggevx API */
            invoke_ggevx(datatype, balanc, jobvl, jobvr, sense, &n_A, A, &lda, B, &ldb, alpha,
                         alphar, alphai, beta, VL, &ldvl, VR, &ldvr, ilo, ihi, lscale, rscale,
                         abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, bwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }

    *time_min_ = t_min;

    free_matrix(A_save);
    free_matrix(B_save);
    if(*sense != 'E')
    {
        free_vector(iwork);
    }
    if(*sense != 'N')
    {
        free_vector(bwork);
    }
    free_vector(rwork);
}

double prepare_lapacke_ggevx_run(integer datatype, int layout, char *balanc, char *jobvl,
                                 char *jobvr, char *sense, integer n_A, void *A, integer lda,
                                 void *B, integer ldb, void *alpha, void *alphar, void *alphai,
                                 void *beta, void *VL, integer ldvl, void *VR, integer ldvr,
                                 integer *ilo, integer *ihi, void *lscale, void *rscale,
                                 void *abnrm, void *bbnrm, void *rconde, void *rcondv,
                                 integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldb_t = ldb;
    integer ldvl_t = ldvl;
    integer ldvr_t = ldvr;
    void *A_t = NULL, *B_t = NULL, *VL_t = NULL, *VR_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggevx_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggevx_ldb, ldb_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggevx_ldvl, ldvl_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n_A, row_major_ggevx_ldvr, ldvr_t);

    A_t = A;
    B_t = B;
    VL_t = VL;
    VR_t = VR;
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n_A, n_A, &A_t, fla_max(n_A, lda_t));
        create_matrix(datatype, layout, n_A, n_A, &B_t, fla_max(n_A, ldb_t));
        if(*jobvl == 'V')
        {
            create_matrix(datatype, layout, n_A, n_A, &VL_t, fla_max(n_A, ldvl_t));
        }
        if(*jobvr == 'V')
        {
            create_matrix(datatype, layout, n_A, n_A, &VR_t, fla_max(n_A, ldvr_t));
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n_A, n_A, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n_A, n_A, B, ldb, B_t, ldb_t);
    }

    exe_time = fla_test_clock();

    /* call LAPACKE ggevx API */
    *info = invoke_lapacke_ggevx(datatype, layout, *balanc, *jobvl, *jobvr, *sense, n_A, A_t, lda_t,
                                 B_t, ldb_t, alpha, alphar, alphai, beta, VL_t, ldvl_t, VR_t,
                                 ldvr_t, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n_A, n_A, A_t, lda_t, A, lda);
        convert_matrix_layout(layout, datatype, n_A, n_A, B_t, ldb_t, B, ldb);

        if(*jobvl == 'V')
        {
            convert_matrix_layout(layout, datatype, n_A, n_A, VL_t, ldvl_t, VL, ldvl);
            free_matrix(VL_t);
        }
        if(*jobvr == 'V')
        {
            convert_matrix_layout(layout, datatype, n_A, n_A, VR_t, ldvr_t, VR, ldvr);
            free_matrix(VR_t);
        }
        /* free temporary buffers */
        free_matrix(A_t);
        free_matrix(B_t);
    }
    return exe_time;
}

void invoke_ggevx(integer datatype, char *balanc, char *jobvl, char *jobvr, char *sense, integer *n,
                  void *a, integer *lda, void *b, integer *ldb, void *alpha, void *alphar,
                  void *alphai, void *beta, void *vl, integer *ldvl, void *vr, integer *ldvr,
                  integer *ilo, integer *ihi, void *lscale, void *rscale, void *abnrm, void *bbnrm,
                  void *rconde, void *rcondv, void *work, integer *lwork, void *rwork,
                  integer *iwork, integer *bwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta,
                              vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde,
                              rcondv, work, lwork, iwork, bwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta,
                              vl, ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde,
                              rcondv, work, lwork, iwork, bwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl,
                              vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                              work, lwork, rwork, iwork, bwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl,
                              vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                              work, lwork, rwork, iwork, bwork, info);
            break;
        }
    }
}
