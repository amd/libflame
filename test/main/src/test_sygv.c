/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#define GET_TRANS_STR(datatype) (((datatype) == FLOAT || (datatype) == DOUBLE) ? "T" : "C")

extern double perf;
extern double time_min;

/* Local prototypes */
void fla_test_sygv_experiment(char *tst_api, test_params_t *params, integer datatype,
                              integer p_cur, integer q_cur, integer pci, integer n_repeats,
                              integer einfo);
void prepare_sygv_run(integer itype, char *jobz, char *uplo, integer n, void *A, integer lda,
                      void *B, integer ldb, void *w, integer datatype, integer *info,
                      test_params_t *params);
void invoke_sygv(integer datatype, integer *itype, char *jobz, char *uplo, integer *n, void *a,
                 integer *lda, void *b, integer *ldb, void *w, void *work, integer *lwork,
                 void *rwork, integer *info);

void fla_test_sygv(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Eigen Decomposition";
    char *front_str = "SYGV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_sygv_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].itype = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].jobz = argv[4][0];
        params->eig_sym_paramslist[0].uplo = argv[5][0];
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ldb = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        if(n_repeats > 0)
        {
            params->eig_sym_paramslist[0].threshold_value = CLI_NORM_THRESH;

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

                fla_test_sygv_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    if(tests_not_run)
    {
        printf("\nIllegal arguments for sygv/hegv\n");
        printf("./<EXE> sygv <precisions - sd> <ITYPE> <JOBZ> <UPLO> <N> <LDA>"
               " <LDB> <LWORK> <repeats>\n");
        printf("./<EXE> hegv <precisions - cz> <ITYPE> <JOBZ> <UPLO> <N> <LDA>"
               " <LDB> <LWORK> <repeats>\n");
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

void fla_test_sygv_experiment(char *tst_api, test_params_t *params, integer datatype,
                              integer p_cur, integer q_cur, integer pci, integer n_repeats,
                              integer einfo)
{
    integer n, lda, ldb, info = 0, itype = 1;
    char jobz, uplo;
    void *A = NULL, *w = NULL, *A_test = NULL;
    void *B = NULL, *B_test = NULL;
    double residual, err_thresh;

    /* Get input matrix dimensions */
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    err_thresh = params->eig_sym_paramslist[pci].threshold_value;

    n = p_cur;
    lda = params->eig_sym_paramslist[pci].lda;
    ldb = params->eig_sym_paramslist[pci].ldb;

    itype = params->eig_sym_paramslist[pci].itype;

    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
        if(ldb == -1)
        {
            ldb = fla_max(1, n);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B, ldb);
    create_realtype_vector(datatype, &w, n);

    if(g_ext_fptr != NULL)
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        init_matrix(datatype, B, n, n, ldb, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        void *EVals = NULL, *U = NULL, *C = NULL, *temp = NULL;

        create_realtype_vector(datatype, &EVals, n);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &U, lda);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &C, lda);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &temp, lda);

        generate_matrix_from_EVs(datatype, 'V', n, U, lda, EVals, 0.1, 1.0,
                                 USE_ABS_EIGEN_VALUES);
        get_triangular_matrix("U", datatype, n, n, U, lda, 1, NON_UNIT_DIAG);

        /* B = U**{T|C} U */
        fla_invoke_gemm(datatype, GET_TRANS_STR(datatype), "N", &n, &n, &n, d_one, U, &lda, U,
                        &lda, d_zero, B, &ldb);

        generate_matrix_from_EVs(datatype, 'R', n, C, lda, EVals, 0.0, 0.0,
                                 USE_ABS_EIGEN_VALUES);

        switch(itype)
        {
            case 1:
                fla_invoke_trmm(datatype, "L", "U", GET_TRANS_STR(datatype), "N", &n, &n, U,
                                &lda, C, &lda);
                fla_invoke_trmm(datatype, "R", "U", "N", "N", &n, &n, U, &lda, C, &lda);
                copy_matrix(datatype, "full", n, n, C, lda, A, lda);
                break;
            case 2:
            case 3:
                copy_matrix(datatype, "full", n, n, C, lda, temp, lda);
                fla_invoke_trsm(datatype, "L", "U", "N", "N", &n, &n, U, &lda, temp, &lda);
                copy_matrix(datatype, "full", n, n, temp, lda, A, lda);
                fla_invoke_trsm(datatype, "R", "U", GET_TRANS_STR(datatype), "N", &n, &n, U,
                                &lda, A, &lda);
                break;
        }

        free_vector(EVals);
        free_matrix(U);
        free_matrix(C);
        free_matrix(temp);
    }

    /* Make a copy of input matrix A and B */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_test, ldb);
    copy_matrix(datatype, "full", n, n, B, ldb, B_test, ldb);

    prepare_sygv_run(itype, &jobz, &uplo, n, A_test, lda, B_test, ldb, w, datatype, &info, params);

    /* performance computation
     * https://support.nag.com/numeric/nl/nagdoc_latest/clhtml/f08/f08sac.html
     * (8/3)n^3 [syev] + (1/3)n^3 [potrf] flops for eigen vectors
     * (4/3)n^3 [syev] + (1/3)n^3 [potrf] flops for eigen values */
    if(same_char(jobz, 'V'))
        perf = (double)(3.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    else
        perf = (double)((5.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);

    /* Free up the buffers */
    free_matrix(A_test);
    free_matrix(B_test);
    free_matrix(A);
    free_matrix(B);
    free_vector(w);
}

void prepare_sygv_run(integer itype, char *jobz, char *uplo, integer n, void *A, integer lda,
                      void *B, integer ldb, void *w, integer datatype, integer *info,
                      test_params_t *params)
{
    void *A_save, *work, *rwork = NULL, *B_save = NULL;
    integer lwork;
    double exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_save, ldb);
    copy_matrix(datatype, "full", n, n, B, ldb, B_save, ldb);

    if(g_lwork <= 0)
    {
        lwork = -1;

        create_vector(datatype, &work, 1);
        invoke_sygv(datatype, &itype, jobz, uplo, &n, NULL, &lda, NULL, &ldb, NULL, work, &lwork,
                    NULL, info);
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
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        copy_matrix(datatype, "full", n, n, B_save, ldb, B, ldb);

        create_vector(datatype, &work, lwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, fla_max(1, 3 * n - 2));
        else
            rwork = NULL;

        exe_time = fla_test_clock();
        invoke_sygv(datatype, &itype, jobz, uplo, &n, A, &lda, B, &ldb, w, work, &lwork, rwork,
                    info);
        exe_time = fla_test_clock() - exe_time;

        FLA_EXEC_LOOP_UPDATE_WITH_INFO

        free_vector(work);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);
    }

    free_matrix(A_save);
    free_matrix(B_save);
}

void invoke_sygv(integer datatype, integer *itype, char *jobz, char *uplo, integer *n, void *a,
                 integer *lda, void *b, integer *ldb, void *w, void *work, integer *lwork,
                 void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info);
            break;
        }
    }
}
