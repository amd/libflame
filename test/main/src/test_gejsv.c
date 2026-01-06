/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

#define GEJSV_VL 0.1
#define GEJSV_VU 10

#define GEJSV_VEC_VL 0.1
#define GEJSV_VEC_VU 10000

/* This is the percentage of small
   singular values that should be considered
   as noise. It is used to test joba = A / R */
#define GEJSV_SMALL_SV_PERCENT 0.25

/* Small singular value range */
#define GEJSV_SMALL_SV_VL(datatype) \
    (get_realtype(datatype) == FLOAT ? FLT_MIN * 1e1 : DBL_MIN * 1e1)
#define GEJSV_SMALL_SV_VU(datatype) \
    (get_realtype(datatype) == FLOAT ? FLT_MIN * 1e2 : DBL_MIN * 1e3)

extern double perf;
extern double time_min;
integer row_major_gejsv_lda;
integer row_major_gejsv_ldu;
integer row_major_gejsv_ldv;

/* Local prototypes */
void fla_test_gejsv_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);

void prepare_gejsv_run(integer datatype, char joba, char jobu, char jobv, char jobr, char jobt,
                       char jobp, integer m, integer n, void *A, integer lda, void *S, void *U,
                       integer ldu, void *V, integer ldv, void *stat, integer *istat,
                       integer interfacetype, integer layout, integer *info, test_params_t *params);

void invoke_gejsv(integer datatype, char *joba, char *jobu, char *jobv, char *jobr, char *jobt,
                  char *jobp, integer *m, integer *n, void *A, integer *lda, void *S, void *U,
                  integer *ldu, void *V, integer *ldv, void *work, integer *lwork, void *rwork,
                  integer *lrwork, integer *iwork, integer *info);

double prepare_lapacke_gejsv_run(integer datatype, int layout, char joba, char jobu, char jobv,
                                 char jobr, char jobt, char jobp, integer m, integer n, void *A,
                                 integer lda, void *S, void *U, integer ldu, void *V, integer ldv,
                                 void *stat, integer *istat, integer *info);

void generate_gejsv_test_matrix(integer datatype, char joba, char jobu, char jobv, integer m,
                                integer n, void *A, integer lda, void *S_test);

integer get_gejsv_lwork(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, integer m,
                        integer n);

/* Helper functions for Bit reproducibility tests */
void store_gejsv_outputs(void *filename, integer datatype, char joba, char jobu, char jobv,
                         char jobr, char jobt, char jobp, integer m, integer n, void *A,
                         integer lda, void *S, void *U, integer ldu, void *V, integer ldv,
                         integer g_lwork, integer g_lrwork, integer g_liwork, void *params);
integer check_bit_reproducibility_gejsv(void *filename, integer datatype, char joba, char jobu,
                                        char jobv, char jobr, char jobt, char jobp, integer m,
                                        integer n, void *A, integer lda, void *S, void *U,
                                        integer ldu, void *V, integer ldv, integer g_lwork,
                                        integer g_lrwork, integer g_liwork, void *params);

void fla_test_gejsv(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Singular Value Decomposition using Jacobi method";
    char *front_str = "GEJSV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    integer num_ranges, i;

    /* Arrays to save original range values */
    integer *orig_m_range_start = NULL;
    integer *orig_m_range_end = NULL;
    integer *orig_n_range_start = NULL;
    integer *orig_n_range_end = NULL;

    if(argc == 1)
    {
        /* Test with parameters from config */
        g_lwork = -1;
        g_lrwork = -1;
        g_liwork = -1;
        g_config_data = 1;

        /* gejsv has different configuration for the matrix sizes.
           Need to map these special sizes to the standard size
           variables used by test op driver. */
        num_ranges = params->svd_paramslist[0].num_ranges;

        /* Save original range values before modifying them */
        orig_m_range_start = (integer *)fla_mem_alloc(num_ranges * sizeof(integer));
        orig_m_range_end = (integer *)fla_mem_alloc(num_ranges * sizeof(integer));
        orig_n_range_start = (integer *)fla_mem_alloc(num_ranges * sizeof(integer));
        orig_n_range_end = (integer *)fla_mem_alloc(num_ranges * sizeof(integer));

        /* Check for allocation failures */
        if(orig_m_range_start == NULL || orig_m_range_end == NULL || orig_n_range_start == NULL
           || orig_n_range_end == NULL)
        {
            printf("\nError: Memory allocation failed for range arrays in GEJSV test\n");

            /* Clean up any successfully allocated memory */
            if(orig_m_range_start != NULL)
                free(orig_m_range_start);
            if(orig_m_range_end != NULL)
                free(orig_m_range_end);
            if(orig_n_range_start != NULL)
                free(orig_n_range_start);
            if(orig_n_range_end != NULL)
                free(orig_n_range_end);

            return;
        }

        for(i = 0; i < num_ranges; i++)
        {
            /* Save original values */
            orig_m_range_start[i] = params->svd_paramslist[i].m_range_start;
            orig_m_range_end[i] = params->svd_paramslist[i].m_range_end;
            orig_n_range_start[i] = params->svd_paramslist[i].n_range_start;
            orig_n_range_end[i] = params->svd_paramslist[i].n_range_end;

            /* Set GEJSV-specific values */
            params->svd_paramslist[i].m_range_start = params->svd_paramslist[i].m_gejsv;
            params->svd_paramslist[i].m_range_end = params->svd_paramslist[i].m_gejsv;
            params->svd_paramslist[i].n_range_start = params->svd_paramslist[i].n_gejsv;
            params->svd_paramslist[i].n_range_end = params->svd_paramslist[i].n_gejsv;
        }

        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, SVD, fla_test_gejsv_experiment);
        tests_not_run = 0;

        /* Restore original range values after GEJSV completes */
        for(i = 0; i < num_ranges; i++)
        {
            params->svd_paramslist[i].m_range_start = orig_m_range_start[i];
            params->svd_paramslist[i].m_range_end = orig_m_range_end[i];
            params->svd_paramslist[i].n_range_start = orig_n_range_start[i];
            params->svd_paramslist[i].n_range_end = orig_n_range_end[i];
        }

        /* Free temporary storage */
        free(orig_m_range_start);
        free(orig_m_range_end);
        free(orig_n_range_start);
        free(orig_n_range_end);
    }
    if(argc == 19)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[18]);
    }
    if(argc >= 18 && argc <= 19)
    {
        /* Test with parameters from commandline */
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->svd_paramslist[0].joba_gejsv = argv[3][0];
        params->svd_paramslist[0].jobu_gejsv = argv[4][0];
        params->svd_paramslist[0].jobv_gejsv = argv[5][0];
        params->svd_paramslist[0].jobr_gejsv = argv[6][0];
        params->svd_paramslist[0].jobt_gejsv = argv[7][0];
        params->svd_paramslist[0].jobp_gejsv = argv[8][0];

        M = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gejsv_lda = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            row_major_gejsv_ldu = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
            row_major_gejsv_ldv = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].lda = M;
            params->svd_paramslist[0].ldu = M;
            params->svd_paramslist[0].ldvt = N;
        }
        else
        {
            params->svd_paramslist[0].lda = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldu = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldvt = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[14], &endptr, CLI_DECIMAL_BASE);
        g_lrwork = strtoimax(argv[15], &endptr, CLI_DECIMAL_BASE);
        g_liwork = strtoimax(argv[16], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[17], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

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
                fla_test_gejsv_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for GEJSV\n");
        printf("./<EXE> gejsv <precisions - sdcz> <joba> <jobu> <jobv> <jobr> <jobt> <jobp> <M> "
               "<N> <LDA> "
               "<LDU> <LDV> <LWORK> <LRWORK> <LIWORK> <repeats>\n");
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

void fla_test_gejsv_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, lda, ldu, ldv;
    integer info = 0;
    void *A = NULL, *A_test = NULL, *U = NULL, *V = NULL, *S = NULL, *S_test = NULL;
    void *stat, *istat, *scal = NULL;
    double residual, err_thresh;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;
    char imatrix = params->imatrix_char;

    /* Test that small noise svds are eliminated
           when joba = A */

    m = p_cur;
    n = q_cur;

    err_thresh = params->svd_paramslist[0].svd_threshold;

    char joba = params->svd_paramslist[0].joba_gejsv;
    char jobu = params->svd_paramslist[0].jobu_gejsv;
    char jobv = params->svd_paramslist[0].jobv_gejsv;
    char jobr = params->svd_paramslist[0].jobr_gejsv;
    char jobt = params->svd_paramslist[0].jobt_gejsv;
    char jobp = params->svd_paramslist[0].jobp_gejsv;

    integer test_eliminated_svds
        = (same_char(joba, 'A') && ((integer)(n * GEJSV_SMALL_SV_PERCENT)) > 1);

    lda = params->svd_paramslist[0].lda;
    ldu = params->svd_paramslist[0].ldu;
    ldv = params->svd_paramslist[0].ldvt;

    lda = lda == -1 ? fla_max(1, m) : lda;
    ldu = ldu == -1 ? fla_max(1, m) : ldu;
    ldv = ldv == -1 ? fla_max(1, n) : ldv;

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    if(!same_char(jobu, 'N'))
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &U, ldu);
    }
    if(!same_char(jobv, 'N'))
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &V, ldv);
    }
    create_vector(get_realtype(datatype), &S, n);
    create_vector(get_realtype(datatype), &S_test, n);
    create_vector(get_realtype(datatype), &stat, 7);
    create_vector(INTEGER, &istat, 4);

    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST) || (FLA_RANDOM_INIT_MODE))
        {
            init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
            init_matrix(get_realtype(datatype), S, n, 1, lda, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            /* Generate matrix A based on the input parameters */
            generate_gejsv_test_matrix(datatype, joba, jobu, jobv, m, n, A, lda, S);
            if(FLA_OVERFLOW_UNDERFLOW_TEST)
            {
                create_vector(get_realtype(datatype), &scal, 1);
                scale_matrix_overflow_underflow_gejsv(datatype, m, n, A, lda, S, imatrix, scal);
            }
        }
    }

    /* This macro is used in the BRT test cases for the following purposes:
     *    - In the Ground truth runs (BRT_char => G, F), the output is stored in a file for future
     * reference
     *    - In the verification runs (BRT_char => V, M), the output is loaded from the file and
     * passed as input to the API
     * */
    FLA_BRT_PROCESS_TWO_INPUT(datatype, m, n, A, lda, get_realtype(datatype), n, 1, S, lda,
                              "ccccccdddddddd", joba, jobu, jobv, jobr, jobt, jobp, m, n, lda, ldu,
                              ldv, g_lwork, g_lrwork, g_liwork)

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_gejsv_run(datatype, joba, jobu, jobv, jobr, jobt, jobp, m, n, A_test, lda, S_test, U,
                      ldu, V, ldv, stat, istat, interfacetype, layout, &info, params);

    /* performance computation */
    /*  6mn^2 + 8n^3 flops */
    perf = (double)((6.0 * m * n * n) + (8.0 * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    IF_FLA_BRT_VALIDATION(
        m, n,
        store_gejsv_outputs(filename, datatype, joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S,
                            U, ldu, V, ldv, g_lwork, g_lrwork, g_liwork, params),
        validate_gejsv(tst_api, joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, S_test, U, ldu,
                       V, ldv, stat, istat, test_eliminated_svds, datatype, residual, scal, imatrix,
                       params),
        check_bit_reproducibility_gejsv(filename, datatype, joba, jobu, jobv, jobr, jobt, jobp, m,
                                        n, A, lda, S, U, ldu, V, ldv, g_lwork, g_lrwork, g_liwork,
                                        params))
    else if(FLA_SKIP_VALIDATION_MODE)
    {
        /* Skip validation for performance modes */
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_gejsv(tst_api, joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, S_test, U, ldu,
                       V, ldv, stat, istat, test_eliminated_svds, datatype, residual, scal, imatrix,
                       params);
    }
    else
    {
        if(!check_extreme_value(datatype, n, 1, S_test, lda, params->imatrix_char))
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
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
            free_vector(scal);
    }
free_buffers:
    FLA_FREE_FILENAME(filename);
    free_matrix(A);
    free_matrix(A_test);
    free_vector(S);
    free_vector(S_test);
    free_vector(stat);
    free_vector(istat);
    if(!same_char(jobu, 'N'))
    {
        free_matrix(U);
    }
    if(!same_char(jobv, 'N'))
    {
        free_matrix(V);
    }
}

void prepare_gejsv_run(integer datatype, char joba, char jobu, char jobv, char jobr, char jobt,
                       char jobp, integer m, integer n, void *A, integer lda, void *S, void *U,
                       integer ldu, void *V, integer ldv, void *stat, integer *istat,
                       integer interfacetype, integer layout, integer *info, test_params_t *params)
{
    void *A_save = NULL, *work = NULL, *rwork = NULL;
    integer lwork, lrwork, liwork, *iwork;
    double exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_save, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_save, lda);

    /*
      Setting workspace size for non-lapacke interfaces
      LAPACKE interface handles workspace query internally
     */
    if(!FLA_IS_LAPACKE_INTERFACE(interfacetype))
    {
        /* gejsv APIs do not implement workspace query for realtype variants
           Hence, we need to set the workspace size based on the internal block size
           and the job parameters */
        if(FLA_IS_REALTYPE(datatype))
        {
            lwork = g_lwork <= 0 ? get_gejsv_lwork(joba, jobu, jobv, jobr, jobt, jobp, m, n)
                                 : g_lwork;
            liwork = g_liwork <= 0 ? fla_max(3, m + 3 * n) : g_liwork;
        }
        else
        {
            /* If any of the workspace length <= 0 then making workspace query */
            if(g_lwork <= 0 || g_lrwork <= 0 || g_liwork <= 0)
            {
                lwork = -1;
                lrwork = -1;
                liwork = -1;
                create_vector(datatype, &work, 2);
                create_vector(get_realtype(datatype), &rwork, 1);
                create_vector(INTEGER, (void **)&iwork, 1);

#if ENABLE_CPP_TEST
                if(interfacetype == LAPACK_CPP_TEST)
                {
                    invoke_cpp_gejsv(datatype, &joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, A,
                                     &lda, S, U, &ldu, V, &ldv, work, &lwork, rwork, &lrwork, iwork,
                                     info);
                }
                else
#endif
                /* call to gejsv API */
                {
                    invoke_gejsv(datatype, &joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, A,
                                 &lda, S, U, &ldu, V, &ldv, work, &lwork, rwork, &lrwork, iwork,
                                 info);
                }
                if(*info == 0)
                {
                    /* If lwork is not provided then set is */
                    if(g_lwork <= 0)
                    {
                        /* The optimal workspace size is saved in the second position */
                        if(datatype == COMPLEX)
                        {
                            lwork = (integer)((scomplex *)work)[1].real;
                        }
                        else if(datatype == DOUBLE_COMPLEX)
                        {
                            lwork = (integer)((dcomplex *)work)[1].real;
                        }
                    }
                    else
                    {
                        lwork = g_lwork;
                    }
                    /* If lrwork is not provided then set recommended value by API */
                    lrwork
                        = g_lrwork <= 0 ? get_work_value(get_realtype(datatype), rwork) : g_lrwork;
                    /* If liwork is not provided then set recommended value by API */
                    liwork = g_liwork <= 0 ? iwork[0] : g_liwork;
                }

                free_vector(work);
                free_vector(rwork);
                free_vector(iwork);
            }
            else
            {
                lwork = g_lwork;
                lrwork = g_lrwork;
                liwork = g_liwork;
            }
        }

        /* Allocating workspace memory */
        create_vector(datatype, &work, lwork);
        if(FLA_IS_COMPLEXTYPE(datatype))
        {
            create_vector(get_realtype(datatype), &rwork, lrwork);
        }
        create_vector(INTEGER, (void **)&iwork, liwork);
    }

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        copy_matrix(datatype, "full", m, n, A_save, lda, A, lda);

        if(FLA_IS_LAPACKE_INTERFACE(interfacetype))
        {
            exe_time
                = prepare_lapacke_gejsv_run(datatype, layout, joba, jobu, jobv, jobr, jobt, jobp, m,
                                            n, A, lda, S, U, ldu, V, ldv, stat, istat, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            invoke_cpp_gejsv(datatype, &joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, A, &lda, S,
                             U, &ldu, V, &ldv, work, &lwork, rwork, &lrwork, iwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            invoke_gejsv(datatype, &joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, A, &lda, S, U,
                         &ldu, V, &ldv, work, &lwork, rwork, &lrwork, iwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /* Copying stat values for Non-lapacke interfaces */
    if(!FLA_IS_LAPACKE_INTERFACE(interfacetype))
    {
        if(FLA_IS_COMPLEXTYPE(datatype))
        {
            copy_realtype_vector(datatype, 7, rwork, 1, stat, 1);
            copy_vector(INTEGER, 4, iwork, 1, istat, 1);
        }
        else
        {
            copy_realtype_vector(datatype, 7, work, 1, stat, 1);
            copy_vector(INTEGER, 3, iwork, 1, istat, 1);
        }
        /* Deallocate work arrays */
        free_vector(work);
        if(FLA_IS_COMPLEXTYPE(datatype))
        {
            free_vector(rwork);
        }
        free_vector(iwork);
    }

    free_matrix(A_save);
}

double prepare_lapacke_gejsv_run(integer datatype, int layout, char joba, char jobu, char jobv,
                                 char jobr, char jobt, char jobp, integer m, integer n, void *A,
                                 integer lda, void *S, void *U, integer ldu, void *V, integer ldv,
                                 void *stat, integer *istat, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldu_t = ldu;
    integer ldv_t = ldv;
    void *A_t = NULL, *U_t = NULL, *V_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_gejsv_lda, lda_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, m, row_major_gejsv_ldu, ldu_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_gejsv_ldv, ldv_t);

    A_t = A;
    U_t = U;
    V_t = V;
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m, n, &A_t, fla_max(n, lda_t));
        if(!same_char(jobu, 'N'))
        {
            create_matrix(datatype, layout, m, m, &U_t, fla_max(m, ldu_t));
        }
        if(!same_char(jobv, 'N'))
        {
            create_matrix(datatype, layout, n, n, &V_t, fla_max(n, ldv_t));
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m, n, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m, m, U, ldu, U_t, ldu_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, V, ldv, V_t, ldv_t);
    }

    exe_time = fla_test_clock();

    /* call LAPACKE gejsv API */
    *info = invoke_lapacke_gejsv(datatype, layout, joba, jobu, jobv, jobr, jobt, jobp, m, n, A_t,
                                 lda_t, S, U_t, ldu_t, V_t, ldv_t, stat, istat);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, m, n, A_t, lda_t, A, lda);
        if(!same_char(jobu, 'N'))
        {
            convert_matrix_layout(layout, datatype, m, m, U_t, ldu_t, U, ldu);
            free_matrix(U_t);
        }
        if(!same_char(jobv, 'N'))
        {
            convert_matrix_layout(layout, datatype, n, n, V_t, ldv_t, V, ldv);
            free_matrix(V_t);
        }
        /* free temporary buffers */
        free_matrix(A_t);
    }
    return exe_time;
}

void invoke_gejsv(integer datatype, char *joba, char *jobu, char *jobv, char *jobr, char *jobt,
                  char *jobp, integer *m, integer *n, void *A, integer *lda, void *S, void *U,
                  integer *ldu, void *V, integer *ldv, void *work, integer *lwork, void *rwork,
                  integer *lrwork, integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U, ldu, V, ldv,
                              work, lwork, iwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U, ldu, V, ldv,
                              work, lwork, iwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U, ldu, V, ldv,
                              work, lwork, rwork, lrwork, iwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, S, U, ldu, V, ldv,
                              work, lwork, rwork, lrwork, iwork, info);
            break;
        }
    }
}

void generate_gejsv_test_matrix(integer datatype, char joba, char jobu, char jobv, integer m,
                                integer n, void *A, integer lda, void *S)
{
    integer info = 0;
    void *V1, *V2;

    /* If either U or V is not calculated then not
       multiplying with V1 or V2 so that known
       singular values are generated which can later
       be validated. */

    integer need_known_singular_values = same_char(jobu, 'N') || same_char(jobu, 'W')
                                         || same_char(jobv, 'N') || same_char(jobv, 'W');
    integer small_values_len, offset;
    integer rdatatype = get_realtype(datatype);

    create_vector(datatype, &V1, m);
    create_vector(datatype, &V2, n);
    switch(joba)
    {
        case 'C':
        case 'E':
            /* generate matrix A with condition number <= 100 */
            create_svd_matrix(datatype, 'U', m, n, A, lda, S, GEJSV_VL, GEJSV_VU, i_zero, i_zero,
                              info);
            /* Generate rand vector of size n */
            if(!need_known_singular_values)
            {
                rand_vector(datatype, n, V2, 1, 0.0, 0.0, 'R');
                multiply_matrix_diag_vector(datatype, 'R', VECTOR_TYPE_COMPLEX, m, n, A, lda, V2,
                                            1);
            }
            break;
        case 'F':
        case 'G':
            /* generate matrix A with condition number <= 100 */
            create_svd_matrix(datatype, 'U', m, n, A, lda, S, GEJSV_VL, GEJSV_VU, i_zero, i_zero,
                              info);
            if(!need_known_singular_values)
            {
                /* Generate rand vector of size m with condition number <= 100000 */
                rand_vector(datatype, m, V1, 1, GEJSV_VEC_VL, GEJSV_VEC_VU, 'V');
                /* Generate rand vector of size n with condition number <= 100000 */
                rand_vector(datatype, n, V2, 1, GEJSV_VEC_VL, GEJSV_VEC_VU, 'V');
                multiply_matrix_diag_vector(datatype, 'L', VECTOR_TYPE_COMPLEX, m, n, A, lda, V1,
                                            1);
                multiply_matrix_diag_vector(datatype, 'R', VECTOR_TYPE_COMPLEX, m, n, A, lda, V2,
                                            1);
            }
            break;
        case 'A':
        case 'R':
            /* generate random matrix with uniform range */
            rand_vector(rdatatype, n, S, i_one, GEJSV_VL, GEJSV_VU, 'U');
            /* populate #small_vec_len values to a small random value */
            small_values_len = (integer)(GEJSV_SMALL_SV_PERCENT * n);
            offset = n - small_values_len;
            /* Populate randomly generated small values */
            rand_vector(rdatatype, small_values_len, get_ptr_at_offset(rdatatype, S, offset), i_one,
                        GEJSV_SMALL_SV_VL(datatype), GEJSV_SMALL_SV_VU(datatype), 'U');
            /* Create svd matrix with given singular values */
            create_svd_matrix(datatype, 'N', m, n, A, lda, S, 0.0, 0.0, i_zero, i_zero, info);
            break;
    }
    free_vector(V1);
    free_vector(V2);
}

integer get_gejsv_lwork(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, integer m,
                        integer n)
{
    integer NBsiz = 128; // Taking 128 as upper bound for NB
    integer lwork = fla_max(7, 2 * m + n);
    integer want_u = same_char(jobu, 'U') || same_char(jobu, 'F');
    integer want_v = same_char(jobv, 'V') || same_char(jobv, 'J');
    integer want_sce = same_char(joba, 'E') || same_char(joba, 'G');
    if(!want_u && !want_v)
        lwork = fla_max(lwork, 3 * n + (n + 1) * NBsiz);
    if(!want_u && !want_v && want_sce)
        lwork = fla_max(lwork, n * n + 4 * n);
    if(!want_u && want_v)
        lwork = fla_max(lwork, 3 * n + (n + 1) * NBsiz);
    if(want_u && !want_v)
        lwork = fla_max(lwork, 3 * n + (n + 1) * NBsiz);
    if(same_char(jobu, 'F') && !want_v)
        lwork = fla_max(lwork, n + m * NBsiz);
    if(want_u && same_char(jobv, 'V'))
        lwork = fla_max(lwork, 6 * n + 2 * n * n);
    if(want_u && same_char(jobv, 'J'))
        lwork = fla_max(lwork, fla_max(4 * n + n * n, 2 * n + n * n + 6));
    if(want_u && want_v)
        lwork = fla_max(lwork, n + m * NBsiz);
    return lwork;
}

void store_gejsv_outputs(void *filename, integer datatype, char joba, char jobu, char jobv,
                         char jobr, char jobt, char jobp, integer m, integer n, void *A,
                         integer lda, void *S, void *U, integer ldu, void *V, integer ldv,
                         integer g_lwork, integer g_lrwork, integer g_liwork, void *params)
{
    /* Create and open a file for storing Ground truth*/
    FLA_OPEN_GT_FILE_STORE

    FLA_STORE_BRT_MATRIX(datatype, m, n, A, lda)
    FLA_STORE_BRT_VECTOR(get_realtype(datatype), n, S)
    if(!same_char(jobu, 'N'))
    {
        FLA_STORE_BRT_MATRIX(datatype, m, n, U, ldu)
    }
    if(!same_char(jobv, 'N'))
    {
        FLA_STORE_BRT_MATRIX(datatype, n, n, V, ldv)
    }

    FLA_CLOSE_GT_FILE_STORE
}

integer check_bit_reproducibility_gejsv(void *filename, integer datatype, char joba, char jobu,
                                        char jobv, char jobr, char jobt, char jobp, integer m,
                                        integer n, void *A, integer lda, void *S, void *U,
                                        integer ldu, void *V, integer ldv, integer g_lwork,
                                        integer g_lrwork, integer g_liwork, void *params)
{
    /* Open the file for reading Ground truth */
    FLA_OPEN_GT_FILE_READ

    FLA_VERIFY_BRT_MATRIX(datatype, m, n, A, lda)
    FLA_VERIFY_BRT_VECTOR(get_realtype(datatype), n, S)
    if(!same_char(jobu, 'N'))
    {
        FLA_VERIFY_BRT_MATRIX(datatype, m, n, U, ldu)
    }
    if(!same_char(jobv, 'N'))
    {
        FLA_VERIFY_BRT_MATRIX(datatype, n, n, V, ldv)
    }

    fclose(gt_file);
    return 1;
}