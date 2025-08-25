/*
    Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/* GESVDX API */

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_gesvdx_lda;
integer row_major_gesvdx_ldu;
integer row_major_gesvdx_ldvt;

/* Local prototypes */
void fla_test_gesvdx_experiment(char *tst_api, test_params_t *params, integer datatype,
                                integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                integer einfo);
void prepare_gesvdx_run(char *jobu, char *jobvt, char *range, integer m_A, integer n_A, void *A,
                        integer lda, void *vl, void *vu, integer il, integer iu, integer *ns,
                        void *s, void *U, integer ldu, void *V, integer ldvt, integer datatype,
                        integer *info, integer interfacetype, int matrix_layout,
                        test_params_t *params);
void invoke_gesvdx(integer datatype, char *jobu, char *jobvt, char *range, integer *m, integer *n,
                   void *a, integer *lda, void *vl, void *vu, integer *il, integer *iu, integer *ns,
                   void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                   integer *lwork, integer *iwork, void *rwork, integer *info);
double prepare_lapacke_gesvdx_run(integer datatype, int matrix_layout, char *jobu, char *jobvt,
                                  char *range, integer m_A, integer n_A, void *A, integer lda,
                                  void *vl, void *vu, integer il, integer iu, integer *ns, void *s,
                                  void *U, integer ldu, void *V, integer ldvt, void *work,
                                  integer *lwork, integer *iwork, void *rwork, integer *info);

/* Helper functions for Bit reproducibility tests */
void store_gesvdx_outputs(void *filename, integer datatype, char jobu, char jobvt, char range,
                          integer m, integer n, void *A, integer lda, integer vl, integer vu,
                          integer il, integer iu, void *s, void *U, integer ldu, void *V,
                          integer ldvt, integer g_lwork, void *params);
integer check_bit_reproducibility_gesvdx(void *filename, integer datatype, char jobu, char jobvt,
                                         char range, integer m, integer n, void *A, integer lda,
                                         integer vl, integer vu, integer il, integer iu, void *s,
                                         void *U, integer ldu, void *V, integer ldvt,
                                         integer g_lwork, void *params);

void fla_test_gesvdx(integer argc, char **argv, test_params_t *params)
{
    srand(10);
    char *op_str = "Singular value decomposition";
    char *front_str = "GESVDX";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        g_lwork = -1;
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, SVD, fla_test_gesvdx_experiment);
        tests_not_run = 0;
    }
    if(argc == 18)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[17]);
    }
    if(argc >= 17 && argc <= 18)
    {
        integer i, num_types, N, M;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;
        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->svd_paramslist[0].jobu_gesvdx = argv[3][0];
        params->svd_paramslist[0].jobvt_gesvdx = argv[4][0];
        params->svd_paramslist[0].range_gesvdx = argv[5][0];
        M = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->svd_paramslist[0].vl = strtod(argv[9], &endptr);
        params->svd_paramslist[0].vu = strtod(argv[10], &endptr);
        params->svd_paramslist[0].il = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        params->svd_paramslist[0].iu = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gesvdx_lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_gesvdx_ldu = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
            row_major_gesvdx_ldvt = strtoimax(argv[14], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].lda = N;
            params->svd_paramslist[0].ldu = N;
            params->svd_paramslist[0].ldvt = N;
        }
        else
        {
            params->svd_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldu = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldvt = strtoimax(argv[14], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[15], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[16], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_gesvdx_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }
    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gesvdx\n");
        printf("./<EXE> gesvdx <precisions - sdcz>  <JOBU> <JOBVT> <RANGE> <M> <N> <LDA> <VL> <VU> "
               "<IL> <IU> <LDU> <LDVT> <LWORK> <repeats>\n");
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

void fla_test_gesvdx_experiment(char *tst_api, test_params_t *params, integer datatype,
                                integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                integer einfo)
{
    char jobu, jobvt, range;
    integer m, n, lda;
    void *vl = NULL, *vu = NULL;
    float d_vl, d_vu;
    integer il, iu, ns, ldu, ldvt;
    integer info = 0;
    void *A = NULL, *U = NULL, *V = NULL, *s = NULL, *A_test = NULL, *s_test = NULL, *scal = NULL;
    double residual, err_thresh;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions. */
    jobu = params->svd_paramslist[pci].jobu_gesvdx;
    jobvt = params->svd_paramslist[pci].jobvt_gesvdx;
    range = params->svd_paramslist[pci].range_gesvdx;
    err_thresh = params->svd_paramslist[pci].svd_threshold;

    m = p_cur;
    n = q_cur;

    lda = params->svd_paramslist[pci].lda;
    d_vl = params->svd_paramslist[pci].vl;
    d_vu = params->svd_paramslist[pci].vu;
    il = params->svd_paramslist[pci].il;
    iu = params->svd_paramslist[pci].iu;
    ldu = params->svd_paramslist[pci].ldu;
    ldvt = params->svd_paramslist[pci].ldvt;
    integer min_m_n = fla_min(m, n);
    /* If leading dimensions = -1, set them to default value
       LDA >= max(1,M). */
    if(g_config_data)
    {

        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
        /* If LDU >= 1; if
          JOBU = 'V', LDU >= M. */
        if(ldu == -1)
        {
            if(same_char(jobu, 'V'))
            {
                ldu = m;
            }
            else
            {
                ldu = 1;
            }
        }
        /* If LDVT >= 1; if
          JOBVT = 'V', LDVT >= NS */
        if(ldvt == -1)
        {
            if(same_char(jobvt, 'V'))
            {
                if(same_char(range, 'I'))
                {
                    ldvt = iu - il + 1;
                }
                else
                {
                    ldvt = min_m_n;
                }
            }
            else
            {
                ldvt = 1;
            }
        }
    }

    /* Create datatype for void */
    create_realtype_vector(datatype, &vl, 1);
    create_realtype_vector(datatype, &vu, 1);

    /* Assign datatype for void */
    assign_value(get_realtype(datatype), vl, d_vl, d_zero);
    assign_value(get_realtype(datatype), vu, d_vu, d_zero);

    /* Create input matrix parameters. */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_realtype_vector(datatype, &s_test, min_m_n);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &U, ldu);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &V, ldvt);
    create_realtype_vector(datatype, &s, min_m_n);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
        create_vector(get_realtype(datatype), &scal, 1);

    /* This code path is run to generate the matrix to be passed to the API. This is the default
     * input generation logic accessed both when BRT is run in Ground truth mode and for non BRT
     * Test cases. For verification runs the input is loaded from the input generated during Ground
     * truth run */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
        {
            /* Initialize input matrix with custom data */
            init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            if(FLA_OVERFLOW_UNDERFLOW_TEST)
            {
                /* Range 'V' treated as 'A' only in case of overflow and underflow */
                if(same_char(range, 'V'))
                    range = 'A';
            }

            /* Generate matrix A by known singular value */
            create_svd_matrix(datatype, range, m, n, A, lda, s_test, d_vl, d_vu, il, iu, info);
            if(FLA_OVERFLOW_UNDERFLOW_TEST)
            {
                /* Initializing matrix with values around overflow underflow */
                init_matrix_overflow_underflow_svdx(datatype, m, n, A, lda, params->imatrix_char,
                                                    scal);
            }
        }
    }

    /* This macro is used in the BRT test cases for the following purposes:
     *    - In the Ground truth runs (BRT_char => G, F), the output is stored in a file for future
     * reference
     *    - In the verification runs (BRT_char => V, M), the output is loaded from the file and
     * passed as input to the API
     * */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, m, n, A, lda, "cccdddffddddd", jobu, jobvt, range, m, n,
                                 lda, d_vl, d_vu, il, iu, ldu, ldvt, g_lwork)

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_gesvdx_run(&jobu, &jobvt, &range, m, n, A_test, lda, vl, vu, il, iu, &ns, s, U, ldu, V,
                       ldvt, datatype, &info, interfacetype, layout, params);

    /* Performance Computation
     * Singular values only, 4mn^2 - 4n^3/3 flops
     * Singular values and some singular vectors U (m x n) and V (n x n), 14mn^2 + 8n^3 flops */

    if(same_char(jobvt, 'N') && same_char(jobu, 'N'))
    {
        if(m >= n)
            perf = (double)((4.0 * m * n * n) - ((4.0 * n * n * n) / 3.0)) / time_min
                   / FLOPS_PER_UNIT_PERF;
        else
            perf = (double)((4.0 * n * m * m) - ((4.0 * m * m * m) / 3.0)) / time_min
                   / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        if(m >= n)
            perf
                = (double)((14.0 * m * n * n) + (8.0 * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
        else
            perf
                = (double)((14.0 * n * m * m) + (8.0 * m * m * m)) / time_min / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* Output Validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    /* Bit reproducibility tests path
     * This path is taken when BRT is enabled.
     *     - In the Ground truth runs (BRT_char => G, F), the output is stored in a file and the
     * default validation function is called
     *     - In the verification runs (BRT_char => V, M), the output is loaded from the file and
     * compared with the generated output
     *  */
    IF_FLA_BRT_VALIDATION(
        m, n,
        store_gesvdx_outputs(filename, datatype, jobu, jobvt, range, m, n, A, lda, d_vl, d_vu, il,
                             iu, s, U, ldu, V, ldvt, g_lwork, params),
        validate_gesvdx(tst_api, &jobu, &jobvt, range, m, n, A, A_test, lda, vl, vu, il, iu, ns, s,
                        s_test, U, ldu, V, ldvt, datatype, residual, g_ext_fptr, scal,
                        params->imatrix_char, params),
        check_bit_reproducibility_gesvdx(filename, datatype, jobu, jobvt, range, m, n, A, lda, d_vl,
                                         d_vu, il, iu, s, U, ldu, V, ldvt, g_lwork, params))
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_gesvdx(tst_api, &jobu, &jobvt, range, m, n, A, A_test, lda, vl, vu, il, iu, ns, s,
                        s_test, U, ldu, V, ldvt, datatype, residual, g_ext_fptr, scal,
                        params->imatrix_char, params);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if((!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char))
           && (!check_extreme_value(datatype, min_m_n, i_one, s_test, i_one, params->imatrix_char))
           && (!check_extreme_value(datatype, m, n, U, ldu, params->imatrix_char))
           && (!check_extreme_value(datatype, m, n, V, ldvt, params->imatrix_char)))
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
free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(A);
    free_matrix(A_test);
    free_vector(s);
    free_vector(vl);
    free_vector(vu);
    free_matrix(U);
    free_matrix(V);
    free_vector(s_test);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
        free_vector(scal);
}

void prepare_gesvdx_run(char *jobu, char *jobvt, char *range, integer m_A, integer n_A, void *A,
                        integer lda, void *vl, void *vu, integer il, integer iu, integer *ns,
                        void *s, void *U, integer ldu, void *V, integer ldvt, integer datatype,
                        integer *info, integer interfacetype, int layout, test_params_t *params)
{
    integer min_m_n, max_m_n;
    void *A_save, *s_test;
    void *work, *rwork;
    void *U_test, *V_test;
    integer lwork, lrwork;
    void *iwork = NULL;
    double exe_time;

    min_m_n = fla_min(m_A, n_A);
    max_m_n = fla_max(m_A, n_A);

    /* Make a copy of the input matrix A. Same input values will be passed in each itertaion */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    /* Get rwork array size since it is not depedent on internal blocks */
    lrwork = fla_max((5 * min_m_n * min_m_n + 5 * min_m_n),
                     (2 * max_m_n * min_m_n + 2 * min_m_n * min_m_n + min_m_n));

    /* Make a workspace query the first time through. This will provide us with and ideal workspace
     * size based on an internal block size. */
    if(g_lwork == -1)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* call gesvdx API */
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_gesvdx(datatype, jobu, jobvt, range, &m_A, &n_A, NULL, &lda, vl, vu, &il,
                              &iu, ns, NULL, NULL, &ldu, NULL, &ldvt, work, &lwork, iwork, NULL,
                              info);
        }
        else
#endif
        {
            invoke_gesvdx(datatype, jobu, jobvt, range, &m_A, &n_A, NULL, &lda, vl, vu, &il, &iu,
                          ns, NULL, NULL, &ldu, NULL, &ldvt, work, &lwork, iwork, NULL, info);
        }
        if(*info == 0)
        {
            /* Get the work size */
            lwork = get_work_value(datatype, work);
        }
        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }
    create_vector(INTEGER, &iwork, (12 * min_m_n));
    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrix A value and allocate memory to output buffers for each iteration */
        copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);
        create_matrix(datatype, LAPACK_COL_MAJOR, m_A, m_A, &U_test, ldu);
        create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &V_test, ldvt);
        create_realtype_vector(datatype, &s_test, min_m_n);
        create_vector(datatype, &work, lwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, lrwork);
        else
            rwork = NULL;

        /* Check if LAPACKE is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gesvdx_run(datatype, layout, jobu, jobvt, range, m_A, n_A, A,
                                                  lda, vl, vu, il, iu, ns, s_test, U_test, ldu,
                                                  V_test, ldvt, work, &lwork, iwork, rwork, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* call CPP gesvdx API  */
            invoke_cpp_gesvdx(datatype, jobu, jobvt, range, &m_A, &n_A, A, &lda, vl, vu, &il, &iu,
                              ns, s_test, U_test, &ldu, V_test, &ldvt, work, &lwork, iwork, rwork,
                              info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();

            /* call to API */
            invoke_gesvdx(datatype, jobu, jobvt, range, &m_A, &n_A, A, &lda, vl, vu, &il, &iu, ns,
                          s_test, U_test, &ldu, V_test, &ldvt, work, &lwork, iwork, rwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO

        /* Make a copy of the output buffers. This is required to validate the API functionality. */
        copy_matrix(datatype, "full", m_A, m_A, U_test, ldu, U, ldu);
        copy_matrix(datatype, "full", min_m_n, n_A, V_test, ldvt, V, ldvt);
        copy_realtype_vector(datatype, min_m_n, s_test, 1, s, 1);

        /* Free up the output buffers */
        free_vector(work);
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            free_vector(rwork);
        }
        free_matrix(U_test);
        free_matrix(V_test);
        free_vector(s_test);
    }

    free_vector(iwork);
    free_matrix(A_save);
}

double prepare_lapacke_gesvdx_run(integer datatype, int layout, char *jobu, char *jobvt,
                                  char *range, integer m_A, integer n_A, void *A, integer lda,
                                  void *vl, void *vu, integer il, integer iu, integer *ns, void *s,
                                  void *U, integer ldu, void *V, integer ldvt, void *work,
                                  integer *lwork, integer *iwork, void *rwork, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldu_t = ldu;
    integer ldvt_t = ldvt;
    void *A_t = NULL, *U_t = NULL, *V_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n_A, row_major_gesvdx_lda, lda_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, m_A, row_major_gesvdx_ldu, ldu_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n_A, row_major_gesvdx_ldvt, ldvt_t);

    A_t = A;
    U_t = U;
    V_t = V;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m_A, n_A, &A_t, fla_max(n_A, lda_t));
        if(!same_char(*jobu, 'N'))
        {
            create_matrix(datatype, layout, m_A, m_A, &U_t, fla_max(m_A, ldu_t));
        }
        if(!same_char(*jobvt, 'N'))
        {
            create_matrix(datatype, layout, n_A, n_A, &V_t, fla_max(n_A, ldvt_t));
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m_A, n_A, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* call LAPACKE gesvd API */
    *info = invoke_lapacke_gesvdx(datatype, layout, *jobu, *jobvt, *range, m_A, n_A, A_t, lda_t, vl,
                                  vu, il, iu, ns, s, U_t, ldu_t, V_t, ldvt_t, work);

    exe_time = fla_test_clock() - exe_time;
    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */
    if(layout == LAPACK_ROW_MAJOR)
    {
        convert_matrix_layout(layout, datatype, m_A, n_A, A_t, lda_t, A, lda);
        if(!same_char(*jobu, 'N') && !same_char(*jobu, 'O'))
        {
            convert_matrix_layout(layout, datatype, m_A, m_A, U_t, ldu_t, U, ldu);
            free_matrix(U_t);
        }
        if(!same_char(*jobvt, 'N') && !same_char(*jobvt, 'O'))
        {
            convert_matrix_layout(layout, datatype, n_A, n_A, V_t, ldvt_t, V, ldvt);
            free_matrix(V_t);
        }
        free_matrix(A_t);
    }

    return exe_time;
}

void invoke_gesvdx(integer datatype, char *jobu, char *jobvt, char *range, integer *m, integer *n,
                   void *a, integer *lda, void *vl, void *vu, integer *il, integer *iu, integer *ns,
                   void *s, void *u, integer *ldu, void *vt, integer *ldvt, void *work,
                   integer *lwork, integer *iwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt,
                               ldvt, work, lwork, iwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt,
                               ldvt, work, lwork, iwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt,
                               ldvt, work, lwork, rwork, iwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt,
                               ldvt, work, lwork, rwork, iwork, info);
            break;
        }
    }
}

void store_gesvdx_outputs(void *filename, integer datatype, char jobu, char jobvt, char range,
                          integer m, integer n, void *A, integer lda, integer vl, integer vu,
                          integer il, integer iu, void *s, void *U, integer ldu, void *V,
                          integer ldvt, integer g_lwork, void *params)
{
    /* Create and open a file for storing Ground truth*/
    FLA_OPEN_GT_FILE_STORE

    integer ns = fla_min(m, n);
    if(same_char(range, 'I'))
    {
        ns = iu - il + 1;
    }

    /* Store the ground truth data */
    FLA_STORE_BRT_MATRIX(datatype, m, n, A, lda)
    FLA_STORE_BRT_VECTOR(get_realtype(datatype), ns, s)
    if(!same_char(jobu, 'N'))
    {
        FLA_STORE_BRT_MATRIX(datatype, m, ns, U, ldu)
    }
    if(!same_char(jobvt, 'N'))
    {
        FLA_STORE_BRT_MATRIX(datatype, ns, n, V, ldvt)
    }

    fclose(gt_file);
}
integer check_bit_reproducibility_gesvdx(void *filename, integer datatype, char jobu, char jobvt,
                                         char range, integer m, integer n, void *A, integer lda,
                                         integer vl, integer vu, integer il, integer iu, void *s,
                                         void *U, integer ldu, void *V, integer ldvt,
                                         integer g_lwork, void *params)
{
    /* Open the file for reading Ground truth */
    FLA_OPEN_GT_FILE_READ

    integer ns = fla_min(m, n);
    if(same_char(range, 'I'))
    {
        ns = iu - il + 1;
    }

    /* Load stored GT and verify with current API outputs */
    FLA_VERIFY_BRT_MATRIX(datatype, m, n, A, lda)
    FLA_VERIFY_BRT_VECTOR(get_realtype(datatype), ns, s)
    if(!same_char(jobu, 'N'))
    {
        FLA_VERIFY_BRT_MATRIX(datatype, m, ns, U, ldu)
    }
    if(!same_char(jobvt, 'N'))
    {
        FLA_VERIFY_BRT_MATRIX(datatype, ns, n, V, ldvt)
    }

    fclose(gt_file);
    return 1;
}