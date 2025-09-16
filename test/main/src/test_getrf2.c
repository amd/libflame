/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#define GETRF2_VL 0.1
#define GETRF2_VU 10
extern double perf;
extern double time_min;

/* Local prototypes */
void fla_test_getrf2_experiment(char *tst_api, test_params_t *params, integer datatype,
                                integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                integer einfo);
void prepare_getrf2_run(integer m_A, integer n_A, void *A, integer lda, integer *ipiv,
                        integer datatype, integer *info, test_params_t *params);
void invoke_getrf2(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv,
                   integer *info);

void fla_test_getrf2(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LU factorization (GETRF2)";
    char *front_str = "GETRF2";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    /* Config mode */
    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_getrf2_experiment);
        tests_not_run = 0;
    }
    /* CLI mode: Parse last arg */
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    /* CLI mode: Parse other args */
    if(argc >= 7 && argc <= 8)
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

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

                /* Check for duplicate datatype presence */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_getrf2_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for getrf2\n");
        printf("./<EXE> getrf2 <precisions - sdcz> <M> <N> <LDA> <repeats>\n");
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

void fla_test_getrf2_experiment(char *tst_api, test_params_t *params, integer datatype,
                                integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                integer einfo)
{
    integer m, n, lda, info = 0;
    void *IPIV = NULL, *A = NULL, *A_test = NULL, *s_test = NULL;
    char range = 'U';
    double residual, err_thresh;
    void *filename = NULL;

    /* Determine the dimensions*/
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_vector(INTEGER, &IPIV, fla_min(m, n));
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);

    /* This code path is run to generate the matrix to be passed to the API. This is the default
     * input generation logic accessed both when BRT is run in Ground truth mode and for non BRT
     * Test cases. For verification runs the input is loaded from the input generated during Ground
     * truth run */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        /* Initialize the test matrices*/
        if(g_ext_fptr != NULL)
        {
            init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            /* Generate input matrix with condition number <= 100 */
            create_realtype_vector(datatype, &s_test, fla_min(m, n));
            create_svd_matrix(datatype, range, m, n, A, lda, s_test, GETRF2_VL, GETRF2_VU, i_zero,
                              i_zero, info);
        }
    }
    /* This macro is used in the BRT test cases for the following purposes:
     *    - In the Ground truth runs (BRT_char => G, F), the output is stored in a file for future
     * reference
     *    - In the verification runs (BRT_char => V, M), the output is loaded from the file and
     * passed as input to the API
     * */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, m, n, A, lda, "ddd", m, n, lda)

    /* Save the original matrix*/
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    /* call to API */
    prepare_getrf2_run(m, n, A_test, lda, IPIV, datatype, &info, params);

    /* performance computation */
    if(m == n)
    {
        perf = (2.0 / 3.0) * n * n * n / time_min / FLOPS_PER_UNIT_PERF;
    }
    else if(m > n)
    {
        perf = (1.0 / 3.0) * n * n * (3.0 * m - n) / time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        perf = (1.0 / 3.0) * m * m * (3.0 * n - m) / time_min / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    /* Bit reproducibility tests path
     * This path is taken when BRT is enabled.
     *     - In the Ground truth runs (BRT_char => G, F), the output is stored in a file and the
     * default validation function is called
     *     - In the verification runs (BRT_char => V, M), the output is loaded from the file and
     * compared with the generated output
     *  */
    IF_FLA_BRT_VALIDATION(m, n,
                          store_outputs_base(filename, params, 1, 1, datatype, m, n, A_test, lda,
                                             INTEGER, fla_min(m, n), IPIV),
                          FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh),
                          check_reproducibility_base(filename, params, 1, 1, datatype, m, n, A_test,
                                                     lda, INTEGER, fla_min(m, n), IPIV))
    /* API functionality validation */
    else
    {
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

    /* Free up the buffers */
free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
    free_vector(s_test);
}

void prepare_getrf2_run(integer m_A, integer n_A, void *A, integer lda, integer *IPIV,
                        integer datatype, integer *info, test_params_t *params)
{
    void *A_save;
    double exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

        exe_time = fla_test_clock();
        /* Call LAPACK getrf2 API */
        invoke_getrf2(datatype, &m_A, &n_A, A_save, &lda, IPIV, info);
        exe_time = fla_test_clock() - exe_time;

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /*  Save the AFACT to matrix A */
    copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);
    free_matrix(A_save);
}

/*
 *  GETRF2_API calls LAPACK interface of
 *  LU factorization with level 2 BLAS - getrf2
 *  */
void invoke_getrf2(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv,
                   integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgetrf2(m, n, a, lda, ipiv, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgetrf2(m, n, a, lda, ipiv, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgetrf2(m, n, a, lda, ipiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgetrf2(m, n, a, lda, ipiv, info);
            break;
        }
    }
}