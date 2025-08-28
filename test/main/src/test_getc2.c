/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#define GETC2_VL 0.1
#define GETC2_VU 10
extern double perf;
extern double time_min;

/* Local prototypes */
void fla_test_getc2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_getc2_run(integer n_A, void *A, integer lda, integer *ipiv, integer *jpiv,
                       integer datatype, integer *info, test_params_t *params);
void invoke_getc2(integer datatype, integer *n, void *a, integer *lda, 
                  integer *ipiv, integer *jpiv, integer *info);

void fla_test_getc2(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LU factorization with complete pivoting (GETC2)";
    char *front_str = "GETC2";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    /* Config mode */
    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_getc2_experiment);
        tests_not_run = 0;
    }
    /* CLI mode: Parse last arg */
    if(argc == 7)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[6]);
    }
    /* CLI mode: Parse other args */
    if(argc >= 6 && argc <= 7)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_getc2_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for getc2\n");
        printf("./<EXE> getc2 <precisions - sdcz> <N> <LDA> <repeats>\n");
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

void fla_test_getc2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, lda, info = 0;
    void *IPIV = NULL, *JPIV = NULL, *A = NULL, *A_test = NULL, *s_test = NULL;
    char range = 'U';
    double residual, err_thresh;

    /* GETC2 only works with square matrices */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_vector(INTEGER, &IPIV, n);      /* Row permutations */
    create_vector(INTEGER, &JPIV, n);      /* Column permutations */
    create_realtype_vector(datatype, &s_test, n);

    /* Initialize the test matrices*/
    if(g_ext_fptr != NULL)
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, n, n, A, lda, s_test, GETC2_VL, GETC2_VU, i_zero, i_zero,
                          info);
    }
    
    /* Save the original matrix*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    /* call to API */
    prepare_getc2_run(n, A_test, lda, IPIV, JPIV, datatype, &info, params);

    /* performance computation*/
    perf = (2.0 / 3.0) * n * n * n / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
    free_vector(JPIV);
    free_vector(s_test);
}

void prepare_getc2_run(integer n_A, void *A, integer lda, integer *IPIV, integer *JPIV,
                       integer datatype, integer *info, test_params_t *params)
{
    void *A_save;
    double exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", n_A, n_A, A, lda, A_save, lda);

        exe_time = fla_test_clock();
        /* Call LAPACK getc2 API */
        invoke_getc2(datatype, &n_A, A_save, &lda, IPIV, JPIV, info);
        exe_time = fla_test_clock() - exe_time;

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /*  Save the AFACT to matrix A */
    copy_matrix(datatype, "full", n_A, n_A, A_save, lda, A, lda);
    free_matrix(A_save);
}

/*
 *  GETC2_API calls LAPACK interface of
 *  LU factorization with complete pivoting - getc2
 *  */
void invoke_getc2(integer datatype, integer *n, void *a, integer *lda, 
                  integer *ipiv, integer *jpiv, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgetc2(n, a, lda, ipiv, jpiv, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgetc2(n, a, lda, ipiv, jpiv, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgetc2(n, a, lda, ipiv, jpiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgetc2(n, a, lda, ipiv, jpiv, info);
            break;
        }
    }
}