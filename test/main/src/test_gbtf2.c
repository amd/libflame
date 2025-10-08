/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

extern double perf;
extern double time_min;

/* Local prototypes */
void fla_test_gbtf2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_gbtf2_run(integer m_A, integer n_A, integer kl, integer ku, void *ab, integer ldab,
                       integer *ipiv, integer datatype, integer *info, test_params_t *params);
void invoke_gbtf2(integer datatype, integer *m, integer *n, integer *kl, integer *ku, void *ab,
                  integer *ldab, integer *ipiv, integer *info);

void fla_test_gbtf2(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LU factorization of banded matrix (GBTF2)";
    char *front_str = "GBTF2";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    /* Config mode */
    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gbtf2_experiment);
        tests_not_run = 0;
    }
    /* CLI mode: Parse last arg */
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }
    /* CLI mode: Parse other args */
    if((argc == 9) || (argc == 10))
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].kl = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ku = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldab = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_gbtf2_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gbtf2\n");
        printf("./<EXE> gbtf2 <precisions - sdcz> <M> <N> <KL> <KU> <LDAB> <repeats>\n");
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

void fla_test_gbtf2_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, kl, ku, ldab;
    integer info = 0;
    void *IPIV;
    void *AB, *AB_test;
    double residual, err_thresh;
    void *filename = NULL;

    /* Determine the dimensions*/
    m = p_cur;
    n = q_cur;
    kl = params->lin_solver_paramslist[pci].kl;
    ku = params->lin_solver_paramslist[pci].ku;
    ldab = params->lin_solver_paramslist[pci].ldab;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(ldab == -1)
        {
            ldab = 2 * kl + ku + 1;
        }
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &AB, ldab);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &AB_test, ldab);
    create_vector(INTEGER, &IPIV, fla_min(m, n));
    reset_vector(INTEGER, IPIV, fla_min(m, n), 1);

    if(!FLA_BRT_VERIFICATION_RUN)
    {
        /* Initialize the test matrices*/
        if(g_ext_fptr != NULL)
        {
            /* Initialize input matrix with custom data from file */
            init_matrix(datatype, AB, m, n, ldab, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            /* Initialize & convert random band matrix into band storage as per API need */
            rand_band_storage_matrix(datatype, m, n, kl, ku, AB, ldab);
        }
    }
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, m, n, AB, ldab, "ddddd", m, n, kl, ku, ldab);

    /* Save the original matrix*/
    copy_matrix(datatype, "full", ldab, n, AB, ldab, AB_test, ldab);

    /* call to API */
    prepare_gbtf2_run(m, n, kl, ku, AB_test, ldab, IPIV, datatype, &info, params);

    /* performance computation */
    perf = (2.0 * n * (kl + ku + 1) * kl) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        perf *= 4.0;
    }

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    IF_FLA_BRT_VALIDATION(m, n,
                          store_outputs_base(filename, params, 1, 1, datatype, m, n, AB_test, ldab,
                                             INTEGER, fla_min(m, n), IPIV),
                          FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh),
                          check_reproducibility_base(filename, params, 1, 1, datatype, m, n,
                                                     AB_test, ldab, INTEGER, fla_min(m, n), IPIV))
    else
    {
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

    /* Free up the buffers */
free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(AB);
    free_matrix(AB_test);
    free_vector(IPIV);
}

void prepare_gbtf2_run(integer m_A, integer n_A, integer kl, integer ku, void *AB, integer ldab,
                       integer *IPIV, integer datatype, integer *info, test_params_t *params)
{
    void *AB_save;
    double exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &AB_save, ldab);
    copy_matrix(datatype, "full", ldab, n_A, AB, ldab, AB_save, ldab);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", ldab, n_A, AB, ldab, AB_save, ldab);

        exe_time = fla_test_clock();
        /* Call LAPACK gbtf2 API */
        invoke_gbtf2(datatype, &m_A, &n_A, &kl, &ku, AB_save, &ldab, IPIV, info);
        exe_time = fla_test_clock() - exe_time;

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /*  Save the ABFACT to matrix AB */
    copy_matrix(datatype, "full", ldab, n_A, AB_save, ldab, AB, ldab);
    free_matrix(AB_save);
}

/*
 *  GBTF2_API calls LAPACK interface of
 *  LU factorization with level 2 BLAS - gbtf2
 *  */
void invoke_gbtf2(integer datatype, integer *m, integer *n, integer *kl, integer *ku, void *ab,
                  integer *ldab, integer *ipiv, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }
    }
}
