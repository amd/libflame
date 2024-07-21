/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_gbtrs_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_gbtrs_run(char trans, integer n_A, integer kl, integer ku, integer nrhs, void *ab,
                       integer ldab, integer *ipiv, void *b, integer ldb, integer datatype,
                       integer n_repeats, double *time_min_, integer *info);
void invoke_gbtrs(integer datatype, char *trans, integer *n, integer *kl, integer *ku,
                  integer *nrhs, void *ab, integer *ldab, integer *ipiv, void *b, integer *ldb,
                  integer *info);
void invoke_gbtrf(integer datatype, integer *m, integer *n, integer *kl, integer *ku, void *ab,
                  integer *ldab, integer *ipiv, integer *info);

void fla_test_gbtrs(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Linear solver of banded matrix";
    char *front_str = "GBTRS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gbtrs_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if((argc == 11) || (argc == 12))
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].transr = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);

        params->lin_solver_paramslist[0].kl = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ku = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldab = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->lin_solver_paramslist[0].solver_threshold = CLI_NORM_THRESH;

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
                fla_test_gbtrs_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, SQUARE_INPUT, N, N, residual,
                                      params->lin_solver_paramslist[0].solver_threshold, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gbtrs\n");
        printf("./<EXE> gbtrs <precisions - sdcz> <TRANS> <N> <KL> <KU> <NRHS> <LDAB> <LDB> "
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
}

void fla_test_gbtrs_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer n, kl, ku, nrhs, ldab, ldb;
    integer info = 0;
    char trans;
    void *IPIV;
    void *AB, *AB_test;
    void *B, *X, *A;
    double time_min = 1e9;

    /* Determine the dimensions*/
    trans = params->lin_solver_paramslist[pci].transr;
    n = p_cur;
    kl = params->lin_solver_paramslist[pci].kl;
    ku = params->lin_solver_paramslist[pci].ku;
    nrhs = params->lin_solver_paramslist[pci].nrhs;
    ldab = params->lin_solver_paramslist[pci].ldab;
    ldb = params->lin_solver_paramslist[pci].ldb;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(ldab == -1)
        {
            ldab = 2 * kl + ku + 1;
        }
        if(ldb == -1)
        {
            ldb = n;
        }
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, matrix_layout, n, n, &AB, ldab);
    create_matrix(datatype, matrix_layout, n, n, &AB_test, ldab);
    create_vector(INTEGER, &IPIV, n);
    create_matrix(datatype, matrix_layout, n, nrhs, &B, ldb);
    create_matrix(datatype, matrix_layout, n, nrhs, &X, ldb);

    /* Initialize the test matrices*/
    if(g_ext_fptr != NULL || FLA_EXTREME_CASE_TEST)
    {
        /* Initialize input matrix with custom data from file */
        init_matrix(datatype, AB, n, n, ldab, g_ext_fptr, params->imatrix_char);
        init_matrix(datatype, IPIV, 1, n, 1, g_ext_fptr, params->imatrix_char);
        init_matrix(datatype, B, n, nrhs, ldb, g_ext_fptr, params->imatrix_char);

        /* Save the original matrix AB */
        copy_matrix(datatype, "full", ldab, n, AB, ldab, AB_test, ldab);
    }
    else
    {
        /* Initialize & convert random band matrix into band storage as per API need */
        rand_band_storage_matrix(datatype, n, n, kl, ku, AB, ldab);

        /* Initialize random B matrix */
        rand_matrix(datatype, B, n, nrhs, ldb);

        /* Oveflow or underflow test initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_gbtrs(datatype, n, nrhs, B, ldb, params->imatrix_char);
        }

        /* Save the original matrix AB */
        copy_matrix(datatype, "full", ldab, n, AB, ldab, AB_test, ldab);

        invoke_gbtrf(datatype, &n, &n, &kl, &ku, AB_test, &ldab, IPIV, &info);
    }
    /* Save the original matrix B */
    copy_matrix(datatype, "full", n, nrhs, B, ldb, X, ldb);

    /* call to API */
    prepare_gbtrs_run(trans, n, kl, ku, nrhs, AB_test, ldab, IPIV, X, ldb, datatype, n_repeats,
                      &time_min, &info);

    /* execution time */
    *t = time_min;

    /* performance computation */
    *perf = (2.0 * n * (ku + 2 * kl)) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        *perf *= 4.0;
    }
    /* output validation */
    if((!FLA_EXTREME_CASE_TEST) && info == 0)
    {
        create_matrix(datatype, matrix_layout, n, n, &A, n);
        reset_matrix(datatype, n, n, A, n);
        /* Get original Band matrix from AB*/
        get_band_matrix_from_band_storage(datatype, n, n, kl, ku, AB, ldab, A, n);
        /* Call validate_getrs() to validate the output*/
        validate_getrs(&trans, n, nrhs, A, n, B, ldb, X, datatype, residual, &info);
        free_matrix(A);
    }

    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(AB);
    free_matrix(AB_test);
    free_vector(IPIV);
    free_matrix(B);
    free_matrix(X);
}

void prepare_gbtrs_run(char trans, integer n_A, integer kl, integer ku, integer nrhs, void *AB,
                       integer ldab, integer *IPIV, void *B, integer ldb, integer datatype,
                       integer n_repeats, double *time_min_, integer *info)
{
    integer i;
    void *B_save;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, matrix_layout, n_A, nrhs, &B_save, ldb);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", n_A, nrhs, B, ldb, B_save, ldb);

        exe_time = fla_test_clock();

        /*  call to API */
        invoke_gbtrs(datatype, &trans, &n_A, &kl, &ku, &nrhs, AB, &ldab, IPIV, B_save, &ldb, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the final result to B matrix*/
    copy_matrix(datatype, "full", n_A, nrhs, B_save, ldb, B, ldb);
    free_matrix(B_save);
}

/*
 *  This function calls LAPACK interface of
 *  LU factorization - GBTRS
 *  */
void invoke_gbtrs(integer datatype, char *trans, integer *n, integer *kl, integer *ku,
                  integer *nrhs, void *ab, integer *ldab, integer *ipiv, void *b, integer *ldb,
                  integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }
    }
}
