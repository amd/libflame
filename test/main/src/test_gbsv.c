/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"


extern double perf;
extern double time_min;

/* Local prototypes */
void fla_test_gbsv_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo);
void prepare_gbsv_run(integer n_A, integer kl, integer ku, integer nrhs, void *AB, integer ldab,
                      void *B, integer ldb, integer *ipiv, integer datatype, integer *info,
                      test_params_t *params);
void invoke_gbsv(integer datatype, integer *n, integer *kl, integer *ku, integer *nrhs, void *ab,
                 integer *ldab, integer *ipiv, void *b, integer *ldb, integer *info);

void fla_test_gbsv(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Linear Solve using LU for Band Matrix";
    char *front_str = "GBSV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gbsv_experiment);
        tests_not_run = 0;
    }
    if(argc == 11)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[10]);
    }
    if(argc >= 10 && argc <= 11)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N, KL, KU;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        KL = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        KU = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldab = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        
        /* Store KL and KU in the structure */
        params->lin_solver_paramslist[0].kl = KL;
        params->lin_solver_paramslist[0].ku = KU;
        n_repeats = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_gbsv_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gbsv\n");
        printf("./<EXE> gbsv <precisions - sdcz> <N> <KL> <KU> <NRHS> <LDAB> <LDB> <repeats>\n");
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

void fla_test_gbsv_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo)
{
    integer n, kl, ku, ldab, ldb, NRHS, info = 0;
    void *IPIV = NULL, *AB = NULL, *AB_save = NULL, *B = NULL, *B_save = NULL;
    double residual, err_thresh;

    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    NRHS = params->lin_solver_paramslist[pci].nrhs;
    kl = params->lin_solver_paramslist[pci].kl;
    ku = params->lin_solver_paramslist[pci].ku;
    /* Determine the dimensions */
    n = p_cur;
    ldab = params->lin_solver_paramslist[pci].ldab;
    ldb = params->lin_solver_paramslist[pci].ldb;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        if(ldab == -1)
        {
            ldab = fla_max(1, 2 * kl + ku + 1);
        }
        if(ldb == -1)
        {
            ldb = fla_max(1, n);
        }
    }

    /* Create the matrices for the current operation */
    create_matrix(datatype, LAPACK_COL_MAJOR, ldab, n, &AB, ldab);
    create_matrix(datatype, LAPACK_COL_MAJOR, ldab, n, &AB_save, ldab);
    create_vector(INTEGER, &IPIV, n);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &B, ldb);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &B_save, ldb);

    /* Initialize the test matrices */
    if(g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data from file */
        init_matrix(datatype, AB, ldab, n, ldab, g_ext_fptr, params->imatrix_char);
        init_matrix(datatype, B, n, NRHS, ldb, g_ext_fptr, params->imatrix_char);

        /* Save the original matrix AB */
        copy_matrix(datatype, "full", ldab, n, AB, ldab, AB_save, ldab);
    }
    else
    {
        /* Initialize & convert random band matrix into band storage as per API need */
        rand_band_storage_matrix(datatype, n, n, kl, ku, AB, ldab);
        /* Initialize random B matrix */
        rand_matrix(datatype, B, n, NRHS, ldb);

        /* Save the original matrix AB */
        copy_matrix(datatype, "full", ldab, n, AB, ldab, AB_save, ldab);
    }
    
    /* Save the original matrix B */
    copy_matrix(datatype, "full", n, NRHS, B, ldb, B_save, ldb);
    
    /* call to API */
    prepare_gbsv_run(n, kl, ku, NRHS, AB_save, ldab, B_save, ldb, IPIV, datatype, &info, params);

    /* performance computation */
    /* For band matrix: O(nkl(kl + ku)) for factorization + O(n(2kl + ku)r) for solution */
    perf = (double)(n * kl * (kl + ku) + n * (2 * kl + ku) * NRHS) / time_min / FLOPS_PER_UNIT_PERF;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);

    /* Free up the buffers */
    free_matrix(AB);
    free_matrix(AB_save);
    free_vector(IPIV);
    free_matrix(B);
    free_matrix(B_save);
}

void prepare_gbsv_run(integer n_A, integer kl, integer ku, integer nrhs, void *AB, integer ldab,
                      void *B, integer ldb, integer *IPIV, integer datatype, integer *info,
                      test_params_t *params)
{
    void *AB_test, *B_test;
    double exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, ldab, n_A, &AB_test, ldab);
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, nrhs, &B_test, ldb);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", ldab, n_A, AB, ldab, AB_test, ldab);
        copy_matrix(datatype, "full", n_A, nrhs, B, ldb, B_test, ldb);
        
        exe_time = fla_test_clock();
        /* call LAPACK gbsv API */
        invoke_gbsv(datatype, &n_A, &kl, &ku, &nrhs, AB_test, &ldab, IPIV, B_test, &ldb, info);
        exe_time = fla_test_clock() - exe_time;
        
        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_WITH_INFO
    }

    /* Save the final result to B matrix */
    copy_matrix(datatype, "full", n_A, nrhs, B_test, ldb, B, ldb);

    free_matrix(AB_test);
    free_matrix(B_test);
}

/*
 * gbsv_API calls LAPACK interface of
 * Band matrix linear system solver - gbsv
 */
void invoke_gbsv(integer datatype, integer *n, integer *kl, integer *ku, integer *nrhs, void *ab,
                 integer *ldab, integer *ipiv, void *b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
            break;
        }
    }
}
