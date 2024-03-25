/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_gtsv_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual);
void prepare_gtsv_run(integer n_A, integer nrhs, void *dl, void *d, void *du, void *B, integer ldb,
                      integer datatype, integer n_repeats, double *time_min_, integer *info);
void invoke_gtsv(integer datatype, integer *nrhs, integer *n, void *dl, void *d, void *du, void *b,
                 integer *ldb, integer *info);

void fla_test_gtsv(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Linear Solve for General Tridiagonal matrix";
    char *front_str = "GTSV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gtsv_experiment);
        tests_not_run = 0;
    }
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if(argc >= 7 && argc <= 8)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_gtsv_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
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
        printf("\nIllegal arguments for gtsv\n");
        printf("./<EXE> gtsv <precisions - sdcz>  <N> <NRHS> <LDB> <repeats>\n");
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

void fla_test_gtsv_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual)
{
    integer n, ldb, NRHS;
    integer info = 0;
    void *dl, *d, *du, *B, *xact = NULL;
    void *dl_save, *d_save, *du_save, *B_save;
    double time_min = 1e9;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    /* Determine the dimensions*/
    n = p_cur;
    NRHS = params->lin_solver_paramslist[pci].nrhs;
    ldb = params->lin_solver_paramslist[pci].ldb;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(ldb == -1)
        {
            ldb = fla_max(1, n);
        }
    }

    /* Create vector and matrix for the current operation */
    create_vector(datatype, &dl, n - 1);
    create_vector(datatype, &d, n);
    create_vector(datatype, &du, n - 1);
    create_matrix(datatype, &B, ldb, NRHS);
    create_matrix(datatype, &xact, ldb, NRHS);

    /* Initializing expected solution xact */
    init_matrix(datatype, xact, n, NRHS, ldb, g_ext_fptr, params->imatrix_char);

    /* Initializing tridiagonal vectors */
    init_vector(datatype, dl, n - 1, i_one, g_ext_fptr, params->imatrix_char);
    init_vector(datatype, d, n, i_one, g_ext_fptr, params->imatrix_char);
    init_vector(datatype, du, n - 1, i_one, g_ext_fptr, params->imatrix_char);

    /* Generating B matrix from xact & tridiagonal vectors */
    tridiag_matrix_multiply(datatype, n, NRHS, dl, d, du, xact, ldb, B, ldb);

    /* Make a copy of input buffers. This is required to validate the API functionality */
    create_vector(datatype, &dl_save, n - 1);
    copy_vector(datatype, n - 1, dl, i_one, dl_save, i_one);

    create_vector(datatype, &d_save, n);
    copy_vector(datatype, n, d, i_one, d_save, i_one);

    create_vector(datatype, &du_save, n - 1);
    copy_vector(datatype, n - 1, du, i_one, du_save, i_one);

    create_matrix(datatype, &B_save, ldb, NRHS);
    copy_matrix(datatype, "full", n, NRHS, B, ldb, B_save, ldb);

    prepare_gtsv_run(n, NRHS, dl_save, d_save, du_save, B_save, ldb, datatype, n_repeats, &time_min,
                     &info);

    /* Execution time */
    *t = time_min;

    /* Performance computation */
    *perf = (double)(4.0 * (n - 1) * NRHS) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Output Validation */
    if(info == 0)
        validate_gtsv(datatype, n, NRHS, B, ldb, B_save, xact, ldb, dl, d, du, dl_save, d_save,
                      du_save, info, residual);

    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_vector(dl);
    free_vector(d);
    free_vector(du);
    free_matrix(B);
    free_vector(dl_save);
    free_vector(d_save);
    free_vector(du_save);
    free_matrix(B_save);
    free_matrix(xact);
}

void prepare_gtsv_run(integer n_A, integer nrhs, void *dl, void *d, void *du, void *B, integer ldb,
                      integer datatype, integer n_repeats, double *time_min_, integer *info)
{
    integer i;
    void *dl_test, *d_test, *du_test, *B_test;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input tridiagonal vectors and matrix. Same input values will be passed in
     * each itertaion.*/
    create_vector(datatype, &dl_test, n_A - 1);
    create_vector(datatype, &d_test, n_A);
    create_vector(datatype, &du_test, n_A - 1);
    create_matrix(datatype, &B_test, ldb, nrhs);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {

        /* Restore input matrix and vector */
        copy_vector(datatype, n_A - 1, dl, i_one, dl_test, i_one);
        copy_vector(datatype, n_A, d, i_one, d_test, i_one);
        copy_vector(datatype, n_A - 1, du, i_one, du_test, i_one);
        copy_matrix(datatype, "full", n_A, nrhs, B, ldb, B_test, ldb);

        exe_time = fla_test_clock();

        /* call  to gtsv API */
        invoke_gtsv(datatype, &n_A, &nrhs, dl_test, d_test, du_test, B_test, &ldb, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /* Make a copy of output buffers. This is required to validate the API functionality */
    copy_vector(datatype, n_A - 1, dl_test, i_one, dl, i_one);
    copy_vector(datatype, n_A, d_test, i_one, d, i_one);
    copy_vector(datatype, n_A - 1, du_test, i_one, du, i_one);
    copy_matrix(datatype, "full", n_A, nrhs, B_test, ldb, B, ldb);

    /* Free up the output buffers */
    free_vector(dl_test);
    free_vector(d_test);
    free_vector(du_test);
    free_matrix(B_test);
}

/* LARFG API call interface - Linear Solve for General Tridiagonal Matrix */
void invoke_gtsv(integer datatype, integer *n, integer *nrhs, void *dl, void *d, void *du, void *b,
                 integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgtsv(n, nrhs, dl, d, du, b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgtsv(n, nrhs, dl, d, du, b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgtsv(n, nrhs, dl, d, du, b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgtsv(n, nrhs, dl, d, du, b, ldb, info);
            break;
        }
    }
}
