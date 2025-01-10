/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;
integer row_major_gtsv_ldb;

/* Local prototypes */
void fla_test_gtsv_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo);
void prepare_gtsv_run(integer n_A, integer nrhs, void *dl, void *d, void *du, void *B, integer ldb,
                      integer datatype, integer n_repeats, double *time_min_, integer *info,
                      integer interfacetype, integer layout);
void invoke_gtsv(integer datatype, integer *nrhs, integer *n, void *dl, void *d, void *du, void *b,
                 integer *ldb, integer *info);
integer invoke_lapacke_gtsv(integer datatype, integer layout, integer n, integer nrhs, void *dl,
                            void *d, void *du, void *B, integer ldb);
double prepare_lapacke_gtsv_run(integer datatype, integer layout, integer n, integer nrhs, void *dl,
                                void *d, void *du, void *B, integer ldb, integer *info);

void fla_test_gtsv(integer argc, char **argv, test_params_t *params)
{
    srand(13); /* Setting the seed for random input genetation values */
    char *op_str = "Linear Solve for General Tridiagonal matrix";
    char *front_str = "GTSV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

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
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gtsv_ldb = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldb = N;
        }
        else
        {
            params->lin_solver_paramslist[0].ldb = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
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
                fla_test_gtsv_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
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

void fla_test_gtsv_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo)
{
    integer n, ldb, ldx, NRHS;
    integer info = 0;
    void *dl, *d, *du, *B, *xact = NULL;
    void *dl_save, *d_save, *du_save, *B_save;
    void *A = NULL, *scal = NULL;
    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;
    double residual, err_thresh;

    /* Determine the dimensions*/
    n = p_cur;
    NRHS = params->lin_solver_paramslist[pci].nrhs;
    ldb = params->lin_solver_paramslist[pci].ldb;

    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(ldb == -1)
        {
            ldb = fla_max(1, n);
        }
    }
    ldx = ldb;

    /* Create vector and matrix for the current operation */
    create_vector(datatype, &dl, n - 1);
    create_vector(datatype, &d, n);
    create_vector(datatype, &du, n - 1);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &B, ldb);

    if(g_ext_fptr != NULL || FLA_EXTREME_CASE_TEST)
    {
        /* Initializing tridiagonal vectors */
        init_vector(datatype, dl, n - 1, i_one, g_ext_fptr, params->imatrix_char);
        init_vector(datatype, d, n, i_one, g_ext_fptr, params->imatrix_char);
        init_vector(datatype, du, n - 1, i_one, g_ext_fptr, params->imatrix_char);

        if(g_ext_fptr != NULL)
        {
            /* Read B from file */
            init_matrix(datatype, B, n, NRHS, ldb, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &xact, ldb);
            /* Initializing expected solution xact without extreme values
             * Extreme values will be introduced in the next step */
            init_matrix(datatype, xact, n, NRHS, ldx, g_ext_fptr, '\0');

            /* Generating B matrix from xact & tridiagonal vectors */
            tridiag_matrix_multiply(datatype, n, NRHS, dl, d, du, xact, ldx, B, ldb);
        }
    }
    else
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &xact, ldb);

        /* Initializing tridiagonal vectors */
        init_vector(datatype, dl, n - 1, i_one, g_ext_fptr, params->imatrix_char);
        init_vector(datatype, d, n, i_one, g_ext_fptr, params->imatrix_char);
        init_vector(datatype, du, n - 1, i_one, g_ext_fptr, params->imatrix_char);

        /* Initializing expected solution xact */
        init_matrix(datatype, xact, n, NRHS, ldx, g_ext_fptr, params->imatrix_char);

        /* Generating B matrix from xact & tridiagonal vectors */
        tridiag_matrix_multiply(datatype, n, NRHS, dl, d, du, xact, ldx, B, ldb);
    }

    /* Scaling the input matrices for overflow or underflow tests.
    Where A is a tridiagonal matrix generated by dl, d, du and again copying back the
    values into dl, d, du from A after the scaling.*/
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, ldb);
        create_vector(get_realtype(datatype), &scal, 1);

        copy_tridiag_matrix(datatype, dl, d, du, n, n, A, ldb);
        init_matrix_overflow_underflow_gtsv(datatype, n, n, A, ldb, params->imatrix_char, scal);
        copy_tridiag_vector(datatype, dl, d, du, n, n, A, ldb);

        /* Generating B matrix from xact & tridiagonal vectors */
        tridiag_matrix_multiply(datatype, n, NRHS, dl, d, du, xact, ldx, B, ldb);
    }

    /* Make a copy of input buffers. This is required to validate the API functionality */
    create_vector(datatype, &dl_save, n - 1);
    copy_vector(datatype, n - 1, dl, i_one, dl_save, i_one);

    create_vector(datatype, &d_save, n);
    copy_vector(datatype, n, d, i_one, d_save, i_one);

    create_vector(datatype, &du_save, n - 1);
    copy_vector(datatype, n - 1, du, i_one, du_save, i_one);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, NRHS, &B_save, ldb);
    copy_matrix(datatype, "full", n, NRHS, B, ldb, B_save, ldb);

    prepare_gtsv_run(n, NRHS, dl_save, d_save, du_save, B_save, ldb, datatype, n_repeats, &time_min,
                     &info, interfacetype, layout);

    /* Performance computation */
    perf = (double)(4.0 * (n - 1) * NRHS) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* Output Validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_gtsv(tst_api, datatype, n, NRHS, B, ldb, B_save, xact, ldb, dl, d, du, dl_save,
                      d_save, du_save, scal, params->imatrix_char, residual);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if((!check_extreme_value(datatype, n, NRHS, dl_save, ldb, params->imatrix_char))
           && (!check_extreme_value(datatype, n, NRHS, d_save, ldb, params->imatrix_char))
           && (!check_extreme_value(datatype, n, NRHS, du_save, ldb, params->imatrix_char))
           && (!check_extreme_value(datatype, n, NRHS, B_save, ldb, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

    /* Free up the buffers */
    free_vector(dl);
    free_vector(d);
    free_vector(du);
    free_matrix(B);
    free_vector(dl_save);
    free_vector(d_save);
    free_vector(du_save);
    free_matrix(B_save);
    if(g_ext_fptr == NULL || FLA_EXTREME_CASE_TEST)
        free_matrix(xact);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_matrix(A);
        free_vector(scal);
    }
}

void prepare_gtsv_run(integer n_A, integer nrhs, void *dl, void *d, void *du, void *B, integer ldb,
                      integer datatype, integer n_repeats, double *time_min_, integer *info,
                      integer interfacetype, integer layout)
{
    integer i;
    void *dl_test, *d_test, *du_test, *B_test;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input tridiagonal vectors and matrix. Same input values will be passed in
     * each itertaion.*/
    create_vector(datatype, &dl_test, n_A - 1);
    create_vector(datatype, &d_test, n_A);
    create_vector(datatype, &du_test, n_A - 1);
    create_matrix(datatype, LAPACK_COL_MAJOR, n_A, nrhs, &B_test, ldb);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {

        /* Restore input matrix and vector */
        copy_vector(datatype, n_A - 1, dl, i_one, dl_test, i_one);
        copy_vector(datatype, n_A, d, i_one, d_test, i_one);
        copy_vector(datatype, n_A - 1, du, i_one, du_test, i_one);
        copy_matrix(datatype, "full", n_A, nrhs, B, ldb, B_test, ldb);

        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gtsv_run(datatype, layout, n_A, nrhs, dl_test, d_test,
                                                du_test, B_test, ldb, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST) /* Call CPP gtsv API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_gtsv(datatype, &n_A, &nrhs, dl_test, d_test, du_test, B_test, &ldb, info);
            ;
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();

            /* call  to gtsv API */
            invoke_gtsv(datatype, &n_A, &nrhs, dl_test, d_test, du_test, B_test, &ldb, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);
    }

    *time_min_ = t_min;
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

double prepare_lapacke_gtsv_run(integer datatype, integer layout, integer n, integer nrhs, void *dl,
                                void *d, void *du, void *B, integer ldb, integer *info)
{
    double exe_time = 0;
    void *dl_t = NULL, *d_t = NULL, *du_t = NULL, *B_t = NULL;
    integer ldb_t = ldb;
    dl_t = dl;
    d_t = d;
    du_t = du;
    B_t = B;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_gtsv_ldb, ldb_t);

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n - 1, i_one, &dl_t, i_one);
        create_matrix(datatype, layout, n, i_one, &d_t, i_one);
        create_matrix(datatype, layout, n - 1, i_one, &du_t, i_one);
        create_matrix(datatype, layout, n, nrhs, &B_t, fla_max(nrhs, ldb_t));

        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n - 1, i_one, dl, n - 1, dl_t, i_one);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, i_one, d, n, d_t, i_one);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n - 1, i_one, du, n - 1, du_t, i_one);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, nrhs, B, ldb, B_t,
                              fla_max(nrhs, ldb_t));
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE gels API */

    *info = invoke_lapacke_gtsv(datatype, layout, n, nrhs, dl_t, d_t, du_t, B_t, ldb_t);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */

    if((layout == LAPACK_ROW_MAJOR))
    {
        convert_matrix_layout(layout, datatype, n - 1, i_one, dl_t, i_one, dl, i_one);
        convert_matrix_layout(layout, datatype, n, i_one, d_t, i_one, d, i_one);
        convert_matrix_layout(layout, datatype, n - 1, i_one, du_t, i_one, du, i_one);
        convert_matrix_layout(layout, datatype, n, nrhs, B_t, ldb_t, B, ldb);

        free_matrix(dl_t);
        free_matrix(d_t);
        free_matrix(du_t);
        free_matrix(B_t);
    }

    return exe_time;
}

/* GTSV API call interface - Linear Solve for General Tridiagonal Matrix */
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

/*
LAPACKE GTSV API invoke function
*/
integer invoke_lapacke_gtsv(integer datatype, integer layout, integer n, integer nrhs, void *dl,
                            void *d, void *du, void *B, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }
    }
    return info;
}
