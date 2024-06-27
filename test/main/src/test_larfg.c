/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_larfg_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_larfg_run(integer datatype, integer n_A, integer incx, void *x, void *tau,
                       integer n_repeats, double *time_min_);
void invoke_larfg(integer datatype, integer *n, void *x, integer *incx, integer *abs_incx,
                  void *tau);

void fla_test_larfg(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Auxilary routines";
    char *front_str = "LARFG";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, AUX, fla_test_larfg_experiment);
        tests_not_run = 0;
    }
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }

    if(argc >= 8 && argc <= 9)
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
        params->aux_paramslist[0].alpha_real = strtod(argv[4], &endptr);
        params->aux_paramslist[0].alpha_imag = strtod(argv[5], &endptr);
        params->aux_paramslist[0].incx_larfg = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->aux_paramslist[0].aux_threshold = CLI_NORM_THRESH;

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
                fla_test_larfg_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, SQUARE_INPUT, N, N, residual,
                                      params->aux_paramslist[0].aux_threshold, time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for larfg \n");
        printf(
            "./<EXE> LARFG <precisions - sdcz>  <N> <alpha_real> <alpha_imag> <incx> <repeats>\n");
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

void fla_test_larfg_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer n, incx, x_length, inc_x, info = 0;
    void *x, *x_test, *tau = NULL;
    double time_min = 1e9;
    double alpha_real, alpha_imag;

    incx = params->aux_paramslist[pci].incx_larfg;
    alpha_real = params->aux_paramslist[pci].alpha_real;
    alpha_imag = params->aux_paramslist[pci].alpha_imag;
    *residual = params->aux_paramslist[pci].aux_threshold;
    /* Determine the dimensions */
    n = p_cur;

    inc_x = fla_i_abs(&incx);
    x_length = (1 + (n - 2) * inc_x) + inc_x;
    create_vector(datatype, &x, x_length);
    create_vector(datatype, &tau, 1);

    /* Initializing input values for vector */
    init_vector(datatype, x, n, inc_x, g_ext_fptr, params->imatrix_char);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_underflow_overflow_larfg(datatype, n, 1, x, incx, params->imatrix_char);
    }

    create_vector(datatype, &x_test, x_length);
    reset_vector(datatype, x_test, x_length, 1);

    /* Assigning alpha value gto X vector in case of command line.
     * In case of config the x[0] is treated as alpha */
    if(!config_data)
        assign_value(datatype, x, alpha_real, alpha_imag);

    /* Make a copy of inputs. This is required to validate the API functionality. */
    copy_vector(datatype, n, x, inc_x, x_test, inc_x);

    /* call to API */
    prepare_larfg_run(datatype, n, incx, x_test, tau, n_repeats, &time_min);
    /* execution time */
    *t = time_min;
    if(time_min == d_zero)
    {
        time_min = 1e-9;
        *t = time_min;
    }
    /* Performance Computation */
    *perf = (double)(2.0 * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        *perf *= 4.0;
    }
    /* Output Validation */
    if((!FLA_EXTREME_CASE_TEST || FLA_OVERFLOW_UNDERFLOW_TEST))
        validate_larfg(datatype, n, incx, x_length, x, x_test, tau, residual);
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if(!check_extreme_value(datatype, n, 1, x_test, incx, params->imatrix_char))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);
    /* Free up the buffers */
    free_vector(x);
    free_vector(x_test);
    free_vector(tau);
}

void prepare_larfg_run(integer datatype, integer n_A, integer incx, void *x, void *tau,
                       integer n_repeats, double *time_min_)
{
    integer i, x_length;
    double time_min = 1e9, exe_time;
    integer inc_x = fla_i_abs(&incx);
    x_length = (1 + (n_A - 2) * inc_x) + inc_x;
    void *x_save;
    create_vector(datatype, &x_save, x_length);

    for(i = 0; i < n_repeats; ++i)
    {
        /* Make a copy of the input vector x. Same input values will be passed in each itertaion. */
        copy_vector(datatype, x_length, x, 1, x_save, 1);
        exe_time = fla_test_clock();

        /*  call  larfg API */
        invoke_larfg(datatype, &n_A, x_save, &incx, &inc_x, tau);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;

    /* Save the final result to x vector */
    copy_vector(datatype, x_length, x_save, 1, x, 1);

    /* Free up the output buffer */
    free_vector(x_save);
}

/* larfg API call interface */
void invoke_larfg(integer datatype, integer *n, void *x, integer *incx, integer *abs_incx,
                  void *tau)
{
    switch(datatype)
    {
        case FLOAT:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            float *x_ptr = x;
            fla_lapack_slarfg(n, &x_ptr[0], &x_ptr[*abs_incx], incx, tau);
            break;
        }
        case DOUBLE:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            double *x_ptr = x;
            fla_lapack_dlarfg(n, &x_ptr[0], &x_ptr[*abs_incx], incx, tau);
            break;
        }
        case COMPLEX:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            scomplex *x_ptr = x;
            fla_lapack_clarfg(n, &x_ptr[0], &x_ptr[*abs_incx], incx, tau);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            dcomplex *x_ptr = x;
            fla_lapack_zlarfg(n, &x_ptr[0], &x_ptr[*abs_incx], incx, tau);
            break;
        }
    }
}