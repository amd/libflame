/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
void fla_test_larf_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual);
void prepare_larf_run(integer datatype, char side, integer m, integer n, void *v, integer incv,
                      void *tau, void *c__, integer ldc__, void *c__out, integer ldc__out,
                      void *work, integer n_repeats, double *time_min_);
void invoke_larf(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                 void *tau, void *c__, integer *ldc, void *work);
void invoke_larfg(integer datatype, integer *n, void *x, integer *incx, integer *abs_incx,
                  void *tau);
void fla_test_larf(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Auxilary routines";
    char *front_str = "LARF";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, AUX, fla_test_larf_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }

    if(argc >= 9 && argc <= 10)
    {
        /* Test with parameters from commandline */
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->aux_paramslist[0].side = argv[3][0];
        M = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].incv = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].ldc = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_larf_experiment(params, datatype, M, N, 0, n_repeats, einfo, &perf,
                                         &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, RECT_INPUT, M, N, residual,
                                      params->aux_paramslist[0].aux_threshold, time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for larf \n");
        printf("./<EXE> LARF <precisions - sdcz>  <SIDE> <M> <N> <INCV> <LDC> <repeats>\n");
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

void fla_test_larf_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual)
{
    integer m, n;
    double time_min = 1e9;
    void *tau = NULL;
    integer v_length;
    void *work = NULL;
    void *v = NULL;
    void *v_tmp = NULL;
    void *c__ = NULL;
    void *c__out = NULL;

    char side = params->aux_paramslist[pci].side;
    integer incv = params->aux_paramslist[pci].incv;
    integer ldc = params->aux_paramslist[pci].ldc;

    m = p_cur;
    n = q_cur;

    integer incv_abs = fla_i_abs(&incv);
    integer v_num_elements;
    integer work_num_elements;

    if(side == 'L')
    {
        v_num_elements = m;
        work_num_elements = n;
    }
    else
    {
        v_num_elements = n;
        work_num_elements = m;
    }

    v_length = 1 + (v_num_elements - 1) * incv_abs;
    create_vector(datatype, &work, work_num_elements);

    create_vector(datatype, &v_tmp, v_length);
    create_vector(datatype, &tau, 1);

    rand_vector(datatype, v_num_elements, v_tmp, incv_abs, d_zero, d_zero, 'R');

    /* Input generation (v_tmp and tau) for larf from larfg
       Increment of v_tmp for larfg must be positive. Hence calling larfg with incv_abs
       Increment of v for larf could be positive or negative. Hence copying
       from v_tmp using incv(which could be positive or negative)
    */
    invoke_larfg(datatype, &v_num_elements, v_tmp, &incv_abs, &incv_abs, tau);
    assign_value(datatype, v_tmp, 1, 0);
    create_vector(datatype, &v, v_length);
    copy_vector(datatype, v_num_elements, v_tmp, incv_abs, v, incv);

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &c__, ldc);
    rand_matrix(datatype, c__, m, n, ldc);

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &c__out, ldc);

    /* call to API */
    prepare_larf_run(datatype, side, m, n, v, incv, tau, c__, ldc, c__out, ldc, work, n_repeats,
                     &time_min);
    /* execution time */
    *t = time_min;
    if(time_min == d_zero)
    {
        time_min = 1e-9;
        *t = time_min;
    }
    /* Performance Computation */
    *perf = (double)(2.0 * m * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        *perf *= 4.0;
    }
    /* Output Validation */
    validate_larf(datatype, side, m, n, v, incv, c__, ldc, c__out, ldc, tau, residual);

    /* Free up the buffers */
    free_matrix(c__);
    free_matrix(c__out);
    free_vector(v);
    free_vector(v_tmp);
    free_vector(work);
    free_vector(tau);
}

void prepare_larf_run(integer datatype, char side, integer m, integer n, void *v, integer incv,
                      void *tau, void *c__, integer ldc__, void *c__out, integer ldc__out,
                      void *work, integer n_repeats, double *time_min_)
{
    integer i;
    double time_min = 1e9, exe_time;

    for(i = 0; i < n_repeats; ++i)
    {
        copy_matrix(datatype, "full", m, n, c__, ldc__, c__out, ldc__out);

        exe_time = fla_test_clock();

        /*  call  larf API */
        invoke_larf(datatype, &side, &m, &n, v, &incv, tau, c__out, &ldc__out, work);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
}

/* larf API call interface */
void invoke_larf(integer datatype, char *side, integer *m, integer *n, void *v, integer *incv,
                 void *tau, void *c__, integer *ldc__, void *work)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dlarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_clarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlarf(side, m, n, v, incv, tau, c__, ldc__, work);
            break;
        }
    }
}
