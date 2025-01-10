/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;
/* Local prototypes */
void fla_test_lartg_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_lartg_run(integer datatype, void *f, void *g, void *r, void *c, void *s,
                       integer n_repeats, double *time_min_, integer interfacetype);
void invoke_lartg(integer datatype, void *f, void *g, void *c, void *s, void *r);

void fla_test_lartg(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Auxilary routines";
    char *front_str = "LARTG";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    integer i, num_types;
    integer datatype, n_repeats;
    char stype, type_flag[4] = {0};
    char *endptr;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        num_types = params->aux_paramslist[0].num_data_types;
        n_repeats = params->aux_paramslist[0].num_repeats;

        if(n_repeats > 0)
        {
            /* Loop over the requested datatypes. */
            for(i = 0; i < num_types; ++i)
            {
                datatype = params->aux_paramslist[0].data_types[i];
                stype = params->aux_paramslist[0].data_types_char[i];

                /* Call the test code */
                fla_test_lartg_experiment(front_str, params, datatype, 2, i_one, 0, n_repeats,
                                          einfo);
                tests_not_run = 0;
            }
        }
    }
    if(argc == 5)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[4]);
    }

    if(argc >= 4 && argc <= 5)
    {
        /* Test with parameters from commandline */
        /* Parse the arguments */
        num_types = strlen(argv[2]);

        n_repeats = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_lartg_experiment(front_str, params, datatype, 2, i_one, 0, n_repeats,
                                          einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\n Illegal arguments for lartg \n");
        printf("./<EXE> lartg <precisions - sdcz> <repeats> [file] \n");
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

void fla_test_lartg_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    void *s = NULL, *c = NULL;
    void *f = NULL, *g = NULL, *r = NULL;
    double err_thresh;
    integer interfacetype = params->interfacetype;

    integer realtype;
    realtype = get_realtype(datatype);

    err_thresh = params->aux_paramslist[pci].aux_threshold;

    create_vector(realtype, &c, 1);
    create_vector(datatype, &s, 1);

    create_vector(datatype, &f, 1);
    create_vector(datatype, &g, 1);
    create_vector(datatype, &r, 1);

    if(g_ext_fptr != NULL)
    {
        init_vector_from_file(datatype, f, 1, 1, g_ext_fptr);
        init_vector_from_file(datatype, g, 1, 1, g_ext_fptr);
    }
    else
    {
        rand_vector(datatype, 1, f, 1, d_zero, d_zero, 'R');
        rand_vector(datatype, 1, g, 1, d_zero, d_zero, 'R');
    }
    /* call to API */
    prepare_lartg_run(datatype, f, g, r, c, s, n_repeats, &time_min, interfacetype);

    /* execution time */
    if(time_min == d_zero)
    {
        time_min = 1e-9;
    }
    /* Compute the performance of the best experiment repeat */
    perf = (double)(6.0) / time_min / FLOPS_PER_UNIT_PERF;

    /* output validation */
    validate_lartg(tst_api, datatype, f, g, r, c, s, err_thresh);

    /* Free up the buffers */
    free_vector(c);
    free_vector(s);
    free_vector(f);
    free_vector(g);
    free_vector(r);
}

void prepare_lartg_run(integer datatype, void *f, void *g, void *r, void *c, void *s,
                       integer n_repeats, double *time_min_, integer interfacetype)
{
    integer i;
    double t_min = 1e9, exe_time;

    for(i = 0; i < n_repeats; ++i)
    {
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call lartg CPP API */
            invoke_cpp_lartg(datatype, f, g, c, s, r);
            exe_time = fla_test_clock() - exe_time;
        }
        else
#endif
        {
            exe_time = fla_test_clock();
            /* Call lartg API */
            invoke_lartg(datatype, f, g, c, s, r);
            exe_time = fla_test_clock() - exe_time;
        }
        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);
    }

    *time_min_ = t_min;
}

void invoke_lartg(integer datatype, void *f, void *g, void *c, void *s, void *r)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slartg(f, g, (float *)c, (float *)s, r);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dlartg(f, g, (double *)c, (double *)s, r);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_clartg(f, g, (float *)c, ((scomplex *)s), r);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlartg(f, g, (double *)c, ((dcomplex *)s), r);
            break;
        }
    }
}
