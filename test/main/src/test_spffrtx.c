/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#if ENABLE_AOCL_EXTENSION_APIS

#include "test_lapack.h"

extern double perf;
extern double time_min;
/* Local prototypes */
void fla_test_spffrtx_experiment(char *tst_api, test_params_t *params, integer datatype,
                                 integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                 integer einfo);
void prepare_spffrtx_run(integer n_A, integer ncolm, integer pn, void *A, integer datatype,
                         test_params_t *params);
void invoke_spffrtx(integer datatype, void *a, integer *n, integer *ncolm, void *work, void *work2);

void fla_test_spffrtx(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Computes LDLT partial factorization";
    char *front_str = "SPFFRTX";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_spffrtx_experiment);
        tests_not_run = 0;
    }
    if(argc == 7)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[6]);
    }
    if(argc >= 6 && argc <= 7)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ncolm = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

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
                fla_test_spffrtx_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for spffrtx \n");
        printf("./<EXE> spffrtx <precisions - sdcz>  <N> <ncolm>  <repeats> [file]\n");
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

void fla_test_spffrtx_experiment(char *tst_api, test_params_t *params, integer datatype,
                                 integer p_cur, integer q_cur, integer pci, integer n_repeats,
                                 integer einfo)
{
    integer n, ncolm, pn;
    void *A, *AP, *AP_save;
    double err_thresh;
    void *filename = NULL;

    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    ncolm = params->lin_solver_paramslist[pci].ncolm;

    /* Determine the dimensions*/
    n = p_cur;
    pn = n * (n + 1) / 2;
    /* Create the matrices for the current operation*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, n);
    create_vector(datatype, &AP, pn);

    /* Initialize input matrix only when not in BRT verification run mode */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST) || (FLA_RANDOM_INIT_MODE))
        {
            /* Initialize input matrix with custom data */
            init_matrix(datatype, A, n, n, n, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            /* Initialize input matrix with random numbers */
            rand_sym_matrix(datatype, A, n, n, n);
        }
        /* Pack a symmetric matrix in column first order */
        pack_matrix_lt(datatype, A, AP, n, n);
    }

    /* BRT macro for processing single packed matrix AP */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, pn, 1, AP, 1, "dd", n, ncolm)
    /* Create a copy of the input */
    create_vector(datatype, &AP_save, pn);
    copy_vector(datatype, pn, AP, i_one, AP_save, i_one);

    /* call to API */
    prepare_spffrtx_run(n, ncolm, pn, AP, datatype, params);

    /* performance computation */
    if(datatype == FLOAT || datatype == DOUBLE)
        perf = (double)ncolm / 6.0f
               * (2.0 * ncolm * ncolm - 6.0 * ncolm * n + 3.0 * ncolm + 6.0 * n * n - 6.0 * n + 7)
               / time_min / FLOPS_PER_UNIT_PERF;
    else if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf
            = (double)ncolm / 3.0f
              * (4.0 * ncolm * ncolm - 12.0 * ncolm * n + 9.0 * ncolm + 12.0 * n * n - 18.0 * n + 8)
              / time_min / FLOPS_PER_UNIT_PERF;

    /* output validation */
    double residual = 0;
    IF_FLA_BRT_VALIDATION(
        pn, 1, store_outputs_base(filename, params, 1, 0, datatype, pn, 1, AP, 1),
        validate_spffrtx(tst_api, n, ncolm, A, AP, datatype, err_thresh, params),
        check_reproducibility_base(filename, params, 1, 0, datatype, pn, 1, AP, 1))
    else if(FLA_RANDOM_INIT_MODE)
    {
        FLA_PRINT_TEST_STATUS(n, n, err_thresh, err_thresh);
    }
    else if(ncolm <= n && n > 0 && ncolm > 0)
    {
        validate_spffrtx(tst_api, n, ncolm, A, AP, datatype, err_thresh, params);
    }

    /* Free up the buffers */
free_buffers:
    FLA_FREE_FILENAME(filename);
    free_vector(AP_save);
    free_matrix(A);
    free_vector(AP);
}

void prepare_spffrtx_run(integer n_A, integer ncolm, integer pn, void *AP, integer datatype,
                         test_params_t *params)
{
    void *AP_save, *work = NULL, *work2 = NULL;
    double exe_time;

    create_vector(datatype, &AP_save, pn);
    create_vector(datatype, &work, 2 * n_A);
    create_vector(datatype, &work2, 2 * n_A);

    FLA_EXEC_LOOP_BEGIN
    {
        /* Copy original input data */
        copy_vector(datatype, pn, AP, i_one, AP_save, i_one);

        exe_time = fla_test_clock();

        /*  call  spffrtx API with AFACT to get A INV */
        invoke_spffrtx(datatype, AP_save, &n_A, &ncolm, work, work2);

        exe_time = fla_test_clock() - exe_time;
        FLA_EXEC_LOOP_UPDATE_NO_INFO
    }

    /*  Save the final result to A matrix*/
    copy_vector(datatype, pn, AP_save, i_one, AP, i_one);
    free_vector(AP_save);
    free_vector(work);
    free_vector(work2);
}

/*
 *  spffrtx_API calls LAPACK interface
 *  */
void invoke_spffrtx(integer datatype, void *ap, integer *n, integer *ncolm, void *work, void *work2)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sspffrtx(ap, n, ncolm, work, work2);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dspffrtx(ap, n, ncolm, work, work2);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cspffrtx(ap, n, ncolm, work, work2);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zspffrtx(ap, n, ncolm, work, work2);
            break;
        }
    }
}

#endif /* ENABLE_AOCL_EXTENSION_APIS */