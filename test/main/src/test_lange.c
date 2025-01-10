/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

/* Local prototypes */
integer i_abs(integer *x);
void fla_test_lange_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_lange_run(integer datatype, char norm_type, void *A, integer m, integer n, integer lda,
                       void *result, integer n_repeats, double *time_min_);
void invoke_lange(integer datatype, char *norm_type, integer *m, integer *n, void *A, integer *lda,
                  void *work, void *result);

void fla_test_lange(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Auxilary routines";
    char front_str[8] = "LANGE  \0";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    integer invalid_normtype = 0;

    if(argc == 1)
    {
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        for(integer i = 0; i < params->aux_paramslist[0].num_ranges; i++)
        {

            invalid_normtype = invalid_normtype
                               || fla_validate_lange_norm_types(
                                   params->aux_paramslist[i].norm_types_str,
                                   params->aux_paramslist[i].norm_types_str, MAX_NUM_NORMTYPES);
        }
        if(!invalid_normtype)
        {
            fla_test_op_driver(front_str, RECT_INPUT, params, AUX, fla_test_lange_experiment);
            tests_not_run = 0;
        }
    }
    if(argc == 9)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[8]);
    }

    if(argc >= 8 && argc <= 9)
    {
        /* Test with parameters from commandline */
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        invalid_normtype = fla_validate_lange_norm_types(
            argv[3], params->aux_paramslist[0].norm_types_str, MAX_NUM_NORMTYPES);

        M = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

        // set the threshold
        params->aux_paramslist[0].aux_threshold = CLI_NORM_THRESH;
        // set front string

        if(n_repeats > 0 && !invalid_normtype)
        {
            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalid dataype */
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
                fla_test_lange_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for lange \n");
        printf("./<EXE> lange <precisions - sdcz> <Norm Types - M1IF> <M> <N> <lda> <repeats> "
               "[file]\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if(invalid_normtype)
    {
        printf("\nInvalid norm types specified, choose valid norm types from 'M1IF' ( Multiple "
               "Norm types can be clubbed together)\n\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
    return;
}

void fla_test_lange_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    void *A, *scal;
    void *result;
    integer m = p_cur;
    integer n = q_cur;
    integer lda = params->aux_paramslist[pci].lda;
    integer i;
    double residual, err_thresh;

    err_thresh = params->aux_paramslist[pci].aux_threshold;
    if(lda == -1)
    {
        lda = fla_max(1, m);
    }

    /* If lda is less than m, then result with invalid param */
    if(lda < m)
    {
        time_min = perf = 0.;
        FLA_PRINT_TEST_STATUS(m, n, DBL_MIN, err_thresh);
        return;
    }

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_vector(get_realtype(datatype), &result, 1);
    create_vector(get_realtype(datatype), &scal, 1);

    for(i = 0; i < MAX_NUM_NORMTYPES; i++)
    {

        char test_norm_type = params->aux_paramslist[pci].norm_types_str[i];
        if(test_norm_type == '\0')
        {
            break;
        }

        residual = err_thresh;
        time_min = 1e9;

        if(g_ext_fptr != NULL)
        {
            /* Initialize input vectors with custom data */
            init_matrix_from_file(datatype, A, m, n, lda, g_ext_fptr);
        }
        else
        {
            /* Initialize input matrix with random numbers */
            rand_matrix(datatype, A, m, n, lda);
        }

        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_lange(datatype, m, n, A, lda, test_norm_type,
                                                  params->imatrix_char, scal);
        }

        prepare_lange_run(datatype, test_norm_type, A, m, n, lda, result, n_repeats, &time_min);

        /* execution time */
        if(time_min == d_zero)
        {
            time_min = 1e-9;
        }
        /* Compute the performance of the best experiment repeat */
        /* 4*n */
        perf = (double)(4.0 * m * n) / time_min / FLOPS_PER_UNIT_PERF;
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        {
            perf *= 2;
        }

        /* output validation */
        validate_lange(tst_api, datatype, test_norm_type, m, n, lda, A, result, residual);
    }

    /* Free up the buffers */
    free_matrix(A);
    free_vector(result);
    free_vector(scal);
}

void prepare_lange_run(integer datatype, char norm_type, void *A, integer m, integer n, integer lda,
                       void *result, integer n_repeats, double *time_min_)
{
    integer i;
    void *work = NULL;
    double exe_time, time_min = 1e9;

    if(norm_type == 'I')
    {
        create_vector(get_realtype(datatype), &work, m);
    }

    for(i = 0; i < n_repeats; ++i)
    {
        exe_time = fla_test_clock();
        /*  call lange API */
        invoke_lange(datatype, &norm_type, &m, &n, A, &lda, work, result);
        exe_time = fla_test_clock() - exe_time;
        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;

    if(norm_type == 'I')
    {
        free_vector(work);
    }
}

/*
 *  lange calls LAPACK interface
 *  */
void invoke_lange(integer datatype, char *norm_type, integer *m, integer *n, void *A, integer *lda,
                  void *work, void *result)
{
    switch(datatype)
    {
        case FLOAT:
        {
            *(float *)result = fla_lapack_slange(norm_type, m, n, A, lda, work);
            break;
        }
        case DOUBLE:
        {
            *(double *)result = fla_lapack_dlange(norm_type, m, n, A, lda, work);
            break;
        }
        case COMPLEX:
        {
            *(float *)result = fla_lapack_clange(norm_type, m, n, A, lda, work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            *(double *)result = fla_lapack_zlange(norm_type, m, n, A, lda, work);
            break;
        }
    }
}
