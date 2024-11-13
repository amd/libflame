/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

/* Local prototypes */
void fla_test_gbtrf_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_gbtrf_run(integer m_A, integer n_A, integer kl, integer ku, void *ab, integer ldab,
                       integer *ipiv, integer datatype, integer n_repeats, double *time_min_,
                       integer *info, integer interfacetype, integer matrix_layout);
void invoke_gbtrf(integer datatype, integer *m, integer *n, integer *kl, integer *ku, void *ab,
                  integer *ldab, integer *ipiv, integer *info);
double prepare_lapacke_gbtrf_run(integer datatype, integer matrix_layout, integer m_A, integer n_A,
                                 integer kl, integer ku, void *ab, integer ldab, integer *ipiv,
                                 integer *info);
integer invoke_lapacke_gbtrf(integer datatype, int matrix_layout, integer m, integer n, integer kl,
                             integer ku, void *ab, integer ldab, integer *ipiv);

void fla_test_gbtrf(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LU factorization of banded matrix";
    char *front_str = "GBTRF";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';
    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gbtrf_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }
    if((argc == 9) || (argc == 10))
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
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
                fla_test_gbtrf_experiment(params, datatype, M, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, RECT_INPUT, M, N, residual,
                                      params->lin_solver_paramslist[0].solver_threshold, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gbtrf\n");
        printf("./<EXE> gbtrf <precisions - sdcz> <M> <N> <KL> <KU> <LDAB> <repeats>\n");
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

void fla_test_gbtrf_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer m, n, kl, ku, ldab;
    integer info = 0, vinfo = 0;
    void *IPIV;
    void *AB, *AB_test, *A;
    double time_min = 1e9;

    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;

    /* Determine the dimensions*/
    m = p_cur;
    n = q_cur;
    kl = params->lin_solver_paramslist[pci].kl;
    ku = params->lin_solver_paramslist[pci].ku;
    ldab = params->lin_solver_paramslist[pci].ldab;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
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

    /* Initialize the test matrices*/
    if(g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data from file or extreme values */
        init_matrix(datatype, AB, m, n, ldab, g_ext_fptr, params->imatrix_char);
    }
    else if(FLA_EXTREME_CASE_TEST)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, m);
        if((params->imatrix_char == 'A') || (params->imatrix_char == 'F'))
        {
            init_matrix_spec_rand_band_matrix_in(datatype, A, m, n, m, kl, ku,
                                                 params->imatrix_char);
        }
        else
        {
            init_matrix_spec_in(datatype, A, m, n, m, params->imatrix_char);
        }

        get_band_storage_matrix(datatype, m, n, kl, ku, A, m, AB, ldab);
        free_matrix(A);
    }
    else
    {
        /* Initialize & convert random band matrix into band storage as per API need */
        rand_band_storage_matrix(datatype, m, n, kl, ku, AB, ldab);
        /* Oveflow or underflow test initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_gbtrf(datatype, m, n, AB, ldab, params->imatrix_char);
        }
    }

    /* Save the original matrix*/
    copy_matrix(datatype, "full", ldab, n, AB, ldab, AB_test, ldab);

    /* call to API */
    prepare_gbtrf_run(m, n, kl, ku, AB_test, ldab, IPIV, datatype, n_repeats, &time_min, &info,
                      interfacetype, layout);

    /* execution time */
    *t = time_min;

    /* performance computation */
    *perf = (2.0 * n * (kl + ku + 1) * kl) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        *perf *= 4.0;
    }

    /* output validation */
    if((!FLA_EXTREME_CASE_TEST) && info == 0)
    {
        validate_gbtrf(m, n, kl, ku, AB, AB_test, ldab, IPIV, datatype, residual, &vinfo);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((info == 0) && !check_extreme_value(datatype, m, n, AB_test, ldab, params->imatrix_char))
        {
            *residual = DBL_MAX;
        }
    }
    else
    {
        FLA_TEST_CHECK_EINFO(residual, info, einfo);
    }

    /* Free up the buffers */
    free_matrix(AB);
    free_matrix(AB_test);
    free_vector(IPIV);
}

void prepare_gbtrf_run(integer m_A, integer n_A, integer kl, integer ku, void *AB, integer ldab,
                       integer *IPIV, integer datatype, integer n_repeats, double *time_min_,
                       integer *info, integer interfacetype, integer layout)
{
    integer i;
    void *AB_save;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &AB_save, ldab);
    copy_matrix(datatype, "full", ldab, n_A, AB, ldab, AB_save, ldab);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Copy original input data */
        copy_matrix(datatype, "full", ldab, n_A, AB, ldab, AB_save, ldab);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gbtrf_run(datatype, layout, m_A, n_A, kl, ku, AB_save, ldab,
                                                 IPIV, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)   /* Call CPP gbtrf API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_gbtrf(datatype, &m_A, &n_A, &kl, &ku, AB_save, &ldab, IPIV, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK gbtrf API */
            invoke_gbtrf(datatype, &m_A, &n_A, &kl, &ku, AB_save, &ldab, IPIV, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the ABFACT to matrix AB */
    copy_matrix(datatype, "full", ldab, n_A, AB_save, ldab, AB, ldab);
    free_matrix(AB_save);
}

double prepare_lapacke_gbtrf_run(integer datatype, integer layout, integer m_A, integer n_A,
                                 integer kl, integer ku, void *ab, integer ldab, integer *ipiv,
                                 integer *info)
{
    double exe_time;
    integer ldab_t = ldab;
    void *ab_t = NULL;

    ab_t = ab;

    if(layout == LAPACK_ROW_MAJOR)
    {
        ldab_t = n_A;
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, ldab, n_A, &ab_t, ldab_t);

        /* Convert column_major matrix layout to row_major matrix layout */
        convert_banded_matrix_layout(LAPACK_COL_MAJOR, datatype, m_A, n_A, ab, ldab, ab_t, ldab_t);
    }

    exe_time = fla_test_clock();

    /*  call LAPACKE gbtrf API */
    *info = invoke_lapacke_gbtrf(datatype, layout, m_A, n_A, kl, ku, ab_t, ldab_t, ipiv);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_banded_matrix_layout(LAPACK_ROW_MAJOR, datatype, m_A, n_A, ab_t, ldab_t, ab, ldab);

        /* free temporary buffers */
        free_matrix(ab_t);
    }

    return exe_time;
}

/*
 *  gbtrf_API calls LAPACK interface of
 *  LU factorization - GBTRF
 *  */
void invoke_gbtrf(integer datatype, integer *m, integer *n, integer *kl, integer *ku, void *ab,
                  integer *ldab, integer *ipiv, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
            break;
        }
    }
}

integer invoke_lapacke_gbtrf(integer datatype, int layout, integer m, integer n, integer kl,
                             integer ku, void *ab, integer ldab, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }
    }
    return info;
}
