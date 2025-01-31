/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"

#define GETRFNP_VL 0.1
#define GETRFNP_VU 10

/* Local prototypes */
void fla_test_getrfnp_experiment(test_params_t *params, integer datatype, integer p_cur,
                                 integer q_cur, integer pci, integer n_repeats, integer einfo,
                                 double *perf, double *t, double *residual);
void prepare_getrfnp_run(integer m_A, integer n_A, void *A, integer lda, integer datatype,
                         integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                         int matrix_layout);
void invoke_getrfnp(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *info);

void fla_test_getrfnp(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "LU factorization";
    char *front_str = "GETRFNP";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_getrfnp_experiment);
        tests_not_run = 0;
    }
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if(argc >= 7 && argc <= 8)
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
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_getrfnp_experiment(params, datatype, M, N, 0, n_repeats, einfo, &perf,
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
        printf("\nIllegal arguments for getrfnp\n");
        printf("./<EXE> getrfnp <precisions - sdcz> <M> <N> <LDA> <repeats>\n");
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

void fla_test_getrfnp_experiment(test_params_t *params, integer datatype, integer p_cur,
                                 integer q_cur, integer pci, integer n_repeats, integer einfo,
                                 double *perf, double *t, double *residual)
{
    integer m, n, lda, info = 0, vinfo = 0, i__, min_mn, max_mn;
    void *IPIV = NULL, *A = NULL, *A_test = NULL, *s_test = NULL;
    void *A_copy;
    char range = 'U';
    double time_min = 1e9;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Determine the dimensions*/
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
    }

    /* Create the matrices for the current operation*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_vector(INTEGER, &IPIV, fla_min(m, n));
    create_realtype_vector(datatype, &s_test, fla_min(m, n));

    /* Initialize the test matrices*/
    if(g_ext_fptr != NULL || FLA_EXTREME_CASE_TEST)
    {
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    }
    /* If the lda is less than m, then do not initilize the input matrix.
        The invalid param error should be reported by the API */
    else if(lda >= m)
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, m, n, A, lda, s_test, GETRFNP_VL, GETRFNP_VU, i_zero,
                          i_zero, info);
        create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_copy, lda);
        copy_matrix(datatype, "full", m, n, A, lda, A_copy, lda);

        /* Invoke getrf to get the optimial permuation vector */
        invoke_getrf(datatype, &m, &n, A_copy, &lda, IPIV, &info);

        /* Swap rows of A as per computed permutation matrix
           to avoid error in LU factorization */
        swap_rows_with_pivot(datatype, m, n, A, lda, IPIV);

        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_getrf(datatype, m, n, A, lda, params->imatrix_char);
        }
        free_matrix(A_copy);
    }

    /* Save the original matrix*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);
    /* Copy data if only lda >= m */
    if(lda >= m)
    {
        copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);
    }

    /* call to API */
    prepare_getrfnp_run(m, n, A_test, lda, datatype, n_repeats, &time_min, &info, interfacetype,
                        layout);

    /* execution time */
    *t = time_min;

    /* performance computation */
    min_mn = fla_min(m, n);
    max_mn = fla_max(m, n);

    *perf = ((1.0 / 3.0) * min_mn * min_mn * (3.0 * max_mn - min_mn)) / time_min
            / FLOPS_PER_UNIT_PERF;

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Fill IPIV specifiying that no permutation has been done */
    for(i__ = 0; i__ < fla_min(m, n); ++i__)
    {
        ((integer *)IPIV)[i__] = i__ + 1;
    }

    /* output validation */
    if((!FLA_EXTREME_CASE_TEST) && info == 0)
    {
        validate_getrf(m, n, A, A_test, lda, IPIV, datatype, residual, &vinfo,
                       params->imatrix_char);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(IPIV);
    free_vector(s_test);
}

void prepare_getrfnp_run(integer m_A, integer n_A, void *A, integer lda, integer datatype,
                         integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                         int layout)
{
    integer i;
    void *A_save;
    double time_min = 1e9, exe_time;

    /* Save the original matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {

        /* Copy original input data */
        copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

        exe_time = fla_test_clock();
        /* Call LAPACK getrf API */
        invoke_getrfnp(datatype, &m_A, &n_A, A_save, &lda, info);
        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);
    }

    *time_min_ = time_min;
    /*  Save the AFACT to matrix A */
    copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);
    free_matrix(A_save);
}

/*
 *  GETRFNP_API calls LAPACK interface of
 *  */
void invoke_getrfnp(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgetrfnp(m, n, a, lda, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgetrfnp(m, n, a, lda, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgetrfnp(m, n, a, lda, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgetrfnp(m, n, a, lda, info);
            break;
        }
    }
}
