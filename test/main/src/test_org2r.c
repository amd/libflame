/*
    Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include "test_common.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_org2r_experiment(test_params_t *params, integer datatype,
                               integer  p_cur, integer  q_cur, integer  pci,
                               integer  n_repeats, integer einfo, double* perf,
                               double* t,double* residual);
void prepare_org2r_run(integer m, integer n, void *A, integer lda, void *T,
                       void* work, integer datatype, integer n_repeats,
                       double* time_min_, integer *info);
void invoke_org2r(integer datatype, integer* m, integer* n, integer *min_A,
                  void* a, integer* lda, void* tau, void* work, integer* info);

void fla_test_org2r(integer argc, char ** argv, test_params_t *params)
{
    char* op_str = "QR factorization";
    char* front_str = "ORG2R";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_org2r_experiment);
        tests_not_run = 0;
    }
    if (argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if (argc >= 7 && argc <= 8)
    {
        integer i, num_types,N,M;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype,type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        g_lwork = -1;

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
                fla_test_org2r_experiment(params, datatype,
                                          M, N,
                                          0,
                                          n_repeats, einfo,
                                          &perf, &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str,
                                      stype,
                                      RECT_INPUT,
                                      M, N,
                                      residual, params->lin_solver_paramslist[0].solver_threshold,
                                      time_min, perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for org2r\n");
        printf("./<EXE> org2r <precisions - sdcz> <M> <N> <lda> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
    if (g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
    return;
}

void fla_test_org2r_experiment(test_params_t *params,
    integer datatype,
    integer p_cur,
    integer q_cur,
    integer pci,
    integer n_repeats,
    integer einfo,
    double* perf,
    double* time_min,
    double* residual)
{
    integer m, n, lda;
    void *A = NULL, *A_test = NULL, *T_test = NULL;
    void *work = NULL, *work_test = NULL;
    void *Q = NULL, *R = NULL;
    integer lwork = -1, info = 0, vinfo = 0;

    /* Get input matrix dimensions.*/
    m = p_cur;
    n = q_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    *time_min = 0.;
    *perf = 0.;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if (config_data)
    {
        if (lda == -1)
        {
            lda = fla_max(1,m);
        }
    }

    if(m >= n)
    {
        /* Create input matrix parameters */
        create_matrix(datatype, &A, lda, n);

        /* create tau vector */
        create_vector(datatype, &T_test, fla_min(m,n));

        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);

        /* Make a copy of input matrix A.
           This is required to validate the API functionality.*/
        create_matrix(datatype, &A_test, lda, n);
        copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

        /* create Q matrix to check orthogonality */
        create_matrix(datatype, &Q, lda, n);
        reset_matrix(datatype, m, n, Q, lda);
 
        /* Make a workspace query the first time. This will provide us with
           and ideal workspace size based on internal block size.*/
        if(g_lwork <= 0)
        {
            lwork = -1;
            create_vector(datatype, &work, 1);

            /* call to  geqrf API */
            invoke_geqrf(datatype, &m, &n, NULL, &lda, NULL, work, &lwork, &info);

            if(info == 0)
            {
                /* Get work size */
                lwork = get_work_value(datatype, work);
            }

            /* Output buffers will be freshly allocated for each iterations, free up
            the current output buffers.*/
            free_vector(work);
        }

        /* create work buffer */
        create_matrix(datatype, &work, lwork, 1);
        create_vector(datatype, &work_test, n);

        /* QR Factorisation on matrix A to generate Q and R */
        invoke_geqrf(datatype, &m, &n, A_test, &lda, T_test, work, &lwork, &info);

        create_matrix(datatype, &R, n, n);
        reset_matrix(datatype, n, n, R, n);
        copy_matrix(datatype, "Upper", n, n, A_test, lda, R, n);

        copy_matrix(datatype, "full", m, n, A_test, lda, Q, lda);

        /*invoke org2r API */
        prepare_org2r_run(m, n, Q, lda, T_test, work_test, datatype, n_repeats, time_min, &info);

        /* performance computation
           (2/3)*n2*(3m - n) */
        *perf = (double)((2.0 * m * n * n) - (( 2.0 / 3.0 ) * n * n * n )) / *time_min / FLOPS_PER_UNIT_PERF;
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            *perf *= 4.0;

        /* output validation */
        if(info == 0)
            validate_orgqr(m, n, A, lda, Q, R, work_test, datatype, residual, &vinfo);

        FLA_TEST_CHECK_EINFO(residual, info, einfo);

        /* Free up the buffers */
        free_matrix(A);
        free_matrix(A_test);
        free_matrix(work);
        free_vector(work_test);
        free_vector(T_test);
        free_matrix(Q);
        free_matrix(R);
    }
}

void prepare_org2r_run(integer m, integer n,
    void* A,
    integer lda,
    void* T,
    void* work,
    integer datatype,
    integer n_repeats,
    double* time_min_,
    integer *info)
{
    integer i;
    void *A_save = NULL;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, &A_save, lda, n);
    copy_matrix(datatype, "full", m, n, A, lda, A_save, lda);

    *info = 0;
    for (i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m, n, A_save, lda, A, lda);

        exe_time = fla_test_clock();

        /* Call to  org2r API */
        invoke_org2r(datatype, &m, &n, &n, A, &lda, T, work, info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

    }

    *time_min_ = time_min;

    free_matrix(A_save);
}


void invoke_org2r(integer datatype, integer* m, integer* n, integer *min_A,
                  void* a, integer* lda, void* tau, void* work, integer* info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sorg2r(m, n, n, a, lda, tau, work, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dorg2r(m, n, n, a, lda, tau, work, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cung2r(m, n, n, a, lda, tau, work, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zung2r(m, n, n, a, lda, tau, work, info);
            break;
        }
    }
}