/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_potrf_lda;

/* Local prototypes.*/
void fla_test_potrf_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_potrf_run(char *uplo, integer m, void *A, integer lda, integer datatype,
                       integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                       int matrix_layout);
void invoke_potrf(char *uplo, integer datatype, integer *m, void *a, integer *lda, integer *info);
double prepare_lapacke_potrf_run(integer datatype, int matrix_layout, char *uplo, integer m,
                                 void *A, integer lda, integer *info);

void fla_test_potrf(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Cholesky factorization";
    char *front_str = "POTRF";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_potrf_experiment);
        tests_not_run = 0;
    }
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if(argc >= 7 && argc <= 8)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].Uplo = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_potrf_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_potrf_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for potrf\n");
        printf("./<EXE> potrf <precisions - sdcz> <Uplo> <N> <LDA> <repeats>\n");
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

void fla_test_potrf_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, lda;
    integer info = 0;
    void *A = NULL, *A_test = NULL;
    char uplo = params->lin_solver_paramslist[pci].Uplo;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;

    /* Get input matrix dimensions */
    m = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A, lda);

    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST && !FLA_OVERFLOW_UNDERFLOW_TEST))
    {
        /* Initialize input matrix with custom data */
        init_matrix(datatype, A, m, m, lda, g_ext_fptr, params->imatrix_char);
        if(params->imatrix_char != '\0')
        {
            char *type = "C";
            if(datatype == FLOAT || datatype == DOUBLE)
            {
                type = "S";
            }
            form_symmetric_matrix(datatype, m, A, lda, type, 'U');
        }
    }
    else
    {
        rand_spd_matrix(datatype, &uplo, A, m, lda);
        /* Oveflow or underflow test initialization */
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_overflow_underflow_potrf(datatype, m, A, lda, params->imatrix_char);
        }
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_test, lda);
    copy_matrix(datatype, "full", m, m, A, lda, A_test, lda);

    prepare_potrf_run(&uplo, m, A_test, lda, datatype, n_repeats, &time_min, &info, interfacetype,
                      layout);

    /* Compute the performance of the best experiment repeat */
    /* (1/3)m^3 for real and (4/3)m^3 for complex*/
    perf = (double)(1.0 / 3.0 * m * m * m) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_potrf(tst_api, &uplo, m, A, A_test, lda, datatype, residual);
    }
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, m, m, A_test, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, m, residual, err_thresh);
    }

    free_matrix(A);
    free_matrix(A_test);
}

void prepare_potrf_run(char *uplo, integer m, void *A, integer lda, integer datatype,
                       integer n_repeats, double *time_min_, integer *info, integer interfacetype,
                       int layout)
{
    void *A_save = NULL;
    double t_min = 1e9, exe_time;
    integer i;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_save, lda);
    copy_matrix(datatype, "full", m, m, A, lda, A_save, lda);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
        for each iteration */
        copy_matrix(datatype, "full", m, m, A_save, lda, A, lda);
        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_potrf_run(datatype, layout, uplo, m, A, lda, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP potrf API */
            invoke_cpp_potrf(uplo, datatype, &m, A, &lda, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK potrf API */
            invoke_potrf(uplo, datatype, &m, A, &lda, info);

            exe_time = fla_test_clock() - exe_time;
        }
        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);
    }

    *time_min_ = t_min;
    free_matrix(A_save);
}

double prepare_lapacke_potrf_run(integer datatype, int layout, char *uplo, integer m, void *A,
                                 integer lda, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, m, row_major_potrf_lda, lda_t);

    A_t = A;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, m, m, &A_t, fla_max(m, lda_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, m, m, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    *info = invoke_lapacke_potrf(datatype, layout, *uplo, m, A_t, lda_t);

    exe_time = fla_test_clock() - exe_time;
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, m, m, A_t, lda_t, A, lda);
        /* free temporary buffers */
        free_matrix(A_t);
    }

    return exe_time;
}

void invoke_potrf(char *uplo, integer datatype, integer *m, void *a, integer *lda, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_spotrf(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dpotrf(uplo, m, a, lda, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cpotrf(uplo, m, a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zpotrf(uplo, m, a, lda, info);
            break;
        }
    }
}

