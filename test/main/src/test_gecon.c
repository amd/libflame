/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

integer row_major_gecon_lda;

void fla_test_gecon_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_gecon_run(integer datatype, char *norm, integer n, void *A, integer lda, void *anorm,
                       void *rcond, void *work, void *lrwork, integer n_repeats, double *time_min,
                       integer *info, integer interfacetype, integer layout);
double prepare_lapacke_gecon_run(integer datatype, integer layout, char norm, integer n, void *A,
                                 integer lda, void *anorm, void *rcond, integer *info);
void invoke_gecon(integer datatype, char *norm, integer *n, void *A, integer *lda, void *anorm,
                  void *rcond, void *work, void *lrwork, integer *info);
integer invoke_lapacke_gecon(integer datatype, integer layout, char norm, integer n, void *A,
                             integer lda, void *anorm, void *rcond);

void fla_test_gecon(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Estimates the reciprocal of the condition number of a general matrix A";
    char *front_str = "GECON";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gecon_experiment);
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
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].norm_gbcon = argv[3][0];
        params->lin_solver_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);

        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gecon_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
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

                /* Check for invalid datatype */
                if(datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }

                /* Check for duplicate datatype */
                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                /* Call the test code */
                fla_test_gecon_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);

                /* Print the result */
                fla_test_print_status(front_str, stype, SQUARE_INPUT, N, N, residual,
                                      params->lin_solver_paramslist[0].solver_threshold, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }
    /* Print error message */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gecon\n");
        printf("./<EXE> gecon <precisions - sdcz> <NORM> <N> <LDA> <repeats>\n");
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

void fla_test_gecon_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual)
{
    integer n, lda, info = 0;
    void *A = NULL, *work = NULL, *rcond = NULL, *anorm = NULL, *lrwork = NULL, *ipiv = NULL,
         *s_test_in = NULL, *A_save = NULL;
    char norm;
    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major, getrfinfo = 0;

    /* Determine the dimensions */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    norm = params->lin_solver_paramslist[pci].norm_gbcon;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
    }

    /* Create the matrices for the current operation */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    create_vector(datatype, &work, n);
    if(datatype == FLOAT || datatype == DOUBLE)
        create_vector(INTEGER, &lrwork, n);
    else if(get_realtype(datatype) == FLOAT)
        create_vector(FLOAT, &lrwork, 2 * n);
    else if(get_realtype(datatype) == DOUBLE)
        create_vector(DOUBLE, &lrwork, 2 * n);
    create_realtype_vector(datatype, &rcond, 1);
    create_realtype_vector(datatype, &anorm, 1);
    create_realtype_vector(datatype, &s_test_in, n);

    /* Initialize the test matrices */
    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        compute_matrix_norm(datatype, norm, n, n, A, lda, anorm, norm, work);
    }
    else
    { /*Generating specific input matrix */
        create_svd_matrix(datatype, 'U', n, n, A, lda, s_test_in, 0.1, 100, i_zero, i_zero,
                          getrfinfo);
        compute_matrix_norm(datatype, norm, n, n, A, lda, anorm, norm, work);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            scale_matrix_underflow_overflow_getrf(datatype, n, n, A, lda, params->imatrix_char);
        }
        create_vector(INTEGER, &ipiv, n);
        getrfinfo = 0;
        invoke_getrf(datatype, &n, &n, A, &lda, ipiv, &getrfinfo);
        copy_matrix(datatype, "Full", n, n, A, lda, A_save, lda);
    }
    /* Save the original matrix */

    /* call to API */
    prepare_gecon_run(datatype, &norm, n, A, lda, anorm, rcond, work, lrwork, n_repeats, t, &info,
                      interfacetype, layout);

    /* Performance computation */

    *perf = (double)(2 * (n * n)) / *t / FLOPS_PER_UNIT_PERF;

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Output validataion */
    if(info == 0 && !FLA_EXTREME_CASE_TEST)
    {
        validate_gecon(datatype, norm, n, A, A_save, lda, residual, params->imatrix_char);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if(!check_extreme_value(datatype, n, n, A, lda, params->imatrix_char))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up buffers */
    free_matrix(A);
    free_matrix(A_save);
    free_vector(ipiv);
    free_vector(work);
    free_vector(s_test_in);
    free_vector(anorm);
    free_vector(rcond);
    free_vector(lrwork);
}

void prepare_gecon_run(integer datatype, char *norm, integer n, void *A, integer lda, void *anorm,
                       void *rcond, void *work, void *lrwork, integer n_repeats, double *time_min,
                       integer *info, integer interfacetype, integer layout)
{
    integer i, lwork = 4 * n;
    void *A_save = NULL;
    double time_min_ = 1e9, exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);

    /* Workspace size calculations for complex datatype */
    if(datatype == COMPLEX && datatype == DOUBLE_COMPLEX)
    {
        lwork = 2 * n;
    }

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; i++)
    {
        /* Copy original input */
        copy_matrix(datatype, "full", lda, n, A, lda, A_save, lda);

        /* Create work buffer */
        create_vector(datatype, &work, lwork);
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            /*call to lapacke gecon*/
            exe_time = prepare_lapacke_gecon_run(datatype, layout, *norm, n, A_save, lda, anorm,
                                                 rcond, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP gecon API */
            invoke_cpp_gecon(datatype, norm, &n, A_save, &lda, anorm, rcond, work, lrwork, info);

            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();

            /*  call to gecon API */
            invoke_gecon(datatype, norm, &n, A_save, &lda, anorm, rcond, work, lrwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min_ = fla_min(time_min_, exe_time);

        free_vector(work);
    }
    *time_min = time_min_;

    /* Save the output to vector A */
    copy_matrix(datatype, "full", lda, n, A_save, lda, A, lda);

    /* Free up buffers */
    free_matrix(A_save);
}

double prepare_lapacke_gecon_run(integer datatype, integer layout, char norm, integer n, void *A,
                                 integer lda, void *anorm, void *rcond, integer *info)
{
    double exe_time = 0;
    void *A_t = NULL;
    integer lda_t = lda;
    A_t = A;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_gecon_lda, lda_t);
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, fla_max(n, lda_t));

        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE gecon API */
    *info = invoke_lapacke_gecon(datatype, layout, norm, n, A_t, lda_t, anorm, rcond);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */

    if((layout == LAPACK_ROW_MAJOR))
    {
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);

        free_matrix(A_t);
    }
    return exe_time;
}

/*
LAPACK gecon API invoke function
*/
void invoke_gecon(integer datatype, char *norm, integer *n, void *A, integer *lda, void *anorm,
                  void *rcond, void *work, void *lrwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgecon(norm, n, A, lda, (float *)anorm, (float *)rcond, work,
                              (integer *)lrwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dgecon(norm, n, A, lda, (double *)anorm, (double *)rcond, work,
                              (integer *)lrwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cgecon(norm, n, A, lda, (float *)anorm, (float *)rcond, work,
                              (float *)lrwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgecon(norm, n, A, lda, (double *)anorm, (double *)rcond, work,
                              (double *)lrwork, info);
            break;
        }
    }
}

/*
LAPACKE gecon API invoke function
*/

integer invoke_lapacke_gecon(integer datatype, integer layout, char norm, integer n, void *A,
                             integer lda, void *anorm, void *rcond)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgecon(layout, norm, n, A, lda, *(float *)anorm, (float *)rcond);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgecon(layout, norm, n, A, lda, *(double *)anorm, (double *)rcond);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgecon(layout, norm, n, A, lda, *(float *)anorm, (float *)rcond);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgecon(layout, norm, n, A, lda, *(double *)anorm, (double *)rcond);
            break;
        }
    }
    return info;
}