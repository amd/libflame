/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"

#define GELS_VL 0.1
#define GELS_VU 10

void invoke_gels(integer datatype, char *trans, integer *m, integer *n, integer *nrhs, void *A,
                 integer *lda, void *B, integer *ldb, void *work, integer *lwork, integer *info);
void fla_test_gels_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual);
void prepare_gels_run(integer datatype, char trans, integer m, integer n, integer nrhs, void *A,
                      integer lda, void *B, integer ldb, void *work, integer lwork,
                      integer n_repeats, double *time_min_, integer *info);
void fla_test_gels(integer argc, char **argv, test_params_t *params)
{
    srand(4);
    char *op_str = "Solves overdetermined or underdetermined systems for GE matrices";
    char *front_str = "GELS";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        config_data = 1;
        g_lwork = -1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, LIN, fla_test_gels_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        integer i, num_types, M, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].nrhs = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].transr = argv[3][0];
        params->lin_solver_paramslist[0].lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ldb = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        g_lwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_gels_experiment(params, datatype, M, N, 0, n_repeats, einfo, &perf,
                                         &time_min, &residual);

                /* Print the result */
                fla_test_print_status(front_str, stype, RECT_INPUT, M, N, residual,
                                      params->lin_solver_paramslist[0].solver_threshold, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }
    /* Print error message */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for gels\n");
        printf("./<EXE> gels <precisions - sdcz> <TRANS> <M> <N> <NRHS> <LDA> <LDB> <LWORK> "
               "<repeats>\n");
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

void fla_test_gels_experiment(test_params_t *params, integer datatype, integer p_cur, integer q_cur,
                              integer pci, integer n_repeats, integer einfo, double *perf,
                              double *t, double *residual)
{
    integer m, n, m_b, nrhs, lda, ldb, lwork = -1, info = 0;
    void *A = NULL, *A_test = NULL, *B = NULL, *B_test = NULL, *work = NULL, *s_test = NULL;
    char trans, range = 'U';

    /* Determine the dimensions */
    m = p_cur;
    n = q_cur;
    nrhs = params->lin_solver_paramslist[pci].nrhs;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;
    *residual = params->lin_solver_paramslist[pci].solver_threshold;
    trans = params->lin_solver_paramslist[pci].transr;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
        if(ldb == -1)
        {
            ldb = fla_max(fla_max(1, m), n);
        }
    }

    /* Based on the value of TRANS the dimension of B changes.
     * Dimension of B is (m, nrhs) if TRANS = "N"
     * Dimension of B is (n, nrhs) if TRANS = "T" */

    m_b = n;
    if(trans == 'N' || trans == 'n')
    {
        trans = 'N';
        m_b = m;
    }
    if(trans == 'T' || trans == 't')
    {
        trans = 'T';
    }

    /* trans for complex number should be equal to 'C' (or 'c') while passing to the GEL api
     */
    if((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && (trans == 'T'))
    {
        trans = 'C';
    }

    /* Create the matrices for the current operation */
    create_matrix(datatype, &A, lda, n);
    create_matrix(datatype, &A_test, lda, n);
    create_matrix(datatype, &B, ldb, nrhs);
    create_matrix(datatype, &B_test, ldb, nrhs);
    create_realtype_vector(datatype, &s_test, fla_min(m, n));

    /* Initialize the test matrices */
    init_matrix(datatype, B, m_b, nrhs, ldb, g_ext_fptr, params->imatrix_char);

    if(params->imatrix_char == NULL && g_ext_fptr == NULL)
    {
        /* Generate input matrix with condition number <= 100 */
        create_svd_matrix(datatype, range, m, n, A, lda, s_test, GELS_VL, GELS_VU, i_zero, i_zero,
                          '\0', NULL, info);
    }
    else
    {
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    }

    /* Save the original matrix */
    copy_matrix(datatype, "full", lda, n, A, lda, A_test, lda);
    copy_matrix(datatype, "full", ldb, nrhs, B, ldb, B_test, ldb);

    /* call to API */
    prepare_gels_run(datatype, trans, m, n, nrhs, A_test, lda, B_test, ldb, work, lwork, n_repeats,
                     t, &info);

    /* Performance computation */
    if(m >= n)
    {
        *perf = (double)((n * n) * (2.0 / 3.0) * ((3 * m) - n)) / *t / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        *perf = (double)((m * m) * (2.0 / 3.0) * ((3 * n) - m)) / *t / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* Output validataion */
    if(!params->imatrix_char && info == 0)
        validate_gels(&trans, m, n, nrhs, A, lda, B, ldb, B_test, datatype, residual, &info);
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char))
           && (!check_extreme_value(datatype, m, n, B_test, ldb, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
        FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Free up buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B);
    free_matrix(B_test);
    free_vector(s_test);
}

void prepare_gels_run(integer datatype, char trans, integer m, integer n, integer nrhs, void *A,
                      integer lda, void *B, integer ldb, void *work, integer lwork,
                      integer n_repeats, double *time_min, integer *info)
{
    integer i;
    void *A_save = NULL, *B_save = NULL;
    double time_min_ = 1e9, exe_time;

    create_matrix(datatype, &A_save, lda, n);
    create_matrix(datatype, &B_save, ldb, nrhs);

    if(g_lwork <= 0)
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

        /* Getting lwork from api by passing lwork = -1 */
        invoke_gels(datatype, &trans, &m, &n, &nrhs, NULL, &lda, NULL, &ldb, work, &lwork, info);
        if(*info == 0)
        {
            lwork = get_work_value(datatype, work);
        }
        free_vector(work);
    }
    else
    {
        lwork = g_lwork;
    }

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; i++)
    {
        /* Copy original input */
        copy_matrix(datatype, "full", lda, n, A, lda, A_save, lda);
        copy_matrix(datatype, "full", ldb, nrhs, B, ldb, B_save, ldb);

        /* Create work buffer */
        create_vector(datatype, &work, lwork);

        exe_time = fla_test_clock();

        /*  call to API */
        invoke_gels(datatype, &trans, &m, &n, &nrhs, A_save, &lda, B_save, &ldb, work, &lwork,
                    info);

        exe_time = fla_test_clock() - exe_time;

        /* Get the best execution time */
        time_min_ = fla_min(time_min_, exe_time);

        free_vector(work);
    }
    *time_min = time_min_;

    /* Save the output to vector A */
    copy_matrix(datatype, "full", lda, n, A_save, lda, A, lda);
    copy_matrix(datatype, "full", ldb, nrhs, B_save, ldb, B, ldb);
    free_matrix(A_save);
    free_matrix(B_save);
}

void invoke_gels(integer datatype, char *trans, integer *m, integer *n, integer *nrhs, void *A,
                 integer *lda, void *B, integer *ldb, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
            break;
        }
    }
}