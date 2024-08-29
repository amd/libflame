/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"

/* Local prototypes.*/
void fla_test_syevd_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_syevd_run(char *jobz, char *uplo, integer n, void *A, integer lda, void *w,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer test_lapacke_interface, int matrix_layout);
void invoke_syevd(integer datatype, char *jobz, char *uplo, integer *n, void *a, integer *lda,
                  void *w, void *work, integer *lwork, void *rwork, integer *lrwork, void *iwork,
                  integer *liwork, integer *info);
double prepare_lapacke_syevd_run(integer datatype, int matrix_layout, char *jobz, char *uplo,
                                 integer n, void *A, integer lda, void *w, integer *info);
integer invoke_lapacke_syevd(integer datatype, int matrix_layout, char jobz, char uplo,
                             integer n, void *a, integer lda, void *w);

#define SYEVD_VL 1
#define SYEVD_VU 5

void fla_test_syevd(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Eigen Decomposition";
    char *front_str = "SYEVD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        g_liwork = -1;
        g_lrwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_syevd_experiment);
        tests_not_run = 0;
    }
    if(argc == 12)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[11]);
    }
    if(argc >= 11 && argc <= 12)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].jobz = argv[3][0];
        params->eig_sym_paramslist[0].uplo = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        g_liwork = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        g_lrwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->eig_sym_paramslist[0].threshold_value = CLI_NORM_THRESH;

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
                fla_test_syevd_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
                                          &time_min, &residual);
                /* Print the results */
                fla_test_print_status(front_str, stype, SQUARE_INPUT, N, N, residual,
                                      params->eig_sym_paramslist[0].threshold_value, time_min,
                                      perf);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for syevd/heevd\n");
        printf("./<EXE> syevd <precisions - sd> <JOBZ> <UPLO> <N> <LDA>"
               " <LWORK> <LIWORK> <LRWORK> <repeats>\n");
        printf("./<EXE> heevd <precisions - cz> <JOBZ> <UPLO> <N> <LDA>"
               " <LWORK> <LIWORK> <LRWORK> <repeats>\n");
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

void fla_test_syevd_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *time_min, double *residual)
{
    integer n, lda, info = 0;
    char jobz, uplo, range = 'R';
    void *A = NULL, *w = NULL, *A_test = NULL, *L = NULL;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    *residual = params->eig_sym_paramslist[pci].threshold_value;

    n = p_cur;
    lda = params->eig_sym_paramslist[pci].lda;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_realtype_vector(datatype, &w, n);

    if(g_ext_fptr != NULL || (params->imatrix_char))
    {
        /* Initialize input matrix with custom data */
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {
        /*  Creating input matrix A by generating random eigen values.
            When range = V, generate EVs in given range (vl,vu)  */
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, range, n, A, lda, L, SYEVD_VL, n * SYEVD_VU);
    }
    /* Make a copy of input matrix A. This is required to validate the API functionality.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_syevd_run(&jobz, &uplo, n, A_test, lda, w, datatype, n_repeats, time_min, &info,
                      test_lapacke_interface, layout);

    /* performance computation
       (8/3)n^3 flops for eigen vectors
       (4/3)n^3 flops for eigen values */
    if(jobz == 'V')
        *perf = (double)((8.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    else
        *perf = (double)((4.0 / 3.0) * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf *= 4.0;

    /* output validation */
    if((info == 0) && (!FLA_EXTREME_CASE_TEST))
    {
        validate_syev(&jobz, &range, n, A, A_test, lda, 0, 0, L, w, NULL, datatype, residual,
                      params->imatrix_char, NULL);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char)))
        {
            *residual = DBL_MAX;
        }
    }
    else
    {
        FLA_TEST_CHECK_EINFO(residual, info, einfo);
    }
    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(w);
    if((g_ext_fptr == NULL) && !(params->imatrix_char))
        free_vector(L);
}

void prepare_syevd_run(char *jobz, char *uplo, integer n, void *A, integer lda, void *w,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer test_lapacke_interface, int layout)
{
    void *A_save, *w_test, *work, *iwork, *rwork = NULL;
    integer lwork, liwork, lrwork;
    integer i;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    if((test_lapacke_interface == 0)
       && (g_lwork <= 0 || ((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && g_lrwork <= 0)
           || g_liwork <= 0))
    {
        lwork = -1;
        liwork = -1;
        lrwork = -1;

        create_vector(datatype, &work, 1);
        create_vector(INTEGER, &iwork, 1);
        create_realtype_vector(datatype, &rwork, 1);
        /* call to  syevd API */
        invoke_syevd(datatype, jobz, uplo, &n, NULL, &lda, NULL, work, &lwork, rwork, &lrwork,
                     iwork, &liwork, info);

        /* Get work size */
        if(*info == 0)
        {
            lwork = get_work_value(datatype, work);
            liwork = get_work_value(INTEGER, iwork);
            lrwork = get_work_value(datatype, rwork);
        }

        free_vector(work);
        free_vector(iwork);
        free_vector(rwork);
    }
    else
    {
        lwork = g_lwork;
        liwork = g_liwork;
        lrwork = g_lrwork;
    }

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);

        create_realtype_vector(datatype, &w_test, n);
        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, liwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, lrwork);
        else
            rwork = NULL;

        /* Check if LAPACKE interface is enabled */
        if(test_lapacke_interface == 1)
        {
            exe_time
                = prepare_lapacke_syevd_run(datatype, layout, jobz, uplo, n, A, lda, w_test, info);
        }
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK API */
            invoke_syevd(datatype, jobz, uplo, &n, A, &lda, w_test, work, &lwork, rwork, &lrwork,
                         iwork, &liwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Make a copy of the output buffers. This is required to validate the API functionality.*/
        copy_realtype_vector(datatype, n, w_test, 1, w, 1);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);

        free_vector(w_test);
    }

    *time_min_ = time_min;

    free_matrix(A_save);
}

double prepare_lapacke_syevd_run(integer datatype, int layout, char *jobz, char *uplo,
                                 integer n, void *A, integer lda, void *w, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;
    A_t = A;

    if(layout == LAPACK_ROW_MAJOR)
    {
        lda_t = fla_max(1, n);
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
    }
    exe_time = fla_test_clock();

    /* call LAPACKE syevd API */
    *info = invoke_lapacke_syevd(datatype, layout, *jobz, *uplo, n, A_t, lda_t, w);

    exe_time = fla_test_clock() - exe_time;
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        /* free temporary buffers */
        free_matrix(A_t);
    }

    return exe_time;
}

void invoke_syevd(integer datatype, char *jobz, char *uplo, integer *n, void *a, integer *lda,
                  void *w, void *work, integer *lwork, void *rwork, integer *lrwork, void *iwork,
                  integer *liwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork,
                              info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork,
                              info);
            break;
        }
    }
}

integer invoke_lapacke_syevd(integer datatype, int layout, char jobz, char uplo, integer n,
                             void *a, integer lda, void *w)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssyevd(layout, jobz, uplo, n, a, lda, w);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsyevd(layout, jobz, uplo, n, a, lda, w);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cheevd(layout, jobz, uplo, n, a, lda, w);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zheevd(layout, jobz, uplo, n, a, lda, w);
            break;
        }
    }
    return info;
}
