/*
    Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_syev_lda;

/* Local prototypes.*/
void fla_test_syev_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo);
void prepare_syev_run(char *jobz, char *uplo, integer n, void *A, integer lda, void *w,
                      integer datatype, integer n_repeats, double *time_min_, integer *info,
                      integer test_lapacke_interface, int matrix_layout);
void invoke_syev(integer datatype, char *jobz, char *uplo, integer *n, void *a, integer *lda,
                 void *w, void *work, integer *lwork, void *rwork, integer *info);
double prepare_lapacke_syev_run(integer datatype, int matrix_layout, char *jobz, char *uplo,
                                integer n, void *A, integer lda, void *w, integer *info);

#define SYEV_VL 1
#define SYEV_VU 5

void fla_test_syev(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Eigen Values and Vectors";
    char *front_str = "SYEV";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_syev_experiment);
        tests_not_run = 0;
    }
    if(argc == 10)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[9]);
    }
    if(argc >= 9 && argc <= 10)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].jobz = argv[3][0];
        params->eig_sym_paramslist[0].uplo = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && params->test_lapacke_interface
           && (params->matrix_major == LAPACK_ROW_MAJOR))
        {
            row_major_syev_lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].lda = N;
        }
        else
        {
            params->eig_sym_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_syev_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for syev/heev\n");
        printf("./<EXE> syev <precisions - sd> <JOBZ> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
        printf("./<EXE> heev <precisions - cz> <JOBZ> <UPLO> <N> <LDA> <LWORK> <repeats>\n");
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
}

void fla_test_syev_experiment(char *tst_api, test_params_t *params, integer datatype, integer p_cur,
                              integer q_cur, integer pci, integer n_repeats, integer einfo)
{
    integer n, lda, info = 0;
    char jobz, uplo, range = 'V';
    void *A = NULL, *w = NULL, *A_test = NULL, *L = NULL, *scal = NULL;
    double residual, err_thresh;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    err_thresh = params->eig_sym_paramslist[pci].threshold_value;

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
    if(g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
    }
    else
    {
        /*  Creating input matrix A by generating random eigen values.
            When range = V, generate EVs in given range (vl,vu)  */
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, range, n, A, lda, L, SYEV_VL, n * SYEV_VU,
                                 USE_ABS_EIGEN_VALUES);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_realtype_vector(get_datatype(datatype), &scal, n);
            scale_matrix_underflow_overflow_syev(datatype, n, A, lda, &params->imatrix_char, scal);
        }
    }
    /* Make a copy of input matrix A. This is required to validate the API functionality.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    prepare_syev_run(&jobz, &uplo, n, A_test, lda, w, datatype, n_repeats, &time_min, &info,
                     test_lapacke_interface, layout);

    /* performance computation
       (8/3)n^3 flops for eigen vectors
       (4/3)n^3 flops for eigen values */
    if(jobz == 'V')
        perf = (double)((8.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    else
        perf = (double)((4.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_syev(tst_api, &jobz, &range, n, A, A_test, lda, 0, 0, L, w, NULL, datatype,
                      residual, params->imatrix_char, scal);
    }
    else
    {
        printf("Extreme Value tests not supported for xSYEV/HEEV APIs\n");
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_vector(w);
    if(g_ext_fptr == NULL)
    {
        free_vector(L);
    }
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal);
    }
}

void prepare_syev_run(char *jobz, char *uplo, integer n, void *A, integer lda, void *w,
                      integer datatype, integer n_repeats, double *time_min_, integer *info,
                      integer test_lapacke_interface, int layout)
{
    void *A_save = NULL, *work = NULL, *rwork = NULL, *w_test = NULL;
    integer i, lwork;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.
       NOTE: LAPACKE interface handles workspace query internally */
    if((test_lapacke_interface == 0) && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* call to  syev API */
        invoke_syev(datatype, jobz, uplo, &n, NULL, &lda, NULL, work, &lwork, rwork, info);
        /* Get work size */
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
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        create_realtype_vector(datatype, &w_test, n);
        create_vector(datatype, &work, lwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, fla_max(1, 3 * n - 2));
        else
            rwork = NULL;

        /* Check if LAPACKE interface is enabled */
        if(test_lapacke_interface == 1)
        {
            exe_time
                = prepare_lapacke_syev_run(datatype, layout, jobz, uplo, n, A, lda, w_test, info);
        }
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK syev API */
            invoke_syev(datatype, jobz, uplo, &n, A, &lda, w_test, work, &lwork, rwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        /* Make a copy of the output buffers. This is required to validate the API functionality.*/
        copy_realtype_vector(datatype, n, w_test, 1, w, 1);

        /* Free up the output buffers */
        free_vector(work);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);

        free_vector(w_test);
    }

    *time_min_ = t_min;

    free_matrix(A_save);
}

double prepare_lapacke_syev_run(integer datatype, int layout, char *jobz, char *uplo, integer n,
                                void *A, integer lda, void *w_test, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_syev_lda, lda_t);

    A_t = A;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, fla_max(n, lda_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
    }
    exe_time = fla_test_clock();

    /* call LAPACKE syev API */
    *info = invoke_lapacke_syev(datatype, layout, *jobz, *uplo, n, A_t, lda_t, w_test);

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

void invoke_syev(integer datatype, char *jobz, char *uplo, integer *n, void *a, integer *lda,
                 void *w, void *work, integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssyev(jobz, uplo, n, a, lda, w, work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
            break;
        }
    }
}
