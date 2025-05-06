/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"

#define GET_TRANS_STR(datatype) (((datatype) == FLOAT || (datatype) == DOUBLE) ? "T" : "C")

extern double perf;
extern double time_min;
/* Local prototypes.*/
void fla_test_sygvd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_sygvd_run(integer itype, char *jobz, char *uplo, integer n, void *A, integer lda,
                       void *B, integer ldb, void *w, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer test_lapacke_interface,
                       int layout);
void invoke_sygvd(integer datatype, integer *itype, char *jobz, char *uplo, integer *n, void *a,
                  integer *lda, void *b, integer *ldb, void *w, void *work, integer *lwork,
                  void *rwork, integer *lrwork, void *iwork, integer *liwork, integer *info);
double prepare_lapacke_sygvd_run(integer datatype, int itype, int layout, char *jobz, char *uplo,
                                 integer n, void *A, integer lda, void *B, integer ldb, void *w,
                                 integer *info);
integer invoke_lapacke_sygvd(integer datatype, int layout, int itype, char jobz, char uplo,
                             integer n, void *a, integer lda, void *b, integer ldb, void *w);

void fla_test_sygvd(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Eigen Decomposition";
    char *front_str = "SYGVD";
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
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_sygvd_experiment);
        tests_not_run = 0;
    }
    if(argc == 13)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[12]);
    }
    if(argc >= 12 && argc <= 13)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].itype = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].jobz = argv[4][0];
        params->eig_sym_paramslist[0].uplo = argv[5][0];
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);

        g_lwork = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        g_liwork = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        g_lrwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_sygvd_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for sygvd/hegvd\n");
        printf("./<EXE> sygvd <precisions - sd> <ITYPE> <JOBZ> <UPLO> <N> <LDA>"
               " <LWORK> <LIWORK> <LRWORK> <repeats>\n");
        printf("./<EXE> hegvd <precisions - cz> <ITYPE> <JOBZ> <UPLO> <N> <LDA>"
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

void fla_test_sygvd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, lda, info = 0, itype = 1;
    char jobz, uplo, range = 'R';
    void *A = NULL, *w = NULL, *A_test = NULL, *EVals = NULL;
    void *B = NULL, *B_test = NULL, *U = NULL, *C = NULL;
    void *scal = NULL, *temp = NULL;
    double residual, err_thresh;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    err_thresh = params->eig_sym_paramslist[pci].threshold_value;

    n = p_cur;
    lda = params->eig_sym_paramslist[pci].lda;

    itype = params->eig_sym_paramslist[pci].itype;

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
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B, lda);

    create_realtype_vector(datatype, &w, n);

    create_realtype_vector(get_realtype(datatype), &scal, 1);

    switch(get_realtype(datatype))
    {
        case FLOAT:
            *(float *)scal = s_one;
            break;
        case DOUBLE:
            *(double *)scal = d_one;
            break;
    }

    if(g_ext_fptr != NULL || (FLA_EXTREME_CASE_TEST))
    {
        /* Initialize input matrix with custom data */
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        init_matrix(datatype, B, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    else
    {

        create_realtype_vector(datatype, &EVals, n);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &U, lda);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &C, lda);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &temp, lda);

        /* Genrating random spd matrix with known chol factor */
        /* Evals is used here as temporary buffer */
        generate_matrix_from_EVs(datatype, 'V', n, U, lda, EVals, 0.1, 1.0, USE_ABS_EIGEN_VALUES);
        get_triangular_matrix("U", datatype, n, n, U, lda, 1, NON_UNIT_DIAG);

        /* B = U**{T|C} U */
        fla_invoke_gemm(datatype, GET_TRANS_STR(datatype), "N", &n, &n, &n, U, &lda, U, &lda, B,
                        &lda);

        /* Genarate matrix C such that C = Q * lambda * Q'
           where L is a diagonal matrix with diagonal values as eigen values
           and Q is orthogonal matrix denoting the eigen vectors
         */
        generate_matrix_from_EVs(datatype, range, n, C, lda, EVals, 0.0, 0.0, USE_ABS_EIGEN_VALUES);

        switch(itype)
        {
            case 1:
                /* For the first type A * X = B * X * lambda
                   Also B = U'U
                   A = U'U * X * lambda
                   inv(U') * A * X = U * X * lambda
                   (inv(U') * A * inv(U)) * (U * X) = (U * X) * lambda

                   C is generated as C * Q = Q * lambda
                   Comparing these two, A = U' * C * U
                */
                fla_invoke_trmm(datatype, "L", "U", GET_TRANS_STR(datatype), "N", &n, &n, U, &lda,
                                C, &lda);
                fla_invoke_trmm(datatype, "R", "U", "N", "N", &n, &n, U, &lda, C, &lda);
                copy_matrix(datatype, "full", n, n, C, lda, A, lda);
                break;
            case 2:
            case 3:
                /* Type 2: A * B * X = X * lambda
                   Type 3: B * A * X = X * lambda

                   Type 2:
                   A * U' * U * X = X * lambda
                   (U * A * U') * (U * X) = (U * X) * lambda

                   Type 3:
                   U' * U * A * X = X * lambda
                   U * A * X = inv(U') * X * lambda
                   (U * A * U') * (inv(U') * X) = (inv(U') * X) * lambda

                   Comparing for both types, A = inv(U) * C * inv(U')

                 */
                copy_matrix(datatype, "full", n, n, C, lda, temp, lda);
                /* temp = inv(U) * C */
                fla_invoke_trsm(datatype, "L", "U", "N", "N", &n, &n, U, &lda, temp, &lda);
                copy_matrix(datatype, "full", n, n, temp, lda, A, lda);
                /* A = inv(U) * C * inv(U') */
                fla_invoke_trsm(datatype, "R", "U", GET_TRANS_STR(datatype), "N", &n, &n, U, &lda,
                                A, &lda);
                break;
        }

        free_matrix(U);
        free_matrix(C);
        free_matrix(temp);
    }

    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_underflow_overflow_sygvd(datatype, n, A, lda, B, lda, itype,
                                              params->imatrix_char, scal);
    }

    /* Make a copy of input matrix A and B. This is required to validate the API functionality.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_test, lda);
    copy_matrix(datatype, "full", n, n, B, lda, B_test, lda);

    prepare_sygvd_run(itype, &jobz, &uplo, n, A_test, lda, B_test, lda, w, datatype, n_repeats,
                      &time_min, &info, test_lapacke_interface, layout);
    /* performance computation
       (8/3)n^3 [syevd] + (1/3)n^3 [potrf] flops for eigen vectors
       (4/3)n^3 [syevd] + (1/3)n^3 [potrf] flops for eigen values */
    if(jobz == 'V')
        perf = (double)(3.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    else
        perf = (double)((5.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_sygvd(tst_api, itype, &jobz, &range, &uplo, n, A, A_test, lda, B, B_test, lda, 0,
                       0, EVals, w, NULL, datatype, residual, params->imatrix_char, scal);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(B);
    free_matrix(B_test);

    free_vector(scal);
    free_vector(w);

    if((g_ext_fptr == NULL) && (!FLA_EXTREME_CASE_TEST))
        free_vector(EVals);
}

void prepare_sygvd_run(integer itype, char *jobz, char *uplo, integer n, void *A, integer lda,
                       void *B, integer ldb, void *w, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer test_lapacke_interface, int layout)
{
    void *A_save, *work, *iwork, *rwork = NULL, *B_save = NULL;
    integer lwork, liwork, lrwork;
    integer i;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix A and B. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_save, ldb);
    copy_matrix(datatype, "full", n, n, B, ldb, B_save, ldb);

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
        invoke_sygvd(datatype, &itype, jobz, uplo, &n, NULL, &lda, NULL, &ldb, NULL, work, &lwork,
                     rwork, &lrwork, iwork, &liwork, info);

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
        /* Restore input matrix A and B value and allocate memory to output buffers
           for each iteration */
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        copy_matrix(datatype, "full", n, n, B_save, ldb, B, ldb);

        create_vector(datatype, &work, lwork);
        create_vector(INTEGER, &iwork, liwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            create_realtype_vector(datatype, &rwork, lrwork);
        else
            rwork = NULL;

        /* Check if LAPACKE interface is enabled */
        if(test_lapacke_interface == 1)
        {
            exe_time = prepare_lapacke_sygvd_run(datatype, itype, layout, jobz, uplo, n, A, lda, B,
                                                 lda, w, info);
        }
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK API */
            invoke_sygvd(datatype, &itype, jobz, uplo, &n, A, &lda, B, &lda, w, work, &lwork, rwork,
                         &lrwork, iwork, &liwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(iwork);

        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
            free_vector(rwork);
    }

    *time_min_ = t_min;

    free_matrix(A_save);
    free_matrix(B_save);
}

double prepare_lapacke_sygvd_run(integer datatype, int itype, int layout, char *jobz, char *uplo,
                                 integer n, void *A, integer lda, void *B, integer ldb, void *w,
                                 integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldb_t = ldb;
    void *A_t = NULL, *B_t = NULL;
    A_t = A;
    B_t = B;

    if(layout == LAPACK_ROW_MAJOR)
    {
        lda_t = fla_max(1, n);
        ldb_t = fla_max(1, n);
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);

        create_matrix(datatype, layout, n, n, &B_t, ldb_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, B, ldb, B_t, ldb_t);
    }
    exe_time = fla_test_clock();

    /* call LAPACKE sygvd API */
    *info
        = invoke_lapacke_sygvd(datatype, layout, itype, *jobz, *uplo, n, A_t, lda_t, B_t, ldb_t, w);

    exe_time = fla_test_clock() - exe_time;
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        convert_matrix_layout(layout, datatype, n, n, B_t, ldb_t, B, ldb);
        /* free temporary buffers */
        free_matrix(A_t);
        free_matrix(B_t);
    }

    return exe_time;
}

void invoke_sygvd(integer datatype, integer *itype, char *jobz, char *uplo, integer *n, void *a,
                  integer *lda, void *b, integer *ldb, void *w, void *work, integer *lwork,
                  void *rwork, integer *lrwork, void *iwork, integer *liwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork,
                              info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork,
                              info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_chegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork,
                              iwork, liwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork,
                              iwork, liwork, info);
            break;
        }
    }
}

integer invoke_lapacke_sygvd(integer datatype, int layout, int itype, char jobz, char uplo,
                             integer n, void *a, integer lda, void *b, integer ldb, void *w)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssygvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsygvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_chegvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zhegvd(layout, itype, jobz, uplo, n, a, lda, b, ldb, w);
            break;
        }
    }
    return info;
}
