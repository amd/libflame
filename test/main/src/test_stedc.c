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
integer row_major_stedc_ldz;

/* Local prototypes. */
void fla_test_stedc_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_stedc_run(char *compz, integer n, void *D, void *E, void *Z, integer ldz,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int matrix_layout);
void invoke_stedc(integer datatype, char *compz, integer *n, void *D, void *E, void *Z,
                  integer *ldz, void *work, integer *lwork, void *rwork, integer *lrwork,
                  integer *iwork, integer *liwork, integer *info);
double prepare_lapacke_stedc_run(integer datatype, int matrix_layout, char *compz, integer n,
                                 void *D, void *E, void *Z, integer ldz, integer *info);

#define STEDC_VL 0.1
#define STEDC_VU 1000

void fla_test_stedc(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Eigenvalues/eigenvectors of symmetric tridiagonal matrix";
    char *front_str = "STEDC";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        g_lwork = -1;
        g_liwork = -1;
        g_lrwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_stedc_experiment);
        tests_not_run = 0;
    }
    if(argc == 11)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[10]);
    }
    if((argc == 10) || (argc == 11))
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].compz = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_stedc_ldz = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldz = N;
        }
        else
        {
            params->eig_sym_paramslist[0].ldz = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }

        g_lwork = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        g_liwork = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        g_lrwork = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_stedc_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("Invalid arguments for stedc\n");
        printf("Usage: ./<EXE> stedc <precisions - sdcz> <COMPZ> <N> <LDZ> <LWORk> <LIWORK> "
               "<LRWORK> <repeats>\n");
    }
    else if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n");
    }
    if(g_ext_fptr != NULL)
    {
        fclose(g_ext_fptr);
        g_ext_fptr = NULL;
    }
}

void fla_test_stedc_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, info = 0, realtype, lda, ldz;
    void *D = NULL, *D_test = NULL, *E = NULL, *E_test = NULL, *Z_test = NULL;
    void *Z = NULL, *A = NULL, *L = NULL, *Q = NULL, *scal = NULL;
    char compz, uplo, range = 'V';
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions. */
    n = q_cur;

    /* Initialize parameter needed for STEDC() call. */
    compz = params->eig_sym_paramslist[pci].compz;
    ldz = params->eig_sym_paramslist[pci].ldz;
    err_thresh = params->eig_sym_paramslist[pci].threshold_value;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(ldz == -1)
        {
            if((same_char(compz, 'N') == 0) && (layout != LAPACK_ROW_MAJOR))
            {
                ldz = 1;
            }
            else
            {
                ldz = fla_max(1, n);
            }
        }
    }

    lda = fla_max(ldz, n);

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z, ldz);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_test, ldz);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q, lda);

    reset_matrix(datatype, lda, n, A, lda);
    reset_matrix(datatype, ldz, n, Z, ldz);
    reset_matrix(datatype, ldz, n, Z_test, ldz);
    reset_matrix(datatype, lda, n, Q, lda);

    realtype = get_realtype(datatype);
    create_vector(realtype, &D, n);
    create_vector(realtype, &E, n - 1);

    if(g_ext_fptr || FLA_EXTREME_CASE_TEST)
    {
        /* Initialize input matrix with random/custom data */
        init_matrix(realtype, D, 1, n, 1, g_ext_fptr, params->imatrix_char);
        init_matrix(realtype, E, 1, n - 1, 1, g_ext_fptr, params->imatrix_char);

        if(same_char(compz, 'V'))
        {
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
            if(g_ext_fptr)
            {
                copy_matrix(datatype, "full", n, n, A, lda, Q, lda);
                /* Input matrix is assumed to be orthogonal matrix for file inputs.
                 * So get the symmetric/hermitian matrix using:
                 * A = Q * T * (Q**T) */
                copy_sym_tridiag_matrix(datatype, D, E, n, n, Z, ldz);
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, Q, &lda, Z, &ldz, Z_test, &ldz);
                fla_invoke_gemm(datatype, "N", "T", &n, &n, &n, Z_test, &ldz, Q, &lda, A, &lda);
            }
            else if(FLA_EXTREME_CASE_TEST)
            {
                /* Get the symmetric/hermitian matrix.*/
                form_symmetric_matrix(datatype, n, A, lda, "C", 'U');
                /* Initialize Q matrix */
                init_matrix(datatype, Q, n, n, lda, g_ext_fptr, params->imatrix_char);
            }
        }
    }
    else
    {
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, range, n, A, lda, L, STEDC_VL, STEDC_VU,
                                 USE_ABS_EIGEN_VALUES);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_realtype_vector(get_datatype(datatype), &scal, n);
            scale_matrix_underflow_overflow_stedc(datatype, n, A, lda, &params->imatrix_char, scal);
        }
        copy_matrix(datatype, "full", n, n, A, lda, Q, lda);
        uplo = 'U';
        invoke_sytrd(datatype, &uplo, compz, n, Q, lda, D, E, &info);
    }
    if(same_char(compz, 'I'))
    {
        set_identity_matrix(datatype, n, n, Z_test, ldz);
        copy_sym_tridiag_matrix(datatype, D, E, n, n, Z, ldz);
    }
    else if(same_char(compz, 'V'))
    {
        copy_matrix(datatype, "full", n, n, Q, lda, Z_test, ldz);
        copy_matrix(datatype, "full", n, n, A, lda, Z, ldz);
    }

    create_vector(realtype, &D_test, n);
    copy_vector(realtype, n, D, 1, D_test, 1);
    create_vector(realtype, &E_test, n - 1);
    copy_vector(realtype, n - 1, E, 1, E_test, 1);

    prepare_stedc_run(&compz, n, D_test, E_test, Z_test, ldz, datatype, n_repeats, &time_min, &info,
                      interfacetype, layout);

    /* Performance computation
       (6)n^3 flops for eigen vectors
       (4/3)n^3 flops for eigen values. */
    perf = (double)((4.0 / 3.0) * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(!same_char(compz, 'N'))
    {
        perf += (double)(6 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    }
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        perf *= 2.0;
    }

    /* Output validation. */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        validate_syev(tst_api, &compz, &range, n, Z, Z_test, ldz, 0, 0, L, D_test, NULL, datatype,
                      residual, params->imatrix_char, scal);
    }
    /* Check for output matrix & vectors when inputs are extreme values */
    else
    {
        residual = err_thresh;
        if(same_char(compz, 'N'))
        {
            if(!check_extreme_value(datatype, n, 1, D_test, 1, params->imatrix_char))
            {
                residual = DBL_MAX;
            }
        }
        else
        {
            if(!check_extreme_value(datatype, n, n, Z_test, lda, params->imatrix_char)
               && !check_extreme_value(datatype, n, 1, D_test, 1, params->imatrix_char))
            {
                residual = DBL_MAX;
            }
        }
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

    /* Free up buffers. */
    free_matrix(Z);
    free_matrix(Z_test);
    free_vector(D_test);
    free_vector(E_test);
    free_matrix(A);
    free_vector(D);
    free_vector(E);
    if(L != NULL)
    {
        free_vector(L);
    }
    free_matrix(Q);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal);
    }
}

void prepare_stedc_run(char *compz, integer n, void *D, void *E, void *Z, integer ldz,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int layout)
{
    integer index, lwork, liwork, lrwork, realtype;
    void *D_save = NULL, *E_save = NULL, *E_test = NULL, *Z_save = NULL;
    void *work = NULL, *iwork = NULL, *rwork = NULL;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrices. Same input values will be passed in
       each itertaion.*/
    if(same_char(*compz, 'V'))
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_save, ldz);
        copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);
    }
    realtype = get_realtype(datatype);
    create_vector(realtype, &D_save, n);
    copy_vector(realtype, n, D, 1, D_save, 1);
    create_vector(realtype, &E_save, n - 1);
    copy_vector(realtype, n - 1, E, 1, E_save, 1);

    /* Call to STEDC() API to get work buffers size.
       NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork == -1 || g_liwork == -1
           || ((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && g_lrwork == -1)))
    {
        /* Make a workspace query the first time. This will provide us with
        and ideal workspace size based on internal block size.*/
        create_vector(datatype, &work, 1);
        create_vector(realtype, &rwork, 1);
        create_vector(INTEGER, &iwork, 1);
        lwork = g_lwork;
        liwork = g_liwork;
        lrwork = g_lrwork;
#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_stedc(datatype, compz, &n, D, E_test, Z, &ldz, work, &lwork, rwork, &lrwork,
                             iwork, &liwork, info);
        }
        else
#endif
        {
            invoke_stedc(datatype, compz, &n, D, E_test, Z, &ldz, work, &lwork, rwork, &lrwork,
                         iwork, &liwork, info);
        }

        /* Get work buffers size. */
        if(*info == 0)
        {
            if(g_lwork == -1)
            {
                lwork = get_work_value(datatype, work);
            }
            if(g_liwork == -1)
            {
                liwork = get_work_value(INTEGER, iwork);
            }
            if((datatype == COMPLEX || datatype == DOUBLE_COMPLEX) && g_lrwork == -1)
            {
                lrwork = get_work_value(realtype, rwork);
            }
        }

        free_vector(work);
        free_vector(rwork);
        free_vector(iwork);
    }
    else
    {
        lwork = g_lwork;
        liwork = g_liwork;
        lrwork = g_lrwork;
    }
    create_vector(datatype, &work, lwork);
    if((datatype == COMPLEX) || (datatype == DOUBLE_COMPLEX))
    {
        create_vector(realtype, &rwork, lrwork);
    }
    create_vector(INTEGER, &iwork, liwork);
    create_vector(realtype, &E_test, n - 1);

    *info = 0;
    for(index = 0; index < n_repeats && *info == 0; ++index)
    {
        /* Restore input matrices and allocate memory to output buffers
           for each iteration. */
        if(same_char(*compz, 'V'))
        {
            copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        }
        copy_vector(realtype, n, D_save, 1, D, 1);
        copy_vector(realtype, n - 1, E_save, 1, E_test, 1);
        reset_vector(datatype, work, lwork, 1);
        if((datatype == COMPLEX) || (datatype == DOUBLE_COMPLEX))
        {
            reset_vector(realtype, rwork, lrwork, 1);
        }
        reset_vector(INTEGER, iwork, liwork, 1);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time
                = prepare_lapacke_stedc_run(datatype, layout, compz, n, D, E_test, Z, ldz, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP stedc API */
            invoke_cpp_stedc(datatype, compz, &n, D, E_test, Z, &ldz, work, &lwork, rwork, &lrwork,
                             iwork, &liwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK stedc API */
            invoke_stedc(datatype, compz, &n, D, E_test, Z, &ldz, work, &lwork, rwork, &lrwork,
                         iwork, &liwork, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time. */
        t_min = fla_min(t_min, exe_time);
    }
    *time_min_ = t_min;

    /* Free up buffers. */
    if(*compz == 'V')
    {
        free_vector(Z_save);
    }
    free_vector(D_save);
    free_vector(E_save);
    free_vector(E_test);
    free_vector(work);
    if((datatype == COMPLEX) || (datatype == DOUBLE_COMPLEX))
    {
        free_vector(rwork);
    }
    free_vector(iwork);
}

double prepare_lapacke_stedc_run(integer datatype, int layout, char *compz, integer n, void *D,
                                 void *E, void *Z, integer ldz, integer *info)
{
    double exe_time;
    integer ldz_t = ldz;
    void *Z_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_stedc_ldz, ldz_t);

    Z_t = Z;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if((!same_char(*compz, 'N')) && (layout == LAPACK_ROW_MAJOR))
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &Z_t, fla_max(n, ldz_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, Z, ldz, Z_t, ldz_t);
    }

    exe_time = fla_test_clock();

    /* Call LAPACKE STEDC() API. */
    *info = invoke_lapacke_stedc(datatype, layout, *compz, n, D, E, Z_t, ldz_t);

    exe_time = fla_test_clock() - exe_time;
    if((!same_char(*compz, 'N')) && (layout == LAPACK_ROW_MAJOR))
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, Z_t, ldz_t, Z, ldz);
        /* free temporary buffers */
        free_matrix(Z_t);
    }

    return exe_time;
}

void invoke_stedc(integer datatype, char *compz, integer *n, void *D, void *E, void *Z,
                  integer *ldz, void *work, integer *lwork, void *rwork, integer *lrwork,
                  integer *iwork, integer *liwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sstedc(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dstedc(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cstedc(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, iwork, liwork,
                              info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zstedc(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, iwork, liwork,
                              info);
            break;
        }
    }
}
