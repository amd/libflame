/*
    Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;
integer row_major_gghrd_lda;
integer row_major_gghrd_ldb;
integer row_major_gghrd_ldq;
integer row_major_gghrd_ldz;

/* Local prototypes */
void fla_test_gghrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_gghrd_run(char *compq, char *compz, integer n, integer *ilo, integer *ihi, void *a,
                       integer lda, void *b, integer ldb, void *q, integer ldq, void *z,
                       integer ldz, integer datatype, integer n_repeats, double *time_min_,
                       integer *info, integer interfacetype, int matrix_layout);
void invoke_gghrd(integer datatype, char *compq, char *compz, integer *n, integer *ilo,
                  integer *ihi, void *a, integer *lda, void *b, integer *ldb, void *q, integer *ldq,
                  void *z, integer *ldz, integer *info);
double prepare_lapacke_gghrd_run(integer datatype, int matrix_layout, char *compq, char *compz,
                                 integer n, integer *ilo, integer *ihi, void *a, integer lda,
                                 void *b, integer ldb, void *q, integer ldq, void *z, integer ldz,
                                 integer *info);
integer invoke_lapacke_gghrd(integer datatype, int matrix_layout, char compq, char compz, integer n,
                             integer ilo, integer ihi, void *a, integer lda, void *b, integer ldb,
                             void *q, integer ldq, void *z, integer ldz);

void fla_test_gghrd(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Reduces a pair matrices (A,B) to generalized upper Hessenberg form";
    char *front_str = "GGHRD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, LIN, fla_test_gghrd_experiment);
        tests_not_run = 0;
    }
    if(argc == 14)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[13]);
    }
    if(argc >= 13 && argc <= 14)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Prase the arguments */
        num_types = strlen(argv[2]);
        params->lin_solver_paramslist[0].compq_gghrd = argv[3][0];
        params->lin_solver_paramslist[0].compz_gghrd = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ilo = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->lin_solver_paramslist[0].ihi = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gghrd_lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_gghrd_ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            row_major_gghrd_ldq = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            row_major_gghrd_ldz = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
            params->lin_solver_paramslist[0].ldb = N;
            params->lin_solver_paramslist[0].ldq = N;
            params->lin_solver_paramslist[0].ldz = N;
        }
        else
        {
            params->lin_solver_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldb = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldq = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].ldz = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        }
        n_repeats = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);

        /* Skip if imatrix is I or F as API doesnot support the INF inputs. */
        if((params->imatrix_char != 'I') && (params->imatrix_char != 'F'))
        {
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
                    fla_test_gghrd_experiment(front_str, params, datatype, N, N, 0, n_repeats,
                                              einfo);
                    tests_not_run = 0;
                }
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for GGHRD\n");
        printf("./<EXE> gghrd <precisions - sdcz> <compq> <compz> <N> <ILO> <IHI> <LDA> <LDB> "
               "<LDQ> <LDZ> <repeats>\n");
        if((params->imatrix_char != 'N') || (params->imatrix_char != 'A'))
        {
            printf("imatrix should be N or A, as API does not support INF as inputs.\n");
        }
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

void fla_test_gghrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, ldz, lda, ldb, ldq;
    integer ilo, ihi, info = 0, ABInitialized = 0;
    void *A = NULL, *Z = NULL, *Q = NULL, *B = NULL, *A_test = NULL, *B_test = NULL, *Q_test = NULL,
         *Z_test = NULL, *temp = NULL, *A_ntest = NULL, *B_ntest = NULL, *Q_ntest = NULL,
         *Z_ntest = NULL;
    char compz, compq;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions. */
    n = p_cur;
    lda = params->lin_solver_paramslist[pci].lda;
    ldb = params->lin_solver_paramslist[pci].ldb;
    ldq = params->lin_solver_paramslist[pci].ldq;
    ldz = params->lin_solver_paramslist[pci].ldz;

    /* Initialize parameter */
    compz = params->lin_solver_paramslist[pci].compz_gghrd;
    compq = params->lin_solver_paramslist[pci].compq_gghrd;
    err_thresh = params->lin_solver_paramslist[pci].solver_threshold;
    ilo = params->lin_solver_paramslist[pci].ilo;
    ihi = params->lin_solver_paramslist[pci].ihi;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
        if(ldb == -1)
        {
            ldb = fla_max(1, n);
        }
        /* LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise */
        if(ldq == -1)
        {
            if((compq == 'V') || (compq == 'I'))
            {
                ldq = n;
            }
            else
            {
                ldq = 1;
            }
        }
        /* LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise */
        if(ldz == -1)
        {
            if((compz == 'V') || (compz == 'I'))
            {
                ldz = n;
            }
            else
            {
                ldz = 1;
            }
        }
    }

    /* Create input matrix parameters*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B, ldb);

    if(compq != 'N')
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q, ldq);
    if(compz != 'N')
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z, ldz);

    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
        init_matrix_from_file(datatype, B, n, n, ldb, g_ext_fptr);
        init_matrix_from_file(datatype, Q, n, n, ldq, g_ext_fptr);
        init_matrix_from_file(datatype, Z, n, n, ldz, g_ext_fptr);
    }
    else
    {
        /* Initialize A, B matrices. */
        if((params->imatrix_char == 'N') || (params->imatrix_char == 'A'))
        {
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
            init_matrix(datatype, B, n, n, ldb, g_ext_fptr, params->imatrix_char);
            ABInitialized = 1;
        }
        else if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            rand_matrix(datatype, A, n, n, lda);
            rand_matrix(datatype, B, n, n, ldb);
            scale_matrix_underflow_overflow_gghrd(datatype, n, A, lda, params->imatrix_char);
            scale_matrix_underflow_overflow_gghrd(datatype, n, B, ldb, params->imatrix_char);
            ABInitialized = 1;
        }
        get_triangular_matrix("U", datatype, n, n, B, ldb, ABInitialized, NON_UNIT_DIAG);
        get_generic_triangular_matrix(datatype, n, A, lda, ilo, ihi, ABInitialized);

        /* Initialize Q matrix. */
        if(compq == 'I')
        {
            set_identity_matrix(datatype, n, n, Q, ldq);
        }
        else if(compq == 'V')
        {
            get_orthogonal_matrix_from_QR(datatype, n, B, ldb, Q, ldq, &info);
        }

        /* Initialize Z matrix. */
        if(compz == 'I')
        {
            set_identity_matrix(datatype, n, n, Z, ldz);
        }
        else if(compz == 'V')
        {
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &temp, n);
            rand_matrix(datatype, temp, n, n, n);
            get_orthogonal_matrix_from_QR(datatype, n, temp, n, Z, ldz, &info);
            free_matrix(temp);
        }
    }

    /* Make copy of matrix A,B,Q and Z. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_test, ldb);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    copy_matrix(datatype, "full", n, n, B, ldb, B_test, ldb);
    if(compq != 'N')
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_test, ldq);
        copy_matrix(datatype, "full", n, n, Q, ldq, Q_test, ldq);
    }
    if(compz != 'N')
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_test, ldz);
        copy_matrix(datatype, "full", n, n, Z, ldz, Z_test, ldz);
    }

    prepare_gghrd_run(&compq, &compz, n, &ilo, &ihi, A_test, lda, B_test, ldb, Q_test, ldq, Z_test,
                      ldz, datatype, n_repeats, &time_min, &info, interfacetype, layout);

    /* If compq=N or/and compz=N, in addition to the first api call, also execute
       second api call with compq=I and compz=I. And validate the H and T matrices from
       the two calls */
    if(compq == 'N' || compz == 'N')
    {
        char compnq = 'I';
        char compnz = 'I';
        double time_minn;
        integer ldqn = fla_max(ldq, n);
        integer ldzn = fla_max(ldz, n);

        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_ntest, lda);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_ntest, ldb);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_ntest, ldqn);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_ntest, ldzn);
        copy_matrix(datatype, "full", n, n, A, lda, A_ntest, lda);
        copy_matrix(datatype, "full", n, n, B, ldb, B_ntest, ldb);
        set_identity_matrix(datatype, n, n, Q_ntest, ldqn);
        set_identity_matrix(datatype, n, n, Z_ntest, ldzn);
        prepare_gghrd_run(&compnq, &compnz, n, &ilo, &ihi, A_ntest, lda, B_ntest, ldb, Q_ntest,
                          ldqn, Z_ntest, ldzn, datatype, n_repeats, &time_minn, &info,
                          interfacetype, layout);
        free_matrix(Q_ntest);
        free_matrix(Z_ntest);
    }
    /* Performance computation
       (7)n^3 flops for eigen vectors for real
       (25)n^3 flops for eigen vectors for complex
       (10)n^3 flops for Schur form is computed for real
       (35)n^3 flops for Schur form is computed for complex
       (20)n^3 flops full Schur factorization is computed for real
       (70)n^3 flops full Schur factorization is computed for complex */

    if(compz == 'N')
    {
        if(datatype == FLOAT || datatype == DOUBLE)
            perf = (double)(7.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
        else
            perf = (double)(25.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    }
    else if(compz == 'I')
    {
        if(datatype == FLOAT || datatype == DOUBLE)
            perf = (double)(10.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
        else
            perf = (double)(35.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    }
    else
    {
        if(datatype == FLOAT || datatype == DOUBLE)
            perf = (double)(20.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
        else
            perf = (double)(70.0 * n * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    }

    /* Output Validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    if(!FLA_EXTREME_CASE_TEST)
    {
        if(compq == 'N' || compz == 'N')
        {
            /* Validation for compq=N or/and compz=N case */
            validate_gghrd(tst_api, &compq, &compz, n, A_test, A_ntest, lda, B_test, B_ntest,
                           ldb, Q, Q_test, ldq, Z, Z_test, ldz, datatype, residual);
        }
        else
        {
            /* Validation for other cases */
            validate_gghrd(tst_api, &compq, &compz, n, A, A_test, lda, B, B_test, ldb, Q,
                           Q_test, ldq, Z, Z_test, ldz, datatype, residual);
        }
    }
    else
    {
        if((params->imatrix_char == 'N') || (params->imatrix_char == 'A'))
        {
            if((!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char))
               && (!check_extreme_value(datatype, n, n, B_test, ldb, params->imatrix_char)))
            {
                residual = DBL_MAX;
            }
            else
            {
                residual = err_thresh;
            }
            FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
        }
        else
        {
            printf("Extreme Value tests not supported for xGGHRD APIs\n");
        }
    }

    /* Free up the buffers */
    free_matrix(A);
    free_matrix(B);
    free_matrix(Q);
    free_matrix(Z);
    free_matrix(A_test);
    free_matrix(B_test);
    if(compq != 'N')
    {
        free_matrix(Q_test);
    }
    if(compz != 'N')
    {
        free_matrix(Z_test);
    }
    if(compq == 'N' || compz == 'N')
    {
        free_matrix(A_ntest);
        free_matrix(B_ntest);
    }
}

void prepare_gghrd_run(char *compq, char *compz, integer n, integer *ilo, integer *ihi, void *A,
                       integer lda, void *B, integer ldb, void *Q, integer ldq, void *Z,
                       integer ldz, integer datatype, integer n_repeats, double *time_min_,
                       integer *info, integer interfacetype, int layout)
{
    void *A_save = NULL, *B_save = NULL, *Q_save = NULL, *Z_save = NULL;
    integer i;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix A,B,Q and Z. Same input values will be passed in each
     * itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_save, ldb);
    if(*compq != 'N')
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_save, ldq);
        copy_matrix(datatype, "full", n, n, Q, ldq, Q_save, ldq);
    }
    if(*compz != 'N')
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_save, ldz);
        copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);
    }
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
    copy_matrix(datatype, "full", n, n, B, ldb, B_save, ldb);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A,B,Q and Z value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        copy_matrix(datatype, "full", n, n, B_save, ldb, B, ldb);
        if(*compq != 'N')
            copy_matrix(datatype, "full", n, n, Q_save, ldq, Q, ldq);
        if(*compz != 'N')
            copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gghrd_run(datatype, layout, compq, compz, n, ilo, ihi, A,
                                                 lda, B, ldb, Q, ldq, Z, ldz, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* Call CPP gghrd API */
            invoke_cpp_gghrd(datatype, compq, compz, &n, ilo, ihi, A, &lda, B, &ldb, Q, &ldq, Z,
                             &ldz, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK gghrd API */
            invoke_gghrd(datatype, compq, compz, &n, ilo, ihi, A, &lda, B, &ldb, Q, &ldq, Z, &ldz,
                         info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);
    }
    *time_min_ = t_min;

    free_matrix(A_save);
    free_matrix(B_save);
    free_matrix(Q_save);
    free_matrix(Z_save);
}

double prepare_lapacke_gghrd_run(integer datatype, int layout, char *compq, char *compz, integer n,
                                 integer *ilo, integer *ihi, void *A, integer lda, void *B,
                                 integer ldb, void *Q, integer ldq, void *Z, integer ldz,
                                 integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldb_t = ldb;
    integer ldq_t = ldq;
    integer ldz_t = ldz;
    void *A_t = NULL, *B_t = NULL, *Q_t = NULL, *Z_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_gghrd_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_gghrd_ldb, ldb_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_gghrd_ldq, ldq_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_gghrd_ldz, ldz_t);

    A_t = A;
    B_t = B;
    Q_t = Q;
    Z_t = Z;
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, fla_max(n, lda_t));
        create_matrix(datatype, layout, n, n, &B_t, fla_max(n, ldb_t));

        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, B, ldb, B_t, ldb_t);
        if(*compq != 'N')
        {
            create_matrix(datatype, layout, n, n, &Q_t, fla_max(n, ldq_t));
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, Q, ldq, Q_t, ldq_t);
        }
        if(*compz != 'N')
        {
            create_matrix(datatype, layout, n, n, &Z_t, fla_max(n, ldz_t));
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, Z, ldz, Z_t, ldz_t);
        }
    }
    exe_time = fla_test_clock();

    /* Call to LAPACKE gghrd API */
    *info = invoke_lapacke_gghrd(datatype, layout, *compq, *compz, n, *ilo, *ihi, A_t, lda_t, B_t,
                                 ldb_t, Q_t, ldq_t, Z_t, ldz_t);
    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        convert_matrix_layout(layout, datatype, n, n, B_t, ldb_t, B, ldb);

        if(*compq != 'N')
        {
            convert_matrix_layout(layout, datatype, n, n, Q_t, ldq_t, Q, ldq);
            free_matrix(Q_t);
        }
        if(*compz != 'N')
        {
            convert_matrix_layout(layout, datatype, n, n, Z_t, ldz_t, Z, ldz);
            free_matrix(Z_t);
        }
        /* free temporary buffers */
        free_matrix(A_t);
        free_matrix(B_t);
    }
    return exe_time;
}

void invoke_gghrd(integer datatype, char *compq, char *compz, integer *n, integer *ilo,
                  integer *ihi, void *a, integer *lda, void *b, integer *ldb, void *q, integer *ldq,
                  void *z, integer *ldz, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info);
            break;
        }
    }
}

integer invoke_lapacke_gghrd(integer datatype, int layout, char compq, char compz, integer n,
                             integer ilo, integer ihi, void *a, integer lda, void *b, integer ldb,
                             void *q, integer ldq, void *z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info
                = LAPACKE_sgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }

        case DOUBLE:
        {
            info
                = LAPACKE_dgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }

        case COMPLEX:
        {
            info
                = LAPACKE_cgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info
                = LAPACKE_zgghrd(layout, compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz);
            break;
        }
    }
    return info;
}
