/*
    Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"

integer row_major_steqr_ldz;

/* Local prototypes.*/
void fla_test_steqr_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *t, double *residual);
void prepare_steqr_run(char *compz, integer n, void *Z, integer ldz, void *D, void *E,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer test_lapacke_interface, int matrix_layout);
void invoke_steqr(integer datatype, char *compz, integer *n, void *z, integer *ldz, void *d,
                  void *e, void *work, integer *info);
double prepare_lapacke_steqr_run(integer datatype, int matrix_layout, char *compz, integer n,
                                 void *Z, integer ldz, void *D, void *E, integer *info);
integer invoke_lapacke_steqr(integer datatype, int matrix_layout, char compz, integer n,
                             void *d, void *e, void *z, integer ldz);

#define STEQR_VL 0.1
#define STEQR_VU 1000

void fla_test_steqr(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Eigen Decomposition of symmetrix tridiagonal matrix";
    char *front_str = "STEQR";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        config_data = 1;
        /* Test with parameters from config */
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_steqr_experiment);
        tests_not_run = 0;
    }
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if(argc >= 7 && argc <= 8)
    {
        /* Test with parameters from commandline */
        integer i, num_types, N;
        integer datatype, n_repeats;
        double perf, time_min, residual;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].compz = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && params->test_lapacke_interface
           && (params->matrix_major == LAPACK_ROW_MAJOR))
        {
            row_major_steqr_ldz = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldz = N;
        }
        else
        {
            params->eig_sym_paramslist[0].ldz = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);

        if(n_repeats > 0)
        {
            params->eig_sym_paramslist[0].threshold_value = CLI_NORM_THRESH;
            params->eig_sym_paramslist[0].uplo = 'L';

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
                fla_test_steqr_experiment(params, datatype, N, N, 0, n_repeats, einfo, &perf,
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
        printf("Invalid arguments for STEQR\n");
        printf("Usage: ./<EXE> steqr <precisions - sdcz> <COMPZ> <N> <LDZ> <repeats>\n");
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
    return;
}

void fla_test_steqr_experiment(test_params_t *params, integer datatype, integer p_cur,
                               integer q_cur, integer pci, integer n_repeats, integer einfo,
                               double *perf, double *time_min, double *residual)
{
    integer n, ldz, lda, info = 0, realtype;
    char compz, uplo, range = 'V';
    void *Z = NULL, *Z_test = NULL, *A = NULL, *Q = NULL;
    void *D = NULL, *D_test = NULL, *E = NULL, *E_test = NULL;
    void *L = NULL, *scal = NULL;
    double resid;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    compz = params->eig_sym_paramslist[pci].compz;
    *residual = params->eig_sym_paramslist[pci].threshold_value;
    uplo = params->eig_sym_paramslist[pci].uplo;

    n = p_cur;

    ldz = params->eig_sym_paramslist[pci].ldz;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(ldz == -1)
        {
            if((compz == 'N') && (layout != LAPACK_ROW_MAJOR))
            {
                ldz = 1;
            }
            else
            {
                ldz = fla_max(1, n);
            }
        }
    }

    lda = fla_max(n, ldz);
    realtype = get_realtype(datatype);

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z, ldz);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_test, ldz);

    reset_matrix(datatype, n, n, Z, ldz);
    reset_matrix(datatype, n, n, A, lda);
    reset_matrix(datatype, n, n, Q, lda);

    create_vector(realtype, &D, n);
    create_vector(realtype, &E, n - 1);

    if(g_ext_fptr != NULL)
    {
        /* Initialize input matrix with custom data */
        init_matrix(realtype, D, 1, n, 1, g_ext_fptr, params->imatrix_char);
        init_matrix(realtype, E, 1, n - 1, 1, g_ext_fptr, params->imatrix_char);
        if(compz == 'V')
        {
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
            copy_matrix(datatype, "full", n, n, A, lda, Q, lda);
            resid = check_orthogonality(datatype, Q, n, n, lda);
            if(resid < *residual)
            {
                /* Input matrix is orthogonal matrix for file inputs.
                    So get the symmetric/hermitian matrix using:
                    A = Q * T * (Q**T) */
                /* Form tridiagonal matrix Z by copying from matrix.*/
                copy_sym_tridiag_matrix(datatype, D, E, n, n, Z, ldz);
                fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, Q, &lda, Z, &ldz, Z_test, &ldz);
                fla_invoke_gemm(datatype, "N", "T", &n, &n, &n, Z_test, &ldz, Q, &lda, A, &lda);
            }
            else
            {
                *residual = DBL_MAX;
                printf("Orthogonality check failed for input matrix Q.\n");
                return;
            }
        }
    }
    else
    {
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, range, n, A, lda, L, STEQR_VL, STEQR_VU);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_realtype_vector(get_datatype(datatype), &scal, n);
            scale_matrix_underflow_overflow_steqr(datatype, n, A, lda, &params->imatrix_char, scal);
        }
        copy_matrix(datatype, "full", n, n, A, lda, Q, lda);
        invoke_sytrd(datatype, &uplo, compz, n, Q, lda, D, E, &info);
    }
    if(compz == 'I')
    {
        set_identity_matrix(datatype, n, n, Z_test, ldz);
        /* Form tridiagonal matrix Z by copying from D, E.*/
        copy_sym_tridiag_matrix(datatype, D, E, n, n, Z, ldz);
    }
    else if(compz == 'V')
    {
        copy_matrix(datatype, "full", n, n, A, lda, Z, ldz);
        copy_matrix(datatype, "full", n, n, Q, lda, Z_test, ldz);
    }

    create_vector(realtype, &D_test, n);
    create_vector(realtype, &E_test, n - 1);
    copy_vector(realtype, n, D, 1, D_test, 1);
    copy_vector(realtype, n - 1, E, 1, E_test, 1);

    prepare_steqr_run(&compz, n, Z_test, ldz, D_test, E_test, datatype, n_repeats, time_min, &info,
                      test_lapacke_interface, layout);

    /* performance computation
       24 n^2 flops for eigen vectors of Z, compz = 'N'
       7 n^3 flops for eigen vectors of Z, compz = 'V' or 'I'
       14 n^3 flops for eigen vectors of Z for complex, compz = 'V' or 'I' */

    if(compz == 'I' || compz == 'V')
        *perf = (double)(7.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    else if(compz == 'N')
        *perf = (double)(24.0 * n * n) / *time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        *perf = (double)(14.0 * n * n * n) / *time_min / FLOPS_PER_UNIT_PERF;

    /* Output validation */
    if(info == 0)
    {
        validate_syev(&compz, &range, n, Z, Z_test, lda, 0, 0, L, D_test, NULL, datatype, residual,
                      params->imatrix_char, scal);
    }
    else
    {
        FLA_TEST_CHECK_EINFO(residual, info, einfo);
    }

    /* Free up the buffers */
    free_matrix(Z);
    free_vector(D);
    free_vector(E);
    free_matrix(A);
    free_matrix(Q);
    free_vector(D_test);
    free_vector(E_test);
    if(L != NULL)
    {
        free_vector(L);
    }
    free_matrix(Z_test);
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal);
    }
}

void prepare_steqr_run(char *compz, integer n, void *Z, integer ldz, void *D, void *E,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer test_lapacke_interface, int layout)
{
    void *Z_save = NULL, *D_save = NULL, *E_save = NULL, *work = NULL;
    integer i, realtype;
    double time_min = 1e9, exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    if(*compz != 'N')
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_save, ldz);
        copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);
    }
    realtype = get_realtype(datatype);
    create_vector(realtype, &D_save, n);
    create_vector(realtype, &E_save, n - 1);
    copy_vector(realtype, n, D, 1, D_save, 1);
    copy_vector(realtype, n - 1, E, 1, E_save, 1);

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        if(*compz != 'N')
        {
            copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        }
        copy_vector(realtype, n, D_save, 1, D, 1);
        copy_vector(realtype, n - 1, E_save, 1, E, 1);

        create_vector(realtype, &work, 2 * (n - 1));
        if(test_lapacke_interface == 1)
        {
            exe_time = prepare_lapacke_steqr_run(datatype, layout, compz, n, Z, ldz, D, E, info);
        }
        else
        {
            exe_time = fla_test_clock();

            /* call to API */
            invoke_steqr(datatype, compz, &n, Z, &ldz, D, E, work, info);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        time_min = fla_min(time_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }

    *time_min_ = time_min;

    if(*compz != 'N')
    {
        free_matrix(Z_save);
    }
    free_matrix(D_save);
    free_matrix(E_save);
}

double prepare_lapacke_steqr_run(integer datatype, int layout, char *compz, integer n, void *Z,
                                 integer ldz, void *D, void *E, integer *info)
{
    double exe_time;
    integer ldz_t = ldz;
    void *Z_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_steqr_ldz, ldz_t);

    Z_t = Z;

    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if((*compz != 'N') && (layout == LAPACK_ROW_MAJOR))
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &Z_t, fla_max(n, ldz_t));
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, Z, ldz, Z_t, ldz_t);
    }

    exe_time = fla_test_clock();

    /* call to LAPACKE steqr API */
    *info = invoke_lapacke_steqr(datatype, layout, *compz, n, D, E, Z_t, ldz_t);

    exe_time = fla_test_clock() - exe_time;
    if((*compz != 'N') && (layout == LAPACK_ROW_MAJOR))
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, Z_t, ldz_t, Z, ldz);
        /* free temporary buffers */
        free_matrix(Z_t);
    }

    return exe_time;
}

void invoke_steqr(integer datatype, char *compz, integer *n, void *z, integer *ldz, void *d,
                  void *e, void *work, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssteqr(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsteqr(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_csteqr(compz, n, d, e, z, ldz, work, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zsteqr(compz, n, d, e, z, ldz, work, info);
            break;
        }
    }
}

integer invoke_lapacke_steqr(integer datatype, int layout, char compz, integer n, void *d,
                             void *e, void *z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_csteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zsteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
    }
    return info;
}