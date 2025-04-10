/*
    Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_hgeqz_ldh;
integer row_major_hgeqz_ldt;
integer row_major_hgeqz_ldq;
integer row_major_hgeqz_ldz;

/* Local prototypes */
void fla_test_hgeqz_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_hgeqz_run(char *job, char *compq, char *compz, integer n, integer *ilo, integer *ihi,
                       void *h, integer ldh, void *t, integer ldt, void *alpha, void *alphar,
                       void *alphai, void *beta, void *q, integer ldq, void *z, integer ldz,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int matrix_layout);
void invoke_hgeqz(integer datatype, char *job, char *compq, char *compz, integer *n, integer *ilo,
                  integer *ihi, void *h, integer *ldh, void *t, integer *ldt, void *alpha,
                  void *alphar, void *alphai, void *beta, void *q, integer *ldq, void *z,
                  integer *ldz, void *work, integer *lwork, void *rwork, integer *info);
double prepare_lapacke_hgeqz_run(integer datatype, int matrix_layout, char *job, char *compq,
                                 char *compz, integer n, integer *ilo, integer *ihi, void *h,
                                 integer ldh, void *t, integer ldt, void *alpha, void *alphar,
                                 void *alphai, void *beta, void *q, integer ldq, void *z,
                                 integer ldz, integer *info);

void fla_test_hgeqz(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Computing Eigen value of a real matrix pair (H,T)";
    char *front_str = "HGEQZ";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_hgeqz_experiment);
        tests_not_run = 0;
    }
    if(argc == 16)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[15]);
    }
    if(argc >= 15 && argc <= 16)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Prase the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].job_seqr = argv[3][0];
        params->eig_sym_paramslist[0].compq_hgeqz = argv[4][0];
        params->eig_sym_paramslist[0].compz_hgeqz = argv[5][0];
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ilo = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ihi = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_hgeqz_ldh = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            row_major_hgeqz_ldt = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            row_major_hgeqz_ldq = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            row_major_hgeqz_ldz = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].lda = N;
            params->eig_sym_paramslist[0].ldb = N;
            params->eig_sym_paramslist[0].ldq = N;
            params->eig_sym_paramslist[0].ldz = N;
        }
        else
        {
            params->eig_sym_paramslist[0].lda = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldb = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldq = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldz = strtoimax(argv[12], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[14], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_hgeqz_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for HGEQZ\n");
        printf("./<EXE> hgeqz <precisions - sdcz> <job> <compq> <compz> <N> <ILO> <IHI> <LDH> "
               "<LDT> <LDQ> <LDZ> <LWORK> <repeats>\n");
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

void fla_test_hgeqz_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, ldz, ldh, ldt, ldq;
    integer ilo, ihi, info = 0;
    void *H = NULL, *Z = NULL, *Q = NULL, *T = NULL, *H_test = NULL, *T_test = NULL, *A = NULL;
    void *B = NULL, *Z_A = NULL, *Q_test = NULL, *Z_test = NULL, *Q_temp = NULL, *tmp = NULL,
         *Q_ntest = NULL, *Z_ntest = NULL;
    void *alpha = NULL, *alphar = NULL, *alphai = NULL, *beta = NULL, *scal_A = NULL,
         *scal_B = NULL, *Q_A = NULL;
    char compz, compq, job;
    void *alphae = NULL, *alphaer = NULL, *alphaei = NULL, *betae = NULL;
    void *H_ntest = NULL, *T_ntest = NULL;
    void *alphan = NULL, *alphanr = NULL, *alphani = NULL, *betan = NULL;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;

    /* Get input matrix dimensions. */
    n = p_cur;
    ldh = params->eig_sym_paramslist[pci].lda;
    ldt = params->eig_sym_paramslist[pci].ldb;
    ldq = params->eig_sym_paramslist[pci].ldq;
    ldz = params->eig_sym_paramslist[pci].ldz;

    /* Initialize parameters */
    job = params->eig_sym_paramslist[pci].job_seqr;
    compz = params->eig_sym_paramslist[pci].compz_hgeqz;
    compq = params->eig_sym_paramslist[pci].compq_hgeqz;
    err_thresh = params->eig_sym_paramslist[pci].threshold_value;
    ilo = params->eig_sym_paramslist[pci].ilo;
    ihi = params->eig_sym_paramslist[pci].ihi;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(ldh == -1)
        {
            ldh = fla_max(1, n);
        }
        if(ldt == -1)
        {
            ldt = fla_max(1, n);
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

    /* Create input matrix */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &H, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &T, ldt);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q, ldq);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z, ldz);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B, ldt);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_A, ldq);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_A, ldz);
    if(datatype == FLOAT || datatype == DOUBLE)
    {
        create_vector(datatype, &alphar, n);
        create_vector(datatype, &alphai, n);
    }
    else
    {
        create_vector(datatype, &alpha, n);
    }
    create_vector(datatype, &beta, n);

    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, H, n, n, ldh, g_ext_fptr);
        init_matrix_from_file(datatype, T, n, n, ldt, g_ext_fptr);
        init_matrix_from_file(datatype, Q, n, n, ldq, g_ext_fptr);
        init_matrix_from_file(datatype, Z, n, n, ldz, g_ext_fptr);
    }
    else
    {
        /* Convert matrix according to ILO and IHI values */
        get_generic_triangular_matrix(datatype, n, A, ldh, ilo, ihi, false);
        /* Initialize matrix with random values */
        rand_matrix(datatype, B, n, n, ldt);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_realtype_vector(datatype, &scal_A, n);
            create_realtype_vector(datatype, &scal_B, n);
            scale_matrix_overflow_underflow_hgeqz_A(datatype, n, A, ldh, params->imatrix_char,
                                                    scal_A);
            scale_matrix_overflow_underflow_hgeqz_B(datatype, n, B, ldt, params->imatrix_char,
                                                    scal_B);
        }
        /* Decompose matrix B in to QR and store orthogonal matrix in Q and R in B */
        if(compq == 'N')
        {
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_temp, ldt);
            get_orthogonal_matrix_from_QR(datatype, n, B, ldt, Q_temp, ldt, &info);
            free_matrix(Q_temp);
        }
        else
        {
            get_orthogonal_matrix_from_QR(datatype, n, B, ldt, Q, ldq, &info);
        }

        /* Create the orthogonal Z matrix */
        if(compz == 'V')
        {
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &tmp, n);
            rand_matrix(datatype, tmp, n, n, n);
            get_orthogonal_matrix_from_QR(datatype, n, tmp, n, Z, ldz, &info);
            free_matrix(tmp);
        }

        info = 0;
        /* Make copy of matrix A and B. This is required to validate the API functionality */
        copy_matrix(datatype, "full", n, n, A, ldh, H, ldh);
        copy_matrix(datatype, "full", n, n, B, ldt, T, ldt);
        if(compq == 'I')
            set_identity_matrix(datatype, n, n, Q, ldq);
        if(compz == 'I')
            set_identity_matrix(datatype, n, n, Z, ldz);
        /* Make copy of matrix Q and Z. This is required to validate the API functionality */
        copy_matrix(datatype, "full", n, n, Q, ldq, Q_A, ldq);
        copy_matrix(datatype, "full", n, n, Z, ldz, Z_A, ldz);
        /* Call to GGHRD API */
        invoke_gghrd(datatype, &compq, &compz, &n, &ilo, &ihi, H, &ldh, T, &ldt, Q, &ldq, Z, &ldz,
                     &info);
        if(compq == 'I')
            set_identity_matrix(datatype, n, n, Q, ldq);
        if(compz == 'I')
            set_identity_matrix(datatype, n, n, Z, ldz);
    }

    /* Make copy of matrix H,T,Q and Z. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &H_test, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &T_test, ldt);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_test, ldq);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_test, ldz);
    copy_matrix(datatype, "full", n, n, H, ldh, H_test, ldh);
    copy_matrix(datatype, "full", n, n, T, ldt, T_test, ldt);
    copy_matrix(datatype, "full", n, n, Q, ldq, Q_test, ldq);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_test, ldz);

    prepare_hgeqz_run(&job, &compq, &compz, n, &ilo, &ihi, H_test, ldh, T_test, ldt, alpha, alphar,
                      alphai, beta, Q_test, ldq, Z_test, ldz, datatype, n_repeats, &time_min, &info,
                      interfacetype, layout);

    /* if job=E, in addition to the first api, also execute the api with jobe=S, compeq=N and
     * compez=N */
    if(job == 'E')
    {
        char jobe = 'S';
        char compeq = 'N';
        char compez = 'N';
        double time_mine;

        if(datatype == FLOAT || datatype == DOUBLE)
        {
            create_vector(datatype, &alphaer, n);
            create_vector(datatype, &alphaei, n);
        }
        else
        {
            create_vector(datatype, &alphae, n);
        }
        create_vector(datatype, &betae, n);

        /* Call the API */
        prepare_hgeqz_run(&jobe, &compeq, &compez, n, &ilo, &ihi, H_test, ldh, T_test, ldt, alphae,
                          alphaer, alphaei, betae, Q_test, ldq, Z_test, ldz, datatype, n_repeats,
                          &time_mine, &info, interfacetype, layout);
    }
    /* If compq=N or/and compz=N, also execute the api with compnq=I and compnz=I */
    else if(compq == 'N' || compz == 'N')
    {
        char compnq = 'I';
        char compnz = 'I';
        double time_minn;

        integer ldqn = fla_max(ldq, n);
        integer ldzn = fla_max(ldz, n);

        /* Create a copy of H and T and send it to the API. H_ntest and T_ntest are used during
         * validation */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &H_ntest, ldh);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &T_ntest, ldt);
        copy_matrix(datatype, "full", n, n, H, ldh, H_ntest, ldh);
        copy_matrix(datatype, "full", n, n, T, ldt, T_ntest, ldt);

        if(datatype == FLOAT || datatype == DOUBLE)
        {
            create_vector(datatype, &alphanr, n);
            create_vector(datatype, &alphani, n);
        }
        else
        {
            create_vector(datatype, &alphan, n);
        }
        create_vector(datatype, &betan, n);

        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_ntest, ldqn);
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_ntest, ldzn);

        /* Set Q_ntest and Z_ntest matrices to identity matrices for compnq=I and compnz=I */
        set_identity_matrix(datatype, n, n, Q_ntest, ldqn);
        set_identity_matrix(datatype, n, n, Z_ntest, ldzn);

        /* Call the API */
        prepare_hgeqz_run(&job, &compnq, &compnz, n, &ilo, &ihi, H_ntest, ldh, T_ntest, ldt, alphan,
                          alphanr, alphani, betan, Q_ntest, ldqn, Z_ntest, ldzn, datatype,
                          n_repeats, &time_minn, &info, interfacetype, layout);
        /* Free the matrices */
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
        /* If job=E, validate eigen values from the first and second api calls */
        if(job == 'E')
        {
            validate_hgeqz_eigen_values(tst_api, datatype, n, alpha, alphar, alphai, beta, alphae,
                                        alphaer, alphaei, betae, residual);
        }
        /* If compq=N or compz=N, validate eigen values from the first and second api
           calls. Then validate H_test and T_test matrices */
        else if(compq == 'N' || compz == 'N')
        {
            validate_hgeqz_comp_n(tst_api, datatype, n, H_test, H_ntest, ldh, T_test, T_ntest, ldt,
                                  alpha, alphar, alphai, beta, alphan, alphanr, alphani, betan,
                                  residual);
        }
        else
        {
            validate_hgeqz(tst_api, &job, &compq, &compz, n, H, H_test, A, ldh, T, T_test, B, ldt,
                           Q, Q_test, Q_A, ldq, Z, Z_test, Z_A, ldz, datatype, residual,
                           params->imatrix_char);
        }
    }
    else
    {
        printf("Extreme Value tests not supported for xHGEQZ APIs\n");
    }

    /* Free up the buffers */
    free_vector(beta);
    free_matrix(H);
    free_matrix(T);
    free_matrix(Q);
    free_matrix(Z);
    free_matrix(H_test);
    free_matrix(T_test);
    free_matrix(Q_test);
    free_matrix(Z_test);
    free_matrix(A);
    free_matrix(B);
    free_matrix(Q_A);
    free_matrix(Z_A);
    if(datatype == FLOAT || datatype == DOUBLE)
    {
        free_vector(alphar);
        free_vector(alphai);
        if(job == 'E')
        {
            free_vector(alphaer);
            free_vector(alphaei);
        }

        if(compq == 'N' || compz == 'N')
        {
            free_vector(alphanr);
            free_vector(alphani);
        }
    }
    else
    {
        free_vector(alpha);
        if(job == 'E')
            free_vector(alphae);
        if(compq == 'N' || compz == 'N')
            free_vector(alphan);
    }
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal_A);
        free_vector(scal_B);
    }

    if(job == 'E')
        free_vector(betae);
    if(compq == 'N' || compz == 'N')
    {
        free_vector(betan);
        free_matrix(H_ntest);
        free_matrix(T_ntest);
    }
}

void prepare_hgeqz_run(char *job, char *compq, char *compz, integer n, integer *ilo, integer *ihi,
                       void *H, integer ldh, void *T, integer ldt, void *alpha, void *alphar,
                       void *alphai, void *beta, void *Q, integer ldq, void *Z, integer ldz,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int layout)
{
    void *H_save = NULL, *T_save = NULL, *Q_save = NULL, *Z_save = NULL, *work = NULL, *rwork;
    integer i, lwork;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix H,T,Q and Z. Same input values will be passed in each
     * itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &H_save, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &T_save, ldt);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_save, ldq);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_save, ldz);
    create_realtype_vector(datatype, &rwork, n);
    copy_matrix(datatype, "full", n, n, H, ldh, H_save, ldh);
    copy_matrix(datatype, "full", n, n, T, ldt, T_save, ldt);
    copy_matrix(datatype, "full", n, n, Q, ldq, Q_save, ldq);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_save, ldz);

    /* Make a workspace query the first time through. This will provide us with
     and ideal workspace size based on an internal block size.
     NOTE: LAPACKE interface handles workspace query internally */
    if((interfacetype != LAPACKE_COLUMN_TEST) && (interfacetype != LAPACKE_ROW_TEST)
       && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);

#if ENABLE_CPP_TEST
        /* call to CPP hgeqz API */
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_hgeqz(datatype, job, compq, compz, &n, ilo, ihi, NULL, &ldh, NULL, &ldt,
                             alpha, alphar, alphai, beta, NULL, &ldq, NULL, &ldz, work, &lwork,
                             rwork, info);
        }
        else
#endif
        {
            /* call to hgeqz API */
            invoke_hgeqz(datatype, job, compq, compz, &n, ilo, ihi, NULL, &ldh, NULL, &ldt, alpha,
                         alphar, alphai, beta, NULL, &ldq, NULL, &ldz, work, &lwork, rwork, info);
        }
        /* Output buffers will be freshly allocated for each iterations, free up
        the current output buffers.*/
        if(*info == 0)
        {
            /* Get work size */
            lwork = get_work_value(datatype, work);
            free_vector(work);
        }
    }
    else
    {
        lwork = g_lwork;
    }

    *info = 0;
    for(i = 0; i < n_repeats && *info == 0; ++i)
    {
        /* Restore input matrix H and Z value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, H_save, ldh, H, ldh);
        copy_matrix(datatype, "full", n, n, T_save, ldt, T, ldt);
        copy_matrix(datatype, "full", n, n, Q_save, ldq, Q, ldq);
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        create_vector(datatype, &work, lwork);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_hgeqz_run(datatype, layout, job, compq, compz, n, ilo, ihi,
                                                 H, ldh, T, ldt, alpha, alphar, alphai, beta, Q,
                                                 ldq, Z, ldz, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST) /* Call CPP hgeqz API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_hgeqz(datatype, job, compq, compz, &n, ilo, ihi, H, &ldh, T, &ldt, alpha,
                             alphar, alphai, beta, Q, &ldq, Z, &ldz, work, &lwork, rwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK hgeqz API */
            invoke_hgeqz(datatype, job, compq, compz, &n, ilo, ihi, H, &ldh, T, &ldt, alpha, alphar,
                         alphai, beta, Q, &ldq, Z, &ldz, work, &lwork, rwork, info);
            exe_time = fla_test_clock() - exe_time;
        }

        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        /* Free up the output buffers */
        free_vector(work);
    }
    *time_min_ = t_min;

    free(H_save);
    free(Z_save);
    free(Q_save);
    free(T_save);
    free_vector(rwork);
}

double prepare_lapacke_hgeqz_run(integer datatype, int layout, char *job, char *compq, char *compz,
                                 integer n, integer *ilo, integer *ihi, void *H, integer ldh,
                                 void *T, integer ldt, void *alpha, void *alphar, void *alphai,
                                 void *beta, void *Q, integer ldq, void *Z, integer ldz,
                                 integer *info)
{
    double exe_time;
    integer ldh_t = ldh;
    integer ldt_t = ldt;
    integer ldq_t = ldq;
    integer ldz_t = ldz;
    void *H_t = NULL, *T_t = NULL, *Q_t = NULL, *Z_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_hgeqz_ldh, ldh_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_hgeqz_ldt, ldt_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_hgeqz_ldq, ldq_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_hgeqz_ldz, ldz_t);

    H_t = H;
    T_t = T;
    Q_t = Q;
    Z_t = Z;
    /* In case of row_major matrix layout,
       convert input matrix to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &H_t, fla_max(n, ldh_t));
        create_matrix(datatype, layout, n, n, &T_t, fla_max(n, ldt_t));

        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, H, ldh, H_t, ldh_t);
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, T, ldt, T_t, ldt_t);
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

    /* Call to hgeqz API */
    *info = invoke_lapacke_hgeqz(datatype, layout, *job, *compq, *compz, n, *ilo, *ihi, H_t, ldh_t,
                                 T_t, ldt_t, alpha, alphar, alphai, beta, Q_t, ldq_t, Z_t, ldz_t);

    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, H_t, ldh_t, H, ldh);
        convert_matrix_layout(layout, datatype, n, n, T_t, ldt_t, T, ldt);

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
        free_matrix(H_t);
        free_matrix(T_t);
    }
    return exe_time;
}

void invoke_hgeqz(integer datatype, char *job, char *compq, char *compz, integer *n, integer *ilo,
                  integer *ihi, void *h, integer *ldh, void *t, integer *ldt, void *alpha,
                  void *alphar, void *alphai, void *beta, void *q, integer *ldq, void *z,
                  integer *ldz, void *work, integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_shgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta,
                              q, ldq, z, ldz, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta,
                              q, ldq, z, ldz, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_chgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq,
                              z, ldz, work, lwork, rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq,
                              z, ldz, work, lwork, rwork, info);
            break;
        }
    }
}
