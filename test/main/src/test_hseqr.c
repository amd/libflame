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
integer row_major_hseqr_ldh;
integer row_major_hseqr_ldz;

/* Local prototypes */
void fla_test_hseqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_hseqr_run(char *job, char *compz, integer n, integer *ilo, integer *ihi, void *h,
                       integer ldh, void *w, void *wr, void *wi, void *z, integer ldz,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int matrix_layout);
void invoke_hseqr(integer datatype, char *job, char *compz, integer *n, integer *ilo, integer *ihi,
                  void *h, integer *ldh, void *w, void *wr, void *wi, void *z, integer *ldz,
                  void *work, integer *lwork, integer *info);
double prepare_lapacke_hseqr_run(integer datatype, int matrix_layout, char *job, char *compz,
                                 integer n, integer *ilo, integer *ihi, void *h, integer ldh,
                                 void *w, void *wr, void *wi, void *z, integer ldz, integer *info);

void fla_test_hseqr(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Computing Eigen value of a Hessenberg matrix";
    char *front_str = "HSEQR";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    if(argc == 1)
    {
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_hseqr_experiment);
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

        /* Prase the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].job_seqr = argv[3][0];
        params->eig_sym_paramslist[0].compz_hseqr = argv[4][0];
        N = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ilo = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->eig_sym_paramslist[0].ihi = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_hseqr_ldh = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_hseqr_ldz = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->lin_solver_paramslist[0].lda = N;
            params->lin_solver_paramslist[0].ldz = N;
        }
        else
        {
            params->eig_sym_paramslist[0].lda = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldz = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        }
        g_lwork = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
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
                fla_test_hseqr_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for HSEQR\n");
        printf("./<EXE> hseqr <precisions - sdcz> <job> <compz> <N> <ILO> <IHI> <LDH> <LDZ> "
               "<LWORK> <repeats>\n");
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

void fla_test_hseqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, ldz, ldh;
    integer ilo, ihi, info = 0;
    void *H = NULL, *w = NULL, *wr = NULL, *wi = NULL, *Z = NULL, *H_test = NULL, *Z_Test = NULL,
         *wr_in = NULL, *wi_in = NULL, *scal_H = NULL;
    char job;
    char compz;
    double residual, err_thresh;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;

    /* Get input matrix dimensions. */
    n = p_cur;

    ldz = params->eig_sym_paramslist[pci].ldz;
    ldh = params->eig_sym_paramslist[pci].lda;

    /* Initialize parameter needed for HSEQR() call. */
    job = params->eig_sym_paramslist[pci].job_seqr;
    compz = params->eig_sym_paramslist[pci].compz_hseqr;
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
        /* if COMPZ = 'I' or COMPZ = 'V', then LDZ >= MAX(1,N)
           Otherwise, LDZ >= 1 */
        if(ldz == -1)
        {
            if((compz == 'I') || (compz == 'V'))
            {
                ldz = fla_max(1, n);
            }
            else
            {
                ldz = 1;
            }
        }
    }

    /* Create input matrix parameters*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &H, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z, ldz);

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        create_vector(datatype, &w, n);
    }
    else
    {
        create_vector(datatype, &wr, n);
        create_vector(datatype, &wi, n);
    }

    if(g_ext_fptr != NULL)
    {
        init_matrix_from_file(datatype, H, n, n, ldh, g_ext_fptr);
        init_matrix_from_file(datatype, Z, n, n, ldz, g_ext_fptr);
    }
    else
    {
        create_vector(datatype, &wr_in, n);
        if(datatype == FLOAT || datatype == DOUBLE)
        {
            create_vector(datatype, &wi_in, n);
            reset_vector(datatype, wi_in, n, 1);
        }

        get_hessenberg_matrix_from_EVs(datatype, n, H, ldh, Z, ldz, &ilo, &ihi, &info, wr_in,
                                       wi_in);
        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_realtype_vector(datatype, &scal_H, n);
            scale_matrix_overflow_underflow_hseqr(datatype, n, H, ldh, params->imatrix_char,
                                                  scal_H);
        }
        if(compz == 'I')
            set_identity_matrix(datatype, n, n, Z, ldz);
    }

    /* Make copy of matrix H and Z. This is required to validate the API functionality */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &H_test, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_Test, ldz);
    copy_matrix(datatype, "full", n, n, H, ldh, H_test, ldh);
    copy_matrix(datatype, "full", n, n, Z, ldz, Z_Test, ldz);

    prepare_hseqr_run(&job, &compz, n, &ilo, &ihi, H_test, ldh, w, wr, wi, Z_Test, ldz, datatype,
                      n_repeats, &time_min, &info, interfacetype, layout);

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
        validate_hseqr(tst_api, &job, &compz, n, H, H_test, ldh, Z, Z_Test, ldz, wr, wr_in, wi,
                       wi_in, w, datatype, residual, &ilo, &ihi, params->imatrix_char, scal_H);
    }
    else
    {
        printf("Extreme Value tests not supported for xHSEQR APIs\n");
    }

    /* Free up the buffers */
    free_matrix(H);
    free_matrix(Z);
    free_matrix(H_test);
    free_matrix(Z_Test);
    free_vector(wr_in);
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
    {
        free_vector(w);
    }
    else
    {
        free_vector(wr);
        free_vector(wi);
        free_vector(wi_in);
    }
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        free_vector(scal_H);
    }
}

void prepare_hseqr_run(char *job, char *compz, integer n, integer *ilo, integer *ihi, void *H,
                       integer ldh, void *w, void *wr, void *wi, void *Z, integer ldz,
                       integer datatype, integer n_repeats, double *time_min_, integer *info,
                       integer interfacetype, int layout)
{
    void *H_save = NULL, *work = NULL, *Z_save = NULL;
    integer i, lwork;
    double t_min = 1e9, exe_time;

    /* Make a copy of the input matrix H and Z. Same input values will be passed in each
     * itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &H_save, ldh);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Z_save, ldz);
    copy_matrix(datatype, "full", n, n, H, ldh, H_save, ldh);
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
        /* call to CPP hseqr API */
        if(interfacetype == LAPACK_CPP_TEST)
        {
            invoke_cpp_hseqr(datatype, job, compz, &n, ilo, ihi, NULL, &ldh, NULL, NULL, NULL, NULL,
                             &ldz, work, &lwork, info);
        }
        else
#endif
        {
            /* call to hseqr API */
            invoke_hseqr(datatype, job, compz, &n, ilo, ihi, NULL, &ldh, NULL, NULL, NULL, NULL,
                         &ldz, work, &lwork, info);
        }

        /* Output buffers will be freshly allocated for each iterations, free up
        the current output buffers.*/
        if(*info == 0)
        {
            /* Get work size */
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
        /* Restore input matrix H and Z value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", n, n, H_save, ldh, H, ldh);
        copy_matrix(datatype, "full", n, n, Z_save, ldz, Z, ldz);
        create_vector(datatype, &work, lwork);

        /* Check if LAPACKE interface is enabled */
        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_hseqr_run(datatype, layout, job, compz, n, ilo, ihi, H, ldh,
                                                 w, wr, wi, Z, ldz, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST) /* Call CPP hseqr API */
        {
            exe_time = fla_test_clock();
            invoke_cpp_hseqr(datatype, job, compz, &n, ilo, ihi, H, &ldh, w, wr, wi, Z, &ldz, work,
                             &lwork, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            /* Call LAPACK hseqr API */
            invoke_hseqr(datatype, job, compz, &n, ilo, ihi, H, &ldh, w, wr, wi, Z, &ldz, work,
                         &lwork, info);

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
}

double prepare_lapacke_hseqr_run(integer datatype, int layout, char *job, char *compz, integer n,
                                 integer *ilo, integer *ihi, void *H, integer ldh, void *w,
                                 void *wr, void *wi, void *Z, integer ldz, integer *info)
{
    double exe_time;
    integer ldh_t = ldh;
    integer ldz_t = ldz;
    void *H_t = NULL, *Z_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_hseqr_ldh, ldh_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_hseqr_ldz, ldz_t);

    H_t = H;
    Z_t = Z;
    /* In case of row_major matrix layout,
       convert input matrices to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &H_t, fla_max(n, ldh_t));

        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, H, ldh, H_t, ldh_t);
        if(*compz != 'N')
        {
            create_matrix(datatype, layout, n, n, &Z_t, fla_max(n, ldz_t));
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, Z, ldz, Z_t, ldz_t);
        }
    }
    exe_time = fla_test_clock();

    /* Call LAPACKE hseqr API */
    *info = invoke_lapacke_hseqr(datatype, layout, *job, *compz, n, *ilo, *ihi, H_t, ldh_t, w, wr,
                                 wi, Z_t, ldz_t);
    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, H_t, ldh_t, H, ldh);
        if(*compz != 'N')
        {
            convert_matrix_layout(layout, datatype, n, n, Z_t, ldz_t, Z, ldz);
            free_matrix(Z_t);
        }
        /* free temporary buffers */
        free_matrix(H_t);
    }
    return exe_time;
}

void invoke_hseqr(integer datatype, char *job, char *compz, integer *n, integer *ilo, integer *ihi,
                  void *h, integer *ldh, void *w, void *wr, void *wi, void *z, integer *ldz,
                  void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_shseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dhseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_chseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zhseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);
            break;
        }
    }
}
