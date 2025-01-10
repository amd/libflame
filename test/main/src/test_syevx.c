/*
    Copyright (C) 2023-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_common.h"
#include "test_lapack.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;
integer row_major_syevx_lda;
integer row_major_syevx_ldz;

/* Local prototypes.*/
void fla_test_syevx_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_syevx_run(char *jobz, char *range, char *uplo, integer n, void *A, integer lda,
                       void *vl, void *vu, integer il, integer iu, void *abstol, void *w,
                       integer ldz, void *ifail, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer test_lapacke_interface,
                       int matrix_layout);
void invoke_syevx(integer datatype, char *jobz, char *range, char *uplo, integer *n, void *a,
                  integer *lda, void *vl, void *vu, integer *il, integer *iu, void *abstol,
                  integer *m, void *w, void *z, integer *ldz, void *work, integer *lwork,
                  void *rwork, void *iwork, void *ifail, integer *info);
double prepare_lapacke_syevx_run(integer datatype, int matrix_layout, char *jobz, char *range,
                                 char *uplo, integer n, void *A, integer lda, void *vl, void *vu,
                                 integer il, integer iu, void *abstol, integer *m, void *w, void *z,
                                 integer ldz, void *ifail, integer *info);
integer invoke_lapacke_syevx(integer datatype, int layout, char jobz, char range, char uplo,
                             integer n, void *a, integer lda, void *vl, void *vu, integer il,
                             integer iu, void *abstol, integer *m, void *w, void *z, integer ldz,
                             integer *ifail);

void fla_test_syevx(integer argc, char **argv, test_params_t *params)
{
    srand(1); /* Setting the seed for random input generation values */
    char *op_str = "Eigen Values and Vectors in specified range";
    char *front_str = "SYEVX";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_lwork = -1;
        config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_SYM, fla_test_syevx_experiment);
        tests_not_run = 0;
    }
    if(argc == 17)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[16]);
    }
    if(argc >= 16 && argc <= 17)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->eig_sym_paramslist[0].jobz = argv[3][0];
        params->eig_sym_paramslist[0].range_x = argv[4][0];
        params->eig_sym_paramslist[0].uplo = argv[5][0];
        N = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && params->test_lapacke_interface
           && (params->matrix_major == LAPACK_ROW_MAJOR))
        {
            row_major_syevx_lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            row_major_syevx_ldz = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldz = N;
            params->eig_sym_paramslist[0].lda = N;
        }
        else
        {
            params->eig_sym_paramslist[0].lda = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
            params->eig_sym_paramslist[0].ldz = strtoimax(argv[13], &endptr, CLI_DECIMAL_BASE);
        }
        params->eig_sym_paramslist[0].VL = atof(argv[8]);
        params->eig_sym_paramslist[0].VU = atof(argv[9]);

        if(params->eig_sym_paramslist[0].range_x == 'I')
        {
            /* 1 <= IL <= IU <= N, if N > 0;
               IL = 1 and IU = 0 if N = 0. */
            if(N == 0)
            {
                params->eig_sym_paramslist[0].IL = 1;
                params->eig_sym_paramslist[0].IU = 0;
            }
            else
            {
                params->eig_sym_paramslist[0].IL = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
                params->eig_sym_paramslist[0].IU = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
            }
        }

        params->eig_sym_paramslist[0].abstol = atof(argv[12]);

        g_lwork = strtoimax(argv[14], &endptr, CLI_DECIMAL_BASE);

        n_repeats = strtoimax(argv[15], &endptr, CLI_DECIMAL_BASE);

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
                fla_test_syevx_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for syevx/heevx\n");
        printf("./<EXE> syevx <precisions - sd> <JOBZ> <RANGE> <UPLO>"
               " <N> <LDA> <VL> <VU> <IL> <IU> <ABSTOL> <LDZ> <LWORK> <repeats>\n");
        printf("./<EXE> heevx <precisions - cz> <JOBZ> <RANGE> <UPLO>"
               " <N> <LDA> <VL> <VU> <IL> <IU> <ABSTOL> <LDZ> <LWORK> <repeats>\n");
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

void fla_test_syevx_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, lda, ldz, il, iu, info = 0;
    char jobz, uplo, range;
    void *A = NULL, *w = NULL, *A_test = NULL, *L = NULL, *ifail = NULL, *scal = NULL;
    void *vl, *vu, *abstol;
    double residual, err_thresh;

    integer test_lapacke_interface = params->test_lapacke_interface;
    int layout = params->matrix_major;

    /* Get input matrix dimensions.*/
    jobz = params->eig_sym_paramslist[pci].jobz;
    uplo = params->eig_sym_paramslist[pci].uplo;
    range = params->eig_sym_paramslist[pci].range_x;
    err_thresh = params->eig_sym_paramslist[pci].threshold_value;

    n = p_cur;
    lda = params->eig_sym_paramslist[pci].lda;
    ldz = params->eig_sym_paramslist[pci].ldz;

    il = params->eig_sym_paramslist[pci].IL;
    iu = params->eig_sym_paramslist[pci].IU;

    create_realtype_vector(datatype, &vl, 1);
    create_realtype_vector(datatype, &vu, 1);
    create_realtype_vector(datatype, &abstol, 1);

    if(datatype == FLOAT || datatype == COMPLEX)
    {
        *(float *)vl = params->eig_sym_paramslist[pci].VL;
        *(float *)vu = params->eig_sym_paramslist[pci].VU;
        *(float *)abstol = params->eig_sym_paramslist[pci].abstol;

        /* When abstol value is set to -1, assign default value.
           NOTE: Eigenvalues will be computed most accurately
                 when ABSTOL is set to twice the underflow
                 threshold 2*SLAMCH('S') */
        if(*(float *)abstol == -1)
            *(float *)abstol = 2 * slamch_("S");
    }
    else
    {
        *(double *)vl = params->eig_sym_paramslist[pci].VL;
        *(double *)vu = params->eig_sym_paramslist[pci].VU;
        *(double *)abstol = params->eig_sym_paramslist[pci].abstol;

        /* When abstol value is set to -1, assign default value.
           NOTE: Eigenvalues will be computed most accurately
                 when ABSTOL is set to twice the underflow
                 threshold 2*DLAMCH('S') */
        if(*(double *)abstol == -1)
            *(double *)abstol = 2 * dlamch_("S");
    }

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
        /* LDZ >= 1;
           if JOBZ = 'V', LDZ >= max(1,N) */
        if(ldz == -1)
        {
            if(jobz == 'V')
            {
                ldz = fla_max(1, n);
            }
            else
            {
                ldz = 1;
            }
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_realtype_vector(datatype, &w, n);
    if((!FLA_EXTREME_CASE_TEST) && (g_ext_fptr == NULL))
    {
        /*  Creating input matrix A by generating random eigen values.
            When range = V, generate EVs in given range (vl,vu)  */
        create_realtype_vector(datatype, &L, n);
        generate_matrix_from_EVs(datatype, range, n, A, lda, L, get_realtype_value(datatype, vl),
                                 get_realtype_value(datatype, vu), USE_ABS_EIGEN_VALUES);

        if(FLA_OVERFLOW_UNDERFLOW_TEST)
        {
            create_vector(get_realtype(datatype), &scal, 1);
            scale_matrix_underflow_overflow_syevx(datatype, n, A, lda, params->imatrix_char, scal);
        }
    }
    else
    {
        /* Initialize input matrix with custom data */
        init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
    }
    /* Make a copy of input matrix A.
       This is required to validate the API functionality.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_test, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_test, lda);
    create_vector(INTEGER, &ifail, n);

    prepare_syevx_run(&jobz, &range, &uplo, n, A_test, lda, vl, vu, il, iu, abstol, w, ldz, ifail,
                      datatype, n_repeats, &time_min, &info, test_lapacke_interface, layout);

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
        validate_syev(tst_api, &jobz, &range, n, A, A_test, lda, il, iu, L, w, ifail, datatype,
                      residual, params->imatrix_char, scal);
    }
    /* check for output matrix when inputs as extreme values */
    else if(FLA_EXTREME_CASE_TEST)
    {
        if((!check_extreme_value(datatype, n, n, A_test, lda, params->imatrix_char))
           && (!check_extreme_value(datatype, n, i_one, w, i_one, params->imatrix_char)))
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
    free_vector(vl);
    free_vector(vu);
    free_vector(abstol);
    free_vector(ifail);
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

void prepare_syevx_run(char *jobz, char *range, char *uplo, integer n, void *A, integer lda,
                       void *vl, void *vu, integer il, integer iu, void *abstol, void *w,
                       integer ldz, void *ifail, integer datatype, integer n_repeats,
                       double *time_min_, integer *info, integer test_lapacke_interface, int layout)
{
    void *A_save = NULL, *work = NULL, *rwork = NULL;
    void *w_test = NULL, *z__ = NULL;
    integer i, m, lwork;
    double t_min = 1e9, exe_time;
    void *iwork = NULL;

    if(*range == 'I')
        m = iu - il + 1;
    else
        m = n;

    /* Make a copy of the input matrix A.
       Same input values will be passed in eaach itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);
    create_vector(INTEGER, &iwork, 5 * n);

    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        create_realtype_vector(datatype, &rwork, (7 * n));
    else
        rwork = NULL;

    /* Make a workspace query the first time through. This will provide us with
       and ideal workspace size based on an internal block size.
       NOTE: LAPACKE interface handles workspace query internally */
    if((test_lapacke_interface == 0) && (g_lwork <= 0))
    {
        lwork = -1;
        create_vector(datatype, &work, 1);
        /* call to  syevx API */
        invoke_syevx(datatype, jobz, range, uplo, &n, NULL, &lda, vl, vu, &il, &iu, abstol, &m,
                     NULL, NULL, &ldz, work, &lwork, rwork, iwork, ifail, info);
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
        create_matrix(datatype, LAPACK_COL_MAJOR, fla_max(1, n), fla_max(1, m), &z__, ldz);

        /* Check if LAPACKE interface is enabled */
        if(test_lapacke_interface == 1)
        {
            exe_time
                = prepare_lapacke_syevx_run(datatype, layout, jobz, range, uplo, n, A, lda, vl, vu,
                                            il, iu, abstol, &m, w_test, z__, ldz, ifail, info);
        }
        else
        {
            exe_time = fla_test_clock();

            /* call to API */
            invoke_syevx(datatype, jobz, range, uplo, &n, A, &lda, vl, vu, &il, &iu, abstol, &m,
                         w_test, z__, &ldz, work, &lwork, rwork, iwork, ifail, info);

            exe_time = fla_test_clock() - exe_time;
        }
        /* Get the best execution time */
        t_min = fla_min(t_min, exe_time);

        /* Make a copy of the output buffers.
           This is required to validate the API functionality.*/
        copy_realtype_vector(datatype, n, w_test, 1, w, 1);

        /* If JOBZ = 'V', the first M columns of Z contain the
           orthonormal eigenvectors of the matrix A corresponding to
           the selected eigenvalues.
           Copy eigen vectors to A to validate API functionality */
        if(*jobz == 'V')
            copy_matrix(datatype, "full", m, m, z__, ldz, A, lda);

        /* Free up the output buffers */
        free_vector(work);
        free_vector(w_test);
        free_matrix(z__);
    }

    *time_min_ = t_min;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        free_vector(rwork);
    free_vector(iwork);
    free_matrix(A_save);
}

double prepare_lapacke_syevx_run(integer datatype, int layout, char *jobz, char *range, char *uplo,
                                 integer n, void *A, integer lda, void *vl, void *vu, integer il,
                                 integer iu, void *abstol, integer *m, void *w, void *z,
                                 integer ldz, void *ifail, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    integer ldz_t = ldz;
    void *A_t = NULL, *Z_t = NULL;

    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_syevx_lda, lda_t);
    SELECT_LDA(g_ext_fptr, config_data, layout, n, row_major_syevx_ldz, ldz_t);

    A_t = A;
    Z_t = z;

    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        create_matrix(datatype, layout, n, n, &A_t, lda_t);
        if(*jobz != 'N')
        {
            create_matrix(datatype, layout, *m, n, &Z_t, ldz_t);
        }
        convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
    }
    exe_time = fla_test_clock();

    /* call LAPACKE syev API */
    *info = invoke_lapacke_syevx(datatype, layout, *jobz, *range, *uplo, n, A_t, lda_t, vl, vu, il,
                                 iu, abstol, m, w, Z_t, ldz_t, ifail);

    exe_time = fla_test_clock() - exe_time;
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* In case of row_major matrix layout, convert output matrices
           to column_major layout */
        convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        if(*jobz != 'N')
        {
            convert_matrix_layout(layout, datatype, *m, n, Z_t, ldz_t, z, ldz);
            free_matrix(Z_t);
        }
        /* free temporary buffers */
        free_matrix(A_t);
    }

    return exe_time;
}

void invoke_syevx(integer datatype, char *jobz, char *range, char *uplo, integer *n, void *a,
                  integer *lda, void *vl, void *vu, integer *il, integer *iu, void *abstol,
                  integer *m, void *w, void *z, integer *ldz, void *work, integer *lwork,
                  void *rwork, void *iwork, void *ifail, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                              work, lwork, iwork, ifail, info);
            break;
        }
        case DOUBLE:
        {
            fla_lapack_dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                              work, lwork, iwork, ifail, info);
            break;
        }
        case COMPLEX:
        {
            fla_lapack_cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                              work, lwork, rwork, iwork, ifail, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            fla_lapack_zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                              work, lwork, rwork, iwork, ifail, info);
            break;
        }
    }
}

integer invoke_lapacke_syevx(integer datatype, int layout, char jobz, char range, char uplo,
                             integer n, void *a, integer lda, void *vl, void *vu, integer il,
                             integer iu, void *abstol, integer *m, void *w, void *z, integer ldz,
                             integer *ifail)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssyevx(layout, jobz, range, uplo, n, a, lda, *(float *)vl, *(float *)vu,
                                  il, iu, *(float *)abstol, m, w, z, ldz, ifail);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsyevx(layout, jobz, range, uplo, n, a, lda, *(double *)vl,
                                  *(double *)vu, il, iu, *(double *)abstol, m, w, z, ldz, ifail);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cheevx(layout, jobz, range, uplo, n, a, lda, *(float *)vl, *(float *)vu,
                                  il, iu, *(float *)abstol, m, w, z, ldz, ifail);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zheevx(layout, jobz, range, uplo, n, a, lda, *(double *)vl,
                                  *(double *)vu, il, iu, *(double *)abstol, m, w, z, ldz, ifail);
            break;
        }
    }
    return info;
}
