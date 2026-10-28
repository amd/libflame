/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_gebal_lda;

/* Local prototypes */
void fla_test_gebal_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_gebal_run(char *job, integer n, void *A, integer lda,
                       integer *ilo, integer *ihi, void *scale, integer datatype, integer *info,
                       integer interfacetype, integer layout, test_params_t *params);
void invoke_gebal(integer datatype, char *job, integer *n, void *A, integer *lda, integer *ilo,
                  integer *ihi, void *scale, integer *info);
double prepare_lapacke_gebal_run(integer datatype, integer layout, char job, integer n, void *A,
                                 integer lda, integer *ilo, integer *ihi,
                                 void *scale, integer *info);
static void gebal_isolate_interior(integer datatype, integer n, void *A, integer lda);
/* Build a matrix with an isolatable interior row and, when possible, an
   isolatable interior column. Diagonal entries stay in place. Those positions
   are not already at the ends GEBAL swaps them to, so JOB P and JOB B must
   exercise the permutation path. */
static void gebal_isolate_interior(integer datatype, integer n, void *A, integer lda)
{
    integer i, row, col;

    if(n < 3 || lda < n || A == NULL)
        return;

    row = n / 2;
    for(i = 0; i < n; i++)
    {
        if(i != row)
            assign_value(datatype, get_m_ptr(datatype, A, row, i, lda), d_zero, d_zero);
    }

    col = (row > 1) ? row - 1 : row + 1;
    if(col <= 0 || col >= n - 1 || col == row)
        return;

    for(i = 0; i < n; i++)
    {
        if(i != col)
            assign_value(datatype, get_m_ptr(datatype, A, i, col, lda), d_zero, d_zero);
    }
}

/* Parse CLI/config parameters and dispatch one experiment per requested precision. */
void fla_test_gebal(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Balance a general matrix";
    char *front_str = "GEBAL";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;

    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, EIG_NSYM, fla_test_gebal_experiment);
        tests_not_run = 0;
    }
    if(argc == 8)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[7]);
    }
    if(argc >= 7 && argc <= 8)
    {
        integer i, num_types, N;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        num_types = strlen(argv[2]);
        params->eig_non_sym_paramslist[0].balance_ggevx = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_gebal_lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
            params->eig_non_sym_paramslist[0].lda = N;
        }
        else
        {
            params->eig_non_sym_paramslist[0].lda = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        }
        n_repeats = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        if(n_repeats > 0)
        {
            params->eig_non_sym_paramslist[0].GenNonSymEigProblem_threshold = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                if(datatype == INVALID_TYPE)
                {
                    invalid_dtype = 1;
                    continue;
                }

                if(type_flag[datatype - FLOAT] == 1)
                    continue;
                type_flag[datatype - FLOAT] = 1;

                fla_test_gebal_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    if(tests_not_run)
    {
        printf("\nIllegal arguments for GEBAL\n");
        printf("./<EXE> gebal <precisions - sdcz> <JOB - N|P|S|B> <N> <LDA> <repeats>\n");
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

/* Allocate test buffers, run GEBAL, and validate the requested JOB mode. */
void fla_test_gebal_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer n, lda;
    integer ilo = 1, ihi = 1, info = 0;
    char job;
    void *A = NULL, *A_Test = NULL, *scale = NULL, *scal = NULL;
    double residual, err_thresh = 0.0;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;
    integer layout = params->matrix_major;

    n = p_cur;
    lda = params->eig_non_sym_paramslist[pci].lda;
    err_thresh = params->eig_non_sym_paramslist[pci].GenNonSymEigProblem_threshold;
    job = params->eig_non_sym_paramslist[pci].balance_ggevx;

    if(g_config_data)
    {
        if(lda == -1)
        {
            lda = fla_max(1, n);
        }
    }

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_Test, lda);
    create_realtype_vector(datatype, &scale, n);
    reset_vector(get_realtype(datatype), scale, n, 1);

    if(!FLA_BRT_VERIFICATION_RUN)
    {
        if(g_ext_fptr != NULL)
        {
            init_matrix_from_file(datatype, A, n, n, lda, g_ext_fptr);
        }
        else if(FLA_EXTREME_CASE_TEST)
        {
            init_matrix(datatype, A, n, n, lda, g_ext_fptr, params->imatrix_char);
        }
        else
        {
            if(FLA_OVERFLOW_UNDERFLOW_TEST)
            {
                rand_matrix(datatype, A, n, n, lda);
                create_vector(get_realtype(datatype), &scal, 1);
                if(n > 0 && lda >= n)
                    init_matrix_overflow_underflow_asym(datatype, n, n, A, lda,
                                                        params->imatrix_char, scal);
            }
            else
            {
                rand_matrix(datatype, A, n, n, lda);
                gebal_isolate_interior(datatype, n, A, lda);
            }
        }
    }

    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, n, n, A, lda, "cdd", job, n, lda)

    copy_matrix(datatype, "full", n, n, A, lda, A_Test, lda);

    prepare_gebal_run(&job, n, A_Test, lda, &ilo, &ihi, scale, datatype, &info,
                      interfacetype, layout, params);

    /* Each entry of A can be read and scaled. */
    perf = (double)n * (double)n / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 2.0;

    FLA_TEST_CHECK_EINFO(residual, info, einfo);
    IF_FLA_BRT_VALIDATION(n, n,
                          store_outputs_base(filename, params, 1, 3, datatype, n, n, A_Test,
                                             lda,
                                             get_realtype(datatype), n, scale, INTEGER, 1, &ilo,
                                             INTEGER, 1, &ihi),
                          validate_gebal(tst_api, job, n, ilo, ihi, A, A_Test, lda, scale,
                                         datatype, residual,
                                         params),
                          check_reproducibility_base(filename, params, 1, 3, datatype, n, n,
                                                     A_Test, lda, get_realtype(datatype),
                                                     n, scale, INTEGER, 1, &ilo, INTEGER, 1, &ihi))
    else if(FLA_SKIP_VALIDATION_MODE)
    {
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_gebal(tst_api, job, n, ilo, ihi, A, A_Test, lda, scale, datatype,
                       residual, params);
    }
    else
    {
        if((!check_extreme_value(datatype, n, n, A_Test, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }

free_buffers:
    FLA_FREE_FILENAME(filename);
    free_matrix(A);
    free_matrix(A_Test);
    free_vector(scale);
    free_vector(scal);
}

/* Run the selected interface in the timing loop and keep the best output. */
void prepare_gebal_run(char *job, integer n, void *A, integer lda,
                       integer *ilo, integer *ihi, void *scale, integer datatype, integer *info,
                       integer interfacetype, integer layout, test_params_t *params)
{
    void *A_save = NULL, *scale_test = NULL;
    double exe_time;

    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_save, lda);
    copy_matrix(datatype, "full", n, n, A, lda, A_save, lda);

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        copy_matrix(datatype, "full", n, n, A_save, lda, A, lda);
        create_realtype_vector(datatype, &scale_test, n);
        reset_vector(get_realtype(datatype), scale_test, n, 1);
        *ilo = 0;
        *ihi = 0;

        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_gebal_run(datatype, layout, *job, n, A, lda,
                                                 ilo, ihi, scale_test, info);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            invoke_cpp_gebal(datatype, job, &n, A, &lda, ilo, ihi, scale_test, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            invoke_gebal(datatype, job, &n, A, &lda, ilo, ihi, scale_test, info);
            exe_time = fla_test_clock() - exe_time;
        }

        FLA_EXEC_LOOP_UPDATE_WITH_INFO

        copy_vector(get_realtype(datatype), n, scale_test, 1, scale, 1);
        free_vector(scale_test);
    }

    free_matrix(A_save);
}

/* Stage row-major LAPACKE inputs when possible, but pass invalid leading
   dimensions directly so LAPACKE performs its own argument validation. */
double prepare_lapacke_gebal_run(integer datatype, integer layout, char job, integer n, void *A,
                                 integer lda, integer *ilo, integer *ihi,
                                 void *scale, integer *info)
{
    double exe_time;
    integer lda_t = lda;
    void *A_t = NULL;

    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_gebal_lda, lda_t);

    A_t = A;

    if(layout == LAPACK_ROW_MAJOR)
    {
        create_matrix(datatype, layout, n, n, &A_t, lda_t);
        /* LAPACKE checks row-major LDA before reading A. Do not stage an
           invalid layout, since its conversion helper writes through LDA. */
        if(lda_t >= fla_max(1, n))
        {
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, n, A, lda, A_t, lda_t);
        }
    }

    exe_time = fla_test_clock();
    *info = invoke_lapacke_gebal(datatype, layout, job, n, A_t, lda_t, ilo, ihi, scale);
    exe_time = fla_test_clock() - exe_time;

    if(layout == LAPACK_ROW_MAJOR)
    {
        if(lda_t >= fla_max(1, n))
        {
            convert_matrix_layout(layout, datatype, n, n, A_t, lda_t, A, lda);
        }
        free_matrix(A_t);
    }
    return exe_time;
}

/* Dispatch to the Fortran-compatible GEBAL symbol for the active datatype. */
void invoke_gebal(integer datatype, char *job, integer *n, void *A, integer *lda, integer *ilo,
                  integer *ihi, void *scale, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_sgebal(job, n, A, lda, ilo, ihi, scale, info);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dgebal(job, n, A, lda, ilo, ihi, scale, info);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_cgebal(job, n, A, lda, ilo, ihi, scale, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zgebal(job, n, A, lda, ilo, ihi, scale, info);
            break;
        }
    }
}
