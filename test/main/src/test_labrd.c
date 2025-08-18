/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif

extern double perf;
extern double time_min;

/* Local prototypes */
void fla_test_labrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_labrd_run(integer m_A, integer n_A, integer nb_A, void *A, integer lda, void *d,
                       void *e, void *tauq, void *taup, void *X, integer ldx, void *Y, integer ldy,
                       integer datatype, integer interfacetype, test_params_t *params);
void invoke_labrd(integer datatype, integer *m, integer *n, integer *nb, void *a, integer *lda,
                  void *d, void *e, void *tauq, void *taup, void *x, integer *ldx, void *y,
                  integer *ldy);

/* Helper functions for Bit reproducibility tests */
void store_labrd_outputs(void *filename, integer datatype, integer m, integer n, integer nb,
                         void *A, integer lda, void *d, void *e, void *tauq, void *taup, void *X,
                         integer ldx, void *Y, integer ldy, void *params);
integer check_bit_reproducibility_labrd(void *filename, integer datatype, integer m, integer n,
                                        integer nb, void *A, integer lda, void *d, void *e,
                                        void *tauq, void *taup, void *X, integer ldx, void *Y,
                                        integer ldy, void *params);

void fla_test_labrd(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Singular value decomposition";
    char *front_str = "LABRD";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, RECT_INPUT, params, AUX, fla_test_labrd_experiment);
        tests_not_run = 0;
    }
    if(argc == 11)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[10]);
    }
    if(argc >= 10 && argc <= 11)
    {
        integer i, num_types, N, M;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        M = strtoimax(argv[3], &endptr, CLI_DECIMAL_BASE);
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        params->aux_paramslist[0].nb = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].lda = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].ldx = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->aux_paramslist[0].ldy = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
        n_repeats = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        if(n_repeats > 0)
        {
            params->aux_paramslist[0].aux_threshold = CLI_NORM_THRESH;

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
                fla_test_labrd_experiment(front_str, params, datatype, M, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run)
    {
        printf("\nIllegal arguments for labrd\n");
        printf("./<EXE> labrd <precisions - sdcz> <M> <N> <NB> <LDA> <LDX> <LDY> <repeats>\n");
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

void fla_test_labrd_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    integer m, n, nb, lda, ldx, ldy;
    void *A = NULL, *X = NULL, *Y = NULL, *d = NULL, *e = NULL, *tauq = NULL, *taup = NULL;
    void *A_test = NULL;
    double residual, err_thresh;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;

    /* Get input matrix dimensions. */
    err_thresh = params->aux_paramslist[pci].aux_threshold;

    m = p_cur;
    n = q_cur;
    nb = params->aux_paramslist[pci].nb;
    /* Adjusting nb */
    if(nb > m || nb > n)
    {
        nb = fla_min(m, n);
    }
    lda = params->aux_paramslist[pci].lda;
    ldx = params->aux_paramslist[pci].ldx;
    ldy = params->aux_paramslist[pci].ldy;

    /* If leading dimensions = -1, set them to default value
       when inputs are from config files */
    if(g_config_data)
    {
        /* LDA >= max(1,M) */
        if(lda == -1)
        {
            lda = fla_max(1, m);
        }
        /* LDX >= max(1,M) */
        if(ldx == -1)
        {
            ldx = fla_max(1, m);
        }
        /* LDY >= max(1,N) */
        if(ldy == -1)
        {
            ldy = fla_max(1, n);
        }
    }

    /* Create input matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_test, lda);

    /* Create output matrix parameters */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, nb, &X, ldx);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, nb, &Y, ldy);

    /* Create output array parameters */
    create_realtype_vector(datatype, &d, nb);
    create_realtype_vector(datatype, &e, nb);

    create_vector(datatype, &tauq, nb);
    create_vector(datatype, &taup, nb);

    /* This code path is run to generate the matrix to be passed to the API. This is the default
     * input generation logic accessed both when BRT is run in Ground truth mode and for non BRT
     * Test cases. For verification runs the input is loaded from the input generated during Ground
     * truth run */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        /* initialize input matrix */
        init_matrix(datatype, A, m, n, lda, g_ext_fptr, params->imatrix_char);
    }

    /* This macro is used in the BRT test cases for the following purposes:
     *    - In the Ground truth runs (BRT_char => G, F), the output is stored in a file for future
     * reference
     *    - In the verification runs (BRT_char => V, M), the output is loaded from the file and
     * passed as input to the API
     * */
    FLA_BRT_PROCESS_SINGLE_INPUT(datatype, m, n, A, lda, "dddddd", m, n, nb, lda, ldx, ldy)

    /* Scaling matrix with values around overflow, underflow for LABRD */
    if(FLA_OVERFLOW_UNDERFLOW_TEST)
    {
        scale_matrix_underflow_overflow_labrd(datatype, m, n, A, lda, params->imatrix_char);
    }

    /* Make a copy of input matrix A. This is required to validate the API functionality. */
    copy_matrix(datatype, "full", m, n, A, lda, A_test, lda);

    prepare_labrd_run(m, n, nb, A_test, lda, d, e, tauq, taup, X, ldx, Y, ldy, datatype,
                      interfacetype, params);

    /* Performance Computation
     * The number of floating point operations in GEBRD is 4n^2(3m - n)/3 if m>=n else 4m^2(3n-m)/3
     * LABRD is a step in GEBRD where we consider only the leading rows and columns for reducing
     * into a bidiagonal matrix.
     * Equation specific to LABRD was not found instead we modify the GEBRD equation.
     * Sum of floating point operations in GEBRD for (m,nb) + (n,nb) - (nb,nb)
     * Link : https://support.nag.com/numeric/nl/nagdoc_latest/clhtml/f08/f08kec.html */

    perf = (double)(((4.0 * nb * nb) * ((3.0 * m) + (3.0 * n) - (4.0 * nb))) / 3.0) / time_min
           / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* Bit reproducibility tests path
     * This path is taken when BRT is enabled.
     *     - In the Ground truth runs (BRT_char => G, F), the output is stored in a file and the
     * default validation function is called
     *     - In the verification runs (BRT_char => V, M), the output is loaded from the file and
     * compared with the generated output
     *  */
    FLA_BRT_OUTPUT_VALIDATION(
        m, n,
        store_labrd_outputs(filename, datatype, m, n, nb, A_test, lda, d, e, tauq, taup, X, ldx, Y,
                            ldy, params),
        validate_labrd(tst_api, m, n, nb, A, A_test, lda, d, e, tauq, taup, X, ldx, Y, ldy,
                       datatype, err_thresh, g_ext_fptr, params->imatrix_char, params),
        check_bit_reproducibility_labrd(filename, datatype, m, n, nb, A_test, lda, d, e, tauq, taup,
                                        X, ldx, Y, ldy, params))
    /* API functionality validation */
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_labrd(tst_api, m, n, nb, A, A_test, lda, d, e, tauq, taup, X, ldx, Y, ldy,
                       datatype, err_thresh, g_ext_fptr, params->imatrix_char, params);
    }
    /* check for output matrix when inputs as extreme values */
    else
    {
        if((!check_extreme_value(datatype, m, n, A_test, lda, params->imatrix_char)))
        {
            residual = DBL_MAX;
        }
        else
        {
            residual = err_thresh;
        }
        FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    }

    /* Free up buffers */
free_buffers:
    FLA_FREE_FILENAME(filename)
    free_matrix(A);
    free_matrix(A_test);
    free_matrix(X);
    free_matrix(Y);
    free_vector(d);
    free_vector(e);
    free_vector(tauq);
    free_vector(taup);
}

void prepare_labrd_run(integer m_A, integer n_A, integer nb_A, void *A, integer lda, void *d,
                       void *e, void *tauq, void *taup, void *X, integer ldx, void *Y, integer ldy,
                       integer datatype, integer interfacetype, test_params_t *params)
{
    void *A_save, *d_test, *e_test, *tauq_test, *taup_test;
    void *X_test, *Y_test;
    double exe_time;

    /* Make a copy of the input matrix A. Same input values will be passed in
       each itertaion.*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m_A, n_A, &A_save, lda);
    copy_matrix(datatype, "full", m_A, n_A, A, lda, A_save, lda);

    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore input matrix A value and allocate memory to output buffers
           for each iteration*/
        copy_matrix(datatype, "full", m_A, n_A, A_save, lda, A, lda);

        /* Create output matrices and vectors */
        create_matrix(datatype, LAPACK_COL_MAJOR, m_A, nb_A, &X_test, ldx);
        create_matrix(datatype, LAPACK_COL_MAJOR, n_A, nb_A, &Y_test, ldy);
        create_realtype_vector(datatype, &d_test, nb_A);
        create_realtype_vector(datatype, &e_test, nb_A);

        create_vector(datatype, &tauq_test, nb_A);
        create_vector(datatype, &taup_test, nb_A);

#if ENABLE_CPP_TEST
        if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            /* call CPP labrd API  */
            invoke_cpp_labrd(datatype, &m_A, &n_A, &nb_A, A, &lda, d_test, e_test, tauq_test,
                             taup_test, X_test, &ldx, Y_test, &ldy);
            exe_time = fla_test_clock() - exe_time;
        }
        else
#endif
        {
            exe_time = fla_test_clock();
            /* Call LAPACK labrd API */
            invoke_labrd(datatype, &m_A, &n_A, &nb_A, A, &lda, d_test, e_test, tauq_test, taup_test,
                         X_test, &ldx, Y_test, &ldy);

            exe_time = fla_test_clock() - exe_time;
        }

        /* Update ctx and loop conditions */
        FLA_EXEC_LOOP_UPDATE_NO_INFO

        /* Make a copy of the output buffers. This is required to validate the API functionality. */
        copy_matrix(datatype, "full", m_A, nb_A, X_test, ldx, X, ldx);
        copy_matrix(datatype, "full", n_A, nb_A, Y_test, ldy, Y, ldy);

        copy_realtype_vector(datatype, nb_A, d_test, 1, d, 1);
        copy_realtype_vector(datatype, nb_A, e_test, 1, e, 1);
        copy_vector(datatype, nb_A, tauq_test, 1, tauq, 1);
        copy_vector(datatype, nb_A, taup_test, 1, taup, 1);

        /* Free up the output buffers */
        free_matrix(X_test);
        free_matrix(Y_test);
        free_vector(d_test);
        free_vector(e_test);
        free_vector(tauq_test);
        free_vector(taup_test);
    }

    free_matrix(A_save);
}

void invoke_labrd(integer datatype, integer *m, integer *n, integer *nb, void *a, integer *lda,
                  void *d, void *e, void *tauq, void *taup, void *x, integer *ldx, void *y,
                  integer *ldy)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_slabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_dlabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
            break;
        }

        case COMPLEX:
        {
            fla_lapack_clabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            fla_lapack_zlabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
            break;
        }
    }
}

void store_labrd_outputs(void *filename, integer datatype, integer m, integer n, integer nb,
                         void *A, integer lda, void *d, void *e, void *tauq, void *taup, void *X,
                         integer ldx, void *Y, integer ldy, void *params)
{
    /* Create and open a file for storing Ground truth*/
    FLA_OPEN_GT_FILE_STORE

    /* Store the ground truth data */
    FLA_STORE_BRT_MATRIX(datatype, m, n, A, lda)
    FLA_STORE_BRT_VECTOR(get_realtype(datatype), nb, d)
    FLA_STORE_BRT_VECTOR(get_realtype(datatype), nb, e)
    FLA_STORE_BRT_VECTOR(datatype, nb, tauq)
    FLA_STORE_BRT_VECTOR(datatype, nb, taup)
    FLA_STORE_BRT_MATRIX_NB_DIAG(datatype, m, nb, nb, X, ldx)
    FLA_STORE_BRT_MATRIX_NB_DIAG(datatype, n, nb, nb, Y, ldy)

    fclose(gt_file);
}

integer check_bit_reproducibility_labrd(void *filename, integer datatype, integer m, integer n,
                                        integer nb, void *A, integer lda, void *d, void *e,
                                        void *tauq, void *taup, void *X, integer ldx, void *Y,
                                        integer ldy, void *params)
{
    /* Open the file for reading Ground truth */
    FLA_OPEN_GT_FILE_READ

    /* Load stored GT and verify with current API outputs */
    FLA_VERIFY_BRT_MATRIX(datatype, m, n, A, lda)
    FLA_VERIFY_BRT_VECTOR(get_realtype(datatype), nb, d)
    FLA_VERIFY_BRT_VECTOR(get_realtype(datatype), nb, e)
    FLA_VERIFY_BRT_VECTOR(datatype, nb, tauq)
    FLA_VERIFY_BRT_VECTOR(datatype, nb, taup)
    FLA_VERIFY_BRT_MATRIX_NB_DIAG(datatype, m, nb, nb, X, ldx)
    FLA_VERIFY_BRT_MATRIX_NB_DIAG(datatype, n, nb, nb, Y, ldy)

    fclose(gt_file);
    return 1;
}