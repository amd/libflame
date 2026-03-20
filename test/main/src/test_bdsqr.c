/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_lapack.h"
#if ENABLE_CPP_TEST
#include <invoke_common.hh>
#endif
#include <invoke_lapacke.h>

extern double perf;
extern double time_min;
integer row_major_bdsqr_ldvt;
integer row_major_bdsqr_ldu;
integer row_major_bdsqr_ldc;

/* Local prototypes */
void fla_test_bdsqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo);
void prepare_bdsqr_run(char uplo, integer n, integer ncvt, integer nru, integer ncc, void *D,
                       void *E, void *VT, integer ldvt, void *U, integer ldu, void *C, integer ldc,
                       integer datatype, integer *info, integer interfacetype, int layout,
                       test_params_t *params);
void invoke_bdsqr(integer datatype, char *uplo, integer *n, integer *ncvt, integer *nru,
                  integer *ncc, void *D, void *E, void *VT, integer *ldvt, void *U, integer *ldu,
                  void *C, integer *ldc, void *work, integer *info);
double prepare_lapacke_bdsqr_run(integer datatype, int layout, char uplo, integer n, integer ncvt,
                                 integer nru, integer ncc, void *D, void *E, void *VT, integer ldvt,
                                 void *U, integer ldu, void *C, integer ldc, integer *info,
                                 void *work);

/* Helper functions for Bit reproducibility tests */
void store_bdsqr_outputs(void *filename, void *params, integer datatype, integer n, integer ncvt,
                         integer nru, integer ncc, void *d, void *e, void *VT, integer ldvt,
                         void *U, integer ldu, void *C, integer ldc);
integer check_bdsqr_reproducibility(void *filename, void *params, integer datatype, integer n,
                                    integer ncvt, integer nru, integer ncc, void *d, void *e,
                                    void *VT, integer ldvt, void *U, integer ldu, void *C,
                                    integer ldc);

void fla_test_bdsqr(integer argc, char **argv, test_params_t *params)
{
    char *op_str = "Singular value decomposition of bidiagonal matrix";
    char *front_str = "BDSQR";
    integer tests_not_run = 1, invalid_dtype = 0, einfo = 0;
    params->imatrix_char = '\0';

    if(argc == 1)
    {
        g_config_data = 1;
        fla_test_output_info("--- %s ---\n", op_str);
        fla_test_output_info("\n");
        fla_test_op_driver(front_str, SQUARE_INPUT, params, SVD, fla_test_bdsqr_experiment);
        tests_not_run = 0;
    }
    if(argc == 13)
    {
        FLA_TEST_PARSE_LAST_ARG(argv[12]);
    }
    if(argc >= 12 && argc <= 13)
    {
        integer i, num_types, N, NCVT, NRU, NCC;
        integer datatype, n_repeats;
        char stype, type_flag[4] = {0};
        char *endptr;

        /* Parse the arguments */
        num_types = strlen(argv[2]);
        params->svd_paramslist[0].uplo = argv[3][0];
        N = strtoimax(argv[4], &endptr, CLI_DECIMAL_BASE);
        NCVT = strtoimax(argv[5], &endptr, CLI_DECIMAL_BASE);
        NRU = strtoimax(argv[6], &endptr, CLI_DECIMAL_BASE);
        NCC = strtoimax(argv[7], &endptr, CLI_DECIMAL_BASE);
        params->svd_paramslist[0].ncvt_bdsqr = NCVT;
        params->svd_paramslist[0].nru_bdsqr = NRU;
        params->svd_paramslist[0].ncc_bdsqr = NCC;

        /* In case of command line inputs for LAPACKE row_major layout save leading dimensions */
        if((g_ext_fptr == NULL) && (params->interfacetype == LAPACKE_ROW_TEST))
        {
            row_major_bdsqr_ldvt = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            row_major_bdsqr_ldu = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            row_major_bdsqr_ldc = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
            /* For internal column-major storage:
               VT: n x ncvt  -> ldvt = n
               U : nru x n   -> ldu  = nru
               C : n x ncc   -> ldc  = n   */
            params->svd_paramslist[0].ldvt = N;
            params->svd_paramslist[0].ldu = NRU;
            params->svd_paramslist[0].ldc = N;
        }
        else
        {
            params->svd_paramslist[0].ldvt = strtoimax(argv[8], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldu = strtoimax(argv[9], &endptr, CLI_DECIMAL_BASE);
            params->svd_paramslist[0].ldc = strtoimax(argv[10], &endptr, CLI_DECIMAL_BASE);
        }
        n_repeats = strtoimax(argv[11], &endptr, CLI_DECIMAL_BASE);
        params->n_repeats = n_repeats;

        if(n_repeats > 0)
        {
            params->svd_paramslist[0].svd_threshold = CLI_NORM_THRESH;

            for(i = 0; i < num_types; i++)
            {
                stype = argv[2][i];
                datatype = get_datatype(stype);

                /* Check for invalid datatype */
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
                fla_test_bdsqr_experiment(front_str, params, datatype, N, N, 0, n_repeats, einfo);
                tests_not_run = 0;
            }
        }
    }

    /* Print error messages */
    if(tests_not_run == 1)
    {
        printf("\nIllegal arguments for bdsqr\n");
        printf("./<EXE> bdsqr <precisions - sdcz> <UPLO> <N> <NCVT> <NRU> <NCC> <LDVT> <LDU> "
               "<LDC> <repeats>\n");
    }
    if(invalid_dtype)
    {
        printf("\nInvalid datatypes specified, choose valid datatypes from 'sdcz'\n\n");
    }
}

void fla_test_bdsqr_experiment(char *tst_api, test_params_t *params, integer datatype,
                               integer p_cur, integer q_cur, integer pci, integer n_repeats,
                               integer einfo)
{
    char uplo;
    integer n, ncvt, nru, ncc, ldvt, ldu, ldc;
    integer info = 0;
    void *d = NULL, *e = NULL, *VT = NULL, *U = NULL, *C = NULL;
    void *d_save = NULL, *e_save = NULL, *U_save = NULL, *VT_save = NULL, *C_save = NULL;
    double residual, err_thresh;
    void *filename = NULL;

    integer interfacetype = params->interfacetype;
    int layout = params->matrix_major;
    integer realtype = get_realtype(datatype);

    /* Get input matrix dimensions */
    uplo = params->svd_paramslist[pci].uplo;
    err_thresh = params->svd_paramslist[pci].svd_threshold;

    n = p_cur;

    /* Get BDSQR-specific dimension from config if available */
    ncvt = params->svd_paramslist[pci].ncvt_bdsqr;
    nru = params->svd_paramslist[pci].nru_bdsqr;
    ncc = params->svd_paramslist[pci].ncc_bdsqr;
    ldvt = params->svd_paramslist[pci].ldvt;
    ldu = params->svd_paramslist[pci].ldu;
    ldc = params->svd_paramslist[pci].ldc;

    /* If coming from config and any LD is -1 assign defaults */
    /* Config mode: scale dimensions proportionally and set default LDs */
    if(g_config_data && n > 0)
    {
        /* Clamp dimensions to n for config-based runs */
        ncvt = (ncvt > 0) ? fla_min(ncvt, n) : 0;
        nru = (nru > 0) ? fla_min(nru, n) : 0;
        ncc = (ncc > 0) ? fla_min(ncc, n) : 0;

        /* Set default leading dimensions if not specified */
        if(ldvt == -1)
            ldvt = (ncvt > 0) ? fla_max(1, n) : 1;
        if(ldu == -1)
            ldu = (nru > 0) ? fla_max(1, nru) : 1;
        if(ldc == -1)
            ldc = (ncc > 0) ? fla_max(1, n) : 1;
    }

    /* Allocate (D,E real; U,VT,C in datatype) */
    create_realtype_vector(realtype, &d, n);
    create_realtype_vector(realtype, &e, (n > 1 ? n - 1 : 1));

    if(ncvt > 0)
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncvt, &VT, ldvt);
    if(nru > 0)
        create_matrix(datatype, LAPACK_COL_MAJOR, nru, n, &U, ldu);
    if(ncc > 0)
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncc, &C, ldc);

    /* This code path is run to generate the matrix to be passed to the API. This is the default
     * input generation logic accessed both when BRT is run in Ground truth mode and for non BRT
     * Test cases. For verification runs the input is loaded from the input generated during Ground
     * truth run */
    if(!FLA_BRT_VERIFICATION_RUN)
    {
        /* Generate U and VT from GEBRD using ORGBR/UNGBR
         * Per Netlib: U and VT should come from bidiagonal reduction.
         * For nru > n or ncvt > n, only the first n rows of U and first n columns
         * of VT come from GEBRD and are used in validation. Extra rows/columns
         * should be initialized to preserve the semi-orthogonal structure.
         */
        if(nru > 0 || ncvt > 0)
        {
            void *A_copy = NULL;
            void *tauq_copy = NULL, *taup_copy = NULL;
            void *work_orgbr = NULL;
            integer lwork_orgbr = -1;
            integer info_orgbr = 0;

            /* Generate GEBRD on n×n matrix - use temporary d_temp, e_temp */
            void *d_temp = NULL, *e_temp = NULL;
            create_realtype_vector(realtype, &d_temp, n);
            create_realtype_vector(realtype, &e_temp, (n > 1 ? n - 1 : 1));

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &A_copy, n);
            rand_matrix(datatype, A_copy, n, n, n);
            create_vector(datatype, &tauq_copy, n);
            create_vector(datatype, &taup_copy, n);

            /* Workspace query */
            create_vector(datatype, &work_orgbr, 1);
            invoke_gebrd(datatype, &n, &n, A_copy, &n, d_temp, e_temp, tauq_copy, taup_copy,
                         work_orgbr, &lwork_orgbr, &info_orgbr);
            lwork_orgbr = get_work_value(datatype, work_orgbr);
            if(lwork_orgbr < 4 * n)
                lwork_orgbr = 4 * n;
            free_vector(work_orgbr);
            create_vector(datatype, &work_orgbr, lwork_orgbr);

            /* Bidiagonal reduction */
            invoke_gebrd(datatype, &n, &n, A_copy, &n, d_temp, e_temp, tauq_copy, taup_copy,
                         work_orgbr, &lwork_orgbr, &info_orgbr);

            /* Copy the bidiagonal elements to d and e for testing */
            copy_realtype_vector(realtype, n, d_temp, 1, d, 1);
            if(n > 1)
                copy_realtype_vector(realtype, n - 1, e_temp, 1, e, 1);

            /* Generate U from Q (left orthogonal matrix from GEBRD) */
            if(nru > 0)
            {
                void *Q_mat = NULL;
                void *work_q = NULL;
                integer lwork_q = -1;
                integer info_q = 0;
                char vect_q = 'Q';

                /* Generate n×n Q matrix */
                create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_mat, n);
                copy_matrix(datatype, "full", n, n, A_copy, n, Q_mat, n);

                create_vector(datatype, &work_q, 1);
                invoke_orgbr(datatype, &vect_q, &n, &n, &n, Q_mat, &n, tauq_copy, work_q, &lwork_q,
                             &info_q);
                lwork_q = get_work_value(datatype, work_q);
                if(lwork_q < n)
                    lwork_q = n;
                free_vector(work_q);
                create_vector(datatype, &work_q, lwork_q);

                /* Generate Q */
                invoke_orgbr(datatype, &vect_q, &n, &n, &n, Q_mat, &n, tauq_copy, work_q, &lwork_q,
                             &info_q);

                /* First, initialize U to zero */
                reset_matrix(datatype, nru, n, U, ldu);

                /* Copy first min(nru,n) rows from Q to U */
                integer rows_to_copy = fla_min(nru, n);
                copy_matrix(datatype, "full", rows_to_copy, n, Q_mat, n, U, ldu);

                /* For nru > n: Set rows n..nru-1 to identity extension
                 * This preserves column orthonormality: U^H * U = I_n */
                if(nru > n)
                {
                    /* Set extra rows to zeros (already done above) */
                    /* This ensures U has n orthonormal columns spanning first n rows */
                }

                free_matrix(Q_mat);
                free_vector(work_q);
            }

            /* Generate VT from P^T (right orthogonal matrix from GEBRD) */
            if(ncvt > 0)
            {
                void *PT_mat = NULL;
                void *work_p = NULL;
                integer lwork_p = -1;
                integer info_p = 0;
                char vect_p = 'P';

                /* Generate n×n P^T matrix */
                create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &PT_mat, n);
                copy_matrix(datatype, "full", n, n, A_copy, n, PT_mat, n);

                create_vector(datatype, &work_p, 1);
                invoke_orgbr(datatype, &vect_p, &n, &n, &n, PT_mat, &n, taup_copy, work_p, &lwork_p,
                             &info_p);
                lwork_p = get_work_value(datatype, work_p);
                if(lwork_p < n)
                    lwork_p = n;
                free_vector(work_p);
                create_vector(datatype, &work_p, lwork_p);

                /* Generate P^T */
                invoke_orgbr(datatype, &vect_p, &n, &n, &n, PT_mat, &n, taup_copy, work_p, &lwork_p,
                             &info_p);

                /* First, initialize VT to zero */
                reset_matrix(datatype, n, ncvt, VT, ldvt);

                /* Copy first min(ncvt,n) columns from P^T to VT */
                integer cols_to_copy = fla_min(ncvt, n);
                copy_matrix(datatype, "full", n, cols_to_copy, PT_mat, n, VT, ldvt);

                /* For ncvt > n: Set columns n..ncvt-1 to identity extension
                 * This preserves row orthonormality: VT * VT^H = I_n */
                if(ncvt > n)
                {
                    /* Set extra columns to zeros (already done above) */
                    /* This ensures VT has n orthonormal rows spanning first n columns */
                }

                free_matrix(PT_mat);
                free_vector(work_p);
            }

            /* Free d_temp and e_temp after both U and VT generation */
            free_vector(d_temp);
            free_vector(e_temp);

            free_matrix(A_copy);
            free_vector(tauq_copy);
            free_vector(taup_copy);
            free_vector(work_orgbr);
        }

        if(g_ext_fptr != NULL || FLA_EXTREME_CASE_TEST)
        {
            init_matrix(realtype, d, n, 1, n, g_ext_fptr, params->imatrix_char);
            if(n > 1)
                init_matrix(realtype, e, n - 1, 1, n - 1, g_ext_fptr, params->imatrix_char);
        }

        if(ncc > 0)
        {
            rand_matrix(datatype, C, n, ncc, ldc);
        }
    }

    /* This macro is used in the BRT test cases for the following purposes:
     *    - In the Ground truth runs (BRT_char => G, F), the inputs are stored in a file for future
     *      reference (so the same inputs can be re-used across runs)
     *    - In the verification runs (BRT_char => V, M), the inputs are loaded from the stored file
     *      instead of re-generating them, ensuring exact same inputs are used for comparison
     *
     * The format string parameter describes how CLI arguments are encoded into the unique filename.
     * We need to call the appropriate macro based on how many inputs we have.
     * BDSQR always has d and e (treated as n×1 and (n-1)×1 matrices), plus optionally VT, U, C
     *
     * IMPORTANT: d and e are ALWAYS real-valued vectors, regardless of whether the main computation
     * is real or complex. For complex BDSQR, d and e are still real (float for complex, double for
     * double complex).
     *
     * To ensure unique filenames for each precision, we explicitly pass the datatype character
     * in the format string, similar to other SVD APIs like GESVD and GEJSV.
     */
    char dtype_char = get_datatype_char(datatype); // Get the actual computation datatype character

    if(ncvt > 0 && nru > 0 && ncc > 0)
    {
        FLA_BRT_PROCESS_FIVE_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                   (n > 1 ? n - 1 : 1), datatype, n, ncvt, VT, ldvt, datatype, nru,
                                   n, U, ldu, datatype, n, ncc, C, ldc, "ccddddddd", dtype_char,
                                   uplo, n, ncvt, nru, ncc, ldvt, ldu, ldc)
    }
    else if(ncvt > 0 && nru > 0)
    {
        FLA_BRT_PROCESS_FOUR_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                   (n > 1 ? n - 1 : 1), datatype, n, ncvt, VT, ldvt, datatype, nru,
                                   n, U, ldu, "ccddddddd", dtype_char, uplo, n, ncvt, nru, ncc,
                                   ldvt, ldu, ldc)
    }
    else if(ncvt > 0 && ncc > 0)
    {
        FLA_BRT_PROCESS_FOUR_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                   (n > 1 ? n - 1 : 1), datatype, n, ncvt, VT, ldvt, datatype, n,
                                   ncc, C, ldc, "ccddddddd", dtype_char, uplo, n, ncvt, nru, ncc,
                                   ldvt, ldu, ldc)
    }
    else if(nru > 0 && ncc > 0)
    {
        FLA_BRT_PROCESS_FOUR_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                   (n > 1 ? n - 1 : 1), datatype, nru, n, U, ldu, datatype, n, ncc,
                                   C, ldc, "ccddddddd", dtype_char, uplo, n, ncvt, nru, ncc, ldvt,
                                   ldu, ldc)
    }
    else if(ncvt > 0)
    {
        FLA_BRT_PROCESS_THREE_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                    (n > 1 ? n - 1 : 1), datatype, n, ncvt, VT, ldvt, "ccddddddd",
                                    dtype_char, uplo, n, ncvt, nru, ncc, ldvt, ldu, ldc)
    }
    else if(nru > 0)
    {
        FLA_BRT_PROCESS_THREE_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                    (n > 1 ? n - 1 : 1), datatype, nru, n, U, ldu, "ccddddddd",
                                    dtype_char, uplo, n, ncvt, nru, ncc, ldvt, ldu, ldc)
    }
    else if(ncc > 0)
    {
        FLA_BRT_PROCESS_THREE_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                    (n > 1 ? n - 1 : 1), datatype, n, ncc, C, ldc, "ccddddddd",
                                    dtype_char, uplo, n, ncvt, nru, ncc, ldvt, ldu, ldc)
    }
    else
    {
        /* 2 inputs: d, e only */
        FLA_BRT_PROCESS_TWO_INPUT(realtype, n, 1, d, n, realtype, (n > 1 ? n - 1 : 1), 1, e,
                                  (n > 1 ? n - 1 : 1), "ccddddddd", dtype_char, uplo, n, ncvt, nru,
                                  ncc, ldvt, ldu, ldc)
    }

    /* Save originals for reconstruction */
    if(n > 0)
    {
        create_realtype_vector(realtype, &d_save, n);
        copy_realtype_vector(realtype, n, d, 1, d_save, 1);
    }
    if(n > 1)
    {
        create_realtype_vector(realtype, &e_save, n - 1);
        copy_realtype_vector(realtype, n - 1, e, 1, e_save, 1);
    }

    /* Save U_in for U reconstruction test */
    if(nru > 0 && n > 0)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, nru, n, &U_save, ldu);
        copy_matrix(datatype, "full", nru, n, U, ldu, U_save, ldu);
    }

    /* Save VT_in for VT reconstruction test */
    if(ncvt > 0 && n > 0)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncvt, &VT_save, ldvt);
        copy_matrix(datatype, "full", n, ncvt, VT, ldvt, VT_save, ldvt);
    }

    /* Save C_in for C reconstruction test */
    if(ncc > 0 && n > 0)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncc, &C_save, ldc);
        copy_matrix(datatype, "full", n, ncc, C, ldc, C_save, ldc);
    }

    /* Call the prepare function */
    prepare_bdsqr_run(uplo, n, ncvt, nru, ncc, d, e, VT, ldvt, U, ldu, C, ldc, datatype, &info,
                      interfacetype, layout, params);

    /* Performance computation :
     * https://support.nag.com/numeric/nl/nagdoc_latest/flhtml/f08/f08mbf.html*/
    perf = (double)(4.0 * n * n) / time_min / FLOPS_PER_UNIT_PERF;
    if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
        perf *= 4.0;

    /* output validation */
    FLA_TEST_CHECK_EINFO(residual, info, einfo);

    /* Bit reproducibility tests path
     * This path is taken when BRT is enabled.
     *  - In Ground-truth runs (BRT char G/F), store_bdsqr_outputs is invoked to write outputs
     *    to the generated filename.
     *  - In Verification runs (BRT char V/M), check_bdsqr_reproducibility loads stored outputs
     *    and verifies current outputs match bit-for-bit.
     */
    IF_FLA_BRT_VALIDATION(n, n,
                          store_bdsqr_outputs(filename, params, datatype, n, ncvt, nru, ncc, d, e,
                                              VT, ldvt, U, ldu, C, ldc),
                          validate_bdsqr(tst_api, n, d, d_save, e_save, U, ldu, VT, ldvt, U_save,
                                         VT_save, ncvt, nru, ncc, C, ldc, C_save, uplo, datatype,
                                         residual, err_thresh, g_ext_fptr, params->imatrix_char,
                                         params),
                          check_bdsqr_reproducibility(filename, params, datatype, n, ncvt, nru, ncc,
                                                      d, e, VT, ldvt, U, ldu, C, ldc))
    else if(FLA_SKIP_VALIDATION_MODE)
    {
        /* perf-only mode */
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    }
    else if(!FLA_EXTREME_CASE_TEST)
    {
        validate_bdsqr(tst_api, n, d, d_save, e_save, U, ldu, VT, ldvt, U_save, VT_save, ncvt, nru,
                       ncc, C, ldc, C_save, uplo, datatype, residual, err_thresh, g_ext_fptr,
                       params->imatrix_char, params);
    }
    /* Extreme-value tests (A/F/N/I from conformance_test.cmake) */
    else
    {
        /* Check that extreme values propagate to output singular values */
        if(!check_extreme_value(realtype, n, 1, d, 1, params->imatrix_char))
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
    FLA_FREE_FILENAME(filename)
    free_vector(d);
    free_vector(e);
    if(C_save)
        free_matrix(C_save);
    if(U_save)
        free_matrix(U_save);
    if(VT)
        free_matrix(VT);
    if(U)
        free_matrix(U);
    if(C)
        free_matrix(C);
    if(d_save)
        free_vector(d_save);
    if(e_save)
        free_vector(e_save);
    if(VT_save)
        free_matrix(VT_save);
}

void prepare_bdsqr_run(char uplo, integer n, integer ncvt, integer nru, integer ncc, void *D,
                       void *E, void *VT, integer ldvt, void *U, integer ldu, void *C, integer ldc,
                       integer datatype, integer *info, integer interfacetype, int layout,
                       test_params_t *params)
{
    void *D_save = NULL, *E_save = NULL;
    void *VT_base = NULL, *U_base = NULL, *C_base = NULL;
    void *work = NULL;
    double exe_time;
    integer realtype = get_realtype(datatype);
    /* LAPACK safe workspace size */
    integer lwork = 4 * fla_max(1, n);

    /* Save original bidiagonal (use realtype) */
    create_realtype_vector(realtype, &D_save, n);
    create_realtype_vector(realtype, &E_save, (n > 1 ? n - 1 : 1));
    copy_realtype_vector(realtype, n, D, 1, D_save, 1);
    if(n > 1)
        copy_realtype_vector(realtype, n - 1, E, 1, E_save, 1);

    /* ------------------------------------------------------------------
       Always build orthogonal bases (U_base, VT_base) using full n x n QR
       then slice. This overrides any imatrix_char for accumulators.
       C_base kept random (not required orthogonal).
       ------------------------------------------------------------------ */
    if(ncvt > 0)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncvt, &VT_base, ldvt);
        copy_matrix(datatype, "full", n, ncvt, VT, ldvt, VT_base, ldvt);
    }
    if(nru > 0)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, nru, n, &U_base, ldu);
        copy_matrix(datatype, "full", nru, n, U, ldu, U_base, ldu);
    }
    if(ncc > 0)
    {
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncc, &C_base, ldc);
        copy_matrix(datatype, "full", n, ncc, C, ldc, C_base, ldc);
    }

    *info = 0;
    FLA_EXEC_LOOP_BEGIN
    {
        /* Restore bidiagonal each repeat */
        copy_realtype_vector(realtype, n, D_save, 1, D, 1);
        if(n > 1)
            copy_realtype_vector(realtype, n - 1, E_save, 1, E, 1);

        /* Copy orthogonal/random bases into user working buffers */
        if(ncvt > 0)
            copy_matrix(datatype, "full", n, ncvt, VT_base, ldvt, VT, ldvt);
        if(nru > 0)
            copy_matrix(datatype, "full", nru, n, U_base, ldu, U, ldu);
        if(ncc > 0)
            copy_matrix(datatype, "full", n, ncc, C_base, ldc, C, ldc);

        create_realtype_vector(realtype, &work, lwork);

        if((interfacetype == LAPACKE_ROW_TEST) || (interfacetype == LAPACKE_COLUMN_TEST))
        {
            exe_time = prepare_lapacke_bdsqr_run(datatype, layout, uplo, n, ncvt, nru, ncc, D, E,
                                                 VT, ldvt, U, ldu, C, ldc, info, work);
        }
#if ENABLE_CPP_TEST
        else if(interfacetype == LAPACK_CPP_TEST)
        {
            exe_time = fla_test_clock();
            invoke_cpp_bdsqr(datatype, &uplo, &n, &ncvt, &nru, &ncc, D, E, VT, &ldvt, U, &ldu, C,
                             &ldc, work, info);
            exe_time = fla_test_clock() - exe_time;
        }
#endif
        else
        {
            exe_time = fla_test_clock();
            invoke_bdsqr(datatype, &uplo, &n, &ncvt, &nru, &ncc, D, E, VT, &ldvt, U, &ldu, C, &ldc,
                         work, info);
            exe_time = fla_test_clock() - exe_time;
        }

        FLA_EXEC_LOOP_UPDATE_WITH_INFO

        free_vector(work);
    }

    /* Free saves */
    free_vector(D_save);
    free_vector(E_save);
    if(VT_base)
        free_matrix(VT_base);
    if(U_base)
        free_matrix(U_base);
    if(C_base)
        free_matrix(C_base);
}

double prepare_lapacke_bdsqr_run(integer datatype, int layout, char uplo, integer n, integer ncvt,
                                 integer nru, integer ncc, void *D, void *E, void *VT, integer ldvt,
                                 void *U, integer ldu, void *C, integer ldc, integer *info,
                                 void *work)
{
    double exe_time;
    integer ldvt_t = ldvt;
    integer ldu_t = ldu;
    integer ldc_t = ldc;
    void *VT_t = NULL, *U_t = NULL, *C_t = NULL;

    /* Configure leading dimensions as per the input matrix layout */
    SELECT_LDA(g_ext_fptr, g_config_data, layout, ncvt, row_major_bdsqr_ldvt, ldvt_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, n, row_major_bdsqr_ldu, ldu_t);
    SELECT_LDA(g_ext_fptr, g_config_data, layout, ncc, row_major_bdsqr_ldc, ldc_t);

    VT_t = VT;
    U_t = U;
    C_t = C;

    /* In case of row_major matrix layout,
       convert input matrices to row_major */
    if(layout == LAPACK_ROW_MAJOR)
    {
        /* Create temporary buffers for converting matrix layout */
        if(ncvt > 0)
        {
            create_matrix(datatype, layout, n, ncvt, &VT_t, fla_max(ncvt, ldvt_t));
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, ncvt, VT, ldvt, VT_t, ldvt_t);
        }
        if(nru > 0)
        {
            create_matrix(datatype, layout, nru, n, &U_t, fla_max(n, ldu_t));
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, nru, n, U, ldu, U_t, ldu_t);
        }
        if(ncc > 0)
        {
            create_matrix(datatype, layout, n, ncc, &C_t, fla_max(ncc, ldc_t));
            convert_matrix_layout(LAPACK_COL_MAJOR, datatype, n, ncc, C, ldc, C_t, ldc_t);
        }
    }

    exe_time = fla_test_clock();

    /* Call LAPACKE bdsqr API */
    *info = invoke_lapacke_bdsqr(datatype, layout, uplo, n, ncvt, nru, ncc, D, E, VT_t, ldvt_t, U_t,
                                 ldu_t, C_t, ldc_t, work);

    exe_time = fla_test_clock() - exe_time;

    /* In case of row_major matrix layout, convert output matrices
       to column_major layout */
    if(layout == LAPACK_ROW_MAJOR)
    {
        if(ncvt > 0)
        {
            convert_matrix_layout(layout, datatype, n, ncvt, VT_t, ldvt_t, VT, ldvt);
            free_matrix(VT_t);
        }
        if(nru > 0)
        {
            convert_matrix_layout(layout, datatype, nru, n, U_t, ldu_t, U, ldu);
            free_matrix(U_t);
        }
        if(ncc > 0)
        {
            convert_matrix_layout(layout, datatype, n, ncc, C_t, ldc_t, C, ldc);
            free_matrix(C_t);
        }
    }

    return exe_time;
}

void invoke_bdsqr(integer datatype, char *uplo, integer *n, integer *ncvt, integer *nru,
                  integer *ncc, void *D, void *E, void *VT, integer *ldvt, void *U, integer *ldu,
                  void *C, integer *ldc, void *work, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
            fla_lapack_sbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, work, info);
            break;

        case DOUBLE:
            fla_lapack_dbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, work, info);
            break;

        case COMPLEX:
            fla_lapack_cbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, work, info);
            break;

        case DOUBLE_COMPLEX:
            fla_lapack_zbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, work, info);
            break;
    }
}

void store_bdsqr_outputs(void *filename, void *params, integer datatype, integer n, integer ncvt,
                         integer nru, integer ncc, void *d, void *e, void *VT, integer ldvt,
                         void *U, integer ldu, void *C, integer ldc)
{
    /* Create and open a file for storing Ground truth*/
    FLA_OPEN_GT_FILE_STORE

    integer realtype = get_realtype(datatype);

    /* Store the ground truth data */
    FLA_STORE_BRT_VECTOR(realtype, n, d)
    if(n > 1)
    {
        FLA_STORE_BRT_VECTOR(realtype, n - 1, e)
    }
    if(ncvt > 0)
    {
        FLA_STORE_BRT_MATRIX(datatype, n, ncvt, VT, ldvt)
    }
    if(nru > 0)
    {
        FLA_STORE_BRT_MATRIX(datatype, nru, n, U, ldu)
    }
    if(ncc > 0)
    {
        FLA_STORE_BRT_MATRIX(datatype, n, ncc, C, ldc)
    }

    FLA_CLOSE_GT_FILE_STORE
}

integer check_bdsqr_reproducibility(void *filename, void *params, integer datatype, integer n,
                                    integer ncvt, integer nru, integer ncc, void *d, void *e,
                                    void *VT, integer ldvt, void *U, integer ldu, void *C,
                                    integer ldc)
{
    /* Open the file for reading Ground truth */
    FLA_OPEN_GT_FILE_READ

    integer realtype = get_realtype(datatype);

    /* Load stored GT and verify with current API outputs */
    FLA_VERIFY_BRT_VECTOR(realtype, n, d)
    if(n > 1)
    {
        FLA_VERIFY_BRT_VECTOR(realtype, n - 1, e)
    }
    if(ncvt > 0)
    {
        FLA_VERIFY_BRT_MATRIX(datatype, n, ncvt, VT, ldvt)
    }
    if(nru > 0)
    {
        FLA_VERIFY_BRT_MATRIX(datatype, nru, n, U, ldu)
    }
    if(ncc > 0)
    {
        FLA_VERIFY_BRT_MATRIX(datatype, n, ncc, C, ldc)
    }

    fclose(gt_file);
    return 1;
}