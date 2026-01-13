/*
    Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_bdsqr.c
 *  @brief Defines validate function of BDSQR() to use in test suite.
 *
 *  NETLIB BDSQR Documentation Reference:
 *  =====================================
 *  BDSQR computes the singular value decomposition (SVD) of a real N-by-N
 *  bidiagonal matrix B:
 *
 *      B = Q * Σ * P^T
 *
 *  where:
 *    - Σ is the diagonal matrix of singular values (stored in D on output)
 *    - Q is an orthogonal matrix of left singular vectors (N × N)
 *    - P is an orthogonal matrix of right singular vectors (N × N)
 *
 *  BDSQR also updates the following matrices:
 *    - U_out  = U_in  * Q      (left accumulator update)
 *    - VT_out = P^T * VT_in    (right accumulator update)
 *    - C_out  = Q^T * C_in     (C matrix update, or Q^H for complex)
 *
 *  Validation Tests:
 *    TEST 01: Singular values ordering - non-negative and non-increasing
 *    TEST 02: B Reconstruction - Verify U_in * B * VT_in = U_out * Σ * VT_out
 *    TEST 03: U/VT Reconstruction with orthogonality check on Q or P^T
 *    TEST 04: C Reconstruction - Verify C_out = Q^T * C_in (for NRU >= N)
 *    TEST 05: Direct B = Q * Σ * P^T verification (when both U and VT are square)
 *    TEST 06: VT Reconstruction - Verify VT_out = P^T * VT_in (for NCVT >= N)
 *    TEST 07: U Reconstruction - Verify U_out = U_in * Q (for NRU >= N)
 */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_bdsqr(char *tst_api, integer n, void *d_out, void *d_in, void *e_in, void *U_out,
                    integer ldu, void *VT_out, integer ldvt, void *U_in, void *VT_in, integer ncvt,
                    integer nru, integer ncc, void *C_out, integer ldc, void *C_in, char uplo,
                    integer datatype, double residual_in, double err_thresh, FILE *g_ext_fptr,
                    char imatrix_char, void *params)
{
    double resid1 = 0., resid2 = 0., resid3 = 0., resid4 = 0., resid5 = 0., resid6 = 0.,
           resid7 = 0.;
    double residual = residual_in;
    double safe_min;
    integer realtype = get_realtype(datatype);
    char trans_char = FLA_IS_COMPLEXTYPE(datatype) ? 'C' : 'T';

    /* Early return: empty problem */
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }

    /* Handle negative / invalid-info tests */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, residual);

    if(residual > err_thresh)
    {
        FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
        return;
    }

    safe_min
        = (realtype == FLOAT) ? (double)fla_lapack_slamch("S") : (double)fla_lapack_dlamch("S");

    /*=========================================================================
     * TEST 01: Check singular values are non-negative and in non-increasing order
     *
     * NETLIB: "On exit, the singular values in non-increasing order"
     *=========================================================================*/
    if(d_out != NULL && n > 0)
    {
        resid1 = svd_check_order(realtype, d_out, n, 1, 0.0);
    }

    /*=========================================================================
     * TEST 02: B Reconstruction Test (works for ALL matrix shapes)
     *
     * Mathematical basis (from NETLIB):
     *   B = Q * Σ * P^T
     *   U_out = U_in * Q
     *   VT_out = P^T * VT_in
     *
     * Therefore: U_in * B * VT_in = U_in * (Q * Σ * P^T) * VT_in
     *                             = (U_in * Q) * Σ * (P^T * VT_in)
     *                             = U_out * Σ * VT_out
     *
     * Test: ||U_in * B * VT_in - U_out * Σ * VT_out||_F / (||B||_F * n * eps)
     *
     * Dimensions:
     *   U_in, U_out:   NRU × N   (can be tall, square, or short)
     *   VT_in, VT_out: N × NCVT  (can be wide, square, or narrow)
     *   B: N × N (bidiagonal), Σ: N × N (diagonal)
     *   Result: NRU × NCVT
     *=========================================================================*/
    if(U_out != NULL && U_in != NULL && VT_out != NULL && VT_in != NULL && nru > 0 && ncvt > 0
       && n > 0 && d_in != NULL && e_in != NULL && d_out != NULL)
    {
        void *B_mat = NULL;
        void *Sigma = NULL;
        void *temp1 = NULL;
        void *A_orig = NULL;
        void *A_rec = NULL;
        integer bidiag_type = (uplo == 'U' || uplo == 'u') ? UPPER_BIDIAG : LOWER_BIDIAG;

        /* Build B matrix (N × N bidiagonal) from d_in and e_in */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_mat, n);
        build_bidiagonal_matrix(datatype, n, n, n, d_in, e_in, B_mat, n, bidiag_type);

        /* Build Σ matrix (N × N diagonal) from d_out */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Sigma, n);
        diagonalize_realtype_vector(datatype, d_out, Sigma, n, n, n);

        /* Compute A_orig = U_in * B * VT_in */
        create_matrix(datatype, LAPACK_COL_MAJOR, nru, n, &temp1, nru);
        create_matrix(datatype, LAPACK_COL_MAJOR, nru, ncvt, &A_orig, nru);

        /* temp1 = U_in * B (NRU × N) */
        fla_invoke_gemm(datatype, "N", "N", &nru, &n, &n, U_in, &ldu, B_mat, &n, temp1, &nru);
        /* A_orig = temp1 * VT_in = U_in * B * VT_in (NRU × NCVT) */
        fla_invoke_gemm(datatype, "N", "N", &nru, &ncvt, &n, temp1, &nru, VT_in, &ldvt, A_orig,
                        &nru);

        /* Compute A_rec = U_out * Σ * VT_out */
        create_matrix(datatype, LAPACK_COL_MAJOR, nru, ncvt, &A_rec, nru);

        /* temp1 = U_out * Σ (NRU × N) */
        fla_invoke_gemm(datatype, "N", "N", &nru, &n, &n, U_out, &ldu, Sigma, &n, temp1, &nru);
        /* A_rec = temp1 * VT_out = U_out * Σ * VT_out (NRU × NCVT) */
        fla_invoke_gemm(datatype, "N", "N", &nru, &ncvt, &n, temp1, &nru, VT_out, &ldvt, A_rec,
                        &nru);

        /* Compute difference: A_rec = A_rec - A_orig */
        matrix_difference(datatype, nru, ncvt, A_rec, nru, A_orig, nru);

        /* Compute norms and residual */
        if(datatype == FLOAT || datatype == COMPLEX)
        {
            float norm_diff_f, norm_B_f, work_norm_f = 0.0f, eps_f = fla_lapack_slamch("P");
            compute_matrix_norm(datatype, 'F', nru, ncvt, A_rec, nru, &norm_diff_f, imatrix_char,
                                &work_norm_f);
            compute_matrix_norm(datatype, 'F', n, n, B_mat, n, &norm_B_f, imatrix_char,
                                &work_norm_f);
            resid2 = (norm_B_f > safe_min) ? (double)(norm_diff_f / (norm_B_f * n * eps_f))
                                           : (double)(norm_diff_f / (n * eps_f));
        }
        else
        {
            double norm_diff_d, norm_B_d, work_norm_d = 0.0, eps_d = fla_lapack_dlamch("P");
            compute_matrix_norm(datatype, 'F', nru, ncvt, A_rec, nru, &norm_diff_d, imatrix_char,
                                &work_norm_d);
            compute_matrix_norm(datatype, 'F', n, n, B_mat, n, &norm_B_d, imatrix_char,
                                &work_norm_d);
            resid2 = (norm_B_d > safe_min) ? norm_diff_d / (norm_B_d * n * eps_d)
                                           : norm_diff_d / (n * eps_d);
        }

        free_matrix(B_mat);
        free_matrix(Sigma);
        free_matrix(temp1);
        free_matrix(A_orig);
        free_matrix(A_rec);
    }

    /*=========================================================================
     * TEST 03: U-only or VT-only Reconstruction with Orthogonality Check
     *
     * Case A: Only U provided (ncvt == 0, nru > 0)
     *   For square U (NRU == N): Extract Q = U_in^{-1} * U_out
     *   Verify Q is orthogonal: Q^T * Q = I (or Q^H * Q = I for complex)
     *
     * Case B: Only VT provided (nru == 0, ncvt > 0)
     *   For square VT (NCVT == N): Extract P^T = VT_out * VT_in^{-1}
     *   Verify P^T is orthogonal: P^T * P = I
     *=========================================================================*/

    /* Case A: U-only with square U_in */
    if(U_out != NULL && U_in != NULL && nru == n && ncvt == 0 && n > 0)
    {
        void *U_in_inv = NULL;
        void *Q_mat = NULL;
        integer info_inv = 0;

        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &U_in_inv, n);
        copy_matrix(datatype, "full", n, n, U_in, ldu, U_in_inv, n);

        info_inv = compute_matrix_inverse(datatype, n, U_in_inv, n);

        if(info_inv == 0)
        {
            /* Q = U_in^{-1} * U_out */
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_mat, n);
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, U_in_inv, &n, U_out, &ldu, Q_mat, &n);

            /* Check orthogonality: Q^T * Q = I */
            resid3 = check_orthogonal_matrix(trans_char, datatype, Q_mat, n, n, n, n, params);

            free_matrix(Q_mat);
        }

        free_matrix(U_in_inv);
    }

    /* Case B: VT-only with square VT_in */
    if(VT_out != NULL && VT_in != NULL && ncvt == n && nru == 0 && n > 0 && resid3 == 0.)
    {
        void *VT_in_inv = NULL;
        void *PT_mat = NULL;
        integer info_inv = 0;

        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &VT_in_inv, n);
        copy_matrix(datatype, "full", n, n, VT_in, ldvt, VT_in_inv, n);

        info_inv = compute_matrix_inverse(datatype, n, VT_in_inv, n);

        if(info_inv == 0)
        {
            /* P^T = VT_out * VT_in^{-1} */
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &PT_mat, n);
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, VT_out, &ldvt, VT_in_inv, &n, PT_mat,
                            &n);

            /* Check orthogonality: P^T * P = I (i.e., 'N' mode: A * A^T = I) */
            resid3 = check_orthogonal_matrix('N', datatype, PT_mat, n, n, n, n, params);

            free_matrix(PT_mat);
        }

        free_matrix(VT_in_inv);
    }

    /*=========================================================================
     * TEST 04: C Reconstruction Test (Extended to NRU >= N)
     *
     * Mathematical basis (from NETLIB):
     *   "Optionally, the subroutine may also compute Q^T*C for a given
     *    real input matrix C" (or Q^H for complex)
     *
     * Equation: C_out = Q^T * C_in (or Q^H * C_in for complex)
     *
     * For NRU >= N, U_in has orthonormal columns (from ORGBR):
     *   U_in^T * U_in = I_N (or U_in^H * U_in = I_N for complex)
     *
     * We can extract Q using pseudo-inverse (which equals transpose):
     *   Q = U_in^T * U_out  (for real)
     *   Q = U_in^H * U_out  (for complex)
     *
     * Then verify: ||Q^T * C_in - C_out||_F / (||C_in||_F * n * eps)
     *=========================================================================*/
    if(C_out != NULL && C_in != NULL && ncc > 0 && U_out != NULL && U_in != NULL && nru >= n
       && n > 0)
    {
        void *Q_mat = NULL;
        void *C_rec = NULL;

        /* Q = U_in^T * U_out (works when U_in has orthonormal columns, i.e., NRU >= N) */
        /* For NRU >= N from ORGBR: U_in^T * U_in = I_N */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_mat, n);

        /* Q = U_in^H * U_out for complex, Q = U_in^T * U_out for real */
        fla_invoke_gemm(datatype, &trans_char, "N", &n, &n, &nru, U_in, &ldu, U_out, &ldu, Q_mat,
                        &n);

        /* C_rec = Q^T * C_in (or Q^H for complex) */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncc, &C_rec, n);
        fla_invoke_gemm(datatype, &trans_char, "N", &n, &ncc, &n, Q_mat, &n, C_in, &ldc, C_rec, &n);

        /* Compute difference and residual */
        matrix_difference(datatype, n, ncc, C_rec, n, C_out, ldc);

        if(datatype == FLOAT || datatype == COMPLEX)
        {
            float norm_diff_f, norm_C_in_f, work_norm_f = 0.0f, eps_f = fla_lapack_slamch("P");
            compute_matrix_norm(datatype, 'F', n, ncc, C_rec, n, &norm_diff_f, imatrix_char,
                                &work_norm_f);
            compute_matrix_norm(datatype, 'F', n, ncc, C_in, ldc, &norm_C_in_f, imatrix_char,
                                &work_norm_f);
            resid4 = (norm_C_in_f > safe_min) ? (double)(norm_diff_f / (norm_C_in_f * n * eps_f))
                                              : (double)(norm_diff_f / (n * eps_f));
        }
        else
        {
            double norm_diff_d, norm_C_in_d, work_norm_d = 0.0, eps_d = fla_lapack_dlamch("P");
            compute_matrix_norm(datatype, 'F', n, ncc, C_rec, n, &norm_diff_d, imatrix_char,
                                &work_norm_d);
            compute_matrix_norm(datatype, 'F', n, ncc, C_in, ldc, &norm_C_in_d, imatrix_char,
                                &work_norm_d);
            resid4 = (norm_C_in_d > safe_min) ? norm_diff_d / (norm_C_in_d * n * eps_d)
                                              : norm_diff_d / (n * eps_d);
        }

        free_matrix(C_rec);
        free_matrix(Q_mat);
    }

    /*=========================================================================
     * TEST 05: Direct B = Q * Σ * P^T Verification
     *
     * When both U and VT are square (NRU == N and NCVT == N):
     *   Q = U_in^{-1} * U_out
     *   P^T = VT_out * VT_in^{-1}
     *
     * Verify: B = Q * Σ * P^T
     * Test: ||B - Q * Σ * P^T||_F / (||B||_F * n * eps)
     *
     * This directly tests the core BDSQR equation.
     *=========================================================================*/
    if(U_out != NULL && U_in != NULL && VT_out != NULL && VT_in != NULL && nru == n && ncvt == n
       && n > 0 && d_in != NULL && e_in != NULL && d_out != NULL)
    {
        void *U_in_inv = NULL;
        void *VT_in_inv = NULL;
        void *Q_mat = NULL;
        void *PT_mat = NULL;
        void *B_mat = NULL;
        void *Sigma = NULL;
        void *temp1 = NULL;
        void *B_rec = NULL;
        integer info_inv_u = 0, info_inv_vt = 0;
        integer bidiag_type = (uplo == 'U' || uplo == 'u') ? UPPER_BIDIAG : LOWER_BIDIAG;

        /* Invert U_in */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &U_in_inv, n);
        copy_matrix(datatype, "full", n, n, U_in, ldu, U_in_inv, n);
        info_inv_u = compute_matrix_inverse(datatype, n, U_in_inv, n);

        /* Invert VT_in */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &VT_in_inv, n);
        copy_matrix(datatype, "full", n, n, VT_in, ldvt, VT_in_inv, n);
        info_inv_vt = compute_matrix_inverse(datatype, n, VT_in_inv, n);

        if(info_inv_u == 0 && info_inv_vt == 0)
        {
            /* Q = U_in^{-1} * U_out */
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_mat, n);
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, U_in_inv, &n, U_out, &ldu, Q_mat, &n);

            /* P^T = VT_out * VT_in^{-1} */
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &PT_mat, n);
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, VT_out, &ldvt, VT_in_inv, &n, PT_mat,
                            &n);

            /* Build B and Σ */
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_mat, n);
            build_bidiagonal_matrix(datatype, n, n, n, d_in, e_in, B_mat, n, bidiag_type);

            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Sigma, n);
            diagonalize_realtype_vector(datatype, d_out, Sigma, n, n, n);

            /* Compute B_rec = Q * Σ * P^T */
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &temp1, n);
            create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &B_rec, n);

            /* temp1 = Q * Σ */
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, Q_mat, &n, Sigma, &n, temp1, &n);
            /* B_rec = temp1 * P^T = Q * Σ * P^T */
            fla_invoke_gemm(datatype, "N", "N", &n, &n, &n, temp1, &n, PT_mat, &n, B_rec, &n);

            /* Compute difference: B_rec = B_rec - B */
            matrix_difference(datatype, n, n, B_rec, n, B_mat, n);

            /* Compute norms and residual */
            if(datatype == FLOAT || datatype == COMPLEX)
            {
                float norm_diff_f, norm_B_f, work_norm_f = 0.0f, eps_f = fla_lapack_slamch("P");
                compute_matrix_norm(datatype, 'F', n, n, B_rec, n, &norm_diff_f, imatrix_char,
                                    &work_norm_f);
                compute_matrix_norm(datatype, 'F', n, n, B_mat, n, &norm_B_f, imatrix_char,
                                    &work_norm_f);
                resid5 = (norm_B_f > safe_min) ? (double)(norm_diff_f / (norm_B_f * n * eps_f))
                                               : (double)(norm_diff_f / (n * eps_f));
            }
            else
            {
                double norm_diff_d, norm_B_d, work_norm_d = 0.0, eps_d = fla_lapack_dlamch("P");
                compute_matrix_norm(datatype, 'F', n, n, B_rec, n, &norm_diff_d, imatrix_char,
                                    &work_norm_d);
                compute_matrix_norm(datatype, 'F', n, n, B_mat, n, &norm_B_d, imatrix_char,
                                    &work_norm_d);
                resid5 = (norm_B_d > safe_min) ? norm_diff_d / (norm_B_d * n * eps_d)
                                               : norm_diff_d / (n * eps_d);
            }

            free_matrix(Q_mat);
            free_matrix(PT_mat);
            free_matrix(B_mat);
            free_matrix(Sigma);
            free_matrix(temp1);
            free_matrix(B_rec);
        }

        free_matrix(U_in_inv);
        free_matrix(VT_in_inv);
    }

    /*=========================================================================
     * TEST 06: VT Reconstruction Test (Extended to NCVT >= N)
     *
     * Mathematical basis (from NETLIB):
     *   VT_out = P^T * VT_in
     *
     * For NCVT >= N, VT_in has orthonormal rows (from ORGBR):
     *   VT_in * VT_in^T = I_N (or VT_in * VT_in^H = I_N for complex)
     *
     * We can extract P^T using right pseudo-inverse (which equals transpose):
     *   P^T = VT_out * VT_in^T  (for real)
     *   P^T = VT_out * VT_in^H  (for complex)
     *
     * Then verify: ||P^T * VT_in - VT_out||_F / (||VT_in||_F * n * eps)
     *
     * This is analogous to TEST 04 but for the right side (VT instead of U/C).
     *=========================================================================*/
    if(VT_out != NULL && VT_in != NULL && ncvt >= n && U_out != NULL && U_in != NULL && nru > 0
       && n > 0)
    {
        void *PT_mat = NULL;
        void *VT_rec = NULL;

        /* P^T = VT_out * VT_in^T (works when VT_in has orthonormal rows, i.e., NCVT >= N) */
        /* For NCVT >= N from ORGBR: VT_in * VT_in^T = I_N */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &PT_mat, n);

        /* P^T = VT_out * VT_in^H for complex, P^T = VT_out * VT_in^T for real */
        fla_invoke_gemm(datatype, "N", &trans_char, &n, &n, &ncvt, VT_out, &ldvt, VT_in, &ldvt,
                        PT_mat, &n);

        /* VT_rec = P^T * VT_in */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, ncvt, &VT_rec, n);
        fla_invoke_gemm(datatype, "N", "N", &n, &ncvt, &n, PT_mat, &n, VT_in, &ldvt, VT_rec, &n);

        /* Compute difference and residual */
        matrix_difference(datatype, n, ncvt, VT_rec, n, VT_out, ldvt);

        if(datatype == FLOAT || datatype == COMPLEX)
        {
            float norm_diff_f, norm_VT_in_f, work_norm_f = 0.0f, eps_f = fla_lapack_slamch("P");
            compute_matrix_norm(datatype, 'F', n, ncvt, VT_rec, n, &norm_diff_f, imatrix_char,
                                &work_norm_f);
            compute_matrix_norm(datatype, 'F', n, ncvt, VT_in, ldvt, &norm_VT_in_f, imatrix_char,
                                &work_norm_f);
            resid6 = (norm_VT_in_f > safe_min) ? (double)(norm_diff_f / (norm_VT_in_f * n * eps_f))
                                               : (double)(norm_diff_f / (n * eps_f));
        }
        else
        {
            double norm_diff_d, norm_VT_in_d, work_norm_d = 0.0, eps_d = fla_lapack_dlamch("P");
            compute_matrix_norm(datatype, 'F', n, ncvt, VT_rec, n, &norm_diff_d, imatrix_char,
                                &work_norm_d);
            compute_matrix_norm(datatype, 'F', n, ncvt, VT_in, ldvt, &norm_VT_in_d, imatrix_char,
                                &work_norm_d);
            resid6 = (norm_VT_in_d > safe_min) ? norm_diff_d / (norm_VT_in_d * n * eps_d)
                                               : norm_diff_d / (n * eps_d);
        }

        free_matrix(VT_rec);
        free_matrix(PT_mat);
    }

    /*=========================================================================
     * TEST 07: U Reconstruction Test (Extended to NRU >= N)
     *
     * Mathematical basis (from NETLIB):
     *   U_out = U_in * Q
     *
     * For NRU >= N, U_in has orthonormal columns (from ORGBR):
     *   U_in^T * U_in = I_N (or U_in^H * U_in = I_N for complex)
     *
     * We can extract Q using left pseudo-inverse (which equals transpose):
     *   Q = U_in^T * U_out  (for real)
     *   Q = U_in^H * U_out  (for complex)
     *
     * Then verify: ||U_in * Q - U_out||_F / (||U_in||_F * n * eps)
     *
     * This is the mirror of TEST 06 for the left side (U instead of VT).
     *=========================================================================*/
    if(U_out != NULL && U_in != NULL && nru >= n && VT_out != NULL && VT_in != NULL && ncvt > 0
       && n > 0)
    {
        void *Q_mat = NULL;
        void *U_rec = NULL;

        /* Q = U_in^T * U_out (works when U_in has orthonormal columns, i.e., NRU >= N) */
        /* For NRU >= N from ORGBR: U_in^T * U_in = I_N */
        create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &Q_mat, n);

        /* Q = U_in^H * U_out for complex, Q = U_in^T * U_out for real */
        fla_invoke_gemm(datatype, &trans_char, "N", &n, &n, &nru, U_in, &ldu, U_out, &ldu, Q_mat,
                        &n);

        /* U_rec = U_in * Q */
        create_matrix(datatype, LAPACK_COL_MAJOR, nru, n, &U_rec, nru);
        fla_invoke_gemm(datatype, "N", "N", &nru, &n, &n, U_in, &ldu, Q_mat, &n, U_rec, &nru);

        /* Compute difference and residual */
        matrix_difference(datatype, nru, n, U_rec, nru, U_out, ldu);

        if(datatype == FLOAT || datatype == COMPLEX)
        {
            float norm_diff_f, norm_U_in_f, work_norm_f = 0.0f, eps_f = fla_lapack_slamch("P");
            compute_matrix_norm(datatype, 'F', nru, n, U_rec, nru, &norm_diff_f, imatrix_char,
                                &work_norm_f);
            compute_matrix_norm(datatype, 'F', nru, n, U_in, ldu, &norm_U_in_f, imatrix_char,
                                &work_norm_f);
            resid7 = (norm_U_in_f > safe_min) ? (double)(norm_diff_f / (norm_U_in_f * n * eps_f))
                                              : (double)(norm_diff_f / (n * eps_f));
        }
        else
        {
            double norm_diff_d, norm_U_in_d, work_norm_d = 0.0, eps_d = fla_lapack_dlamch("P");
            compute_matrix_norm(datatype, 'F', nru, n, U_rec, nru, &norm_diff_d, imatrix_char,
                                &work_norm_d);
            compute_matrix_norm(datatype, 'F', nru, n, U_in, ldu, &norm_U_in_d, imatrix_char,
                                &work_norm_d);
            resid7 = (norm_U_in_d > safe_min) ? norm_diff_d / (norm_U_in_d * n * eps_d)
                                              : norm_diff_d / (n * eps_d);
        }

        free_matrix(U_rec);
        free_matrix(Q_mat);
    }

    /* Combine residuals - take the maximum */
    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(residual, resid3);
    residual = fla_test_max(residual, resid4);
    residual = fla_test_max(residual, resid5);
    residual = fla_test_max(residual, resid6);
    residual = fla_test_max(residual, resid7);

    /* Print results */
    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
    FLA_PRINT_SUBTEST_STATUS(resid5, err_thresh, "05");
    FLA_PRINT_SUBTEST_STATUS(resid6, err_thresh, "06");
    FLA_PRINT_SUBTEST_STATUS(resid7, err_thresh, "07");
}