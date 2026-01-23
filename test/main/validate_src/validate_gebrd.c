/*
    Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gebrd.c
 *  @brief Defines validate function of GEBRD() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

/*
 * Helper: Generate Q or P orthogonal matrix from reflectors in A, tau, for GEBRD.
 * If gen_q is true, generates Q (m x m), else generates P (n x n).
 */
static void gebrd_generate_orthogonal(integer datatype, integer m, integer n, integer k, void *QorP,
                                      integer ldqorp, void *tau, int gen_q)
{
    integer lwork = -1, info = 0;
    void *work = NULL;
    double dwork;
    float swork;
    scomplex cwork;
    dcomplex zwork;
    char *vect = gen_q ? "Q" : "P";

    /* Generate Q or P using the appropriate LAPACK routine and workspace query */
    switch(datatype)
    {
        case FLOAT:
        {
            swork = 0.f;
            /* Workspace query for SORGBR */
            fla_lapack_sorgbr(vect, &m, &n, &k, NULL, &ldqorp, NULL, &swork, &lwork, &info);
            lwork = get_work_value(datatype, &swork);
            create_vector(datatype, &work, lwork);
            /* Generate Q or P */
            fla_lapack_sorgbr(vect, &m, &n, &k, QorP, &ldqorp, tau, work, &lwork, &info);
            break;
        }
        case DOUBLE:
        {
            dwork = 0.0;
            /* Workspace query for DORGBR */
            fla_lapack_dorgbr(vect, &m, &n, &k, NULL, &ldqorp, NULL, &dwork, &lwork, &info);
            lwork = get_work_value(datatype, &dwork);
            create_vector(datatype, &work, lwork);
            /* Generate Q or P */
            fla_lapack_dorgbr(vect, &m, &n, &k, QorP, &ldqorp, tau, work, &lwork, &info);
            break;
        }
        case COMPLEX:
        {
            cwork.real = 0.f;
            cwork.imag = 0.f;
            /* Workspace query for CUNGBR */
            fla_lapack_cungbr(vect, &m, &n, &k, NULL, &ldqorp, NULL, &cwork, &lwork, &info);
            lwork = get_work_value(datatype, &cwork);
            create_vector(datatype, &work, lwork);
            /* Generate Q or P */
            fla_lapack_cungbr(vect, &m, &n, &k, QorP, &ldqorp, tau, work, &lwork, &info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            zwork.real = 0.0;
            zwork.imag = 0.0;
            /* Workspace query for ZUNGBR */
            fla_lapack_zungbr(vect, &m, &n, &k, NULL, &ldqorp, NULL, &zwork, &lwork, &info);
            lwork = get_work_value(datatype, &zwork);
            create_vector(datatype, &work, lwork);
            /* Generate Q or P */
            fla_lapack_zungbr(vect, &m, &n, &k, QorP, &ldqorp, tau, work, &lwork, &info);
            break;
        }
        default:
            /* Handle unexpected datatype */
            break;
    }
    free_vector(work);
}

void validate_gebrd(integer datatype, char *tst_api, integer m, integer n, void *A, integer lda,
                    void *A_test, integer ldat, void *d, void *e, void *tauq, void *taup,
                    double err_thresh, void *params)
{
    double residual = err_thresh;
    double resid1 = 0., resid2, resid3;
    integer k = fla_min(m, n);
    integer kq, kp;
    void *Q = NULL, *P = NULL, *B = NULL, *A_recon = NULL;

    /* Early return conditions */
    if(m == 0 || n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    /* Determine number of householder vectors kq and kp for Q and P */
    kq = n;
    kp = m;

    /* Allocate Q (m x m), P (n x n), B (m x n), A_recon (m x n) */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, k, &Q, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, k, n, &P, k);
    create_matrix(datatype, LAPACK_COL_MAJOR, k, k, &B, k);

    /* Generate Q and P from reflectors in A_test, tauq, taup */
    copy_matrix(datatype, "full", m, k, A_test, ldat, Q, m); /* Copy v vectors on to Q */
    gebrd_generate_orthogonal(datatype, m, k, kq, Q, m, tauq, 1); /* Q: m x m, kq = n */
    copy_matrix(datatype, "full", k, n, A_test, ldat, P, k);
    gebrd_generate_orthogonal(datatype, k, n, kp, P, k, taup, 0); /* P: n x n, kp = m */

    /* Form bidiagonal matrix B from d, e */
    if(m >= n)
    {
        build_bidiagonal_matrix(datatype, k, k, k, d, e, B, k, UPPER_BIDIAG);
    }
    else
    {
        build_bidiagonal_matrix(datatype, k, k, k, d, e, B, k, LOWER_BIDIAG);
    }

    /* Reconstruct A: A_recon = Q * B * P^T */
    void *temp1 = NULL;
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &A_recon, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, k, &temp1, m);
    fla_invoke_gemm(datatype, "N", "N", &m, &k, &k, d_one, Q, &m, B, &k, d_zero, temp1, &m);
    fla_invoke_gemm(datatype, "N", "N", &m, &n, &k, d_one, temp1, &m, P, &k, d_zero, A_recon, &m);
    free_matrix(temp1);

    /* Compute residual: norm(A_recon - A) / (eps * norm(A) * n) */
    double norm, norm_A;
    matrix_difference(datatype, m, n, A_recon, m, A, lda);
    switch(datatype)
    {
        case FLOAT:
        {
            norm_A = fla_lapack_slange("1", &m, &n, A, &lda, NULL);
            norm = fla_lapack_slange("1", &m, &n, A_recon, &m, NULL);
            break;
        }
        case DOUBLE:
        {
            norm_A = fla_lapack_dlange("1", &m, &n, A, &lda, NULL);
            norm = fla_lapack_dlange("1", &m, &n, A_recon, &m, NULL);
            break;
        }
        case COMPLEX:
        {
            norm_A = fla_lapack_clange("1", &m, &n, A, &lda, NULL);
            norm = fla_lapack_clange("1", &m, &n, A_recon, &m, NULL);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            norm_A = fla_lapack_zlange("1", &m, &n, A, &lda, NULL);
            norm = fla_lapack_zlange("1", &m, &n, A_recon, &m, NULL);
            break;
        }
        default:
        {
            norm_A = 1.0;
            norm = 0.0;
            break;
        }
    }
    resid1 = fla_compute_residual(datatype, 'P', norm, norm_A, n, params);

    /* Check orthogonality of Q and P */
    resid2 = check_orthogonality(datatype, Q, m, k, m, params);
    resid3 = (double)check_orthogonal_matrix('N', datatype, P, k, n, k, k, params);

    residual = fla_max(resid1, resid2);
    residual = fla_max(residual, resid3);
    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");

    free_matrix(Q);
    free_matrix(P);
    free_matrix(B);
    free_matrix(A_recon);
}
