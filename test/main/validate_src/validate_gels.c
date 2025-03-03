/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gels.c
 *  @brief Defines validate function of GELS() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_gels(char *tst_api, char *trans, integer m, integer n, integer nrhs, void *A,
                   integer lda, void *B, integer ldb, void *x, integer datatype, double err_thresh,
                   char imatrix)
{
    integer m1 = m, n1 = n, ldc = fla_max(m, fla_max(n, nrhs));
    char NORM = '1';
    void *work = NULL;
    void *C = NULL;
    double residual;
    double resid1 = 0., resid2 = 0., resid3 = 0., resid4 = 0.;

    /* Early return conditions */
    if(m == 0 || n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
        return;
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    if(same_char(*trans, 'T') || same_char(*trans, 'C'))
    {
        m1 = n;
        n1 = m;
        NORM = 'I';
    }

    create_matrix(datatype, LAPACK_COL_MAJOR, nrhs, n1, &C, ldc);
    create_vector(datatype, &work, fla_max(m, nrhs));

    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            eps = fla_lapack_slamch("E");
            if((same_char(*trans, 'N') && m > n) || (same_char(*trans, 'T') && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, &resid1);
            }
            else
            {
                /* Test - 2
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m1, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n1, nrhs, x, i_one, &norm_x, imatrix, work);
                sgemm_(trans, "N", &m1, &nrhs, &n1, &s_n_one, A, &lda, x, &ldb, &s_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m1, i_one, B, ldb, &norm, imatrix, work);
                resid2 = (norm / norm_a) / (fla_max(m1, n1) * norm_x * eps);

                /* Test - 3
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                sgemm_("T", trans, &nrhs, &n1, &m1, &s_one, B, &ldb, A, &lda, &s_zero, C, &ldc);
                compute_matrix_norm(datatype, NORM, nrhs, n1, C, ldc, &norm, imatrix, work);
                resid3 = (norm / norm_a) / (fla_max(m1, fla_max(n1, nrhs)) * norm_b * eps);

                /* Test - 4
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                // check_vector_in_rowspace(datatype, trans, m, n, nrhs, A, lda, x, ldb, &resid4);
            }
            break;
        }
        case DOUBLE:
        {
            double eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            eps = fla_lapack_dlamch("E");
            if((same_char(*trans, 'N') && m > n) || (same_char(*trans, 'T') && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, &resid1);
            }
            else
            {
                /* Test - 2
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m1, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n1, nrhs, x, i_one, &norm_x, imatrix, work);
                dgemm_(trans, "N", &m1, &nrhs, &n1, &d_n_one, A, &lda, x, &ldb, &d_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m1, i_one, B, ldb, &norm, imatrix, work);
                resid2 = (norm / norm_a) / (fla_max(m1, n1) * norm_x * eps);

                /* Test - 3
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                dgemm_("T", trans, &nrhs, &n1, &m1, &d_one, B, &ldb, A, &lda, &d_zero, C, &ldc);
                compute_matrix_norm(datatype, NORM, nrhs, n1, C, ldc, &norm, imatrix, work);
                resid3 = (norm / norm_a) / (fla_max(m1, fla_max(n1, nrhs)) * norm_b * eps);

                /* Test - 4
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                // check_vector_in_rowspace(datatype, trans, m, n, nrhs, A, lda, x, ldb, &resid4);
            }
            break;
        }
        case COMPLEX:
        {
            float eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            eps = fla_lapack_slamch("E");
            if((same_char(*trans, 'N') && m > n) || (same_char(*trans, 'C') && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, &resid1);
            }
            else
            {
                /* Test - 2
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m1, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n1, nrhs, x, i_one, &norm_x, imatrix, work);
                cgemm_(trans, "N", &m1, &nrhs, &n1, &c_n_one, A, &lda, x, &ldb, &c_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m1, i_one, B, ldb, &norm, imatrix, work);
                resid2 = (norm / norm_a) / (fla_max(m1, n1) * norm_x * eps);

                /* Test - 3
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                cgemm_("T", trans, &nrhs, &n1, &m1, &c_one, B, &ldb, A, &lda, &c_zero, C, &ldc);
                compute_matrix_norm(datatype, NORM, nrhs, n1, C, ldc, &norm, imatrix, work);
                resid3 = (norm / norm_a) / (fla_max(m1, fla_max(n1, nrhs)) * norm_b * eps);

                /* Test - 4
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                // check_vector_in_rowspace(datatype, trans, m, n, nrhs, A, lda, x, ldb, &resid4);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            eps = fla_lapack_dlamch("E");
            if((same_char(*trans, 'N') && m > n) || (same_char(*trans, 'C') && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */

                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, &resid1);
            }
            else
            {
                /* Test - 2
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m1, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n1, nrhs, x, i_one, &norm_x, imatrix, work);
                zgemm_(trans, "N", &m1, &nrhs, &n1, &z_n_one, A, &lda, x, &ldb, &z_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m1, i_one, B, ldb, &norm, imatrix, work);
                resid2 = (norm / norm_a) / (fla_max(m1, n1) * norm_x * eps);

                /* Test - 3
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                zgemm_("T", trans, &nrhs, &n1, &m1, &z_one, B, &ldb, A, &lda, &z_zero, C, &ldc);
                compute_matrix_norm(datatype, NORM, nrhs, n1, C, ldc, &norm, imatrix, work);
                resid3 = (norm / norm_a) / (fla_max(m1, fla_max(n1, nrhs)) * norm_b * eps);

                /* Test - 4
                 * checks whether X is in the row space of A or A'.  It does so
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing a QR factorization of [A,X]
                 * (if TRANS = 'T') or an LQ factorization of [A',X]' (if TRANS = 'N'),
                 * and returning the norm of the trailing triangle, scaled by
                 * MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                // check_vector_in_rowspace(datatype, trans, m, n, nrhs, A, lda, x, ldb, &resid4);
            }
            break;
        }
    }
    free_vector(work);
    free_matrix(C);

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);

    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
}
