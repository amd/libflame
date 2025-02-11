/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/* > \brief \b validate_gelsd.c                                              */
/* =========== DOCUMENTATION ===========                                     */
/* Definition:                                                               */
/* ===========                                                               */
/* SUBROUTINE validate_gelsd(integer m, integer n, integer NRHS, void *A,....*/
/*                           integer lda, void *B, integer ldb, void *S,     */
/*                           void *X, void *rcond, integer *rank,            */
/*                           integer datatype, double* residual);            */
/* > \par Purpose:                                                           */
/* =============                                                             */
/* >                                                                         */
/* > \verbatim                                                               */
/* >                                                                         */
/* > Defines validate function of GELSD() to use in test suite               */
/* >                                                                         */
/* > This API validates GELSD() functionality by performing                  */
/* > the folllowing tests:                                                   */
/* >                                                                         */
/* > TEST 1: If m >= n and RANK = n, the residual sum-of-squares             */
/* >         for the solution in the i-th column is given by the sum of      */
/* >         squares of elements n+1:m in that column                        */
/* >                                                                         */
/* > TEST 2: Check AX = B.....................................               */
/* >                                                                         */
/* > TEST 3: checks whether X is in the row space of A or A'                 */
/* >                                                                         */
/* > Test SVD: Check if Singular values generated are in order               */
/* >                                                                         */
/* > \endverbatim                                                            */
/* ========================================================================= */

#include "test_common.h"

extern double perf;
extern double time_min;

void validate_gelsd(char *tst_api, integer m, integer n, integer nrhs, void *A, integer lda,
                    void *B, integer ldb, void *S, void *X, void *rcond, integer *rank,
                    integer datatype, double err_thresh, char imatrix)
{
    void *work = NULL;
    char NORM = '1';
    integer ldx;
    ldx = ldb;
    double residual, resid1 = 0., resid2 = 0.;
    double resid3 = 0., resid4 = 0.;

    /* Early return conditions */
    if(m == 0 || n == 0 || nrhs == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    /* Test SVD
       Check the order of Singular values generated */
    resid1 = svd_check_order(datatype, S, m, n, err_thresh);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_slamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column */
            if((m >= n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, &resid2);
            }
            else
            {
                /* Test 2 */
                /* Compute |B-AX| = 0 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);

                /* Compute B-AX */
                sgemm_("N", "N", &m, &nrhs, &n, &s_n_one, A, &lda, X, &ldx, &s_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm, imatrix, work);

                resid3 = (double)(norm / norm_a)
                         / (norm_x * fla_max(m, fla_max(n, nrhs)) * norm_b * eps);

                /* Test 3
                 * checks whether X is in the row space of A or A'
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing an LQ factorization
                 * of [A',X]', and returning the norm of the trailing triangle,
                 * scaled by MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                /*if(!((imatrix == 'U') || (imatrix == 'O')))
                {
                    check_vector_in_rowspace(datatype, "N", m, n, nrhs, A, lda, X, ldb, resid4);
                }
                */
            }
            break;
        }
        case DOUBLE:
        {
            double norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_dlamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column*/
            if((m >= n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, &resid2);
            }
            else
            {
                /* Compute |B-AX| = 0 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);

                /* Compute B-AX */
                dgemm_("N", "N", &m, &nrhs, &n, &d_n_one, A, &lda, X, &ldx, &d_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm, imatrix, work);

                resid3 = (double)(norm / norm_a)
                         / (norm_x * fla_max(m, fla_max(n, nrhs)) * norm_b * eps);

                /* Test 3
                 * checks whether X is in the row space of A or A'
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing an LQ factorization
                 * of [A',X]', and returning the norm of the trailing triangle,
                 * scaled by MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                /*if(!((imatrix == 'U') || (imatrix == 'O')))
                {
                    check_vector_in_rowspace(datatype, "N", m, n, nrhs, A, lda, X, ldb, resid4);
                }*/
            }
            break;
        }
        case COMPLEX:
        {
            float norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_slamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column*/
            if((m >= n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, &resid2);
            }
            else
            {
                /* Compute |B-AX| = 0 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
                eps = fla_lapack_slamch("E");

                /* Compute B-AX */
                cgemm_("N", "N", &m, &nrhs, &n, &c_n_one, A, &lda, X, &ldx, &c_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm, imatrix, work);

                resid3 = (double)(norm / norm_a)
                         / (norm_x * fla_max(m, fla_max(n, nrhs)) * norm_b * eps);

                /* Test 3
                 * checks whether X is in the row space of A or A'
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing an LQ factorization
                 * of [A',X]', and returning the norm of the trailing triangle,
                 * scaled by MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                /*if(!((imatrix == 'U') || (imatrix == 'O')))
                {
                    check_vector_in_rowspace(datatype, "N", m, n, nrhs, A, lda, X, ldb, resid4);
                }*/
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_a, norm_b, norm_x, eps, norm = 0;
            eps = fla_lapack_dlamch("E");

            /* Test 1 */
            /* If m >= n and RANK = n, the residual
               sum-of-squares for the solution in the i-th column is given
               by the sum of squares of elements n+1:m in that column*/
            if((m >= n) && (*rank == n))
            {
                residual_sum_of_squares(datatype, m, n, nrhs, X, ldb, &resid2);
            }
            else
            {
                /* Compute |B-AX| = 0 */
                compute_matrix_norm(datatype, NORM, m, n, A, lda, &norm_a, imatrix, work);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm_b, imatrix, work);
                compute_matrix_norm(datatype, NORM, n, nrhs, X, ldx, &norm_x, imatrix, work);
                eps = fla_lapack_dlamch("E");

                /* Compute B-AX */
                zgemm_("N", "N", &m, &nrhs, &n, &z_n_one, A, &lda, X, &ldx, &z_one, B, &ldb);
                compute_matrix_norm(datatype, NORM, m, nrhs, B, ldb, &norm, imatrix, work);

                resid3 = (double)(norm / norm_a)
                         / (norm_x * fla_max(m, fla_max(n, nrhs)) * norm_b * eps);

                /* Test 3
                 * checks whether X is in the row space of A or A'
                 * by scaling both X and A such that their norms are in the range
                 * [sqrt(eps), 1/sqrt(eps)], then computing an LQ factorization
                 * of [A',X]', and returning the norm of the trailing triangle,
                 * scaled by MAX(M,N,NRHS)*eps.
                 * Currently disabled because of random failures: TODO.
                 */
                /*if(!((imatrix == 'U') || (imatrix == 'O')))
                {
                    check_vector_in_rowspace(datatype, "N", m, n, nrhs, A, lda, X, ldb, resid4);
                }*/
            }
            break;
        }
    }
    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);

    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
}
