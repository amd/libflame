/******************************************************************************
 * Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gels.c
 *  @brief Defines validate function of GELS() to use in test suite.
 *  */

#include "test_common.h"

void validate_gels(char *trans, integer m, integer n, integer nrhs, void *A, integer lda, void *B,
                   integer ldb, void *x, integer datatype, double *residual, integer *info)
{
    integer m1 = m, n1 = n, ldc = fla_max(m, fla_max(n, nrhs));
    char NORM = '1';
    void *work = NULL;
    void *C = NULL;

    if(*trans == 'T' || *trans == 'C')
    {
        m1 = n;
        n1 = m;
        NORM = 'I';
    }

    create_matrix(datatype, &C, ldc, n1);
    create_vector(datatype, &work, m);

    if(m == 0 || n == 0)
    {
        return;
    }
    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            float resid1 = 0, resid2 = 0;

            eps = fla_lapack_slamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'T' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, residual);
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_slange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_slange(&NORM, &m1, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_slange(&NORM, &n1, &nrhs, x, &i_one, work);
                sgemm_(trans, "N", &m1, &nrhs, &n1, &s_n_one, A, &lda, x, &ldb, &s_one, B, &ldb);
                norm = fla_lapack_slange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1, n1) * norm_a * norm_x * eps);

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                sgemm_("T", trans, &nrhs, &n1, &m1, &s_one, B, &ldb, A, &lda, &s_zero, C, &ldc);
                norm = fla_lapack_slange("1", &nrhs, &n1, C, &ldc, work);
                resid2 = norm / (fla_max(m1, fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                *residual = (double)fla_max(resid1, resid2);
            }
            break;
        }
        case DOUBLE:
        {
            double eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            double resid1 = 0, resid2 = 0;

            eps = fla_lapack_dlamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'T' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, residual);
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_dlange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_dlange(&NORM, &m1, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_dlange(&NORM, &n1, &nrhs, x, &i_one, work);
                dgemm_(trans, "N", &m1, &nrhs, &n1, &d_n_one, A, &lda, x, &ldb, &d_one, B, &ldb);
                norm = fla_lapack_dlange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1, n1) * norm_a * norm_x * eps);

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                dgemm_("T", trans, &nrhs, &n1, &m1, &d_one, B, &ldb, A, &lda, &d_zero, C, &ldc);
                norm = fla_lapack_dlange("1", &nrhs, &n1, C, &ldc, work);
                resid2 = norm / (fla_max(m1, fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                *residual = fla_max(resid1, resid2);
            }
            break;
        }
        case COMPLEX:
        {
            float eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            float resid1 = 0, resid2 = 0;

            eps = fla_lapack_slamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'C' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */
                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, residual);
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_clange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_clange(&NORM, &m1, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_clange(&NORM, &n1, &nrhs, x, &i_one, work);
                cgemm_(trans, "N", &m1, &nrhs, &n1, &c_n_one, A, &lda, x, &ldb, &c_one, B, &ldb);
                norm = fla_lapack_clange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1, n1) * norm_a * norm_x * eps);
                *residual = (double)resid1;

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                cgemm_("T", trans, &nrhs, &n1, &m1, &c_one, B, &ldb, A, &lda, &c_zero, C, &ldc);
                norm = fla_lapack_clange("1", &nrhs, &n1, C, &ldc, work);
                resid2 = norm / (fla_max(m1, fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                *residual = (double)fla_max(resid1, resid2);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm = 0, norm_a = 0, norm_b = 0, norm_x = 0;
            double resid1 = 0, resid2 = 0;

            eps = fla_lapack_dlamch("E");
            if((*trans == 'N' && m > n) || (*trans == 'C' && m < n))
            {
                /* Test - 1
                 * If m > n and Trans == 'N' or m < n and Trans = 'T'
                 * then the residual sum of squares for the solution in each column
                 * is given by the sum of squares of elements m+1(n+1) to n(m) in that column.
                 */

                residual_sum_of_squares(datatype, m1, n1, nrhs, x, ldb, residual);
            }
            else
            {
                /* Test - 1
                 * Compute norm(B - A*x) / (max(m, n) * norm(A) * norm(x) * eps)
                 */
                norm_a = fla_lapack_zlange(&NORM, &m, &n, A, &lda, work);
                norm_b = fla_lapack_zlange(&NORM, &m1, &nrhs, B, &ldb, work);
                norm_x = fla_lapack_zlange(&NORM, &n1, &nrhs, x, &i_one, work);
                zgemm_(trans, "N", &m1, &nrhs, &n1, &z_n_one, A, &lda, x, &ldb, &z_one, B, &ldb);
                norm = fla_lapack_zlange(&NORM, &m1, &i_one, B, &ldb, work);
                resid1 = norm / (fla_max(m1, n1) * norm_a * norm_x * eps);
                *residual = (double)resid1;

                /* Test - 2
                 * Compute norm(B - A*x)**T * A / (max(m, n, nrhs) * norm(A) * norm(B) * eps)
                 */
                zgemm_("T", trans, &nrhs, &n1, &m1, &z_one, B, &ldb, A, &lda, &z_zero, C, &ldc);
                norm = fla_lapack_zlange("1", &nrhs, &n1, C, &ldc, work);
                resid2 = norm / (fla_max(m1, fla_max(n1, nrhs)) * norm_a * norm_b * eps);

                *residual = fla_max(resid1, resid2);
            }
            break;
        }
    }
    free_vector(work);
    free_matrix(C);
}
