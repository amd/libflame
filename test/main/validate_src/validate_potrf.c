/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_potrf.c
 *  @brief Defines validate function of POTRF() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_potrf(char *tst_api, char *uplo, integer m, void *A, void *A_test, integer lda,
                    integer datatype, double err_thresh)
{
    void *b = NULL, *x = NULL;
    void *x_test = NULL, *b_test = NULL;
    void *work = NULL;
    void *A_save = NULL;
    void *buff_A = NULL, *buff_B = NULL;
    integer nrhs = 1, incx = 1, incy = 1;
    char trans_A, trans_B;
    integer info = 0;
    double residual, resid1 = 0., resid2 = 0.;

    /* Early return conditions */
    if(m == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, m, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, m, err_thresh);

    /* Create Matrix buff_A and buff_B for computing LL'*/
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &buff_A, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &buff_B, m);

    /* Reset the matrix to avoid junk entries */
    reset_matrix(datatype, m, m, buff_A, m);
    reset_matrix(datatype, m, m, buff_B, m);

    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &A_save, lda);

    /* Create vector to compute Ax-b */
    create_vector(datatype, &b, m);
    create_vector(datatype, &x, m);
    create_vector(datatype, &b_test, m);
    create_vector(datatype, &x_test, m);

    /* Generate random vector b */
    rand_vector(datatype, m, b, 1, d_zero, d_zero, 'R');

    /* Copy lower or upper triangular matrix based on uplo to buff_A, buff_B */
    copy_matrix(datatype, uplo, m, m, A_test, lda, buff_A, m);
    copy_matrix(datatype, uplo, m, m, A_test, lda, buff_B, m);
    copy_matrix(datatype, "full", m, m, A, lda, A_save, lda);

    copy_vector(datatype, m, b, 1, b_test, 1);

    /* set transpose flag based on uplo */
    set_transpose(datatype, uplo, &trans_A, &trans_B);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_b, norm, eps;
            float norm_A;

            /* Test 1 */
            norm_A = fla_lapack_slange("1", &m, &m, A, &lda, work);

            /* Compute LL'-A */
            sgemm_(&trans_A, &trans_B, &m, &m, &m, &s_one, buff_A, &m, buff_B, &m, &s_n_one, A,
                   &lda);

            norm = fla_lapack_slange("1", &m, &m, A, &lda, work);
            eps = fla_lapack_slamch("P");

            resid1 = norm / (eps * norm_A * (float)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, lda, A, lda);
            norm_b = snrm2_(&m, b, &incx);

            /* Find x to compute Ax-b */
            fla_lapack_spotrs(uplo, &m, &nrhs, A_test, &lda, b_test, &m, &info);
            if(info < 0)
                break;

            copy_vector(datatype, m, b_test, 1, x, 1);

            /* Compute Ax-b */
            sgemv_("N", &m, &m, &s_one, A, &lda, x, &incx, &s_n_one, b, &incy);
            norm = snrm2_(&m, b, &incx);

            resid2 = norm / (eps * norm_b * (float)m);
            break;
        }
        case DOUBLE:
        {
            double norm_b, norm, eps;
            double norm_A;

            /* Test 1 */
            norm_A = fla_lapack_dlange("1", &m, &m, A, &lda, work);

            /* compute L*L'-A */
            dgemm_(&trans_A, &trans_B, &m, &m, &m, &d_one, buff_A, &m, buff_B, &m, &d_n_one, A,
                   &lda);

            norm = fla_lapack_dlange("1", &m, &m, A, &lda, work);
            eps = fla_lapack_dlamch("P");

            resid1 = norm / (eps * norm_A * (double)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, lda, A, lda);
            norm_b = dnrm2_(&m, b, &incx);

            /*Compute Ax-b Linear equations and find x */
            fla_lapack_dpotrs(uplo, &m, &nrhs, A_test, &lda, b_test, &m, &info);
            if(info < 0)
                break;

            copy_vector(datatype, m, b_test, 1, x, 1);

            /* Compute Ax-b */
            dgemv_("N", &m, &m, &d_one, A, &lda, x, &incx, &d_n_one, b, &incy);
            norm = dnrm2_(&m, b, &incx);

            resid2 = norm / (eps * norm_b * (double)m);
            break;
        }
        case COMPLEX:
        {
            float norm_b, norm, eps;
            float norm_A;

            /* Test 1 */
            norm_A = fla_lapack_clange("1", &m, &m, A, &lda, work);

            /* compute L*L'-A */
            cgemm_(&trans_A, &trans_B, &m, &m, &m, &c_one, buff_A, &m, buff_B, &m, &c_n_one, A,
                   &lda);

            norm = fla_lapack_clange("1", &m, &m, A, &lda, work);
            eps = fla_lapack_slamch("P");

            resid1 = norm / (eps * norm_A * (float)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, lda, A, lda);
            norm_b = scnrm2_(&m, b, &incx);

            /*Find x to compute Ax-b */
            fla_lapack_cpotrs(uplo, &m, &nrhs, A_test, &lda, b_test, &m, &info);
            if(info < 0)
                break;

            copy_vector(datatype, m, b_test, 1, x, 1);

            /* Compute Ax-b */
            cgemv_("N", &m, &m, &c_one, A, &lda, x, &incx, &c_n_one, b, &incy);
            norm = scnrm2_(&m, b, &incx);

            resid2 = norm / (eps * norm_b * (float)m);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_b, norm, eps;
            double norm_A;

            /* Test 1 */
            norm_A = fla_lapack_zlange("1", &m, &m, A, &lda, work);

            /* compute L*L'-A */
            zgemm_(&trans_A, &trans_B, &m, &m, &m, &z_one, buff_A, &m, buff_B, &m, &z_n_one, A,
                   &lda);

            norm = fla_lapack_zlange("1", &m, &m, A, &lda, work);
            eps = fla_lapack_dlamch("P");

            resid1 = norm / (eps * norm_A * (double)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, lda, A, lda);
            norm_b = dznrm2_(&m, b, &incx);

            /* Find x to compute Ax-b */
            fla_lapack_zpotrs(uplo, &m, &nrhs, A_test, &lda, b_test, &m, &info);
            if(info < 0)
                break;

            copy_vector(datatype, m, b_test, 1, x, 1);

            /* Compute Ax-b */
            zgemv_("N", &m, &m, &z_one, A, &lda, x, &incx, &z_n_one, b, &incy);

            norm = dznrm2_(&m, b, &incx);
            eps = fla_lapack_dlamch("P");

            resid2 = norm / (eps * norm_b * (double)m);
            break;
        }
    }
    /* Free up buffers */
    free_matrix(b);
    free_matrix(x);
    free_matrix(x_test);
    free_matrix(b_test);
    free_matrix(A_save);
    free_matrix(buff_A);
    free_matrix(buff_B);

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(m, m, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
