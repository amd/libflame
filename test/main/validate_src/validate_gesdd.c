/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gesdd.c
 *  @brief Defines validate function of GESDD() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_gesdd(char *tst_api, char *jobz, integer m, integer n, void *A, void *A_test,
                    integer lda, void *s, void *s_in, void *U, integer ldu, void *V, integer ldvt,
                    integer datatype, double err_thresh, char imatrix, void *scal)
{
    void *sigma = NULL, *Usigma = NULL;
    void *work = NULL, *U_temp = NULL, *V_temp = NULL;
    integer ns = fla_min(m, n);
    integer n_U, m_V, ldu_t = ldu, ldvt_t = ldvt;
    double residual, resid1 = 0., resid2 = 0.;
    double resid3 = 0., resid4 = 0., resid5 = 0.;

    /* Early return conditions */
    if(m == 0 || n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    n_U = (!same_char(*jobz,'A')) ? ns : m;
    m_V = (!same_char(*jobz,'A')) ? ns : n;

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &sigma, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &Usigma, m);

    reset_matrix(datatype, m, n, Usigma, m);

    diagonalize_realtype_vector(datatype, s, sigma, m, n, m);

    /* In case of JOBZ = O, modify ldu, ldvt to make use of the same U buffer
       for further validation similar to JOBZ=A/S */
    if(same_char(*jobz , 'O'))
    {
        if(m >= n)
        {
            ldu = m;
        }
        else
        {
            ldvt = n;
        }
    }
    /* Create temporary buffers(U_temp, V_temp) and copy U,V or A contents
       to use them commonly across all cases JOBZ=A/S/O */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &U_temp, ldu);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &V_temp, ldvt);
    if(same_char(*jobz , 'A') || same_char(*jobz , 'S'))
    {
        copy_matrix(datatype, "FULL", m, n_U, U, ldu_t, U_temp, ldu);
        copy_matrix(datatype, "FULL", m_V, n, V, ldvt_t, V_temp, ldvt);
    }
    /* if JOBZ = 'O',
       if M >= N, A is overwritten with the first N columns
                  of U (the left singular vectors, stored columnwise);
       if M < N, A is overwritten with the first M rows
                 of V**T (the right singular vectors, stored rowwise) */
    else if(same_char(*jobz , 'O'))
    {
        if(m >= n)
        {
            copy_matrix(datatype, "FULL", m, n, A_test, lda, U_temp, ldu);
            copy_matrix(datatype, "FULL", m_V, n, V, ldvt_t, V_temp, ldvt);
        }
        else
        {
            copy_matrix(datatype, "FULL", m, n_U, U, ldu_t, U_temp, ldu);
            copy_matrix(datatype, "FULL", m, n, A_test, lda, V_temp, ldvt);
        }
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps;
            norm = norm_A = 0.f;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(!same_char(*jobz,'N'))
            {
                norm_A = fla_lapack_slange("1", &m, &n, A, &lda, work);
                sgemm_("N", "N", &m, &n, &n_U, &s_one, U_temp, &ldu, sigma, &m, &s_zero, Usigma,
                       &m);
                sgemm_("N", "N", &m, &n, &m_V, &s_one, Usigma, &m, V_temp, &ldvt, &s_n_one, A,
                       &lda);
                norm = fla_lapack_slange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)fla_max(m, n));

                /* Test 2
                   compute norm(I - U'*U) / (N * EPS)*/
                resid2 = (float)check_orthogonal_matrix('T', datatype, U_temp, ns, m, ns, ldu);
                /* Test 3
                   compute norm(I - V*V') / (N * EPS)*/
                resid3 = (float)check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);
            }

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, err_thresh);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((same_char(imatrix, 'O') || same_char(imatrix, 'U')) && (scal != NULL))
                {
                    sscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps;
            norm = norm_A = 0.;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(!same_char(*jobz, 'N'))
            {
                norm_A = fla_lapack_dlange("1", &m, &n, A, &lda, work);
                dgemm_("N", "N", &m, &n, &n_U, &d_one, U_temp, &ldu, sigma, &m, &d_zero, Usigma,
                       &m);
                dgemm_("N", "N", &m, &n, &m_V, &d_one, Usigma, &m, V_temp, &ldvt, &d_n_one, A,
                       &lda);
                norm = fla_lapack_dlange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (double)fla_max(m, n));

                /* Test 2
                   compute norm(I - U'*U) / (N * EPS)*/
                resid2 = check_orthogonal_matrix('T', datatype, U_temp, ns, m, ns, ldu);
                /* Test 3
                   compute norm(I - V*V') / (N * EPS)*/
                resid3 = check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);
            }

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, err_thresh);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((same_char(imatrix, 'O') || same_char(imatrix, 'U')) && (scal != NULL))
                {
                    dscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps;
            norm = norm_A = 0.f;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(!same_char(*jobz,'N'))
            {
                norm_A = fla_lapack_clange("1", &m, &n, A, &lda, work);
                cgemm_("N", "N", &m, &n, &n_U, &c_one, U_temp, &ldu, sigma, &m, &c_zero, Usigma,
                       &m);
                cgemm_("N", "N", &m, &n, &m_V, &c_one, Usigma, &m, V_temp, &ldvt, &c_n_one, A,
                       &lda);
                norm = fla_lapack_clange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)fla_max(m, n));

                /* Test 2
                   compute norm(I - U'*U) / (N * EPS)*/
                resid2 = (float)check_orthogonal_matrix('C', datatype, U_temp, ns, m, ns, ldu);

                /* Test 3
                   compute norm(I - V*V') / (N * EPS)*/
                resid3 = (float)check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);
            }

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, err_thresh);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((same_char(imatrix, 'O') || same_char(imatrix, 'U')) && (scal != NULL))
                {
                    sscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps;
            norm = norm_A = 0.;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(!same_char(*jobz,'N'))
            {
                norm_A = fla_lapack_zlange("1", &m, &n, A, &lda, work);
                zgemm_("N", "N", &m, &n, &n_U, &z_one, U_temp, &ldu, sigma, &m, &z_zero, Usigma,
                       &m);
                zgemm_("N", "N", &m, &n, &m_V, &z_one, Usigma, &m, V_temp, &ldvt, &z_n_one, A,
                       &lda);
                norm = fla_lapack_zlange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (double)fla_max(m, n));

                /* Test 2
                   compute norm(I - U'*U) / (N * EPS)*/
                resid2 = check_orthogonal_matrix('C', datatype, U_temp, ns, m, ns, ldu);
                /* Test 3
                   compute norm(I - V*V') / (N * EPS)*/
                resid3 = check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);
            }

            /* Test 4
               Test to Check order of Singular values of SVD  (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, err_thresh);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((same_char(imatrix, 'O') || same_char(imatrix, 'U')) && (scal != NULL))
                {
                    dscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            break;
        }
    }
    free_matrix(U_temp);
    free_matrix(V_temp);
    free_matrix(sigma);
    free_matrix(Usigma);

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);
    residual = fla_test_max(resid4, residual);
    residual = fla_test_max(resid5, residual);

    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
    FLA_PRINT_SUBTEST_STATUS(resid4, err_thresh, "04");
    FLA_PRINT_SUBTEST_STATUS(resid5, err_thresh, "05");
}
