/******************************************************************************
 * Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gesdd.c
 *  @brief Defines validate function of GESDD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesdd(char *jobz, integer m, integer n, void *A, void *A_test, integer lda, void *s,
                    void *s_in, void *U, integer ldu, void *V, integer ldvt, integer datatype,
                    double *residual, integer *info, char imatrix, void *scal)
{
    if(m == 0 || n == 0)
        return;
    void *sigma = NULL, *Usigma = NULL;
    void *work = NULL;
    integer ns = fla_min(m, n);
    integer n_U, m_V;

    n_U = (*jobz != 'A') ? ns : m;
    m_V = (*jobz != 'A') ? ns : n;
    *info = 0;
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &sigma, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &Usigma, m);

    reset_matrix(datatype, m, n, Usigma, m);

    diagonalize_realtype_vector(datatype, s, sigma, m, n, m);
    /* if JOBZ = 'O',
       if M >= N, A is overwritten with the first N columns
                  of U (the left singular vectors, stored columnwise);
       if M < N, A is overwritten with the first M rows
                 of V**T (the right singular vectors, stored rowwise) */
    if(*jobz == 'O' && m >= n)
    {
        copy_matrix(datatype, "FULL", m, m, A_test, lda, U, ldu);
    }
    else if(*jobz == 'O' && m < n)
    {
        copy_matrix(datatype, "FULL", m, n, A_test, lda, V, ldvt);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2, resid3, resid4, resid5;
            norm = norm_A = resid1 = resid2 = resid3 = resid4 = resid5 = FLT_MIN;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(*jobz != 'N')
            {
                norm_A = fla_lapack_slange("1", &m, &n, A, &lda, work);
                sgemm_("N", "N", &m, &n, &n_U, &s_one, U, &ldu, sigma, &m, &s_zero, Usigma, &m);
                sgemm_("N", "N", &m, &n, &m_V, &s_one, Usigma, &m, V, &ldvt, &s_n_one, A, &lda);
                norm = fla_lapack_slange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid2 = (float)check_orthogonal_matrix('T', datatype, U, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid3 = (float)check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, *residual);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    sscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            *residual = (double)fla_max(fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4),
                                        resid5);
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2, resid3, resid4, resid5;
            norm = norm_A = resid1 = resid2 = resid3 = resid4 = resid5 = DBL_MIN;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(*jobz != 'N')
            {
                norm_A = fla_lapack_dlange("1", &m, &n, A, &lda, work);
                dgemm_("N", "N", &m, &n, &n_U, &d_one, U, &ldu, sigma, &m, &d_zero, Usigma, &m);
                dgemm_("N", "N", &m, &n, &m_V, &d_one, Usigma, &m, V, &ldvt, &d_n_one, A, &lda);
                norm = fla_lapack_dlange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (double)fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid2 = check_orthogonal_matrix('T', datatype, U, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid3 = check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, *residual);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    dscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            *residual = fla_max(fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4), resid5);
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2, resid3, resid4, resid5;
            norm = norm_A = resid1 = resid2 = resid3 = resid4 = resid5 = FLT_MIN;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(*jobz != 'N')
            {
                norm_A = fla_lapack_clange("1", &m, &n, A, &lda, work);
                cgemm_("N", "N", &m, &n, &n_U, &c_one, U, &ldu, sigma, &m, &c_zero, Usigma, &m);
                cgemm_("N", "N", &m, &n, &m_V, &c_one, Usigma, &m, V, &ldvt, &c_n_one, A, &lda);
                norm = fla_lapack_clange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (float)fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid2 = (float)check_orthogonal_matrix('C', datatype, U, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid3 = (float)check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, *residual);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    sscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_slange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            *residual = (double)fla_max(fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4),
                                        resid5);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2, resid3, resid4, resid5;
            norm = norm_A = resid1 = resid2 = resid3 = resid4 = resid5 = DBL_MIN;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if(*jobz != 'N')
            {
                norm_A = fla_lapack_zlange("1", &m, &n, A, &lda, work);
                zgemm_("N", "N", &m, &n, &n_U, &z_one, U, &ldu, sigma, &m, &z_zero, Usigma, &m);
                zgemm_("N", "N", &m, &n, &m_V, &z_one, Usigma, &m, V, &ldvt, &z_n_one, A, &lda);
                norm = fla_lapack_zlange("1", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * (double)fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid2 = check_orthogonal_matrix('C', datatype, U, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobz == 'A' || *jobz == 'S')
                resid3 = check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular values of SVD  (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, *residual);

            /* Test 5: In case of specific input generation, compare input and
               output singular values */
            if(s_in != NULL)
            {
                if((imatrix == 'O' || imatrix == 'U') && (scal != NULL))
                {
                    dscal_(&ns, scal, s_in, &i_one);
                }
                norm_A = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_in, &i_one);
                norm = fla_lapack_dlange("1", &ns, &i_one, s_in, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            *residual = fla_max(fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4), resid5);
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(Usigma);
}
