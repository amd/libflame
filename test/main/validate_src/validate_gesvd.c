/******************************************************************************
 * Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_gesvd.c
 *  @brief Defines validate function of GESVD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesvd(char *jobu, char *jobvt, integer m, integer n, void *A, void *A_test,
                    integer lda, void *s, void *s_test, void *U, integer ldu, void *V, integer ldvt,
                    integer datatype, double *residual, integer *info, FILE *g_ext_fptr,
                    char imatrix, void *scal)
{
    if(m == 0 || n == 0)
        return;
    void *sigma = NULL, *Usigma = NULL;
    void *work = NULL;
    *info = 0;
    integer n_U, m_V, ns = fla_min(m, n);
    n_U = (*jobu != 'A') ? ns : m;
    m_V = (*jobvt != 'A') ? ns : n;

    create_matrix(datatype, &sigma, m, n);
    create_matrix(datatype, &Usigma, m, n);
    reset_matrix(datatype, m, n, Usigma, m);

    diagonalize_realtype_vector(datatype, s, sigma, m, n, m);
    /* If jobu or jobvt is 'O' .The first min(m,n) columns/rows of singular vectors
       are overwritten on A output matrix (A_test).*/
    if(*jobu == 'O')
        copy_matrix(datatype, "FULL", m, fla_min(m, n), A_test, lda, U, ldu);
    else if(*jobvt == 'O')
        copy_matrix(datatype, "FULL", fla_min(m, n), n, A_test, lda, V, ldvt);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2, resid3, resid4, resid5;
            norm = norm_A = resid1 = resid2 = resid3 = resid4 = resid5 = FLT_MIN;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if((*jobu != 'N' && *jobvt != 'N'))
            {
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        float *vector = (float *)A + i * lda;
                        norm_A = fla_max(norm_A, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_slange("F", &m, &n, A, &lda, work);
                sgemm_("N", "N", &m, &n, &n_U, &s_one, U, &ldu, sigma, &m, &s_zero, Usigma, &m);
                sgemm_("N", "N", &m, &n, &m_V, &s_one, Usigma, &m, V, &ldvt, &s_n_one, A, &lda);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        float *vector = (float *)A + i * lda;
                        norm = fla_max(norm, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_slange("F", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = (float)check_orthogonal_matrix('T', datatype, U, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = (float)check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, *residual);

            /*
             * Test 5: To check the functionality: Calculate the difference between
             * the known singular value (sigma 's_test') from input generation and the
             * output (s) from the API
             * */
            if(g_ext_fptr == NULL)
            {
                /* Test 2: To check functionality compute (s_test - s) */

                if(imatrix == 'O' || imatrix == 'U')
                {
                    *(float *)scal = 1.00 / *(float *)scal;
                    sscal_(&ns, scal, s, &i_one);
                }
                norm_A = fla_lapack_slange("F", &ns, &i_one, s_test, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_test, &i_one);
                norm = fla_lapack_slange("F", &ns, &i_one, s_test, &i_one, work);
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
            if((*jobu != 'N' && *jobvt != 'N'))
            {
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        double *vector = (double *)A + i * lda;
                        norm_A = fla_max(norm_A, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_dlange("F", &m, &n, A, &lda, work);
                dgemm_("N", "N", &m, &n, &n_U, &d_one, U, &ldu, sigma, &m, &d_zero, Usigma, &m);
                dgemm_("N", "N", &m, &n, &m_V, &d_one, Usigma, &m, V, &ldvt, &d_n_one, A, &lda);
                norm = fla_lapack_dlange("F", &m, &n, A, &lda, work);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        double *vector = (double *)A + i * lda;
                        norm = fla_max(norm, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_dlange("F", &m, &n, A, &lda, work);

                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = check_orthogonal_matrix('T', datatype, U, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, *residual);

            /*
             * Test 5: To check the functionality: Calculate the difference between
             * the known singular value (sigma 's_test') from input generation and the
             * output (s) from the API
             * */
            if(g_ext_fptr == NULL)
            {
                if(imatrix == 'O' || imatrix == 'U')
                {
                    *(double *)scal = 1.00 / *(double *)scal;
                    dscal_(&ns, scal, s, &i_one);
                }
                norm_A = fla_lapack_dlange("F", &ns, &i_one, s_test, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_test, &i_one);
                norm = fla_lapack_dlange("F", &ns, &i_one, s_test, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            *residual = (double)fla_max(fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4),
                                        resid5);
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2, resid3, resid4, resid5;
            norm = norm_A = resid1 = resid2 = resid3 = resid4 = resid5 = FLT_MIN;

            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            if((*jobu != 'N' && *jobvt != 'N'))
            {
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        scomplex *vector = (scomplex *)A + i * lda;
                        norm_A = fla_max(norm_A, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_clange("F", &m, &n, A, &lda, work);
                cgemm_("N", "N", &m, &n, &n_U, &c_one, U, &ldu, sigma, &m, &c_zero, Usigma, &m);
                cgemm_("N", "N", &m, &n, &m_V, &c_one, Usigma, &m, V, &ldvt, &c_n_one, A, &lda);

                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        scomplex *vector = (scomplex *)A + i * lda;
                        norm = fla_max(norm, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_clange("F", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = (float)check_orthogonal_matrix('C', datatype, U, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = (float)check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, *residual);

            /*
             * Test 5: To check the functionality: Calculate the difference between
             * the known singular value (sigma 's_test') from input generation and the
             * output (s) from the API
             * */
            if(g_ext_fptr == NULL)
            {
                /* To do : Validate for input generation for overflow and underflow. */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    /* Scaledown the sigma during overflow */
                    *(float *)scal = s_one / *(float *)scal;
                    sscal_(&ns, scal, s, &i_one);
                }
                norm_A = fla_lapack_slange("F", &ns, &i_one, s_test, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_test, &i_one);
                norm = fla_lapack_slange("F", &ns, &i_one, s_test, &i_one, work);
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
            if((*jobu != 'N' && *jobvt != 'N'))
            {
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        dcomplex *vector = (dcomplex *)A + i * lda;
                        norm_A = fla_max(norm_A, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_zlange("F", &m, &n, A, &lda, work);
                zgemm_("N", "N", &m, &n, &n_U, &z_one, U, &ldu, sigma, &m, &z_zero, Usigma, &m);
                zgemm_("N", "N", &m, &n, &m_V, &z_one, Usigma, &m, V, &ldvt, &z_n_one, A, &lda);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        dcomplex *vector = (dcomplex *)A + i * lda;
                        norm = fla_max(norm, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_zlange("F", &m, &n, A, &lda, work);

                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = check_orthogonal_matrix('C', datatype, U, ns, m, ns, ldu);
            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, *residual);

            /*
             * Test 5: To check the functionality: Calculate the difference between
             * the known singular value (sigma 's_test') from input generation and the
             * output (s) from the API
             * */
            if(g_ext_fptr == NULL)
            {
                /*To do : Validate for input generation for overflow and underflow.*/

                if(imatrix == 'O' || imatrix == 'U')
                {
                    /* Scaledown the sigma during overflow */
                    *(double *)scal = d_one / *(double *)scal;
                    dscal_(&ns, scal, s, &i_one);
                }
                norm_A = fla_lapack_dlange("F", &ns, &i_one, s_test, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_test, &i_one);
                norm = fla_lapack_dlange("F", &ns, &i_one, s_test, &i_one, work);
                resid5 = norm / (eps * norm_A * ns);
            }
            *residual = (double)fla_max(fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4),
                                        resid5);
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(Usigma);
}
