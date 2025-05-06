/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gesvd.c
 *  @brief Defines validate function of GESVD() to use in test suite.
 *  */

#include "test_common.h"

extern double perf;
extern double time_min;

void validate_gesvd(char *tst_api, char *jobu, char *jobvt, integer m, integer n, void *A,
                    void *A_test, integer lda, void *s, void *s_test, void *U, integer ldu, void *V,
                    integer ldvt, integer datatype, double err_thresh, FILE *g_ext_fptr,
                    char imatrix, void *scal)
{
    void *sigma = NULL, *Usigma = NULL;
    void *work = NULL, *U_temp = NULL, *V_temp = NULL;
    integer n_U, m_V, ns = fla_min(m, n), ldu_t = ldu, ldvt_t = ldvt;
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

    n_U = (*jobu != 'A') ? ns : m;
    m_V = (*jobvt != 'A') ? ns : n;

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &sigma, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &Usigma, m);
    reset_matrix(datatype, m, n, Usigma, m);

    diagonalize_realtype_vector(datatype, s, sigma, m, n, m);
    /* In case of JOBU/JOBVT = O, modify ldu, ldvt to make use of the same U,V
       buffers for further validation similar to JOBZ=A/S */
    if(*jobu == 'O')
    {
        ldu = m;
    }
    else if(*jobvt == 'O')
    {
        ldvt = n;
    }
    /* Create temporary buffers(U_temp, V_temp) and copy U,V or A contents
       to use them commonly across all cases JOBZ=A/S/O */
    create_matrix(datatype, LAPACK_COL_MAJOR, m, m, &U_temp, ldu);
    create_matrix(datatype, LAPACK_COL_MAJOR, n, n, &V_temp, ldvt);

    if(*jobu == 'A' || *jobu == 'S')
    {
        copy_matrix(datatype, "FULL", m, n_U, U, ldu_t, U_temp, ldu);
    }
    if(*jobvt == 'A' || *jobvt == 'S')
    {
        copy_matrix(datatype, "FULL", m_V, n, V, ldvt_t, V_temp, ldvt);
    }
    /* If jobu or jobvt is 'O' .The first min(m,n) columns/rows of singular vectors
       are overwritten on A output matrix (A_test).*/
    if(*jobu == 'O')
        copy_matrix(datatype, "FULL", m, n_U, A_test, lda, U_temp, ldu);
    else if(*jobvt == 'O')
        copy_matrix(datatype, "FULL", m_V, n, A_test, lda, V_temp, ldvt);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps;
            norm = norm_A = 0.f;
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
                        norm_A = fla_test_max(norm_A, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_slange("F", &m, &n, A, &lda, work);
                sgemm_("N", "N", &m, &n, &n_U, &s_one, U_temp, &ldu, sigma, &m, &s_zero, Usigma,
                       &m);
                sgemm_("N", "N", &m, &n, &m_V, &s_one, Usigma, &m, V_temp, &ldvt, &s_n_one, A,
                       &lda);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        float *vector = (float *)A + i * lda;
                        norm = fla_test_max(norm, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_slange("F", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = (float)check_orthogonal_matrix('T', datatype, U_temp, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = (float)check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, err_thresh);

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
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps;
            norm = norm_A = 0.;
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
                        norm_A = fla_test_max(norm_A, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_dlange("F", &m, &n, A, &lda, work);
                dgemm_("N", "N", &m, &n, &n_U, &d_one, U_temp, &ldu, sigma, &m, &d_zero, Usigma,
                       &m);
                dgemm_("N", "N", &m, &n, &m_V, &d_one, Usigma, &m, V_temp, &ldvt, &d_n_one, A,
                       &lda);
                norm = fla_lapack_dlange("F", &m, &n, A, &lda, work);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        double *vector = (double *)A + i * lda;
                        norm = fla_test_max(norm, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_dlange("F", &m, &n, A, &lda, work);

                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = check_orthogonal_matrix('T', datatype, U_temp, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, err_thresh);

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
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps;
            norm = norm_A = 0.f;

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
                        norm_A = fla_test_max(norm_A, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_clange("F", &m, &n, A, &lda, work);
                cgemm_("N", "N", &m, &n, &n_U, &c_one, U_temp, &ldu, sigma, &m, &c_zero, Usigma,
                       &m);
                cgemm_("N", "N", &m, &n, &m_V, &c_one, Usigma, &m, V_temp, &ldvt, &c_n_one, A,
                       &lda);

                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        scomplex *vector = (scomplex *)A + i * lda;
                        norm = fla_test_max(norm, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_clange("F", &m, &n, A, &lda, work);
                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = (float)check_orthogonal_matrix('C', datatype, U_temp, ns, m, ns, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = (float)check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = (float)svd_check_order(datatype, s, m, n, err_thresh);

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
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps;
            norm = norm_A = 0.;
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
                        norm_A = fla_test_max(norm_A, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_zlange("F", &m, &n, A, &lda, work);
                zgemm_("N", "N", &m, &n, &n_U, &z_one, U_temp, &ldu, sigma, &m, &z_zero, Usigma,
                       &m);
                zgemm_("N", "N", &m, &n, &m_V, &z_one, Usigma, &m, V_temp, &ldvt, &z_n_one, A,
                       &lda);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < n; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        dcomplex *vector = (dcomplex *)A + i * lda;
                        norm = fla_test_max(norm, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_zlange("F", &m, &n, A, &lda, work);

                resid1 = norm / (eps * norm_A * fla_max(m, n));
            }

            /* Test 2
               compute norm(I - U'*U) / (N * EPS)*/
            if(*jobu != 'N')
                resid2 = check_orthogonal_matrix('C', datatype, U_temp, ns, m, ns, ldu);
            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            if(*jobvt != 'N')
                resid3 = check_orthogonal_matrix('N', datatype, V_temp, ns, n, ns, ldvt);

            /* Test 4
               Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid4 = svd_check_order(datatype, s, m, n, err_thresh);

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
