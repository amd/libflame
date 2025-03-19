/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_gesvdx.c
 *  @brief Defines validate function of GESVDX() to use in test suite.
 *  */
#include "test_common.h"

extern double perf;
extern double time_min;

void validate_gesvdx(char *tst_api, char *jobu, char *jobvt, char range, integer m, integer n,
                     void *A, void *A_test, integer lda, void *vl, void *vu, integer il, integer iu,
                     integer ns, void *s, void *s_test, void *U, integer ldu, void *V, integer ldvt,
                     integer datatype, double err_thresh, FILE *g_ext_fptr, void *scal,
                     char imatrix)
{
    void *sigma = NULL, *U_A = NULL;
    void *work = NULL;
    integer min_m_n = fla_min(m, n);
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

    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &sigma, m);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &U_A, m);
    reset_matrix(datatype, m, n, U_A, m);

    diagonalize_realtype_vector(datatype, s, sigma, min_m_n, n, min_m_n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm1, norm2, norm_s, norm_sigma, eps;
            norm1 = norm2 = norm_s = norm_sigma = 0.f;
            eps = fla_lapack_slamch("P");
            if((*jobu == 'V' && *jobvt == 'V') || (*jobu == 'v' && *jobvt == 'v'))
            {
                /* Test 1: Compute (sigma_test - (U'A VT')) */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        float *vector = (float *)sigma + i * ns;
                        norm_sigma = fla_test_max(norm_sigma, snrm2_(&ns, vector, &i_one));
                    }
                }
                else
                    norm_sigma = fla_lapack_slange("1", &ns, &ns, sigma, &min_m_n, work);
                sgemm_("T", "N", &m, &n, &m, &s_one, U, &ldu, A, &lda, &s_zero, U_A, &m);
                sgemm_("N", "T", &min_m_n, &min_m_n, &n, &s_one, U_A, &m, V, &ldvt, &s_n_one, sigma,
                       &min_m_n);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        float *vector = (float *)sigma + i * ns;
                        norm1 = fla_test_max(norm1, snrm2_(&ns, vector, &i_one));
                    }
                }
                else
                    norm1 = fla_lapack_slange("F", &ns, &ns, sigma, &min_m_n, work);
                resid1 = norm1 / (eps * norm_sigma * fla_max(m, n));
            }
            /* In case of g_ext_fptr the svd matrix is not called*/
            if(g_ext_fptr == NULL)
            {
                if((imatrix == 'O' || imatrix == 'U') && (range != 'V' && range != 'v'))
                {
                    *(float *)scal = 1.00 / *(float *)scal;
                    sscal_(&ns, scal, s, &i_one);
                }
                /* Test 2: To check functionality compute (s_test - s) */
                norm_s = fla_lapack_slange("F", &ns, &i_one, s_test, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_test, &i_one);
                norm2 = fla_lapack_slange("F", &ns, &i_one, s_test, &i_one, work);
                resid2 = norm2 / (eps * norm_s * ns);
            }

            /* Test 3: Checking Orthogonal property: U' * U = I */
            if(*jobu == 'V' || *jobu == 'v')
                resid3 = (float)check_orthogonal_matrix('T', datatype, U, ns, m, ns, ldu);

            /* Test 4: Checking Orthogonal property: V * V' = I */
            if(*jobvt == 'V' || *jobvt == 'v')
                resid4 = (float)check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 5: Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid5 = (float)svd_check_order(datatype, s, m, n, err_thresh);
            break;
        }

        case DOUBLE:
        {
            double norm1, norm2, norm_s, norm_sigma, eps;
            norm1 = norm2 = norm_s = norm_sigma = 0.;
            eps = fla_lapack_dlamch("P");
            if((*jobu == 'V' && *jobvt == 'V') || (*jobu == 'v' && *jobvt == 'v'))
            {
                /* Test 1: Compute (sigma_test - (U'A VT')) */

                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        double *vector = (double *)sigma + i * ns;
                        norm_sigma = fla_test_max(norm_sigma, dnrm2_(&ns, vector, &i_one));
                    }
                }
                else
                    norm_sigma = fla_lapack_dlange("F", &ns, &ns, sigma, &min_m_n, work);
                dgemm_("T", "N", &m, &n, &m, &d_one, U, &ldu, A, &lda, &d_zero, U_A, &m);
                dgemm_("N", "T", &min_m_n, &min_m_n, &n, &d_one, U_A, &m, V, &ldvt, &d_n_one, sigma,
                       &min_m_n);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        double *vector = (double *)sigma + i * ns;
                        norm1 = fla_test_max(norm1, dnrm2_(&ns, vector, &i_one));
                    }
                }
                else
                    norm1 = fla_lapack_dlange("F", &ns, &ns, sigma, &min_m_n, work);
                resid1 = norm1 / (eps * norm_sigma * fla_max(m, n));
            }
            if(g_ext_fptr == NULL)
            {
                /* Test 2: To check functionality compute (s_test - s) */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    *(double *)scal = 1.00 / *(double *)scal;
                    dscal_(&ns, scal, s, &i_one);
                }
                norm_s = fla_lapack_dlange("F", &ns, &i_one, s_test, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_test, &i_one);
                norm2 = fla_lapack_dlange("F", &ns, &i_one, s_test, &i_one, work);
                resid2 = norm2 / (eps * norm_s * min_m_n);
            }

            /* Test 3: Checking Orthogonal property: U' * U = I */
            if(*jobu == 'V' || *jobu == 'v')
                resid3 = check_orthogonal_matrix('T', datatype, U, ns, m, ns, ldu);

            /* Test 4: Checking Orthogonal property: V * V' = I */
            if(*jobvt == 'V' || *jobvt == 'v')
                resid4 = check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 5: Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid5 = svd_check_order(datatype, s, m, n, err_thresh);
            break;
        }

        case COMPLEX:
        {
            float norm1, norm2, norm_s, norm_sigma, eps;
            norm1 = norm2 = norm_s = norm_sigma = 0.f;
            eps = fla_lapack_slamch("P");
            if((*jobu == 'V' && *jobvt == 'V') || (*jobu == 'v' && *jobvt == 'v'))
            {
                /* Test 1: Compute (sigma_test - (U'A VT')) */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        scomplex *vector = (scomplex *)sigma + i * ns;
                        norm_sigma = fla_test_max(norm_sigma, scnrm2_(&ns, vector, &i_one));
                    }
                }
                else
                    norm_sigma = fla_lapack_clange("F", &ns, &ns, sigma, &min_m_n, work);
                cgemm_("C", "N", &m, &n, &m, &c_one, U, &ldu, A, &lda, &c_zero, U_A, &m);
                cgemm_("N", "C", &min_m_n, &min_m_n, &n, &c_one, U_A, &m, V, &ldvt, &c_n_one, sigma,
                       &min_m_n);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        scomplex *vector = (scomplex *)sigma + i * ns;
                        norm1 = fla_test_max(norm1, scnrm2_(&ns, vector, &i_one));
                    }
                }
                else
                    norm1 = fla_lapack_clange("1", &ns, &ns, sigma, &min_m_n, work);
                resid1 = norm1 / (eps * norm_sigma * fla_max(m, n));
            }
            if(g_ext_fptr == NULL)
            {
                /* Test 2: To check functionality compute (s_test - s) */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    /* Scaledown the sigma during overflow */
                    *(float *)scal = s_one / *(float *)scal;
                    sscal_(&ns, scal, s, &i_one);
                }
                norm_s = fla_lapack_slange("1", &ns, &i_one, s, &i_one, work);
                saxpy_(&ns, &s_n_one, s, &i_one, s_test, &i_one);
                norm2 = fla_lapack_slange("1", &ns, &i_one, s_test, &i_one, work);
                resid2 = norm2 / (eps * norm_s * min_m_n);
            }

            /* Test 3: Checking Orthogonal property: U' * U = I */
            if(*jobu == 'V' || *jobu == 'v')
                resid3 = (float)check_orthogonal_matrix('C', datatype, U, ns, m, ns, ldu);

            /* Test 4: Checking Orthogonal property: V * V' = I */
            if(*jobvt == 'V' || *jobvt == 'v')
                resid4 = (float)check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 5: Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid5 = (float)svd_check_order(datatype, s, m, n, err_thresh);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm1, norm2, norm_s, norm_sigma, eps;
            norm1 = norm2 = norm_s = norm_sigma = 0.;
            eps = fla_lapack_dlamch("P");
            if((*jobu == 'V' && *jobvt == 'V') || (*jobu == 'v' && *jobvt == 'v'))
            {
                /* Test 1: Compute (sigma_test - (U'A VT')) */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        dcomplex *vector = (dcomplex *)sigma + i * ns;
                        norm_sigma = fla_test_max(norm_sigma, dznrm2_(&ns, vector, &i_one));
                    }
                }
                else
                    norm_sigma = fla_lapack_zlange("F", &ns, &ns, sigma, &min_m_n, work);
                zgemm_("C", "N", &m, &n, &m, &z_one, U, &ldu, A, &lda, &z_zero, U_A, &m);
                zgemm_("N", "C", &min_m_n, &min_m_n, &n, &z_one, U_A, &m, V, &ldvt, &z_n_one, sigma,
                       &min_m_n);
                if(imatrix == 'O')
                {
                    for(int i = 0; i < ns; i++)
                    {
                        /* To handle large size values nrm2 is used */
                        dcomplex *vector = (dcomplex *)sigma + i * ns;
                        norm1 = fla_test_max(norm1, dznrm2_(&ns, vector, &i_one));
                    }
                }
                norm1 = fla_lapack_zlange("F", &ns, &ns, sigma, &min_m_n, work);
                resid1 = norm1 / (eps * norm_sigma * fla_max(m, n));
            }
            if(g_ext_fptr == NULL)
            {

                /* Test 2: To check functionality compute (s_test - s) */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    /* Scaledown the sigma during overflow and underflow */
                    *(double *)scal = d_one / *(double *)scal;
                    dscal_(&ns, scal, s, &i_one);
                }
                norm_s = fla_lapack_dlange("1", &ns, &i_one, s_test, &i_one, work);
                daxpy_(&ns, &d_n_one, s, &i_one, s_test, &i_one);
                norm2 = fla_lapack_dlange("1", &ns, &i_one, s_test, &i_one, work);
                resid2 = norm2 / (eps * norm_s * min_m_n);
            }

            /* Test 3: Checking Orthogonal property: U' * U = I */
            if(*jobu == 'V' || *jobu == 'v')
                resid3 = check_orthogonal_matrix('C', datatype, U, ns, m, ns, ldu);

            /* Test 4: Checking Orthogonal property: V * V' = I */
            if(*jobvt == 'V' || *jobvt == 'v')
                resid4 = check_orthogonal_matrix('N', datatype, V, ns, n, ns, ldvt);

            /* Test 5: Test to Check order of Singular SVD values (positive and non-decreasing) */
            resid5 = svd_check_order(datatype, s, m, n, err_thresh);
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(U_A);

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
