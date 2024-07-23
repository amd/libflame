/******************************************************************************
 * Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_getrf.c
 *  @brief Defines validate function of GETRF() to use in test suite.
 *  */

#include "test_common.h"

void validate_getrf(integer m_A, integer n_A, void *A, void *A_test, /*AFACT*/
                    integer lda, integer *IPIV, integer datatype, double *residual, integer *info,
                    char imatrix)
{
    if(m_A == 0 || n_A == 0)
        return;
    /* System generated locals */
    integer m_n_vector, min_A;
    integer m_L, n_L, m_U, n_U, k;
    void *L, *U, *T, *work, *X, *B, *A_save;
    integer nrhs = 1;
    *info = 0;

    m_n_vector = m_A * n_A;
    min_A = fla_min(m_A, n_A);
    create_vector(datatype, &X, m_A);
    create_vector(datatype, &B, m_A);
    if(m_A > n_A)
    {
        m_L = m_A;
        n_L = m_U = n_U = n_A;
        k = n_A;
    }
    else
    {
        m_L = n_L = m_U = m_A;
        n_U = n_A;
        k = m_A;
    }
    create_matrix(datatype, matrix_layout, m_L, n_L, &L, m_L);
    create_matrix(datatype, matrix_layout, m_U, n_U, &U, m_L);
    create_matrix(datatype, matrix_layout, m_A, n_A, &T, m_A);
    // Create matrix of A(MxN) for subtracting in Test 2 (T-A).
    create_matrix(datatype, matrix_layout, m_A, n_A, &A_save, m_A);

    reset_matrix(datatype, m_U, n_U, U, m_U);
    create_vector(datatype, &work, 2 * m_A);

    rand_vector(datatype, m_A, X, 1, d_zero, d_zero, 'R');

    /* Lower triangular matrix should be sqare matrix m x m */
    /* For m==i OR  m < n OR m > n -->  A(mxn) = L(mxm) * U(mxn) */
    copy_matrix(datatype, "Lower", m_L, n_L, A_test, lda, L, m_L);
    copy_matrix(datatype, "Upper", m_U, n_U, A_test, lda, U, m_U);

    copy_matrix(datatype, "Full", m_A, n_A, A, lda, A_save, m_A);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm = 0, norm_A = 0, norm_X = 0, eps, resid1, resid2;
            eps = fla_lapack_slamch("Epsilon");
            /* Test 1 */
            if(m_A == n_A)
            {
                norm_X = snrm2_(&m_A, X, &i_one);
                /* B = A * X */
                sgemv_("N", &m_A, &m_A, &s_one, A, &lda, X, &i_one, &s_zero, B, &i_one);
                /* Compute X' by passing A_test and B */
                fla_lapack_sgetrs("N", &m_A, &nrhs, A_test, &lda, IPIV, B, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute X - X' */
                saxpy_(&m_A, &s_n_one, B, &i_one, X, &i_one);
                norm = snrm2_(&m_A, X, &i_one);
                resid1 = (norm / norm_X) / (m_A * eps);
            }
            else
            {
                resid1 = 0.0;
            }

            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            fla_lapack_slaset("U", &m_L, &n_L, &s_zero, &s_one, L, &m_L);
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    float *vector = (float *)A + i * lda;
                    norm_A = fla_max(norm_A, snrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm_A = fla_lapack_slange("F", &m_A, &n_A, A, &lda, work);
            }

            /* T = L * U  */
            sgemm_("N", "N", &m_A, &n_A, &k, &s_one, L, &m_L, U, &m_U, &s_zero, T, &m_A);
            /*  Row interchanges based on IPIV values */
            fla_lapack_slaswp(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            saxpy_(&m_n_vector, &s_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    float *vector = (float *)T + i * m_A;
                    norm = fla_max(norm, snrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm = fla_lapack_slange("F", &m_A, &n_A, T, &m_A, work);
            }
            resid2 = (norm / norm_A) / (n_A * eps);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double norm = 0, norm_A = 0, norm_X, eps, resid1, resid2;

            eps = fla_lapack_dlamch("Epsilon");
            /* Test 1 */
            if(m_A == n_A)
            {
                norm_X = dnrm2_(&m_A, X, &i_one);
                /* B = A * X */
                dgemv_("N", &m_A, &m_A, &d_one, A, &lda, X, &i_one, &d_zero, B, &i_one);
                /* Compute X' by passing A and X */
                fla_lapack_dgetrs("N", &m_A, &nrhs, A_test, &lda, IPIV, B, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute X - X' */
                daxpy_(&m_A, &d_n_one, B, &i_one, X, &i_one);
                norm = dnrm2_(&m_A, X, &i_one);
                resid1 = (norm / norm_X) / (m_A * eps);
            }
            else
            {
                resid1 = 0.0;
            }

            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            fla_lapack_dlaset("U", &m_L, &n_L, &d_zero, &d_one, L, &m_L);
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    double *vector = (double *)A + i * lda;
                    norm_A = fla_max(norm_A, dnrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm_A = fla_lapack_dlange("F", &m_A, &n_A, A, &lda, work);
            }
            /* T = L * U  */
            dgemm_("N", "N", &m_A, &n_A, &k, &d_one, L, &m_L, U, &m_U, &d_zero, T, &m_A);
            /*  Row interchanges based on IPIV values*/
            fla_lapack_dlaswp(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            daxpy_(&m_n_vector, &d_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    double *vector = (double *)T + i * m_A;
                    norm = fla_max(norm, dnrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm = fla_lapack_dlange("F", &m_A, &n_A, T, &m_A, work);
            }
            resid2 = (norm / norm_A) / (n_A * eps);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm = 0, norm_A = 0, norm_X = 0, eps, resid1, resid2;

            eps = fla_lapack_slamch("Epsilon");
            /* Test 1 */
            if(m_A == n_A)
            {
                norm_X = scnrm2_(&m_A, X, &i_one);
                /* B = A * X */
                cgemv_("N", &m_A, &m_A, &c_one, A, &lda, X, &i_one, &c_zero, B, &i_one);
                /* Compute X' by passing A and X */
                fla_lapack_cgetrs("N", &m_A, &nrhs, A_test, &lda, IPIV, B, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute X - X' */
                caxpy_(&m_A, &c_n_one, B, &i_one, X, &i_one);
                norm = scnrm2_(&m_A, X, &i_one);
                resid1 = (norm / norm_X) / (m_A * eps);
            }
            else
            {
                resid1 = 0.0;
            }

            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            fla_lapack_claset("U", &m_L, &n_L, &c_zero, &c_one, L, &m_L);
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    scomplex *vector = (scomplex *)A + i * lda;
                    norm_A = fla_max(norm_A, scnrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm_A = fla_lapack_clange("F", &m_A, &n_A, A, &lda, work);
            }
            /* T = L * U  */
            cgemm_("N", "N", &m_A, &n_A, &k, &c_one, L, &m_L, U, &m_U, &c_zero, T, &m_A);
            /*  Row interchanges based on IPIV values*/
            fla_lapack_claswp(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            caxpy_(&m_n_vector, &c_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    scomplex *vector = (scomplex *)T + i * m_A;
                    norm = fla_max(norm, scnrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm = fla_lapack_clange("F", &m_A, &n_A, T, &m_A, work);
            }
            resid2 = (norm / norm_A) / (n_A * eps);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm = 0, norm_A = 0, norm_X = 0, eps, resid1, resid2;

            eps = fla_lapack_dlamch("Epsilon");
            /* Test 1 */
            if(m_A == n_A)
            {
                norm_X = dznrm2_(&m_A, X, &i_one);
                /* B = A * X */
                zgemv_("N", &m_A, &m_A, &z_one, A, &lda, X, &i_one, &z_zero, B, &i_one);
                /* Compute X' by passing A and X */
                fla_lapack_zgetrs("N", &m_A, &nrhs, A_test, &lda, IPIV, B, &m_A, info);
                if(*info < 0)
                    break;
                /* Compute X - X' */
                zaxpy_(&m_A, &z_n_one, B, &i_one, X, &i_one);
                norm = dznrm2_(&m_A, X, &i_one);
                resid1 = (norm / norm_X) / (m_A * eps);
            }
            else
            {
                resid1 = 0.0;
            }

            /* Test 2 */
            /* Unity diagonal elements to Lower triangular matrix */
            fla_lapack_zlaset("U", &m_L, &n_L, &z_zero, &z_one, L, &m_L);
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    dcomplex *vector = (dcomplex *)A + i * lda;
                    norm_A = fla_max(norm_A, dznrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm_A = fla_lapack_zlange("F", &m_A, &n_A, A, &lda, work);
            }
            /* T = L * U  */
            zgemm_("N", "N", &m_A, &n_A, &k, &z_one, L, &m_L, U, &m_U, &z_zero, T, &m_A);
            /*  Row interchanges based on IPIV values*/
            fla_lapack_zlaswp(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
            /* T - A --> L*U - A */
            zaxpy_(&m_n_vector, &z_n_one, A_save, &i_one, T, &i_one);
            /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
            if(imatrix == 'O')
            {
                for(int i = 0; i < n_A; i++)
                {
                    /* To handle large size values nrm2 is used */
                    dcomplex *vector = (dcomplex *)T + i * m_A;
                    norm = fla_max(norm, dznrm2_(&m_A, vector, &i_one));
                }
            }
            else
            {
                norm = fla_lapack_zlange("F", &m_A, &n_A, T, &m_A, work);
            }
            resid2 = (norm / norm_A) / (n_A * eps);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }

    // Free up buffers
    free_matrix(L);
    free_matrix(U);
    free_matrix(T);
    free_vector(work);
    free_vector(X);
    free_vector(B);
    free_vector(A_save);
}