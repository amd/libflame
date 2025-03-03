/******************************************************************************
 * Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_ormqr.c
 *  @brief Defines validate function of ormqr() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

void validate_ormqr(char *tst_api, char side, char trans, integer m, integer n, integer k, void *A,
                    integer lda, void *C, void *Tau, integer ldc, void *C_test, integer datatype,
                    double err_thresh, char imatrix)
{
    double residual = 0.;
    char norm = '1';
    char trans_g;
    integer lwork = -1;
    void *Q = NULL, *CC = NULL, *work = NULL;
    integer m_mod = m;
    integer info = 0;
    /* Early return conditions */
    if(m <= 0 || n <= 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, n, err_thresh);

    /* Modify m based on side */
    if(same_char(side, 'R'))
    {
        m_mod = n;
    }

    /* Allocate memory for Q and CC */
    create_matrix(datatype, LAPACK_COL_MAJOR, m_mod, m_mod, &Q, lda);
    create_matrix(datatype, LAPACK_COL_MAJOR, m, n, &CC, ldc);

    /* Initialize Q and CC to zero */
    reset_matrix(datatype, m_mod, m_mod, Q, lda);
    reset_matrix(datatype, m, n, CC, ldc);

    /* Assign trans_g based on trans */
    if(same_char(trans, 'N'))
    {
        trans_g = 'N';
    }
    else
    {
        trans_g = (datatype == COMPLEX || datatype == DOUBLE_COMPLEX) ? 'C' : 'T';
    }
    /* Copy the first k columns of the factorization to the array Q */
    copy_matrix(datatype, "full", m_mod, k, A, lda, Q, lda);

    /* Copy C to CC */
    copy_matrix(datatype, "full", m, n, C, ldc, CC, ldc);

    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm_C, norm_diff, twork;
            eps = fla_lapack_slamch("P");

            /* Generating the Q martrix using the elementary reflectors and scalar
               factor values */
            fla_lapack_sorgqr(&m_mod, &m_mod, &k, NULL, &lda, NULL, &twork, &lwork, &info);
            if(info < 0)
            {
                residual = DBL_MIN;
                break;
            }

            lwork = twork;
            create_vector(datatype, &work, lwork);

            fla_lapack_sorgqr(&m_mod, &m_mod, &k, Q, &lda, Tau, work, &lwork, &info);

            /* Form Q * C and subtract from output C from ORMQR  */
            if(same_char(side, 'L'))
            {
                sgemm_(&trans_g, "N", &m, &n, &m_mod, &s_n_one, Q, &lda, C_test, &ldc, &s_one, CC,
                       &ldc);
            }
            else
            {
                sgemm_("N", &trans_g, &m, &n, &m_mod, &s_n_one, C_test, &ldc, Q, &lda, &s_one, CC,
                       &ldc);
            }

            /* Compute error in the difference */
            norm_C = fla_lapack_slange(&norm, &m, &n, C, &ldc, work);
            norm_diff = fla_lapack_slange(&norm, &m, &n, CC, &ldc, work);
            residual = (double)((norm_diff / norm_C) / (eps * (float)m));

            break;
        }
        case DOUBLE:
        {
            double eps, norm_C, norm_diff, twork;
            eps = fla_lapack_dlamch("P");

            /* Generating the Q martrix using the elementary reflectors and scalar
               factor values */
            fla_lapack_dorgqr(&m_mod, &m_mod, &k, NULL, &lda, NULL, &twork, &lwork, &info);
            if(info < 0)
            {
                residual = DBL_MIN;
                break;
            }

            lwork = twork;
            create_vector(datatype, &work, lwork);

            fla_lapack_dorgqr(&m_mod, &m_mod, &k, Q, &lda, Tau, work, &lwork, &info);

            /* Form Q * C and subtract from output C from ORMQR  */
            if(same_char(side, 'L'))
            {
                dgemm_(&trans_g, "N", &m, &n, &m_mod, &d_n_one, Q, &lda, C_test, &ldc, &d_one, CC,
                       &ldc);
            }
            else
            {
                dgemm_("N", &trans_g, &m, &n, &m_mod, &d_n_one, C_test, &ldc, Q, &lda, &d_one, CC,
                       &ldc);
            }

            /* Compute error in the difference */
            norm_C = fla_lapack_dlange(&norm, &m, &n, C, &ldc, work);
            norm_diff = fla_lapack_dlange(&norm, &m, &n, CC, &ldc, work);
            residual = (norm_diff / norm_C) / (eps * (double)m);

            break;
        }
        case COMPLEX:
        {
            float eps, norm_C, norm_diff;
            scomplex twork;
            eps = fla_lapack_slamch("P");

            /* Generating the Q martrix using the elementary reflectors and scalar
               factor values */
            fla_lapack_cungqr(&m_mod, &m_mod, &k, NULL, &lda, NULL, &twork, &lwork, &info);
            if(info < 0)
            {
                residual = DBL_MIN;
                break;
            }

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_cungqr(&m_mod, &m_mod, &k, Q, &lda, Tau, work, &lwork, &info);

            /* Form Q * C and subtract from output C from ORMQR  */
            if(same_char(side, 'L'))
            {
                cgemm_(&trans_g, "N", &m, &n, &m_mod, &c_n_one, Q, &lda, C_test, &ldc, &c_one, CC,
                       &ldc);
            }
            else
            {
                cgemm_("N", &trans_g, &m, &n, &m_mod, &c_n_one, C_test, &ldc, Q, &lda, &c_one, CC,
                       &ldc);
            }

            /* Compute error in the difference */
            norm_C = fla_lapack_clange(&norm, &m, &n, C, &ldc, work);
            norm_diff = fla_lapack_clange(&norm, &m, &n, CC, &ldc, work);
            residual = (double)((norm_diff / norm_C) / (eps * (float)m));

            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm_C, norm_diff;
            dcomplex twork;
            eps = fla_lapack_dlamch("P");

            /* Generating the Q martrix using the elementary reflectors and scalar
               factor values */
            fla_lapack_zungqr(&m_mod, &m_mod, &k, NULL, &lda, NULL, &twork, &lwork, &info);
            if(info < 0)
            {
                residual = DBL_MIN;
                break;
            }

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_zungqr(&m_mod, &m_mod, &k, Q, &lda, Tau, work, &lwork, &info);

            /* Form Q * C and subtract from output C from ORMQR  */
            if(same_char(side, 'L'))
            {
                zgemm_(&trans_g, "N", &m, &n, &m_mod, &z_n_one, Q, &lda, C_test, &ldc, &z_one, CC,
                       &ldc);
            }
            else
            {
                zgemm_("N", &trans_g, &m, &n, &m_mod, &z_n_one, C_test, &ldc, Q, &lda, &z_one, CC,
                       &ldc);
            }

            /* Compute error in the difference */
            norm_C = fla_lapack_zlange(&norm, &m, &n, C, &ldc, work);
            norm_diff = fla_lapack_zlange(&norm, &m, &n, CC, &ldc, work);
            residual = (norm_diff / norm_C) / (eps * (double)m);

            break;
        }
    }

    /* Free allocated memory */
    free_matrix(Q);
    free_matrix(CC);

    FLA_PRINT_TEST_STATUS(m, n, residual, err_thresh);
}