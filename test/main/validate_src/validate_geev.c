/******************************************************************************
 * Copyright (C) 2022-2024, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file validate_geev.c
 *  @brief Defines validate function of GEEV() to use in test suite.
 *  */

#include "test_common.h"

void validate_geev(char *jobvl, char *jobvr, integer m, void *A, void *A_test, integer lda,
                   void *VL, integer ldvl, void *VR, integer ldvr, void *w, void *wr, void *wi,
                   integer datatype, char imatrix, void *scal, double *residual, integer *info,
                   void *wr_in, void *wi_in)
{
    if(m == 0)
        return;
    void *work = NULL;
    void *lambda = NULL, *Vlambda = NULL;
    *info = 0;
    integer incr = m + 1;

    create_matrix(datatype, matrix_layout, m, m, &lambda, m);
    create_matrix(datatype, matrix_layout, m, m, &Vlambda, m);

    reset_matrix(datatype, m, m, lambda, m);
    reset_matrix(datatype, m, m, Vlambda, m);
    sort_vector(datatype, "A", m, wr_in, 1);

    if(datatype == FLOAT || datatype == DOUBLE)
    {
        create_block_diagonal_matrix(datatype, wr, wi, lambda, m, m, m);
        sort_vector(datatype, "A", m, wr, 1);
        sort_vector(datatype, "A", m, wi, 1);

        add_negative_values(datatype, wi_in, m);
        sort_vector(datatype, "A", m, wi_in, 1);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, norm_W, resid1, resid2, eps;
            norm = norm_A = norm_W = resid1 = resid2 = FLT_MIN;
            eps = fla_lapack_slamch("P");
            if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                sgemm_("N", "N", &m, &m, &m, &s_one, A, &lda, VR, &ldvr, &s_zero, Vlambda, &m);
                 /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 1; i < m; i++)
                    {
                        float *vector = (float *)Vlambda + i * m;
                        norm_A = fla_max(norm_A, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_slange("F", &m, &m, Vlambda, &m, work);
                sgemm_("N", "N", &m, &m, &m, &s_one, VR, &ldvr, lambda, &m, &s_n_one, Vlambda, &m);
                /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 1; i < m; i++)
                    {
                        float *vector = (float *)Vlambda + i * m;
                        norm = fla_max(norm, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_slange("F", &m, &m, Vlambda, &m, work);
                resid1 = norm / (eps * norm_A * (float)m);
            }
            if(*jobvl == 'V')
            {
                /* Test 2
                 * compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)
                 */
                sgemm_("C", "N", &m, &m, &m, &s_one, A, &lda, VL, &ldvl, &s_zero, Vlambda, &m);
                /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 1; i < m; i++)
                    {
                        /* To handle large size values (3.40E+38) nrm2 is used */
                        float *vector = (float *)Vlambda + i * m;
                        norm_A = fla_max(norm_A, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_slange("F", &m, &m, Vlambda, &m, work);
                sgemm_("N", "C", &m, &m, &m, &s_one, VL, &ldvl, lambda, &m, &s_n_one, Vlambda, &m);
                /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 1; i < m; i++)
                    {
                        /* To handle large size values (3.40E+38) nrm2 is used */
                        float *vector = (float *)Vlambda + i * m;
                        norm = fla_max(norm, snrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_slange("F", &m, &m, Vlambda, &m, work);
                resid2 = norm / (eps * norm_A * (float)m);
            }
            *residual = (double)fla_max(resid1, resid2);
            if(wr_in != NULL && wi_in != NULL)
            {
                /* Test 3: In case of specific input generation, compare input and
                           output eigen values */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    *(float *)scal = 1.00 / *(float *)scal;
                    sscal_(&m, scal, wr, &i_one);
                    sscal_(&m, scal, wi, &i_one);
                }
                norm_W = fla_lapack_slange("1", &m, &i_one, wr_in, &i_one, work);
                saxpy_(&m, &s_n_one, wr, &i_one, wr_in, &i_one);
                norm = fla_lapack_slange("1", &m, &i_one, wr_in, &i_one, work);
                resid1 = norm / (eps * norm_W * m);

                norm_W = fla_lapack_slange("1", &m, &i_one, wi_in, &i_one, work);
                saxpy_(&m, &s_n_one, wi, &i_one, wi_in, &i_one);
                norm = fla_lapack_slange("1", &m, &i_one, wi_in, &i_one, work);
                resid2 = norm / (eps * norm_W * m);
                *residual = (double)fla_max(*residual, (double)fla_max(resid1, resid2));
            }
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, norm_W, eps, resid1, resid2;
            norm = norm_A = norm_W = resid1 = resid2 = DBL_MIN;
            eps = fla_lapack_dlamch("P");

            if(*jobvr == 'V')
            {
                /* Test 1
                 * compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)
                 */
                dgemm_("N", "N", &m, &m, &m, &d_one, A, &lda, VR, &ldvr, &d_zero, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        /* To handle large size values (1.79E+308) nrm2 is used */
                        double *vector = (double *)Vlambda + i * m;
                        norm_A = fla_max(norm_A, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_dlange("F", &m, &m, Vlambda, &m, work);
                dgemm_("N", "N", &m, &m, &m, &d_one, VR, &ldvr, lambda, &m, &d_n_one, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        double *vector = (double *)Vlambda + i * m;
                        norm = fla_max(norm, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_dlange("F", &m, &m, Vlambda, &m, work);
                resid1 = norm / (eps * norm_A * (double)m);
            }
            if(*jobvl == 'V')
            {
                /* Test 2
                 * compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)
                 */
                dgemm_("C", "N", &m, &m, &m, &d_one, A, &lda, VL, &ldvl, &d_zero, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        /* To handle large size values (1.79E+308) nrm2 is used */
                        double *vector = (double *)Vlambda + i * m;
                        norm_A = fla_max(norm_A, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_dlange("F", &m, &m, Vlambda, &m, work);
                dgemm_("N", "C", &m, &m, &m, &d_one, VL, &ldvl, lambda, &m, &d_n_one, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        /* To handle large size values (1.79E+308) nrm2 is used */
                        double *vector = (double *)Vlambda + i * m;
                        norm = fla_max(norm, dnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_dlange("F", &m, &m, Vlambda, &m, work);
                resid2 = norm / (eps * norm_A * (double)m);
            }
            *residual = (double)fla_max(resid1, resid2);
            if(wr_in != NULL && wi_in != NULL)
            {
                /* Test 3: In case of specific input generation, compare input and
                           output eigen values */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    *(double *)scal = 1.00 / *(double *)scal;
                    dscal_(&m, scal, wr, &i_one);
                    dscal_(&m, scal, wi, &i_one);
                }
                norm_W = fla_lapack_dlange("1", &m, &i_one, wr_in, &i_one, work);
                daxpy_(&m, &d_n_one, wr, &i_one, wr_in, &i_one);
                norm = fla_lapack_dlange("1", &m, &i_one, wr_in, &i_one, work);
                resid1 = norm / (eps * norm_W * m);

                norm_W = fla_lapack_dlange("1", &m, &i_one, wi_in, &i_one, work);
                daxpy_(&m, &d_n_one, wi, &i_one, wi_in, &i_one);
                norm = fla_lapack_dlange("1", &m, &i_one, wi_in, &i_one, work);
                resid2 = norm / (eps * norm_W * m);
                *residual = fla_max(*residual, fla_max(resid1, resid2));
            }
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, norm_W, eps, resid1, resid2, resid3;
            norm = norm_A = norm_W = resid1 = resid2 = resid3 = FLT_MIN;
            eps = fla_lapack_slamch("P");
            /* Scaleup the output during underflow to avoid
             the very least values during validation*/

            reset_matrix(datatype, m, m, lambda, m);
            ccopy_(&m, w, &i_one, lambda, &incr);

            if(*jobvr == 'V')
            {
                /* Test 1
                 * compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)
                 */
                cgemm_("N", "N", &m, &m, &m, &c_one, A, &lda, VR, &ldvr, &c_zero, Vlambda, &m);
                /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        scomplex *vector = (scomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm_A = fla_max(norm_A, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_clange("F", &m, &m, Vlambda, &m, work);
                cgemm_("N", "N", &m, &m, &m, &c_one, VR, &ldvr, lambda, &m, &c_n_one, Vlambda, &m);
                /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        scomplex *vector = (scomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm = fla_max(norm, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_clange("F", &m, &m, Vlambda, &m, work);
                resid1 = norm / (eps * norm_A * (float)m);
            }
            if(*jobvl == 'V')
            {
                /* Test 2
                 * compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)
                 */
                cgemm_("C", "N", &m, &m, &m, &c_one, A, &lda, VL, &ldvl, &c_zero, Vlambda, &m);
                /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        scomplex *vector = (scomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm_A = fla_max(norm_A, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_clange("F", &m, &m, Vlambda, &m, work);
                cgemm_("N", "C", &m, &m, &m, &c_one, VL, &ldvl, lambda, &m, &c_n_one, Vlambda, &m);
                /* To handle large size values (3.40E+38) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        scomplex *vector = (scomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm = fla_max(norm, scnrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_clange("F", &m, &m, Vlambda, &m, work);
                resid2 = norm / (eps * norm_A * (float)m);
            }
            *residual = (double)fla_max(resid1, resid2);
            if(wr_in != NULL)
            {
                /* Test 3: In case of specific input generation, compare input and
                           output eigen values (A-B = 0) */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    *(float *)scal = 1.00 / *(float *)scal;
                    csscal_(&m, scal, w, &i_one);
                }
                sort_vector(datatype, "A", m, w, 1);
                norm_W = fla_lapack_clange("1", &m, &i_one, wr_in, &i_one, work);
                caxpy_(&m, &c_n_one, w, &i_one, wr_in, &i_one);
                norm = fla_lapack_clange("1", &m, &i_one, wr_in, &i_one, work);
                resid3 = norm / (eps * norm_W * m);
                *residual = (double)fla_max(*residual, resid3);
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, norm_W, eps, resid1, resid2, resid3;
            norm = norm_A = norm_W = resid1 = resid2 = DBL_MIN;
            eps = fla_lapack_dlamch("P");
            /* Scaleup the output during underflow to avoid
             the very least values during validation*/

            reset_matrix(datatype, m, m, lambda, m);
            zcopy_(&m, w, &i_one, lambda, &incr);
            if(*jobvr == 'V')
            {
                /* Test 1
                 * compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)
                 */
                zgemm_("N", "N", &m, &m, &m, &z_one, A, &lda, VR, &ldvr, &z_zero, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        dcomplex *vector = (dcomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm_A = fla_max(norm_A, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_zlange("F", &m, &m, Vlambda, &m, work);
                zgemm_("N", "N", &m, &m, &m, &z_one, VR, &ldvr, lambda, &m, &z_n_one, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        dcomplex *vector = (dcomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm = fla_max(norm, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_zlange("F", &m, &m, Vlambda, &m, work);
                resid1 = norm / (eps * norm_A * (double)m);
            }
            if(*jobvl == 'V')
            {
                /* Test 2
                 * compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)
                 */
                zgemm_("C", "N", &m, &m, &m, &z_one, A, &lda, VL, &ldvl, &z_zero, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        dcomplex *vector = (dcomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm_A = fla_max(norm_A, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm_A = fla_lapack_zlange("F", &m, &m, Vlambda, &m, work);
                zgemm_("N", "C", &m, &m, &m, &z_one, VL, &ldvl, lambda, &m, &z_n_one, Vlambda, &m);
                /* To handle large size values (1.79E+308) nrm2 is used */
                if(imatrix == 'O')
                {
                    for(int i = 0; i < m; i++)
                    {
                        dcomplex *vector = (dcomplex *)Vlambda + i * m;
                        /* To handle large values using snrm2 for norm calc */
                        norm = fla_max(norm, dznrm2_(&m, vector, &i_one));
                    }
                }
                else
                    norm = fla_lapack_zlange("F", &m, &m, Vlambda, &m, work);
                resid2 = norm / (eps * norm_A * (double)m);
            }
            *residual = fla_max(resid1, resid2);
            if(wr_in != NULL)
            {
                /* Test 3: In case of specific input generation, compare input and
                           output eigen values (A-B = 0) */
                if(imatrix == 'O' || imatrix == 'U')
                {
                    *(double *)scal = 1.00 / *(double *)scal;
                    zdscal_(&m, scal, w, &i_one);
                }
                sort_vector(datatype, "A", m, w, 1);
                norm_W = fla_lapack_zlange("1", &m, &i_one, wr_in, &i_one, work);
                zaxpy_(&m, &z_n_one, w, &i_one, wr_in, &i_one);
                norm = fla_lapack_zlange("1", &m, &i_one, wr_in, &i_one, work);
                resid3 = norm / (eps * norm_W * m);
                *residual = fla_max(*residual, resid3);
            }
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(Vlambda);
}