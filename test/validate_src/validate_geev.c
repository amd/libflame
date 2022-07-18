/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_geev.c
 *  @brief Defines validate function of GEEV() to use in test suite.
 *  */

#include "test_common.h"

void validate_geev(char* jobvl, char* jobvr,
    integer m,
    void* A,
    void* A_test,
    void* VL,
    void* VR,
    void* w,
    void* wr,
    void* wi,
    integer datatype,
    double* residual)
{
    void *work = NULL;
    void *lambda = NULL, *Vlambda = NULL;
    integer lda;

    lda = m;

    create_matrix(datatype, &lambda, m, m);
    create_matrix(datatype, &Vlambda, m, m);

    reset_matrix(datatype, m, m, lambda, m);
    reset_matrix(datatype, m, m, Vlambda, m);

    if (datatype == FLOAT || datatype == DOUBLE)
    {
        create_block_diagonal_matrix(datatype, wr, wi, lambda, m, m, lda);
    }

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                sgemm_("N", "N", &m, &m, &m, &s_one, A, &m, VR, &m, &s_zero, Vlambda, &m);
                norm_A = slange_("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "N", &m, &m, &m, &s_one, VR, &m, lambda, &m, &s_n_one, Vlambda, &m);
                norm = slange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);

                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                sgemm_("C", "N", &m, &m, &m, &s_one, A, &m, VL, &m, &s_zero, Vlambda, &m);
                norm_A = slange_("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "C", &m, &m, &m, &s_one, VL, &m, lambda, &m, &s_n_one, Vlambda, &m);
                norm = slange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (float)m);
                *residual = (double)max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                sgemm_("N", "N", &m, &m, &m, &s_one, A, &m, VR, &m, &s_zero, Vlambda, &m);
                norm_A = slange_("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "N", &m, &m, &m, &s_one, VR, &m, lambda, &m, &s_n_one, Vlambda, &m);
                norm = slange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid1;
            }
            else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                sgemm_("C", "N", &m, &m, &m, &s_one, A, &m, VL, &m, &s_zero, Vlambda, &m);
                norm_A = slange_("1", &m, &m, Vlambda, &m, work);
                sgemm_("N", "C", &m, &m, &m, &s_one, VL, &m, lambda, &m, &s_n_one, Vlambda, &m);
                norm = slange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid2;
            }
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = dlamch_("P");

            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                dgemm_("N", "N", &m, &m, &m, &d_one, A, &m, VR, &m, &d_zero, Vlambda, &m);
                norm_A = dlange_("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "N", &m, &m, &m, &d_one, VR, &m, lambda, &m, &d_n_one, Vlambda, &m);
                norm = dlange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);

                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                dgemm_("C", "N", &m, &m, &m, &d_one, A, &m, VL, &m, &d_zero, Vlambda, &m);
                norm_A = dlange_("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "C", &m, &m, &m, &d_one, VL, &m, lambda, &m, &d_n_one, Vlambda, &m);
                norm = dlange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                dgemm_("N", "N", &m, &m, &m, &d_one, A, &m, VR, &m, &d_zero, Vlambda, &m);
                norm_A = dlange_("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "N", &m, &m, &m, &d_one, VR, &m, lambda, &m, &d_n_one, Vlambda, &m);
                norm = dlange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid1;
            }
            else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                dgemm_("C", "N", &m, &m, &m, &d_one, A, &m, VL, &m, &d_zero, Vlambda, &m);
                norm_A = dlange_("1", &m, &m, Vlambda, &m, work);
                dgemm_("N", "C", &m, &m, &m, &d_one, VL, &m, lambda, &m, &d_n_one, Vlambda, &m);
                norm = dlange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid2;
            }
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            integer incr;
            incr = m+1;
            eps = slamch_("P");
            ccopy_(&m, w, &i_one, lambda, &incr);

            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                cgemm_("N", "N", &m, &m, &m, &c_one, A, &m, VR, &m, &c_zero, Vlambda, &m);
                norm_A = clange_("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "N", &m, &m, &m, &c_one, VR, &m, lambda, &m, &c_n_one, Vlambda, &m);
                norm = clange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);

                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                cgemm_("C", "N", &m, &m, &m, &c_one, A, &m, VL, &m, &c_zero, Vlambda, &m);
                norm_A = clange_("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "C", &m, &m, &m, &c_one, VL, &m, lambda, &m, &c_n_one, Vlambda, &m);
                norm = clange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (float)m);
                *residual = (double)max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                cgemm_("N", "N", &m, &m, &m, &c_one, A, &m, VR, &m, &c_zero, Vlambda, &m);
                norm_A = clange_("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "N", &m, &m, &m, &c_one, VR, &m, lambda, &m, &c_n_one, Vlambda, &m);
                norm = clange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid1;
            }
	    else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                cgemm_("C", "N", &m, &m, &m, &c_one, A, &m, VL, &m, &c_zero, Vlambda, &m);
                norm_A = clange_("1", &m, &m, Vlambda, &m, work);
                cgemm_("N", "C", &m, &m, &m, &c_one, VL, &m, lambda, &m, &c_n_one, Vlambda, &m);
                norm = clange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (float)m);
                *residual = (double)resid2;
            }
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            integer incr;
            incr = m+1;
            eps = dlamch_("P");
            zcopy_(&m, w, &i_one, lambda, &incr);

            if(*jobvl == 'V' && *jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                zgemm_("N", "N", &m, &m, &m, &z_one, A, &m, VR, &m, &z_zero, Vlambda, &m);
                norm_A = zlange_("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "N", &m, &m, &m, &z_one, VR, &m, lambda, &m, &z_n_one, Vlambda, &m);
                norm = zlange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);
                /* Test 2
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                zgemm_("C", "N", &m, &m, &m, &z_one, A, &m, VL, &m, &z_zero, Vlambda, &m);
                norm_A = zlange_("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "C", &m, &m, &m, &z_one, VL, &m, lambda, &m, &z_n_one, Vlambda, &m);
                norm = zlange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)max(resid1, resid2);
            }
            else if(*jobvr == 'V')
            {
                /* Test 1
                   compute norm((A*V = V*lambda)) / (V * norm(A) * EPS)*/
                zgemm_("N", "N", &m, &m, &m, &z_one, A, &m, VR, &m, &z_zero, Vlambda, &m);
                norm_A = zlange_("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "N", &m, &m, &m, &z_one, VR, &m, lambda, &m, &z_n_one, Vlambda, &m);
                norm = zlange_("1", &m, &m, Vlambda, &m, work);
                resid1 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid1;
            }
            else if(*jobvl == 'V')
            {
                /* Test 1
                   compute norm (A**H * VL - VL * W**H) / (V * norm(A) * EPS)*/
                zgemm_("C", "N", &m, &m, &m, &z_one, A, &m, VL, &m, &z_zero, Vlambda, &m);
                norm_A = zlange_("1", &m, &m, Vlambda, &m, work);
                zgemm_("N", "C", &m, &m, &m, &z_one, VL, &m, lambda, &m, &z_n_one, Vlambda, &m);
                norm = zlange_("1", &m, &m, Vlambda, &m, work);
                resid2 = norm/(eps * norm_A * (double)m);
                *residual = (double)resid2;
            }
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(Vlambda);
}