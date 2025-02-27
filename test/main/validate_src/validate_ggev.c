/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_ggev.c
 *  @brief Defines validate function of GGEV() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_ggev(char *tst_api, char *jobvl, char *jobvr, integer n, void *A, integer lda,
                   void *B, integer ldb, void *alpha, void *alphar, void *alphai, void *beta,
                   void *VL, integer ldvl, void *VR, integer ldvr, integer datatype,
                   double err_thresh)
{
    integer i, j;
    void *work = NULL;
    double residual;

    /* Early return conditions */
    if(n == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm_A, alphar_t, max_val = 0.0;
            void *YC = NULL, *Y = NULL, *VRTemp = NULL, *VLTemp = NULL, *VRC = NULL, *VLC = NULL;
            scomplex alphac;

            create_vector(datatype, &Y, n);
            create_vector(COMPLEX, &YC, n);
            eps = fla_lapack_slamch("P");
            norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            if(*jobvr == 'V')
            {
                create_vector(COMPLEX, &VRC, n);
                create_vector(datatype, &VRTemp, n);

                for(i = 0; i < n; i++)
                {
                    if(((float *)alphai)[i] != 0.0)
                    {
                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        for(j = 0; j < n; j++)
                        {
                            ((scomplex *)VRC)[j].real = ((float *)VR)[i * ldvr + j];
                            ((scomplex *)VRC)[j].imag = ((float *)VR)[(i + 1) * ldvr + j];
                        }
                        alphac.real = ((float *)beta)[i];
                        alphac.imag = 0;
                        /* beta * A * VR */
                        scgemv('N', 0, n, n, &alphac, A, lda, VRC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float *)alphar)[i];
                        alphac.imag = ((float *)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('N', 0, n, n, &alphac, B, ldb, VRC, i_one, s_n_one, YC, i_one);
                        max_val = fla_test_max(max_val, scnrm2_(&n, YC, &i_one));

                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        alphac.real = ((float *)beta)[i + 1];
                        alphac.imag = 0;
                        for(j = 0; j < n; j++)
                        {
                            ((scomplex *)VRC)[j].real = ((float *)VR)[i * ldvr + j];
                            ((scomplex *)VRC)[j].imag = -((float *)VR)[(i + 1) * ldvr + j];
                        }
                        /* beta * A * VR */
                        scgemv('N', 0, n, n, &alphac, A, lda, VRC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float *)alphar)[i + 1];
                        alphac.imag = ((float *)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('N', 0, n, n, &alphac, B, ldb, VRC, i_one, s_n_one, YC, i_one);
                        i++;
                        max_val = fla_test_max(max_val, scnrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(FLOAT, n, 1, Y, lda);
                        for(j = 0; j < n; j++)
                        {
                            ((float *)VRTemp)[j] = ((float *)VR)[i * ldvr + j];
                        }
                        /* beta*A*vr(i) */
                        alphar_t = ((float *)beta)[i];
                        sgemv_("N", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &s_zero, Y, &i_one);
                        /* alpha*B*vr(i) - beta*A*vr(i) */
                        alphar_t = ((float *)alphar)[i];
                        sgemv_("N", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &s_n_one, Y,
                               &i_one);
                        max_val = fla_test_max(max_val, snrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VRC);
                free_vector(VRTemp);
            }
            if(*jobvl == 'V')
            {
                create_vector(COMPLEX, &VLC, n);
                create_vector(FLOAT, &VLTemp, n);

                for(i = 0; i < n; i++)
                {
                    if(((float *)alphai)[i] != 0.0)
                    {
                        alphac.real = ((float *)beta)[i];
                        alphac.imag = 0;
                        for(j = 0; j < n; j++)
                        {
                            ((scomplex *)VLC)[j].real = ((float *)VL)[i * ldvl + j];
                            ((scomplex *)VLC)[j].imag = ((float *)VL)[(i + 1) * ldvl + j];
                        }
                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        /* beta * A * VL */
                        scgemv('T', 0, n, n, &alphac, A, lda, VLC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float *)alphar)[i];
                        alphac.imag = -((float *)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('T', 0, n, n, &alphac, B, ldb, VLC, i_one, s_n_one, YC, i_one);
                        max_val = fla_test_max(max_val, scnrm2_(&n, YC, &i_one));

                        reset_matrix(COMPLEX, n, 1, YC, lda);
                        alphac.real = ((float *)beta)[i + 1];
                        alphac.imag = 0;
                        for(j = 0; j < n; j++)
                        {
                            ((scomplex *)VLC)[j].real = ((float *)VL)[i * ldvl + j];
                            ((scomplex *)VLC)[j].imag = -((float *)VL)[(i + 1) * ldvl + j];
                        }
                        /* beta * A * VR */
                        scgemv('T', 0, n, n, &alphac, A, lda, VLC, i_one, s_zero, YC, i_one);
                        alphac.real = ((float *)alphar)[i + 1];
                        alphac.imag = -((float *)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        scgemv('T', 0, n, n, &alphac, B, ldb, VLC, i_one, s_n_one, YC, i_one);
                        i++;
                        max_val = fla_test_max(max_val, scnrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(FLOAT, n, 1, Y, lda);
                        for(j = 0; j < n; j++)
                        {
                            ((float *)VLTemp)[j] = ((float *)VL)[i * ldvl + j];
                        }
                        /* beta*A*vl(i) */
                        alphar_t = ((float *)beta)[i];
                        sgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &s_zero, Y, &i_one);
                        /* alpha*B*vl**H(i) - beta*A*vl(i) */
                        alphar_t = ((float *)alphar)[i];
                        sgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &s_n_one, Y,
                               &i_one);
                        max_val = fla_test_max(max_val, snrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VLC);
                free_vector(VLTemp);
            }

            residual = (double)max_val / (eps * norm_A * (float)n);
            free_vector(Y);
            free_vector(YC);
            break;
        }
        case DOUBLE:
        {
            double norm_A, eps, alphar_t, max_val = 0.0;
            void *YC = NULL, *Y = NULL, *VRTemp = NULL, *VLTemp = NULL, *VRC = NULL, *VLC = NULL;
            dcomplex alphac;

            create_vector(datatype, &Y, n);
            create_vector(DOUBLE_COMPLEX, &YC, n);
            eps = fla_lapack_dlamch("P");
            norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            if(*jobvr == 'V')
            {
                create_vector(DOUBLE_COMPLEX, &VRC, n);
                create_vector(datatype, &VRTemp, n);

                for(i = 0; i < n; i++)
                {
                    if(((double *)alphai)[i] != 0.0)
                    {
                        alphac.real = ((double *)beta)[i];
                        alphac.imag = 0;
                        for(j = 0; j < n; j++)
                        {
                            ((dcomplex *)VRC)[j].real = ((double *)VR)[i * ldvr + j];
                            ((dcomplex *)VRC)[j].imag = ((double *)VR)[(i + 1) * ldvr + j];
                        }
                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda);
                        /* beta * A * VR */
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, A, &lda, VRC, &ldvr, &z_zero, YC,
                                &n);
                        alphac.real = ((double *)alphar)[i];
                        alphac.imag = ((double *)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, B, &ldb, VRC, &ldvr, &z_n_one,
                                YC, &n);
                        max_val = fla_test_max(max_val, dznrm2_(&n, YC, &i_one));

                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda);
                        alphac.real = ((double *)beta)[i + 1];
                        alphac.imag = 0;
                        for(j = 0; j < n; j++)
                        {
                            ((dcomplex *)VRC)[j].real = ((double *)VR)[i * ldvr + j];
                            ((dcomplex *)VRC)[j].imag = -((double *)VR)[(i + 1) * ldvr + j];
                        }
                        /* beta * A * VR */
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, A, &lda, VRC, &ldvr, &z_zero, YC,
                                &n);
                        alphac.real = ((double *)alphar)[i + 1];
                        alphac.imag = ((double *)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("N", "N", &n, &i_one, &n, &alphac, B, &ldb, VRC, &ldvr, &z_n_one,
                                YC, &n);
                        i++;
                        max_val = fla_test_max(max_val, dznrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(datatype, n, 1, Y, lda);
                        for(j = 0; j < n; j++)
                        {
                            ((double *)VRTemp)[j] = ((double *)VR)[i * ldvr + j];
                        }
                        /* beta*A*vr(i) */
                        alphar_t = ((double *)beta)[i];
                        dgemv_("N", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &d_zero, Y, &i_one);
                        /* alpha*B*vr(i) - beta*A*vr(i) */
                        alphar_t = ((double *)alphar)[i];
                        dgemv_("N", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &d_n_one, Y,
                               &i_one);
                        max_val = fla_test_max(max_val, dnrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VRC);
                free_vector(VRTemp);
            }
            if(*jobvl == 'V')
            {
                create_vector(DOUBLE_COMPLEX, &VLC, n);
                create_vector(datatype, &VLTemp, n);

                for(i = 0; i < n; i++)
                {
                    if(((double *)alphai)[i] != 0.0)
                    {
                        alphac.real = ((double *)beta)[i];
                        alphac.imag = 0;
                        for(j = 0; j < n; j++)
                        {
                            ((dcomplex *)VLC)[j].real = ((double *)VL)[i * ldvl + j];
                            ((dcomplex *)VLC)[j].imag = ((double *)VL)[(i + 1) * ldvl + j];
                        }
                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda);
                        /* beta * VL**H * A --> beta A**H * VL */

                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, A, &lda, VLC, &ldvl, &z_zero, YC,
                                &n);
                        alphac.real = ((double *)alphar)[i];
                        alphac.imag = -((double *)alphai)[i];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, B, &ldb, VLC, &ldvl, &z_n_one,
                                YC, &n);
                        max_val = fla_test_max(max_val, dznrm2_(&n, YC, &i_one));

                        reset_matrix(DOUBLE_COMPLEX, n, 1, YC, lda);
                        alphac.real = ((double *)beta)[i + 1];
                        alphac.imag = 0;
                        for(j = 0; j < n; j++)
                        {
                            ((dcomplex *)VLC)[j].real = ((double *)VL)[i * ldvl + j];
                            ((dcomplex *)VLC)[j].imag = -((double *)VL)[(i + 1) * ldvl + j];
                        }
                        /* beta * A * VR */
                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, A, &lda, VLC, &ldvl, &z_zero, YC,
                                &n);
                        alphac.real = ((double *)alphar)[i + 1];
                        alphac.imag = -((double *)alphai)[i + 1];
                        /* alpha*B*vr(i) - beta*a*vr(i) */
                        dzgemm_("C", "N", &n, &i_one, &n, &alphac, B, &ldb, VLC, &ldvl, &z_n_one,
                                YC, &n);
                        i++;
                        max_val = fla_test_max(max_val, dznrm2_(&n, YC, &i_one));
                    }
                    else
                    {
                        reset_matrix(datatype, n, 1, Y, lda);
                        for(j = 0; j < n; j++)
                        {
                            ((double *)VLTemp)[j] = ((double *)VL)[i * ldvl + j];
                        }
                        /* beta*A*vl(i) */
                        alphar_t = ((double *)beta)[i];
                        dgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &d_zero, Y, &i_one);
                        /* alpha*B*vl(i) - beta*A*vl(i) */
                        alphar_t = ((double *)alphar)[i];
                        dgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &d_n_one, Y,
                               &i_one);
                        max_val = fla_test_max(max_val, dnrm2_(&n, Y, &i_one));
                    }
                }
                free_vector(VLC);
                free_vector(VLTemp);
            }

            residual = (double)max_val / (eps * norm_A * (double)n);
            free_vector(Y);
            free_vector(YC);
            break;
        }
        case COMPLEX:
        {
            float eps, norm_A, max_val = 0.0;
            void *VRTemp = NULL, *VLTemp = NULL;
            void *Y = NULL;
            scomplex alphar_t;
            eps = fla_lapack_slamch("P");
            norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            create_vector(datatype, &VRTemp, n);
            create_vector(datatype, &VLTemp, n);
            create_vector(datatype, &Y, n);

            if(*jobvr == 'V')
            {
                for(i = 0; i < n; i++)
                {
                    reset_matrix(datatype, n, 1, Y, lda);
                    for(j = 0; j < n; j++)
                    {
                        ((scomplex *)VRTemp)[j] = ((scomplex *)VR)[i * ldvr + j];
                    }
                    /* beta * A * VR */
                    alphar_t = ((scomplex *)beta)[i];
                    cgemv_("N", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &c_zero, Y, &i_one);
                    /* alpha * B * VR - beta * A * VR */
                    alphar_t = ((scomplex *)alpha)[i];
                    cgemv_("N", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &c_n_one, Y, &i_one);
                    max_val = fla_test_max(max_val, scnrm2_(&n, Y, &i_one));
                }
            }

            if(*jobvl == 'V')
            {

                for(i = 0; i < n; i++)
                {
                    reset_matrix(datatype, n, 1, Y, lda);
                    for(j = 0; j < n; j++)
                    {
                        ((scomplex *)VLTemp)[j] = ((scomplex *)VL)[i * ldvl + j];
                    }
                    /* betaC * a**H * vl */
                    ((scomplex *)beta)[i].imag = -((scomplex *)beta)[i].imag;
                    alphar_t = ((scomplex *)beta)[i];
                    cgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &c_zero, Y, &i_one);
                    /* alphaC * b**H * vl - beta * a**H * vr */
                    ((scomplex *)alpha)[i].imag = -((scomplex *)alpha)[i].imag;
                    alphar_t = ((scomplex *)alpha)[i];
                    cgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &c_n_one, Y, &i_one);
                    max_val = fla_test_max(max_val, scnrm2_(&n, Y, &i_one));
                }
            }

            residual = (double)max_val / (eps * norm_A * (float)n);
            free_vector(VLTemp);
            free_vector(Y);
            free_vector(VRTemp);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps, norm_A, max_val = 0.0;
            void *VRTemp = NULL, *VLTemp = NULL;
            void *Y = NULL;
            dcomplex alphar_t;
            eps = fla_lapack_dlamch("P");
            norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
            /* Test 1 */
            /* Validation for 'V' 'V' combination */
            create_vector(datatype, &VRTemp, n);
            create_vector(datatype, &VLTemp, n);
            create_vector(datatype, &Y, n);

            if(*jobvr == 'V')
            {
                for(i = 0; i < n; i++)
                {

                    reset_matrix(datatype, n, 1, Y, lda);
                    for(j = 0; j < n; j++)
                    {
                        ((dcomplex *)VRTemp)[j] = ((dcomplex *)VR)[i * ldvr + j];
                    }
                    /* beta * A * VR */
                    alphar_t = ((dcomplex *)beta)[i];
                    zgemv_("NoTran", &n, &n, &alphar_t, A, &lda, VRTemp, &i_one, &z_zero, Y,
                           &i_one);
                    /* alpha * B * VR - beta * A * VR */
                    alphar_t = ((dcomplex *)alpha)[i];
                    zgemv_("NoTran", &n, &n, &alphar_t, B, &ldb, VRTemp, &i_one, &z_n_one, Y,
                           &i_one);
                    max_val = fla_test_max(max_val, dznrm2_(&n, Y, &i_one));
                }
            }

            if(*jobvl == 'V')
            {

                for(i = 0; i < n; i++)
                {
                    reset_matrix(datatype, n, 1, Y, lda);
                    for(j = 0; j < n; j++)
                    {
                        ((dcomplex *)VLTemp)[j] = ((dcomplex *)VL)[i * ldvl + j];
                    }
                    /* betaC * a**H * vl */
                    ((dcomplex *)beta)[i].imag = -((dcomplex *)beta)[i].imag;
                    alphar_t = ((dcomplex *)beta)[i];
                    zgemv_("C", &n, &n, &alphar_t, A, &lda, VLTemp, &i_one, &z_zero, Y, &i_one);
                    /* alphaC * b**H * vl - beta * a**H * vr */
                    ((dcomplex *)alpha)[i].imag = -((dcomplex *)alpha)[i].imag;
                    alphar_t = ((dcomplex *)alpha)[i];
                    zgemv_("C", &n, &n, &alphar_t, B, &ldb, VLTemp, &i_one, &z_n_one, Y, &i_one);
                    max_val = fla_test_max(max_val, dznrm2_(&n, Y, &i_one));
                }
            }

            residual = (double)max_val / (eps * norm_A * (double)n);
            free_vector(VLTemp);
            free_vector(Y);
            free_vector(VRTemp);
            break;
        }
        default:
            residual = err_thresh;
            break;
    }

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}

/* Test case to check if the given set of alpha, beta and alpha_copy, beta_copy are same or not */
void validate_ggev_EVs(char *tst_api, integer m, void *alpha, void *alphar, void *alphai,
                       void *beta, void *alpha_copy, void *alphar_copy, void *alphai_copy,
                       void *beta_copy, integer datatype, double err_thresh)
{
    double residual;
    double resid1 = 0., resid2 = 0., resid3 = 0.;

    /* Early return conditions */
    if(m == 0)
    {
        FLA_TEST_PRINT_STATUS_AND_RETURN(m, m, err_thresh);
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(m, m, err_thresh);

    void *work = NULL;
    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_a, norm_b, eps;
            eps = fla_lapack_slamch("P");

            norm_a = fla_lapack_slange("1", &m, &i_one, alphar_copy, &i_one, work);
            saxpy_(&m, &s_n_one, alphar, &i_one, alphar_copy, &i_one);
            norm = fla_lapack_slange("1", &m, &i_one, alphar_copy, &i_one, work);
            resid1 = norm / (eps * norm_a * m);

            norm_a = fla_lapack_slange("1", &m, &i_one, alphai_copy, &i_one, work);
            saxpy_(&m, &s_n_one, alphai, &i_one, alphai_copy, &i_one);
            norm = fla_lapack_slange("1", &m, &i_one, alphai_copy, &i_one, work);
            resid2 = norm / (eps * norm_a * m);

            norm_b = fla_lapack_slange("1", &m, &i_one, beta_copy, &i_one, work);
            saxpy_(&m, &s_n_one, beta, &i_one, beta_copy, &i_one);
            norm = fla_lapack_slange("1", &m, &i_one, beta_copy, &i_one, work);
            resid3 = norm / (eps * norm_b * m);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_a, norm_b, eps;
            eps = fla_lapack_dlamch("P");

            norm_a = fla_lapack_dlange("1", &m, &i_one, alphar_copy, &i_one, work);
            daxpy_(&m, &d_n_one, alphar, &i_one, alphar_copy, &i_one);
            norm = fla_lapack_dlange("1", &m, &i_one, alphar_copy, &i_one, work);
            resid1 = norm / (eps * norm_a * m);

            norm_a = fla_lapack_dlange("1", &m, &i_one, alphai_copy, &i_one, work);
            daxpy_(&m, &d_n_one, alphai, &i_one, alphai_copy, &i_one);
            norm = fla_lapack_dlange("1", &m, &i_one, alphai_copy, &i_one, work);
            resid2 = norm / (eps * norm_a * m);

            norm_b = fla_lapack_dlange("1", &m, &i_one, beta_copy, &i_one, work);
            daxpy_(&m, &d_n_one, beta, &i_one, beta_copy, &i_one);
            norm = fla_lapack_dlange("1", &m, &i_one, beta_copy, &i_one, work);
            resid3 = norm / (eps * norm_b * m);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_a, norm_b, eps;
            eps = fla_lapack_slamch("P");

            norm_a = fla_lapack_slange("1", &m, &i_one, alpha_copy, &i_one, work);
            saxpy_(&m, &s_n_one, alpha, &i_one, alpha_copy, &i_one);
            norm = fla_lapack_slange("1", &m, &i_one, alpha_copy, &i_one, work);
            resid1 = norm / (eps * norm_a * m);

            norm_b = fla_lapack_slange("1", &m, &i_one, beta_copy, &i_one, work);
            saxpy_(&m, &s_n_one, beta, &i_one, beta_copy, &i_one);
            norm = fla_lapack_slange("1", &m, &i_one, beta_copy, &i_one, work);
            resid2 = norm / (eps * norm_b * m);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_a, norm_b, eps;
            eps = fla_lapack_dlamch("P");

            norm_a = fla_lapack_dlange("1", &m, &i_one, alpha_copy, &i_one, work);
            daxpy_(&m, &d_n_one, alpha, &i_one, alpha_copy, &i_one);
            norm = fla_lapack_dlange("1", &m, &i_one, alpha_copy, &i_one, work);
            resid1 = norm / (eps * norm_a * m);

            norm_b = fla_lapack_dlange("1", &m, &i_one, beta_copy, &i_one, work);
            daxpy_(&m, &d_n_one, beta, &i_one, beta_copy, &i_one);
            norm = fla_lapack_dlange("1", &m, &i_one, beta_copy, &i_one, work);
            resid2 = norm / (eps * norm_b * m);
            break;
        }
    }

    residual = fla_test_max(resid1, resid2);
    residual = fla_test_max(resid3, residual);

    FLA_PRINT_TEST_STATUS(m, m, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
    FLA_PRINT_SUBTEST_STATUS(resid3, err_thresh, "03");
}
