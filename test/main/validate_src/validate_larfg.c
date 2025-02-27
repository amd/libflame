/*
    Copyright (C) 2024-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_larfg.c
 *  @brief Defines validate function of LARFG() to use in test suite.
 *  */
#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_larfg(char *tst_api, integer datatype, integer n, integer incx, integer x_length,
                    void *x, void *v, void *tau, double err_thresh)
{
    void *beta, *x_temp, *v_temp, *work;
    double residual;

    /* Early return for n <= 0, tau = 0.0f */
    if(n <= 0)
    {
        residual = is_value_zero(datatype, tau, err_thresh);
        FLA_TEST_PRINT_STATUS_AND_RETURN(n, n, err_thresh);
    }
    /* return invalid param for incx < 1 */
    if(incx < 1)
    {
        residual = DBL_MIN;
        return;
    }
    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(n, n, err_thresh);

    create_vector(datatype, &work, i_one);
    create_vector(datatype, &x_temp, n);
    create_vector(datatype, &v_temp, n);
    create_vector(datatype, &beta, n);

    reset_vector(datatype, beta, n, 1);
    /* x_temp consists of [alpha(input), x(input)] elements */
    copy_vector(datatype, n, x, incx, x_temp, 1);
    /* First element of V is the beta value */
    copy_vector(datatype, 1, v, 1, beta, 1);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_beta, eps;
            ((float *)v)[0] = s_one;
            /* v_temp consists of [1, v(output)] elements */
            copy_vector(datatype, n, v, incx, v_temp, 1);

            /* Test1 : Compute norm(beta - (H * x_temp)) / (EPS * norm(beta) * n) */
            eps = fla_lapack_slamch("P");
            norm_beta = fla_lapack_slange("1", &n, &i_one, beta, &i_one, work);

            /* By using larf API we apply elemntary reflector H to x_temp
             * H = I - tau * v_temp * v_temp**T
             * v_temp is the representation of H
             * x_temp is overwriiten by H * x_temp , it returns x_temp as (beta, 0,0,...n)
             */
            slarf_("L", &n, &i_one, v_temp, &i_one, tau, x_temp, &n, work);
            saxpy_(&n, &s_n_one, x_temp, &i_one, beta, &i_one);
            norm = fla_lapack_slange("1", &n, &i_one, beta, &i_one, work);
            residual = (double)(norm / (eps * norm_beta * n));
            break;
        }
        case DOUBLE:
        {
            double norm, norm_beta, eps;
            ((double *)v)[0] = d_one;
            /* v_temp consists of [1, v(output)] elements */
            copy_vector(datatype, n, v, incx, v_temp, 1);

            /* Test1 : Compute norm(beta - (H * x_temp)) / (EPS * norm(beta) * n) */
            eps = fla_lapack_dlamch("P");
            norm_beta = fla_lapack_dlange("1", &n, &i_one, beta, &i_one, work);

            /* By using larf API we apply elemntary reflector H to x_temp
             * H = I - tau * v_temp * v_temp**T
             * v_temp is the representation of H
             * x_temp is overwriiten by H * x_temp , it returns x_temp as (beta, 0,0,...n)
             */
            dlarf_("L", &n, &i_one, v_temp, &i_one, tau, x_temp, &n, work);
            daxpy_(&n, &d_n_one, x_temp, &i_one, beta, &i_one);
            norm = fla_lapack_dlange("1", &n, &i_one, beta, &i_one, work);
            residual = norm / (eps * norm_beta * n);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_beta, eps;
            ((scomplex *)v)[0] = c_one;
            /* v_temp consists of [1, v(output)] elements */
            copy_vector(datatype, n, v, incx, v_temp, 1);

            /* Get Conjugate of tau */
            ((scomplex *)tau)[0].imag = (-1.0) * ((scomplex *)tau)[0].imag;

            /* Test1 : Compute norm(beta - (H * x_temp)) / (EPS * norm(beta) * n) */
            eps = fla_lapack_slamch("P");
            norm_beta = fla_lapack_clange("1", &n, &i_one, beta, &i_one, work);

            /* By using larf API we apply elemntary reflector H to x_temp
             * H = I - tau * v_temp * v_temp**H
             * v_temp is the representation of H
             * x_temp is overwriiten by H * x_temp , it returns x_temp as (beta, 0,0,...n)
             */
            clarf_("L", &n, &i_one, v_temp, &i_one, tau, x_temp, &n, work);
            caxpy_(&n, &c_n_one, x_temp, &i_one, beta, &i_one);
            norm = fla_lapack_clange("1", &n, &i_one, beta, &i_one, work);
            residual = (double)(norm / (eps * norm_beta * n));
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_beta, eps;
            ((dcomplex *)v)[0] = z_one;
            /* v_temp consists of [1, v(output)] elements */
            copy_vector(datatype, n, v, incx, v_temp, 1);

            /* Get conjugate of tau */
            ((dcomplex *)tau)[0].imag = (-1.0) * ((dcomplex *)tau)[0].imag;

            /* Test1 : Compute norm(beta - (H_matrix * x_temp)) / (EPS * norm(Beta) * n) */
            eps = fla_lapack_dlamch("P");
            norm_beta = fla_lapack_zlange("1", &n, &i_one, beta, &i_one, work);

            /* By using larf API we apply elemntary reflector H to x_temp
             * H = I - tau * v_temp * v_temp**H
             * v_temp is the representation of H
             * x_temp is overwriiten by H * x_temp , it returns x_temp as (beta, 0,0,...n)
             */
            zlarf_("L", &n, &i_one, v_temp, &i_one, tau, x_temp, &n, work);
            zaxpy_(&n, &z_n_one, x_temp, &i_one, beta, &i_one);
            norm = fla_lapack_zlange("1", &n, &i_one, beta, &i_one, work);
            residual = norm / (eps * norm_beta * n);
            break;
        }
        default:
            residual = err_thresh;
            break;
    }
    free_vector(beta);
    free_vector(v_temp);
    free_vector(x_temp);
    free_vector(work);

    FLA_PRINT_TEST_STATUS(n, n, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(residual, err_thresh, "01");
}
