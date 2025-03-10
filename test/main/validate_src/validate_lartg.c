/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_rot.c
 *  @brief Defines validate function of LARTG() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_lartg(char *tst_api, integer datatype, void *f, void *g, void *r, void *c, void *s,
                    double err_thresh)
{
    void *out_zero = NULL;
    double residual, resid1 = 0., resid2 = 0.;

    /* print overall status if incoming threshold is
     * an extreme value indicating that API returned
     * unexpected info value */
    FLA_TEST_PRINT_INVALID_STATUS(2, 1, err_thresh);

    create_vector(datatype, &out_zero, 1);
    reset_vector(datatype, out_zero, 1, 1);

    switch(datatype)
    {
        case FLOAT:
        {
            float eps, norm_r, norm_f, norm_res, norm_1;
            float res = 0.0;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = snrm2_(&i_one, &s_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((float *)c)[0] * ((float *)c)[0]) + (((float *)s)[0] * ((float *)s)[0]));
            /*res->res-1*/
            saxpy_(&i_one, &s_n_one, &s_one, &i_one, &res, &i_one);
            norm_res = snrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res / norm_1 / eps);

            /*Test 2 Validating R*/
            norm_r = snrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((float *)s)[0] = -((float *)s)[0];
            /*Rotating output vectors*/
            fla_lapack_srot(&i_one, r, &i_one, out_zero, &i_one, ((float *)c), ((float *)s));

            /*r-f*/
            saxpy_(&i_one, &s_n_one, f, &i_one, r, &i_one);
            saxpy_(&i_one, &s_n_one, g, &i_one, out_zero, &i_one);
            norm_f = snrm2_(&i_one, r, &i_one);
            resid2 = (norm_f / norm_r / eps);
            break;
        }

        case DOUBLE:
        {
            double eps, norm_r, norm_f, norm_res, norm_1;
            double res = 0.0;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = dnrm2_(&i_one, &d_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((double *)c)[0] * ((double *)c)[0])
                       + (((double *)s)[0] * ((double *)s)[0]));
            /*res->res -1*/
            daxpy_(&i_one, &d_n_one, &d_one, &i_one, &res, &i_one);
            norm_res = dnrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res / norm_1 / eps);

            /*Test 2 Validating R*/
            norm_r = dnrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((double *)s)[0] = -((double *)s)[0];
            /*Rotating output vectors*/
            fla_lapack_drot(&i_one, r, &i_one, out_zero, &i_one, ((double *)c), ((double *)s));

            /*r-f*/
            daxpy_(&i_one, &d_n_one, f, &i_one, r, &i_one);
            daxpy_(&i_one, &d_n_one, g, &i_one, out_zero, &i_one);
            norm_f = dnrm2_(&i_one, r, &i_one);
            resid2 = (norm_f / norm_r / eps);
            break;
        }

        case COMPLEX:
        {
            float eps, norm_r, norm_f, norm_res, norm_1;
            float res = 0.0;
            ;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = snrm2_(&i_one, &s_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((float *)c)[0] * ((float *)c)[0])
                       + (((scomplex *)s)[0].real * ((scomplex *)s)[0].real)
                       + (((scomplex *)s)[0].imag * ((scomplex *)s)[0].imag));
            /*res->res - 1*/
            saxpy_(&i_one, &s_n_one, &s_one, &i_one, &res, &i_one);
            norm_res = snrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res / norm_1 / eps);

            /*Test 2 Validating R*/
            norm_r = scnrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((scomplex *)s)[0].real = -((scomplex *)s)[0].real;
            ((scomplex *)s)[0].imag = -((scomplex *)s)[0].imag;
            /*Rotating output vectors*/
            fla_lapack_crot(&i_one, r, &i_one, out_zero, &i_one, ((float *)c), ((scomplex *)s));

            /*r-f*/
            caxpy_(&i_one, &c_n_one, f, &i_one, r, &i_one);
            caxpy_(&i_one, &c_n_one, g, &i_one, out_zero, &i_one);
            norm_f = scnrm2_(&i_one, r, &i_one);
            resid2 = (norm_f / norm_r / eps);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double eps, norm_r, norm_f, norm_res, norm_1;
            double res = 0.0;

            eps = fla_lapack_slamch("P");

            /*Test 1 validating C and S */
            norm_1 = snrm2_(&i_one, &s_one, &i_one);
            /*square root of ((c*c) + (s*s)) = 1*/
            res = sqrt((((double *)c)[0] * ((double *)c)[0])
                       + (((dcomplex *)s)[0].real * ((dcomplex *)s)[0].real)
                       + (((dcomplex *)s)[0].imag * ((dcomplex *)s)[0].imag));
            /*res-> res - 1*/
            daxpy_(&i_one, &d_n_one, &d_one, &i_one, &res, &i_one);
            norm_res = dnrm2_(&i_one, &res, &i_one);
            resid1 = (norm_res / norm_1 / eps);

            /*Test 2 Validating R*/
            norm_r = dznrm2_(&i_one, r, &i_one);
            /*Hermitian Transpose of original rotation vector*/
            ((dcomplex *)s)[0].real = -((dcomplex *)s)[0].real;
            ((dcomplex *)s)[0].imag = -((dcomplex *)s)[0].imag;
            /*Rotating output vectors*/
            fla_lapack_zrot(&i_one, r, &i_one, out_zero, &i_one, ((double *)c), ((dcomplex *)s));

            /*r-f*/
            zaxpy_(&i_one, &z_n_one, f, &i_one, r, &i_one);
            zaxpy_(&i_one, &z_n_one, g, &i_one, out_zero, &i_one);
            norm_f = dznrm2_(&i_one, r, &i_one);
            resid2 = (norm_f / norm_r / eps);
            break;
        }
    }
    free_vector(out_zero);

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(2, 1, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
