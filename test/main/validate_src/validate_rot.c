/*
    Copyright (C) 2022-2025, Advanced Micro Devices, Inc. All rights reserved.
*/

/*! @file validate_rot.c
 *  @brief Defines validate function of ROT() to use in test suite.
 *  */

#include "test_common.h"
#include "test_prototype.h"

extern double perf;
extern double time_min;

void validate_rot(char *tst_api, integer datatype, integer n, void *cx, void *cx_test, integer incx,
                  void *cy, void *cy_test, integer incy, void *c, void *s, double err_thresh)
{
    double residual, resid1 = 0., resid2 = 0.;
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
            float eps = fla_lapack_slamch("P");
            float norm_cx, norm_cy;

            norm_cx = snrm2_(&n, cx, &incx);
            norm_cy = snrm2_(&n, cy, &incy);
            /*Hermitian Transpose of original rotation vector*/
            ((float *)s)[0] = -((float *)s)[0];

            /*Rotating output vectors*/
            fla_lapack_srot(&n, cx, &incx, cy, &incy, ((float *)c), ((float *)s));

            /*cx-cx_test*/
            saxpy_(&n, &s_n_one, cx_test, &incx, cx, &incx);
            saxpy_(&n, &s_n_one, cy_test, &incy, cy, &incy);

            resid1 = snrm2_(&n, cx, &incx);
            resid1 = resid1 / eps / norm_cx / ((float)n);
            resid2 = snrm2_(&n, cy, &incy);
            resid2 = resid2 / eps / norm_cy / ((float)n);
            break;
        }

        case DOUBLE:
        {
            double eps = fla_lapack_slamch("P");
            double norm_cx, norm_cy;

            norm_cx = dnrm2_(&n, cx, &incx);
            norm_cy = dnrm2_(&n, cy, &incy);
            /*Hermitian Transpose of original rotation vector*/
            ((double *)s)[0] = -((double *)s)[0];

            /*Rotating output vectors*/
            fla_lapack_drot(&n, cx, &incx, cy, &incy, ((double *)c), ((double *)s));

            /*cx-cx_test*/
            daxpy_(&n, &d_n_one, cx_test, &incx, cx, &incx);
            daxpy_(&n, &d_n_one, cy_test, &incy, cy, &incy);

            resid1 = dnrm2_(&n, cx, &incx);
            resid1 = resid1 / eps / norm_cx / ((double)n);
            resid2 = dnrm2_(&n, cy, &incy);
            resid2 = resid2 / eps / norm_cy / ((double)n);
            break;
        }
        case COMPLEX:
        {
            float eps = fla_lapack_slamch("P");
            float norm_cx, norm_cy;

            norm_cx = scnrm2_(&n, cx, &incx);
            norm_cy = scnrm2_(&n, cy, &incy);

            /*Hermitian Transpose of original rotation vector*/
            ((scomplex *)s)[0].real = -((scomplex *)s)[0].real;
            ((scomplex *)s)[0].imag = -((scomplex *)s)[0].imag;

            /*Rotating output vectors*/
            fla_lapack_crot(&n, cx, &incx, cy, &incy, ((float *)c), ((scomplex *)s));

            /*cx-cx_test*/
            caxpy_(&n, &c_n_one, cx_test, &incx, cx, &incx);
            caxpy_(&n, &c_n_one, cy_test, &incy, cy, &incy);

            resid1 = scnrm2_(&n, cx, &incx);
            resid1 = resid1 / eps / norm_cx / ((float)n);
            resid2 = scnrm2_(&n, cy, &incy);
            resid2 = resid2 / eps / norm_cy / ((float)n);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double eps = fla_lapack_slamch("P");
            double norm_cx, norm_cy;

            norm_cx = dznrm2_(&n, cx, &incx);
            norm_cy = dznrm2_(&n, cy, &incy);

            /*Hermitian Transpose of original rotation vector*/
            ((dcomplex *)s)[0].real = -((dcomplex *)s)[0].real;
            ((dcomplex *)s)[0].imag = -((dcomplex *)s)[0].imag;

            /*Rotating output vectors*/
            fla_lapack_zrot(&n, cx, &incx, cy, &incy, ((double *)c), ((dcomplex *)s));

            /*cx-cx_test*/
            zaxpy_(&n, &z_n_one, cx_test, &incx, cx, &incx);
            zaxpy_(&n, &z_n_one, cy_test, &incy, cy, &incy);

            resid1 = dznrm2_(&n, cx, &incx);
            resid1 = resid1 / eps / norm_cx / ((double)n);
            resid2 = dznrm2_(&n, cy, &incy);
            resid2 = resid2 / eps / norm_cy / ((double)n);
            break;
        }
    }

    residual = fla_test_max(resid1, resid2);
    FLA_PRINT_TEST_STATUS(n, 2, residual, err_thresh);
    FLA_PRINT_SUBTEST_STATUS(resid1, err_thresh, "01");
    FLA_PRINT_SUBTEST_STATUS(resid2, err_thresh, "02");
}
