/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_larfg(integer datatype, integer *n, void *x, integer *incx,
                             integer *abs_incx, void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            float *x_ptr = x;
            info = LAPACKE_slarfg(*n, x, &x_ptr[*abs_incx], *incx, tau);
            break;
        }
        case DOUBLE:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            double *x_ptr = x;
            info = LAPACKE_dlarfg(*n, x, &x_ptr[*abs_incx], *incx, tau);
            break;
        }
        case COMPLEX:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            scomplex *x_ptr = x;
            info = LAPACKE_clarfg(*n, x, (void *)(&x_ptr[*abs_incx]), *incx, tau);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            dcomplex *x_ptr = x;
            info = LAPACKE_zlarfg(*n, x, (void *)(&x_ptr[*abs_incx]), *incx, tau);
            break;
        }
    }
    return info;
}
