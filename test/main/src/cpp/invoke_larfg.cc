/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_larfg(integer datatype, integer *n, void *x, integer *incx, integer *abs_incx,
                      void *tau)
{
    switch(datatype)
    {
        case FLOAT:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            float *x_ptr = (float *)x;
            libflame::larfg<float>(n, &x_ptr[0], &x_ptr[*abs_incx], incx, (float *)tau);
            break;
        }
        case DOUBLE:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            double *x_ptr = (double *)x;
            libflame::larfg<double>(n, &x_ptr[0], &x_ptr[*abs_incx], incx, (double *)tau);
            break;
        }
        case COMPLEX:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            scomplex *x_ptr = (scomplex *)x;
            libflame::larfg<scomplex>(n, &x_ptr[0], &x_ptr[*abs_incx], incx, (scomplex *)tau);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            /* First value of the x vector is treated as Alpha -> x[0] */
            dcomplex *x_ptr = (dcomplex *)x;
            libflame::larfg<dcomplex>(n, &x_ptr[0], &x_ptr[*abs_incx], incx, (dcomplex *)tau);
            break;
        }
    }
}
