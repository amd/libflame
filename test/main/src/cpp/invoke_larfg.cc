/*
    Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>

#include "invoke_common.hh"

void invoke_cpp_larfg(integer datatype, integer *n, void *alpha, void *x, integer *incx, void *tau)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::larfg<float>(n, (float *)alpha, (float *)x, incx, (float *)tau);
            break;
        }
        case DOUBLE:
        {
            libflame::larfg<double>(n, (double *)alpha, (double *)x, incx, (double *)tau);
            break;
        }
        case COMPLEX:
        {
            libflame::larfg<scomplex>(n, (scomplex *)alpha, (scomplex *)x, incx, (scomplex *)tau);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::larfg<dcomplex>(n, (dcomplex *)alpha, (dcomplex *)x, incx, (dcomplex *)tau);
            break;
        }
    }
}
