/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_rot(integer datatype, integer *n, void *cx, integer *incx, void *cy, integer *incy,
                    void *c, void *s)
{
    switch(datatype)
    {
        case FLOAT:
        {
            fla_lapack_srot(n, (float *)cx, incx, (float *)cy, incy, (float *)c, (float *)s);
            break;
        }

        case DOUBLE:
        {
            fla_lapack_drot(n, (double *)cx, incx, (double *)cy, incy, (double *)c, (double *)s);
            break;
        }
        case COMPLEX:
        {
            libflame::rot<scomplex, float>(n, (scomplex *)cx, incx, (scomplex *)cy, incy,
                                           (float *)c, (scomplex *)s);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::rot<dcomplex, double>(n, (dcomplex *)cx, incx, (dcomplex *)cy, incy,
                                            (double *)c, (dcomplex *)s);
            break;
        }
    }
}
