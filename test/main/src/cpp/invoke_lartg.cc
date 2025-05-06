/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_lartg(integer datatype, void *f, void *g, void *c, void *s, void *r)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::lartg<float>((float *)f, (float *)g, (float *)c, (float *)s, (float *)r);
            break;
        }

        case DOUBLE:
        {
            libflame::lartg<double>((double *)f, (double *)g, (double *)c, (double *)s, (double *)r);
            break;
        }
        case COMPLEX:
        {
            libflame::lartg<scomplex, float> ((scomplex *)f, (scomplex *)g, (float *)c, ((scomplex *)s), (scomplex *)r);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::lartg<dcomplex, double> ((dcomplex *)f, (dcomplex *)g, (double *)c, ((dcomplex *)s), (dcomplex *)r);
            break;
        }
    }
}
