/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_steqr(integer datatype, char *compz, integer *n, void *z, integer *ldz, void *d,
                      void *e, void *work, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::steqr<float>(compz, n, (float *)d, (float *)e, (float *)z, ldz, (float *)work,
                                   info);
            break;
        }
        case DOUBLE:
        {
            libflame::steqr<double>(compz, n, (double *)d, (double *)e, (double *)z, ldz,
                                    (double *)work, info);
            break;
        }
        case COMPLEX:
        {
            libflame::steqr<scomplex, float>(compz, n, (float *)d, (float *)e, (scomplex *)z, ldz,
                                             (float *)work, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::steqr<dcomplex, double>(compz, n, (double *)d, (double *)e, (dcomplex *)z,
                                              ldz, (double *)work, info);
            break;
        }
    }
}
