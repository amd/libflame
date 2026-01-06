/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>

#include "invoke_common.hh"

void invoke_cpp_bdsqr(integer datatype, char *uplo, integer *n, integer *ncvt, integer *nru,
                      integer *ncc, void *d, void *e, void *vt, integer *ldvt, void *u,
                      integer *ldu, void *c, integer *ldc, void *work, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::bdsqr<float>(uplo, n, ncvt, nru, ncc, (float *)d, (float *)e, (float *)vt,
                                   ldvt, (float *)u, ldu, (float *)c, ldc, (float *)work, info);
            break;
        }

        case DOUBLE:
        {
            libflame::bdsqr<double>(uplo, n, ncvt, nru, ncc, (double *)d, (double *)e, (double *)vt,
                                    ldvt, (double *)u, ldu, (double *)c, ldc, (double *)work, info);
            break;
        }

        case COMPLEX:
        {
            libflame::bdsqr<scomplex, float>(uplo, n, ncvt, nru, ncc, (float *)d, (float *)e,
                                             (scomplex *)vt, ldvt, (scomplex *)u, ldu,
                                             (scomplex *)c, ldc, (float *)work, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::bdsqr<dcomplex, double>(uplo, n, ncvt, nru, ncc, (double *)d, (double *)e,
                                              (dcomplex *)vt, ldvt, (dcomplex *)u, ldu,
                                              (dcomplex *)c, ldc, (double *)work, info);
            break;
        }
    }
}