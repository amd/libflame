/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_stedc(integer datatype, char *compz, integer *n, void *D, void *E, void *Z,
                      integer *ldz, void *work, integer *lwork, void *rwork, integer *lrwork,
                      integer *iwork, integer *liwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::stedc<float>(compz, n, (float *)D, (float *)E, (float *)Z, ldz, (float *)work,
                                   lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::stedc<double>(compz, n, (double *)D, (double *)E, (double *)Z, ldz,
                                    (double *)work, lwork, iwork, liwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::stedc<scomplex, float>(compz, n, (float *)D, (float *)E, (scomplex *)Z, ldz,
                                             (scomplex *)work, lwork, (float *)rwork, lrwork, iwork,
                                             liwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::stedc<dcomplex, double>(compz, n, (double *)D, (double *)E, (dcomplex *)Z,
                                              ldz, (dcomplex *)work, lwork, (double *)rwork, lrwork,
                                              iwork, liwork, info);
            break;
        }
    }
}
