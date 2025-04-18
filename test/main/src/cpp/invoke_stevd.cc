/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_stevd(integer datatype, char *jobz, integer *n, void *z, integer *ldz, void *d,
                      void *e, void *work, integer *lwork, integer *iwork, integer *liwork,
                      integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::stevd<float>(jobz, n, (float *)d, (float *)e, (float *)z, ldz, (float *)work,
                                   lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::stevd<double>(jobz, n, (double *)d, (double *)e, (double *)z, ldz,
                                    (double *)work, lwork, iwork, liwork, info);
            break;
        }
    }
}
