/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gtsv(integer datatype, integer *n, integer *nrhs, void *dl, void *d, void *du, void *b,
                     integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gtsv<float>(n, nrhs, (float *)dl, (float *)d, (float *)du, (float *)b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gtsv<double>(n, nrhs, (double *)dl, (double *)d, (double *)du, (double *)b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gtsv<scomplex>(n, nrhs, (scomplex *)dl, (scomplex *)d, (scomplex *)du, (scomplex *)b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gtsv<dcomplex>(n, nrhs, (dcomplex *)dl, (dcomplex *)d, (dcomplex *)du, (dcomplex *)b, ldb, info);
            break;
        }
    }
}
