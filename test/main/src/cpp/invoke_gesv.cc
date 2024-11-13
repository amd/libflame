/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gesv(integer datatype, integer *n, integer *nrhs, void *a, integer *lda, integer *ipiv,
                 void *b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gesv<float>(n, nrhs, (float *)a, lda, ipiv, (float *)b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gesv<double>(n, nrhs, (double *)a, lda, ipiv, (double *)b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gesv<scomplex>(n, nrhs, (scomplex *)a, lda, ipiv, (scomplex *)b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gesv<dcomplex>(n, nrhs, (dcomplex *)a, lda, ipiv, (dcomplex *)b, ldb, info);
            break;
        }
    }
}
