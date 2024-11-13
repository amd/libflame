/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_getrs(integer datatype, char *trans, integer *n, integer *nrhs, void *a, integer *lda,
                      integer *ipiv, void *b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::getrs<float>(trans, n, nrhs, (float *)a, lda, ipiv, (float *)b, ldb, info);
            break;
        }

        case DOUBLE:
        {
            libflame::getrs<double>(trans, n, nrhs, (double *)a, lda, ipiv, (double *)b, ldb, info);
            break;
        }

        case COMPLEX:
        {
            libflame::getrs<scomplex>(trans, n, nrhs, (scomplex *)a, lda, ipiv, (scomplex *)b, ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::getrs<dcomplex>(trans, n, nrhs, (dcomplex *)a, lda, ipiv, (dcomplex *)b, ldb, info);
            break;
        }
    }
}
