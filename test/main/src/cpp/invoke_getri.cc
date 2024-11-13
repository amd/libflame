/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_getri(integer datatype, integer *n, void *a, integer *lda, integer *ipiv, void *work,
                      integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::getri<float>(n, (float *)a, lda, ipiv, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::getri<double>(n, (double *)a, lda, ipiv, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::getri<scomplex>(n, (scomplex *)a, lda, ipiv, (scomplex *)work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::getri<dcomplex>(n, (dcomplex *)a, lda, ipiv, (dcomplex *)work, lwork, info);
            break;
        }
    }
}
