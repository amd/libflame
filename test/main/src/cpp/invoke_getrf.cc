/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv,
                    integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::getrf<float>(m, n, (float *)a, lda, ipiv, info);
            break;
        }

        case DOUBLE:
        {
            libflame::getrf<double>(m, n, (double *)a, lda, ipiv, info);
            break;
        }

        case COMPLEX:
        {
            libflame::getrf<scomplex>(m, n, (scomplex *)a, lda, ipiv, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::getrf<dcomplex>(m, n, (dcomplex *)a, lda, ipiv, info);
            break;
        }
    }
}
