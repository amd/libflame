/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_potri(char *uplo, integer datatype, integer *n, void *a, integer *lda,
                      integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::potri<float>(uplo, n, (float *)a, lda, info);
            break;
        }
        case DOUBLE:
        {
            libflame::potri<double>(uplo, n, (double *)a, lda, info);
            break;
        }
        case COMPLEX:
        {
            libflame::potri<scomplex>(uplo, n, (scomplex *)a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::potri<dcomplex>(uplo, n, (dcomplex *)a, lda, info);
            break;
        }
    }
}
