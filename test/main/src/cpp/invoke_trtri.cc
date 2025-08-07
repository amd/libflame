/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_trtri(integer datatype, char *uplo, char *diag, integer *n, void *a, integer *lda,
                      integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::trtri<float>(uplo, diag, n, (float *)a, lda, info);
            break;
        }

        case DOUBLE:
        {
            libflame::trtri<double>(uplo, diag, n, (double *)a, lda, info);
            break;
        }

        case COMPLEX:
        {
            libflame::trtri<scomplex>(uplo, diag, n, (scomplex *)a, lda, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::trtri<dcomplex>(uplo, diag, n, (dcomplex *)a, lda, info);
            break;
        }
    }
}