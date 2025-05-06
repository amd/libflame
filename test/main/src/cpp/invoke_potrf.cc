/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_potrf(char *uplo, integer datatype, integer *m, void *a, integer *lda,
                      integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::potrf<float>(uplo, m, (float *)a, lda, info);
            break;
        }
        case DOUBLE:
        {
            libflame::potrf<double>(uplo, m, (double *)a, lda, info);
            break;
        }
        case COMPLEX:
        {
            libflame::potrf<scomplex>(uplo, m, (scomplex *)a, lda, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::potrf<dcomplex>(uplo, m, (dcomplex *)a, lda, info);
            break;
        }
    }
}
