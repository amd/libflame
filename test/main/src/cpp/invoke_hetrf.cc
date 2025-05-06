/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/
#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_hetrf(integer datatype, char *uplo, integer *n, void *A, integer *lda,
                      integer *ipiv, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case COMPLEX:
        {
            libflame::hetrf<scomplex>(uplo, n, (scomplex *)A, lda, ipiv, (scomplex *)work, lwork,
                                      info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::hetrf<dcomplex>(uplo, n, (dcomplex *)A, lda, ipiv, (dcomplex *)work, lwork,
                                      info);
            break;
        }
    }
}
