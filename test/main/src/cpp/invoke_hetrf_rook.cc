/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_hetrf_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                           integer *ipiv, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case COMPLEX:
        {
            libflame::hetrf_rook<scomplex>(uplo, n, (scomplex *)a, lda, ipiv, (scomplex *)work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::hetrf_rook<dcomplex>(uplo, n, (dcomplex *)a, lda, ipiv, (dcomplex *)work, lwork, info);
            break;
        }
    }
}
