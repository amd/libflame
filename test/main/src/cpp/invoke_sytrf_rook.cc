/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_sytrf_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                           integer *ipiv, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::sytrf_rook<float>(uplo, n, (float *)a, lda, ipiv, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::sytrf_rook<double>(uplo, n, (double *)a, lda, ipiv, (double *)work, lwork,
                                         info);
            break;
        }

        case COMPLEX:
        {
            libflame::sytrf_rook<scomplex>(uplo, n, (scomplex *)a, lda, ipiv, (scomplex *)work,
                                           lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::sytrf_rook<dcomplex>(uplo, n, (dcomplex *)a, lda, ipiv, (dcomplex *)work,
                                           lwork, info);
            break;
        }
    }
}
