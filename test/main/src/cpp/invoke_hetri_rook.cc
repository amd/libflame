/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/
#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_hetri_rook(integer datatype, char *uplo, integer *n, void *a, integer *lda,
                           integer *ipiv, void *work, integer *info)
{
    switch(datatype)
    {
        case COMPLEX:
        {
            libflame::hetri_rook<scomplex>(uplo, n, (scomplex *)a, lda, ipiv, (scomplex *)work,
                                           info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::hetri_rook<dcomplex>(uplo, n, (dcomplex *)a, lda, ipiv, (dcomplex *)work,
                                           info);
            break;
        }
    }
}
