/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_potrs(char *uplo, integer datatype, integer *n, void *A, integer *lda,
                      integer *nrhs, void *B, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::potrs<float>(uplo, n, nrhs, (float *)A, lda, (float *)B, ldb, info);
            break;
        }
        case DOUBLE:
        {
            libflame::potrs<double>(uplo, n, nrhs, (double *)A, lda, (double *)B, ldb, info);
            break;
        }
        case COMPLEX:
        {
            libflame::potrs<scomplex>(uplo, n, nrhs, (scomplex *)A, lda, (scomplex *)B, ldb, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::potrs<dcomplex>(uplo, n, nrhs, (dcomplex *)A, lda, (dcomplex *)B, ldb, info);
            break;
        }
    }
}
