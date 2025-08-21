/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/
#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_trtrs(char *uplo, char *trans, char *diag, integer datatype, integer *n, void *a,
                      integer *lda, integer *nrhs, void *b, integer *ldb, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::trtrs<float>(uplo, trans, diag, n, nrhs, (float *)a, lda, (float *)b, ldb,
                                   info);
            break;
        }

        case DOUBLE:
        {
            libflame::trtrs<double>(uplo, trans, diag, n, nrhs, (double *)a, lda, (double *)b, ldb,
                                    info);
            break;
        }

        case COMPLEX:
        {
            libflame::trtrs<scomplex>(uplo, trans, diag, n, nrhs, (scomplex *)a, lda, (scomplex *)b,
                                      ldb, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::trtrs<dcomplex>(uplo, trans, diag, n, nrhs, (dcomplex *)a, lda, (dcomplex *)b,
                                      ldb, info);
            break;
        }
    }
}
