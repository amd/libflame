/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_sytrd(integer datatype, char *uplo, integer *n, void *A, integer *lda, void *D,
                      void *E, void *tau, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::sytrd<float>(uplo, n, (float *)A, lda, (float *)D, (float *)E, (float *)tau,
                                   (float *)work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::sytrd<double>(uplo, n, (double *)A, lda, (double *)D, (double *)E,
                                    (double *)tau, (double *)work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::hetrd<scomplex, float>(uplo, n, (scomplex *)A, lda, (float *)D, (float *)E,
                                             (scomplex *)tau, (scomplex *)work, lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::hetrd<dcomplex, double>(uplo, n, (dcomplex *)A, lda, (double *)D, (double *)E,
                                              (dcomplex *)tau, (dcomplex *)work, lwork, info);
            break;
        }
    }
}
