/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_syev(integer datatype, char *jobz, char *uplo, integer *n, void *a, integer *lda,
                     void *w, void *work, integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::syev<float>(jobz, uplo, n, (float *)a, lda, (float *)w, (float *)work, lwork,
                                  info);
            break;
        }
        case DOUBLE:
        {
            libflame::syev<double>(jobz, uplo, n, (double *)a, lda, (double *)w, (double *)work,
                                   lwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::heev<scomplex, float>(jobz, uplo, n, (scomplex *)a, lda, (float *)w,
                                            (scomplex *)work, lwork, (float *)rwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::heev<dcomplex, double>(jobz, uplo, n, (dcomplex *)a, lda, (double *)w,
                                             (dcomplex *)work, lwork, (double *)rwork, info);
            break;
        }
    }
}
