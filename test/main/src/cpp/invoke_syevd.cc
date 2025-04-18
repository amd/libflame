/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_syevd(integer datatype, char *jobz, char *uplo, integer *n, void *a, integer *lda,
                      void *w, void *work, integer *lwork, void *rwork, integer *lrwork,
                      integer *iwork, integer *liwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::syevd<float>(jobz, uplo, n, (float *)a, lda, (float *)w, (float *)work, lwork,
                                   iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::syevd<double>(jobz, uplo, n, (double *)a, lda, (double *)w, (double *)work,
                                    lwork, iwork, liwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::heevd<scomplex, float>(jobz, uplo, n, (scomplex *)a, lda, (float *)w,
                                             (scomplex *)work, lwork, (float *)rwork, lrwork, iwork,
                                             liwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::heevd<dcomplex, double>(jobz, uplo, n, (dcomplex *)a, lda, (double *)w,
                                              (dcomplex *)work, lwork, (double *)rwork, lrwork,
                                              iwork, liwork, info);
            break;
        }
    }
}
