/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_sygvd(integer datatype, integer *itype, char *jobz, char *uplo, integer *n, void *a,
                      integer *lda, void *b, integer *ldb, void *w, void *work, integer *lwork,
                      void *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::sygvd<float>(itype, jobz, uplo, n, (float *)a, lda, (float *)b, ldb,
                                   (float *)w, (float *)work, lwork, iwork, liwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::sygvd<double>(itype, jobz, uplo, n, (double *)a, lda, (double *)b, ldb,
                                    (double *)w, (double *)work, lwork, iwork, liwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::hegvd<scomplex, float>(itype, jobz, uplo, n, (scomplex *)a, lda,
                                             (scomplex *)b, ldb, (float *)w, (scomplex *)work,
                                             lwork, (float *)rwork, lrwork, iwork, liwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::hegvd<dcomplex, double>(itype, jobz, uplo, n, (dcomplex *)a, lda,
                                              (dcomplex *)b, ldb, (double *)w, (dcomplex *)work,
                                              lwork, (double *)rwork, lrwork, iwork, liwork, info);
            break;
        }
    }
}
