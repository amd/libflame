/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gelsd(integer datatype, integer *m, integer *n, integer *nrhs, void *a, integer *lda,
                  void *b, integer *ldb, void *s, void *rcond, integer *rank, void *work,
                  integer *lwork, void *rwork, integer *iwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gelsd<float>(m, n, nrhs, (float *)a, lda, (float *)b, ldb, (float *)s, (float *)rcond, rank, (float *)work, lwork, iwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gelsd<double>(m, n, nrhs, (double *)a, lda, (double *)b, ldb, (double *)s, (double *)rcond, rank, (double *)work, lwork, iwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gelsd<scomplex, float>(m, n, nrhs, (scomplex *)a, lda, (scomplex *)b, ldb, (float *)s, (float *)rcond, rank, (scomplex *)work, lwork, (float *)rwork, iwork,
                              info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gelsd<dcomplex, double>(m, n, nrhs, (dcomplex *)a, lda, (dcomplex *)b, ldb, (double *)s, (double *)rcond, rank, (dcomplex *)work, lwork, (double *)rwork, iwork,
                              info);
            break;
        }
    }
}
