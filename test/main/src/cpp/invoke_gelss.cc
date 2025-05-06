/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gelss(integer datatype, integer *m, integer *n, integer *nrhs, void *A, integer *lda,
                  void *B, integer *ldb, void *s, void *rcond, integer *rank, void *work,
                  integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gelss<float>(m, n, nrhs, (float *)A, lda, (float *)B, ldb, (float *)s, (float *)rcond, rank, (float *)work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::gelss<double>(m, n, nrhs, (double *)A, lda, (double *)B, ldb, (double *)s, (double *)rcond, rank, (double *)work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::gelss<scomplex, float>(m, n, nrhs, (scomplex *)A, lda, (scomplex *)B, ldb, (float *)s, (float *)rcond, rank, (scomplex *)work, lwork, (float *)rwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::gelss<dcomplex, double>(m, n, nrhs, (dcomplex *)A, lda, (dcomplex *)B, ldb, (double *)s, (double *)rcond, rank, (dcomplex *)work, lwork, (double *)rwork, info);
            break;
        }
    }
}
