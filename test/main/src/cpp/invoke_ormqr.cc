/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/
#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"
void invoke_cpp_ormqr(integer datatype, char *side, char *trans, integer *m, integer *n, integer *k,
                      void *A, integer *lda, void *tau, void *C, integer *ldc, void *work,
                      integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::ormqr<float>(side, trans, m, n, k, (float *)A, lda, (float *)tau, (float *)C,
                                   ldc, (float *)work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::ormqr<double>(side, trans, m, n, k, (double *)A, lda, (double *)tau,
                                    (double *)C, ldc, (double *)work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::unmqr<scomplex>(side, trans, m, n, k, (scomplex *)A, lda, (scomplex *)tau,
                                      (scomplex *)C, ldc, (scomplex *)work, lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::unmqr<dcomplex>(side, trans, m, n, k, (dcomplex *)A, lda, (dcomplex *)tau,
                                      (dcomplex *)C, ldc, (dcomplex *)work, lwork, info);
            break;
        }
    }
}