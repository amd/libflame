/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_orgqr(integer datatype, integer *m, integer *n, integer *min_A, void *a,
                      integer *lda, void *tau, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::orgqr<float>(m, n, n, (float *)a, lda, (float *)tau, (float *)work, lwork,
                                   info);
            break;
        }

        case DOUBLE:
        {
            libflame::orgqr<double>(m, n, n, (double *)a, lda, (double *)tau, (double *)work, lwork,
                                    info);
            break;
        }

        case COMPLEX:
        {
            libflame::ungqr<scomplex>(m, n, n, (scomplex *)a, lda, (scomplex *)tau,
                                      (scomplex *)work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::ungqr<dcomplex>(m, n, n, (dcomplex *)a, lda, (dcomplex *)tau,
                                      (dcomplex *)work, lwork, info);
            break;
        }
    }
}
