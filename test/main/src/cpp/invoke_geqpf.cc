/*
    Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_geqpf(integer datatype, integer *m, integer *n, void *a, integer *lda,
                      integer *jpvt, void *tau, void *work, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::geqpf<float>(m, n, (float *)a, lda, jpvt, (float *)tau, (float *)work, info);
            break;
        }

        case DOUBLE:
        {
            libflame::geqpf<double>(m, n, (double *)a, lda, jpvt, (double *)tau, (double *)work,
                                    info);
            break;
        }

        case COMPLEX:
        {
            libflame::geqpf<scomplex, float>(m, n, (scomplex *)a, lda, jpvt, (scomplex *)tau,
                                             (scomplex *)work, (float *)rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::geqpf<dcomplex, double>(m, n, (dcomplex *)a, lda, jpvt, (dcomplex *)tau,
                                              (dcomplex *)work, (double *)rwork, info);
            break;
        }
    }
}