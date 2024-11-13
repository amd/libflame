/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_geqp3(integer datatype, integer *m, integer *n, void *a, integer *lda, integer* jpvt, void *tau,
                      void *work, integer *lwork, void *rwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::geqp3<float>(m, n, (float *)a, lda, jpvt, (float *)tau, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::geqp3<double>(m, n, (double *)a, lda, jpvt, (double *)tau, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::geqp3<scomplex, float>(m, n, (scomplex *)a, lda, jpvt, (scomplex *)tau, (scomplex *)work, lwork, (float *)rwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::geqp3<dcomplex, double>(m, n, (dcomplex *)a, lda, jpvt, (dcomplex *)tau, (dcomplex *)work, lwork, (double *)rwork, info);
            break;
        }
    }
}
