/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gerqf(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                      void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gerqf<float>(m, n, (float *)a, lda, (float *)tau, (float *)work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::gerqf<double>(m, n, (double *)a, lda, (double *)tau, (double *)work, lwork,
                                    info);
            break;
        }
        case COMPLEX:
        {
            libflame::gerqf<scomplex>(m, n, (scomplex *)a, lda, (scomplex *)tau, (scomplex *)work,
                                      lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::gerqf<dcomplex>(m, n, (dcomplex *)a, lda, (dcomplex *)tau, (dcomplex *)work,
                                      lwork, info);
            break;
        }
    }
}
