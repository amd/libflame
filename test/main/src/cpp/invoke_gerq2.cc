/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gerq2(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau,
                      void *work, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gerq2<float>(m, n, (float *)a, lda, (float *)tau, (float *)work, info);
            break;
        }
        case DOUBLE:
        {
            libflame::gerq2<double>(m, n, (double *)a, lda, (double *)tau, (double *)work, info);
            break;
        }
        case COMPLEX:
        {
            libflame::gerq2<scomplex>(m, n, (scomplex *)a, lda, (scomplex *)tau, (scomplex *)work,
                                      info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::gerq2<dcomplex>(m, n, (dcomplex *)a, lda, (dcomplex *)tau, (dcomplex *)work,
                                      info);
            break;
        }
    }
}
