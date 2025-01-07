/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_org2r(integer datatype, integer *m, integer *n, integer *min_A, void *a,
                      integer *lda, void *tau, void *work, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::org2r<float>(m, n, n, (float *)a, lda, (float *)tau, (float *)work, info);
            break;
        }

        case DOUBLE:
        {
            libflame::org2r<double>(m, n, n, (double *)a, lda, (double *)tau, (double *)work, info);
            break;
        }

        case COMPLEX:
        {
            libflame::ung2r<scomplex>(m, n, n, (scomplex *)a, lda, (scomplex *)tau,
                                      (scomplex *)work, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::ung2r<dcomplex>(m, n, n, (dcomplex *)a, lda, (dcomplex *)tau,
                                      (dcomplex *)work, info);
            break;
        }
    }
}
