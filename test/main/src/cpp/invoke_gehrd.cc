/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gehrd(integer datatype, integer *n, integer *ilo, integer *ihi, void *A, integer *lda,
                  void *tau, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gehrd<float>(n, ilo, ihi, (float *)A, lda, (float *)tau, (float *)work, lwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gehrd<double>(n, ilo, ihi, (double *)A, lda, (double *)tau, (double *)work, lwork, info);
            break;
        }

        case COMPLEX:
        {
            libflame::gehrd<scomplex>(n, ilo, ihi, (scomplex *)A, lda, (scomplex *)tau, (scomplex *)work, lwork, info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gehrd<dcomplex>(n, ilo, ihi, (dcomplex *)A, lda, (dcomplex *)tau, (dcomplex *)work, lwork, info);
            break;
        }
    }
}
