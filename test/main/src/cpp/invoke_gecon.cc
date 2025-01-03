/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gecon(integer datatype, char *norm, integer *n, void *A, integer *lda, void *anorm,
                      void *rcond, void *work, void *lrwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gecon<float>(norm, n, (float *)A, lda, (float *)anorm, (float *)rcond,
                                   (float *)work, (integer *)lrwork, info);
            break;
        }

        case DOUBLE:
        {
            libflame::gecon<double>(norm, n, (double *)A, lda, (double *)anorm, (double *)rcond,
                                    (double *)work, (integer *)lrwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::gecon<scomplex, float>(norm, n, (scomplex *)A, lda, (float *)anorm,
                                             (float *)rcond, (scomplex *)work, (float *)lrwork,
                                             info);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            libflame::gecon<dcomplex, double>(norm, n, (dcomplex *)A, lda, (double *)anorm,
                                              (double *)rcond, (dcomplex *)work, (double *)lrwork,
                                              info);
            break;
        }
    }
}
