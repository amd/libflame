/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_lange(integer datatype, char *norm_type, integer *m, integer *n, void *A,
                      integer *lda, void *work, void *result)
{
    switch(datatype)
    {
        case FLOAT:
        {
            *(float *)result
                = libflame::lange<float>(norm_type, m, n, (float *)A, lda, (float *)work);
            break;
        }
        case DOUBLE:
        {
            *(double *)result
                = libflame::lange<double>(norm_type, m, n, (double *)A, lda, (double *)work);
            break;
        }
        case COMPLEX:
        {
            *(float *)result = libflame::lange<scomplex, float>(norm_type, m, n, (scomplex *)A, lda,
                                                                (float *)work);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            *(double *)result = libflame::lange<dcomplex, double>(norm_type, m, n, (dcomplex *)A,
                                                                  lda, (double *)work);
            break;
        }
    }
}
