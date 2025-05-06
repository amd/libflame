/*
    Copyright (C) 2024, Advanced Micro Devices, Inc. All rights reserved.
*/

#include <test_lapack.h>
#include <libflame_interface.hh>
#include "invoke_common.hh"

void invoke_cpp_gels(integer datatype, char *trans, integer *m, integer *n, integer *nrhs, void *A,
                 integer *lda, void *B, integer *ldb, void *work, integer *lwork, integer *info)
{
    switch(datatype)
    {
        case FLOAT:
        {
            libflame::gels<float>(trans, m, n, nrhs, (float *)A, lda, (float *)B, ldb, (float *)work, lwork, info);
            break;
        }
        case DOUBLE:
        {
            libflame::gels<double>(trans, m, n, nrhs, (double *)A, lda, (double *)B, ldb, (double *)work, lwork, info);
            break;
        }
        case COMPLEX:
        {
            libflame::gels<scomplex>(trans, m, n, nrhs, (scomplex *)A, lda, (scomplex *)B, ldb, (scomplex *)work, lwork, info);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            libflame::gels<dcomplex>(trans, m, n, nrhs, (dcomplex *)A, lda, (dcomplex *)B, ldb, (dcomplex *)work, lwork, info);
            break;
        }
    }
}
