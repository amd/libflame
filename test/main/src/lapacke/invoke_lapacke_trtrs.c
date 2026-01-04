/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/
#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_trtrs(integer datatype, int layout, char *uplo, char *trans, char *diag,
                             integer n, integer nrhs, void *A, integer lda, void *b, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_strtrs(layout, *uplo, *trans, *diag, n, nrhs, A, lda, b, ldb);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dtrtrs(layout, *uplo, *trans, *diag, n, nrhs, A, lda, b, ldb);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_ctrtrs(layout, *uplo, *trans, *diag, n, nrhs, A, lda, b, ldb);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_ztrtrs(layout, *uplo, *trans, *diag, n, nrhs, A, lda, b, ldb);
            break;
        }
    }
    return info;
}