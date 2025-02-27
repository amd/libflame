/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_potrs(integer datatype, int layout, char uplo, integer n, integer nrhs,
                             const void *A, integer lda, void *B, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_spotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dpotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cpotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zpotrs(layout, uplo, n, nrhs, A, lda, B, ldb);
            break;
        }
    }
    return info;
}
