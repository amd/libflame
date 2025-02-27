/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_potrf(integer datatype, int layout, char uplo, integer n, void *a,
                             integer lda)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_spotrf(layout, uplo, n, a, lda);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dpotrf(layout, uplo, n, a, lda);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cpotrf(layout, uplo, n, a, lda);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zpotrf(layout, uplo, n, a, lda);
            break;
        }
    }
    return info;
}
