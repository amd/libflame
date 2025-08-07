/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_trtri(integer datatype, int layout, char uplo, char diag, integer n, void *a,
                             integer lda)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_strtri(layout, uplo, diag, n, a, lda);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dtrtri(layout, uplo, diag, n, a, lda);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_ctrtri(layout, uplo, diag, n, a, lda);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_ztrtri(layout, uplo, diag, n, a, lda);
            break;
        }
    }
    return info;
}
