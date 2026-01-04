/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE HETRF_ROOK API invoke function
*/
integer invoke_lapacke_hetrf_rook(integer datatype, integer layout, char uplo, integer n, void *a,
                                  integer lda, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case COMPLEX:
        {
            info = LAPACKE_chetrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zhetrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }
    }
    return info;
}
