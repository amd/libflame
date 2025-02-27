/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE SYTRF_ROOK API invoke function
*/
integer invoke_lapacke_sytrf_rook(integer datatype, integer layout, char uplo, integer n, void *a,
                                  integer lda, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssytrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dsytrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_csytrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zsytrf_rook(layout, uplo, n, a, lda, ipiv);
            break;
        }
    }
    return info;
}
