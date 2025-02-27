/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE SYTRF API invoke function
*/
integer invoke_lapacke_sytrf(integer datatype, integer layout, char uplo, integer n, void *a,
                             integer lda, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssytrf(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dsytrf(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_csytrf(layout, uplo, n, a, lda, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zsytrf(layout, uplo, n, a, lda, ipiv);
            break;
        }
    }
    return info;
}
