/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE SYTRD API invoke function
*/
integer invoke_lapacke_sytrd(integer datatype, integer layout, char uplo, integer n, void *A,
                             integer lda, void *D, void *E, void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssytrd(layout, uplo, n, A, lda, D, E, tau);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dsytrd(layout, uplo, n, A, lda, D, E, tau);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_chetrd(layout, uplo, n, A, lda, D, E, tau);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zhetrd(layout, uplo, n, A, lda, D, E, tau);
            break;
        }
    }
    return info;
}
