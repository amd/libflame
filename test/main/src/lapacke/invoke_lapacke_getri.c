/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_getri(integer datatype, int layout, integer n, void *a, integer lda,
                             const integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgetri(layout, n, a, lda, ipiv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgetri(layout, n, a, lda, ipiv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgetri(layout, n, a, lda, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgetri(layout, n, a, lda, ipiv);
            break;
        }
    }
    return info;
}
