/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_getrf(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgetrf(layout, m, n, a, lda, ipiv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgetrf(layout, m, n, a, lda, ipiv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgetrf(layout, m, n, a, lda, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgetrf(layout, m, n, a, lda, ipiv);
            break;
        }
    }
    return info;
}
