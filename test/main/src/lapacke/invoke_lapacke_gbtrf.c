/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/
#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gbtrf(integer datatype, int layout, integer m, integer n, integer kl,
                             integer ku, void *ab, integer ldab, integer *ipiv)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgbtrf(layout, m, n, kl, ku, ab, ldab, ipiv);
            break;
        }
    }
    return info;
}
