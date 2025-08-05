/*
    Copyright (c) 2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gebrd(integer datatype, int layout, integer m, integer n, void *a,
                             integer lda, void *d, void *e, void *tauq, void *taup)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgebrd(layout, m, n, a, lda, d, e, tauq, taup);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dgebrd(layout, m, n, a, lda, d, e, tauq, taup);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cgebrd(layout, m, n, a, lda, d, e, tauq, taup);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgebrd(layout, m, n, a, lda, d, e, tauq, taup);
            break;
        }
    }
    return info;
}

