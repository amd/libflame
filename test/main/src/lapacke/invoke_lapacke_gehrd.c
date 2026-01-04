/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_gehrd(integer datatype, int layout, integer n, integer ilo, integer ihi,
                             void *a, integer lda, void *tau)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgehrd(layout, n, ilo, ihi, a, lda, tau);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgehrd(layout, n, ilo, ihi, a, lda, tau);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgehrd(layout, n, ilo, ihi, a, lda, tau);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgehrd(layout, n, ilo, ihi, a, lda, tau);
            break;
        }
    }
    return info;
}
