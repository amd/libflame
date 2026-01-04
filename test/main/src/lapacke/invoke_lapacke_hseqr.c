/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_hseqr(integer datatype, int layout, char job, char compz, integer n,
                             integer ilo, integer ihi, void *h, integer ldh, void *w, void *wr,
                             void *wi, void *z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_shseqr(layout, job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dhseqr(layout, job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_chseqr(layout, job, compz, n, ilo, ihi, h, ldh, w, z, ldz);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zhseqr(layout, job, compz, n, ilo, ihi, h, ldh, w, z, ldz);
            break;
        }
    }
    return info;
}
