/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_stevd(integer datatype, int layout, char jobz, integer n, void *d, void *e,
                             void *z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sstevd(layout, jobz, n, d, e, z, ldz);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dstevd(layout, jobz, n, d, e, z, ldz);
            break;
        }
    }
    return info;
}
