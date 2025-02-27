/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_steqr(integer datatype, int layout, char compz, integer n, void *d, void *e,
                             void *z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_ssteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dsteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_csteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zsteqr(layout, compz, n, d, e, z, ldz);
            break;
        }
    }
    return info;
}
