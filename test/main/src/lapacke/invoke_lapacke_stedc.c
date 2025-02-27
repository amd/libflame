/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

integer invoke_lapacke_stedc(integer datatype, int layout, char compz, integer n, void *D, void *E,
                             void *Z, integer ldz)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sstedc(layout, compz, n, D, E, Z, ldz);
            break;
        }
        case DOUBLE:
        {
            info = LAPACKE_dstedc(layout, compz, n, D, E, Z, ldz);
            break;
        }
        case COMPLEX:
        {
            info = LAPACKE_cstedc(layout, compz, n, D, E, Z, ldz);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zstedc(layout, compz, n, D, E, Z, ldz);
            break;
        }
    }
    return info;
}
