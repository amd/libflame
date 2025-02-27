/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE GTSV API invoke function
*/
integer invoke_lapacke_gtsv(integer datatype, integer layout, integer n, integer nrhs, void *dl,
                            void *d, void *du, void *B, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgtsv(layout, n, nrhs, dl, d, du, B, ldb);
            break;
        }
    }
    return info;
}
