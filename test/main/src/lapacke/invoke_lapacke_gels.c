/*
    Copyright (C) 2025, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "lapacke.h"
#include "test_common.h"

/*
LAPACKE GELS API invoke function
*/
integer invoke_lapacke_gels(integer datatype, integer layout, char trans, integer m, integer n,
                            integer nrhs, void *A, integer lda, void *B, integer ldb)
{
    integer info = 0;
    switch(datatype)
    {
        case FLOAT:
        {
            info = LAPACKE_sgels(layout, trans, m, n, nrhs, A, lda, B, ldb);
            break;
        }

        case DOUBLE:
        {
            info = LAPACKE_dgels(layout, trans, m, n, nrhs, A, lda, B, ldb);
            break;
        }

        case COMPLEX:
        {
            info = LAPACKE_cgels(layout, trans, m, n, nrhs, A, lda, B, ldb);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            info = LAPACKE_zgels(layout, trans, m, n, nrhs, A, lda, B, ldb);
            break;
        }
    }
    return info;
}
