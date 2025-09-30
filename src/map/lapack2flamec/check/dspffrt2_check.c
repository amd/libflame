/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Nov 06, 2020
*/

#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"

int dspffrt2_check(double *ap, aocl_int64_t *n, aocl_int64_t *ncolm, double *work, double *work2)
{
    aocl_int64_t ret_val = LAPACK_SUCCESS;

    if(*n < 0)
    {
        ret_val = LAPACK_FAILURE;
    }
    else if(*ncolm < 0 || *ncolm > *n)
    {
        ret_val = LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*n == 0 || *ncolm == 0)
    {
        ret_val = LAPACK_QUICK_RETURN;
    }
    return ret_val;
}
