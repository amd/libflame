/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_idamax_small.c
 *  @brief finds index of first element having max abs value
 *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT

/* IDAMAX for small sizes */
fla_dim_t fla_idamax_small(fla_dim_t *n, doublereal *dx, fla_dim_t *incx)
{
    /* System generated locals */
    fla_dim_t ret_val, i__1;
    --dx;
    /* Builtin functions */
    /* Function Body */
    ret_val = 0;
    if (*n < 1 || *incx <= 0)
    {
        return ret_val;
    }
    ret_val = 1;
    if (*n == 1)
    {
        return ret_val;
    }
    if (*incx == 1)
    {
        doublereal temp;
        i__1 = *n;
        doublereal dmax = f2c_dabs(dx[1]);

        /* index of the first element having maximum absolute value */
        for(fla_dim_t i = 2; i<= i__1; i++ )
        {
            temp = f2c_dabs(dx[i]);
            if(temp > dmax)
            {
                dmax = temp;
                ret_val = i;
            }
        }
    }
    else
    {
        aocl_blas_idamax(n, dx, incx);
    }
    
    return ret_val;
}
#endif
