/*
 * Copyright (c) 2020 Advanced Micro Devices, Inc.
 */

#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"

int sgetrfnp_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "dpotri inputs: uplo %c, n %d, lda %d\n", *uplo, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SGETRFNP", &i__1, (ftnlen)8);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
