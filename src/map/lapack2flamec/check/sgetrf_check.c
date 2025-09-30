#include "FLA_f2c.h" /* Table of constant values */
#include "FLA_lapack2flame_return_defs.h"

int sgetrf_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int_t *ipiv, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
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
        aocl_blas_xerbla("SGETRF", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
