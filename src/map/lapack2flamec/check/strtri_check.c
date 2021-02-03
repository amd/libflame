#include "FLA_f2c.h" /* Table of constant values */
#include "FLA_lapack2flame_return_defs.h"

int strtri_check(char *uplo, char *diag, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;

    /* Local variables */
    logical upper;
    logical nounit;
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "strtri inputs: uplo %c, diag %c, n %d, lda %d\n", *uplo, *diag, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(!nounit && !lsame_(diag, "U", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("STRTRI", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    /* Check for singularity if non-unit. */
    if(nounit)
    {
        i__1 = *n;
        for(*info = 1; *info <= i__1; ++(*info))
        {
            if(a[*info + *info * a_dim1] == 0.f)
            {
                return LAPACK_FAILURE;
            }
        }
        *info = 0;
    }
    return LAPACK_SUCCESS;
}
