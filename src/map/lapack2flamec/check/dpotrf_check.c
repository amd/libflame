#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"

int dpotrf_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, aocl_int64_t *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1;    
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "dpotrf inputs: uplo %c, n %d, lda %d\n", *uplo, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* Local variables */
    logical upper;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DPOTRF", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
