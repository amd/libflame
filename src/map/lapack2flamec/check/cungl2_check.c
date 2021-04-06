#include "FLA_f2c.h"
int cungl2_check(integer *m, integer *n, integer *k, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "cungl2 inputs: m %d, n %d, k %d, lda %d\n", *m, *n, *k, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif

    /* Builtin functions */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    if(*m < 0)
    {
        *info = -1;
    }
    else if(*n < *m)
    {
        *info = -2;
    }
    else if(*k < 0 || *k > *m)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -5;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CUNGL2", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*m <= 0)
    {
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
