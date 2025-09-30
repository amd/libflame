#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
int cgeqpf_check(aocl_int64_t *m, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt, scomplex *tau,
                 scomplex *work, float *rwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;
    --rwork;
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
        aocl_blas_xerbla("CGEQPF", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
