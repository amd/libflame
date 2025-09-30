#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;

int ssytrd_check(char *uplo, aocl_int64_t *n, float *a, aocl_int64_t *lda, float *d__, float *e, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;
    /* Local variables */
    aocl_int64_t nb;
    logical upper;
    aocl_int64_t lwkopt;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    lquery = *lwork == -1;
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
    else if(*lwork < 1 && !lquery)
    {
        *info = -9;
    }
    if(*info == 0)
    {
        /* Determine the block size. */
        nb = aocl_lapack_ilaenv(&c__1, "SSYTRD", uplo, n, &c_n1, &c_n1, &c_n1);
        lwkopt = *n * nb;
        work[1] = (float)lwkopt;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SSYTRD", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        work[1] = 1.f;
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
