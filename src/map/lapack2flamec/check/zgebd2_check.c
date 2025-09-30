#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"

int zgebd2_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, doublereal *d__, doublereal *e,
                 dcomplex *tauq, dcomplex *taup, dcomplex *work, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --d__;
    --e;
    --tauq;
    --taup;
    --work;
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
    if(*info < 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZGEBD2", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
