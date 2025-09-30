#include "FLA_f2c.h" /* Table of constant values */
#include "FLA_lapack2flame_return_defs.h"
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;

int zgeqp3_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, aocl_int64_t *jpvt, dcomplex *tau,
                 dcomplex *work, aocl_int64_t *lwork, double *rwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;

    /* Local variables */
    aocl_int64_t nb, iws, minmn;
    aocl_int64_t lwkopt;
    logical lquery;

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
    lquery = *lwork == -1;
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
    if(*info == 0)
    {
        minmn = fla_min(*m, *n);
        if(minmn == 0)
        {
            iws = 1;
            lwkopt = 1;
        }
        else
        {
            iws = *n + 1;
            nb = aocl_lapack_ilaenv(&c__1, "ZGEQRF", " ", m, n, &c_n1, &c_n1);
            lwkopt = (*n + 1) * nb;
        }
        work[1].real = (double)lwkopt;
        work[1].imag = 0.; // , expr subst
        if(*lwork < iws && !lquery)
        {
            *info = -8;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZGEQP3", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible. */
    if(minmn == 0)
    {
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
