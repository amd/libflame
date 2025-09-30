#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"

static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;

int zgebrd_check(aocl_int64_t *m, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda, double *d__, double *e,
                 dcomplex *tauq, dcomplex *taup, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2;
    double d__1;

    /* Local variables */
    aocl_int64_t nb;
    aocl_int64_t minmn;
    aocl_int64_t lwkopt;
    logical lquery;

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
    /* Computing MAX */
    i__1 = 1;
    i__2 = aocl_lapack_ilaenv(&c__1, "ZGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
    nb = fla_max(i__1, i__2);
    lwkopt = (*m + *n) * nb;
    d__1 = (double)lwkopt;
    work[1].real = d__1;
    work[1].imag = 0.; // , expr subst
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
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = fla_max(1, *m);
        if(*lwork < fla_max(i__1, *n) && !lquery)
        {
            *info = -10;
        }
    }
    if(*info < 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZGEBRD", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    minmn = fla_min(*m, *n);
    if(minmn == 0)
    {
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
