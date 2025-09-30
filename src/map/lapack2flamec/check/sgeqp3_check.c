#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;

int sgeqp3_check(aocl_int64_t *m, aocl_int64_t *n, float *a, aocl_int64_t *lda, aocl_int64_t *jpvt, float *tau,
                 float *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1;
    /* Local variables */
    aocl_int64_t iws, nb;
    aocl_int64_t minmn;
    aocl_int64_t lwkopt;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --jpvt;
    --tau;
    --work;
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
            iws = *n * 3 + 1;
            nb = aocl_lapack_ilaenv(&c__1, "SGEQRF", " ", m, n, &c_n1, &c_n1);
            lwkopt = (*n << 1) + (*n + 1) * nb;
        }
        work[1] = (float)lwkopt;
        if(*lwork < iws && !lquery)
        {
            *info = -8;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("SGEQP3", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }

    return LAPACK_SUCCESS;
}
