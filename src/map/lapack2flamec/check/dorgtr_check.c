#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;

int dorgtr_check(char *uplo, aocl_int64_t *n, double *a, aocl_int64_t *lda, double *tau, double *work,
                 aocl_int64_t *lwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    aocl_int64_t nb;
    logical upper;
    aocl_int64_t lwkopt;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    lquery = *lwork == -1;
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
    else /* if(complicated condition) */
    {
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n - 1; // , expr subst
        if(*lwork < fla_max(i__1, i__2) && !lquery)
        {
            *info = -7;
        }
    }
    if(*info == 0)
    {
        if(upper)
        {
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            nb = aocl_lapack_ilaenv(&c__1, "DORGQL", " ", &i__1, &i__2, &i__3, &c_n1);
        }
        else
        {
            i__1 = *n - 1;
            i__2 = *n - 1;
            i__3 = *n - 1;
            nb = aocl_lapack_ilaenv(&c__1, "DORGQR", " ", &i__1, &i__2, &i__3, &c_n1);
        }
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n - 1; // , expr subst
        lwkopt = fla_max(i__1, i__2) * nb;
        work[1] = (double)lwkopt;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DORGTR", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        work[1] = 1.;
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
