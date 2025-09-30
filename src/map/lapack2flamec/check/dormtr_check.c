#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;

int dormtr_check(char *side, char *uplo, char *trans, aocl_int64_t *m, aocl_int64_t *n, double *a,
                 aocl_int64_t *lda, double *tau, double *c__, aocl_int64_t *ldc, double *work, aocl_int64_t *lwork,
                 aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, c_dim1, c_offset, i__2, i__3;
    char ch__1[2];
    /* Builtin functions */
    /* Subroutine */

    /* Local variables */
    aocl_int64_t nb, nq, nw;
    logical left;
    logical upper;
    aocl_int64_t lwkopt;
    logical lquery;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --work;
    /* Function Body */
    *info = 0;
    left = lsame_(side, "L", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    lquery = *lwork == -1;
    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if(left)
    {
        nq = *m;
        nw = *n;
    }
    else
    {
        nq = *n;
        nw = *m;
    }
    if(!left && !lsame_(side, "R", 1, 1))
    {
        *info = -1;
    }
    else if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -2;
    }
    else if(!lsame_(trans, "N", 1, 1) && !lsame_(trans, "T", 1, 1))
    {
        *info = -3;
    }
    else if(*m < 0)
    {
        *info = -4;
    }
    else if(*n < 0)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, nq))
    {
        *info = -7;
    }
    else if(*ldc < fla_max(1, *m))
    {
        *info = -10;
    }
    else if(*lwork < fla_max(1, nw) && !lquery)
    {
        *info = -12;
    }
    if(*info == 0)
    {
        if(upper)
        {
            if(left)
            {
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = aocl_lapack_ilaenv(&c__1, "DORMQL", ch__1, &i__2, n, &i__3, &c_n1);
            }
            else
            {
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = aocl_lapack_ilaenv(&c__1, "DORMQL", ch__1, m, &i__2, &i__3, &c_n1);
            }
        }
        else
        {
            if(left)
            {
                i__2 = *m - 1;
                i__3 = *m - 1;
                nb = aocl_lapack_ilaenv(&c__1, "DORMQR", ch__1, &i__2, n, &i__3, &c_n1);
            }
            else
            {
                i__2 = *n - 1;
                i__3 = *n - 1;
                nb = aocl_lapack_ilaenv(&c__1, "DORMQR", ch__1, m, &i__2, &i__3, &c_n1);
            }
        }
        lwkopt = fla_max(1, nw) * nb;
        work[1] = (double)lwkopt;
    }
    if(*info != 0)
    {
        i__2 = -(*info);
        aocl_blas_xerbla("DORMTR", &i__2, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0 || nq == 1)
    {
        work[1] = 1.;
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
