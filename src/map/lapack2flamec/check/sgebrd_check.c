#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static integer c__1 = 1;
static integer c_n1 = -1;

int sgebrd_check(integer *m, integer *n, float *a, integer *lda, float *d__, float *e, float *tauq,
                 float *taup, float *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer nb, minmn;
    integer lwkopt;
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
    i__2 = ilaenv_(&c__1, "SGEBRD", " ", m, n, &c_n1, &c_n1); // , expr subst
    nb = fla_max(i__1, i__2);
    lwkopt = (*m + *n) * nb;
    work[1] = (float)lwkopt;
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
        xerbla_("SGEBRD", &i__1, (ftnlen)6);
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
        work[1] = 1.f;
        return LAPACK_QUICK_RETURN;
    }

    return LAPACK_SUCCESS;
}
