#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
int cgeqr2p_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work,
                  integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
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
        xerbla_("CGEQR2P", &i__1, (ftnlen)7);
        return LAPACK_FAILURE;
    }
    return LAPACK_SUCCESS;
}
