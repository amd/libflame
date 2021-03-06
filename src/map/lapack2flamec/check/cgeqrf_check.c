#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"
static integer c__1 = 1;
static integer c_n1 = -1;
int cgeqrf_check(integer *m, integer *n, scomplex *a, integer *lda, scomplex *tau, scomplex *work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    integer k, nb, lwkopt;
    logical lquery;
#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "cgeqrf inputs: m %d, n %d, lda %d\n", *m, *n, *lda);
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;
    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "CGEQRF", " ", m, n, &c_n1, &c_n1);
    lwkopt = *n * nb;
    work[1].real = (float) lwkopt;
    work[1].imag = 0.f; // , expr subst
    lquery = *lwork == -1;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < max(1,*m))
    {
        *info = -4;
    }
    else if (*lwork < max(1,*n) && ! lquery)
    {
        *info = -7;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEQRF", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    k = min(*m,*n);
    if (k == 0)
    {
        work[1].real = 1.f;
        work[1].imag = 0.f; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
