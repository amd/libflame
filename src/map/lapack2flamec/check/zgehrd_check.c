#include "FLA_lapack2flame_return_defs.h"
#include "FLA_f2c.h"

static integer c__1 = 1;
static integer c_n1 = -1;

int zgehrd_check(integer *n, integer *ilo, integer *ihi, dcomplex *a, integer *lda, dcomplex *tau, dcomplex * work, integer *lwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer i__;
    integer nb, nh;
    integer lwkopt;
    logical lquery;

#if AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "zgehrd inputs: n %d, ilo %d, ihi %d, lda %d\n", *n, *ilo, *ihi, *lda);
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
    /* Computing MIN */
    i__1 = 64;
    i__2 = ilaenv_(&c__1, "ZGEHRD", " ", n, ilo, ihi, &c_n1); // , expr subst
    nb = min(i__1,i__2);
    lwkopt = *n * nb;
    work[1].real = (double) lwkopt;
    work[1].imag = 0.; // , expr subst
    lquery = *lwork == -1;
    if (*n < 0)
    {
        *info = -1;
    }
    else if (*ilo < 1 || *ilo > max(1,*n))
    {
        *info = -2;
    }
    else if (*ihi < min(*ilo,*n) || *ihi > *n)
    {
        *info = -3;
    }
    else if (*lda < max(1,*n))
    {
        *info = -5;
    }
    else if (*lwork < max(1,*n) && ! lquery)
    {
        *info = -8;
    }
    if (*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZGEHRD", &i__1);
        return LAPACK_FAILURE;
    }
    else if (lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */
    i__1 = *ilo - 1;
    for (i__ = 1;
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        tau[i__2].real = 0.;
        tau[i__2].imag = 0.; // , expr subst
        /* L10: */
    }
    i__1 = *n - 1;
    for (i__ = max(1,*ihi);
            i__ <= i__1;
            ++i__)
    {
        i__2 = i__;
        tau[i__2].real = 0.;
        tau[i__2].imag = 0.; // , expr subst
        /* L20: */
    }
    /* Quick return if possible */
    nh = *ihi - *ilo + 1;
    if (nh <= 1)
    {
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
