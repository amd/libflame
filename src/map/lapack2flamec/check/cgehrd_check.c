#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static aocl_int64_t c__1 = 1;
static aocl_int64_t c_n1 = -1;

int cgehrd_check(aocl_int64_t *n, aocl_int64_t *ilo, aocl_int64_t *ihi, scomplex *a, aocl_int64_t *lda, scomplex *tau,
                 scomplex *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    aocl_int64_t i__;
    aocl_int64_t nb, nh;
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
    /* Computing MIN */
    i__1 = 64;
    i__2 = aocl_lapack_ilaenv(&c__1, "CGEHRD", " ", n, ilo, ihi, &c_n1); // , expr subst
    nb = fla_min(i__1, i__2);
    lwkopt = *n * nb;
    work[1].real = (float)lwkopt;
    work[1].imag = 0.f; // , expr subst
    lquery = *lwork == -1;
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*ilo < 1 || *ilo > fla_max(1, *n))
    {
        *info = -2;
    }
    else if(*ihi < fla_min(*ilo, *n) || *ihi > *n)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*lwork < fla_max(1, *n) && !lquery)
    {
        *info = -8;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGEHRD", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        return LAPACK_QUERY_RETURN;
    }
    /* Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */
    i__1 = *ilo - 1;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = i__;
        tau[i__2].real = 0.f;
        tau[i__2].imag = 0.f; // , expr subst
        /* L10: */
    }
    i__1 = *n - 1;
    for(i__ = fla_max(1, *ihi); i__ <= i__1; ++i__)
    {
        i__2 = i__;
        tau[i__2].real = 0.f;
        tau[i__2].imag = 0.f; // , expr subst
        /* L20: */
    }
    /* Quick return if possible */
    nh = *ihi - *ilo + 1;
    if(nh <= 1)
    {
        work[1].real = 1.f;
        work[1].imag = 0.f; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
