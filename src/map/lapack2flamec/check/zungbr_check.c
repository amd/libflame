#include "FLA_f2c.h"
#include "FLA_lapack2flame_return_defs.h"
static aocl_int64_t c_n1 = -1;

int zungbr_check(char *vect, aocl_int64_t *m, aocl_int64_t *n, aocl_int64_t *k, dcomplex *a, aocl_int64_t *lda,
                 dcomplex *tau, dcomplex *work, aocl_int64_t *lwork, aocl_int64_t *info)
{
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    aocl_int64_t mn;
    aocl_int64_t iinfo;
    logical wantq;
    aocl_int64_t lwkopt;
    logical lquery;
    extern integer zunglq_check(aocl_int64_t *, aocl_int64_t *, aocl_int64_t *, dcomplex *, aocl_int64_t *, dcomplex *,
                                dcomplex *, aocl_int64_t *, aocl_int64_t *),
        zungqr_check(aocl_int64_t *, aocl_int64_t *, aocl_int64_t *, dcomplex *, aocl_int64_t *, dcomplex *, dcomplex *,
                     aocl_int64_t *, aocl_int64_t *);

#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
    sprintf(buffer, "zungbr inputs: vect %c, m %d, n %d, k %d, lda %d\n", *vect, *m, *n, *k, *lda);
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
    wantq = lsame_(vect, "Q", 1, 1);
    mn = fla_min(*m, *n);
    lquery = *lwork == -1;
    if(!wantq && !lsame_(vect, "P", 1, 1))
    {
        *info = -1;
    }
    else if(*m < 0)
    {
        *info = -2;
    }
    else if(*n < 0 || wantq && (*n > *m || *n < fla_min(*m, *k))
            || !wantq && (*m > *n || *m < fla_min(*n, *k)))
    {
        *info = -3;
    }
    else if(*k < 0)
    {
        *info = -4;
    }
    else if(*lda < fla_max(1, *m))
    {
        *info = -6;
    }
    else if(*lwork < fla_max(1, mn) && !lquery)
    {
        *info = -9;
    }
    if(*info == 0)
    {
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        if(wantq)
        {
            if(*m >= *k)
            {
                zungqr_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if(*m > 1)
                {
                    i__1 = *m - 1;
                    i__2 = *m - 1;
                    i__3 = *m - 1;
                    zungqr_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], &work[1],
                                 &c_n1, &iinfo);
                }
            }
        }
        else
        {
            if(*k < *n)
            {
                zunglq_check(m, n, k, &a[a_offset], lda, &tau[1], &work[1], &c_n1, &iinfo);
            }
            else
            {
                if(*n > 1)
                {
                    i__1 = *n - 1;
                    i__2 = *n - 1;
                    i__3 = *n - 1;
                    zunglq_check(&i__1, &i__2, &i__3, &a[(a_dim1 << 1) + 2], lda, &tau[1], &work[1],
                                 &c_n1, &iinfo);
                }
            }
        }
        lwkopt = (integer)work[1].real;
        lwkopt = fla_max(lwkopt, mn);
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZUNGBR", &i__1, (ftnlen)6);
        return LAPACK_FAILURE;
    }
    else if(lquery)
    {
        work[1].real = (double)lwkopt;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUERY_RETURN;
    }
    /* Quick return if possible */
    if(*m == 0 || *n == 0)
    {
        work[1].real = 1.;
        work[1].imag = 0.; // , expr subst
        return LAPACK_QUICK_RETURN;
    }
    return LAPACK_SUCCESS;
}
