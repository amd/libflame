/* dtpsv.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int dtpsv_(char *uplo, char *trans, char *diag, integer *n, doublereal *ap, doublereal *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    integer info;
    doublereal temp;
    integer i__, j, k;
    extern logical lsame_(char *, char *, integer, integer);
    integer kk, ix, jx, kx;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical nounit;
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* DTPSV solves one of the systems of equations */
    /* A*x = b, or A'*x = b, */
    /* where b and x are n element vectors and A is an n by n unit, or */
    /* non-unit, upper or lower triangular matrix, supplied in packed form. */
    /* No test for singularity or near-singularity is included in this */
    /* routine. Such tests must be performed before calling this routine. */
    /* Parameters */
    /* ========== */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the matrix is an upper or */
    /* lower triangular matrix as follows: */
    /* UPLO = 'U' or 'u' A is an upper triangular matrix. */
    /* UPLO = 'L' or 'l' A is a lower triangular matrix. */
    /* Unchanged on exit. */
    /* TRANS - CHARACTER*1. */
    /* On entry, TRANS specifies the equations to be solved as */
    /* follows: */
    /* TRANS = 'N' or 'n' A*x = b. */
    /* TRANS = 'T' or 't' A'*x = b. */
    /* TRANS = 'C' or 'c' A'*x = b. */
    /* Unchanged on exit. */
    /* DIAG - CHARACTER*1. */
    /* On entry, DIAG specifies whether or not A is unit */
    /* triangular as follows: */
    /* DIAG = 'U' or 'u' A is assumed to be unit triangular. */
    /* DIAG = 'N' or 'n' A is not assumed to be unit */
    /* triangular. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the order of the matrix A. */
    /* N must be at least zero. */
    /* Unchanged on exit. */
    /* AP - DOUBLE PRECISION array of DIMENSION at least */
    /* ( ( n*( n + 1 ) )/2 ). */
    /* Before entry with UPLO = 'U' or 'u', the array AP must */
    /* contain the upper triangular matrix packed sequentially, */
    /* column by column, so that AP( 1 ) contains a( 1, 1 ), */
    /* AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 ) */
    /* respectively, and so on. */
    /* Before entry with UPLO = 'L' or 'l', the array AP must */
    /* contain the lower triangular matrix packed sequentially, */
    /* column by column, so that AP( 1 ) contains a( 1, 1 ), */
    /* AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 ) */
    /* respectively, and so on. */
    /* Note that when DIAG = 'U' or 'u', the diagonal elements of */
    /* A are not referenced, but are assumed to be unity. */
    /* Unchanged on exit. */
    /* X - DOUBLE PRECISION array of dimension at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCX ) ). */
    /* Before entry, the incremented array X must contain the n */
    /* element right-hand side vector b. On exit, X is overwritten */
    /* with the solution vector x. */
    /* INCX - INTEGER. */
    /* On entry, INCX specifies the increment for the elements of */
    /* X. INCX must not be zero. */
    /* Unchanged on exit. */
    /* Level 2 Blas routine. */
    /* -- Written on 22-October-1986. */
    /* Jack Dongarra, Argonne National Lab. */
    /* Jeremy Du Croz, Nag Central Office. */
    /* Sven Hammarling, Nag Central Office. */
    /* Richard Hanson, Sandia National Labs. */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. External Functions .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --x;
    --ap;
    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", 1, 1) && ! lsame_(uplo, "L", 1, 1))
    {
        info = 1;
    }
    else if (! lsame_(trans, "N", 1, 1) && ! lsame_(trans, "T", 1, 1) && ! lsame_(trans, "C", 1, 1))
    {
        info = 2;
    }
    else if (! lsame_(diag, "U", 1, 1) && ! lsame_(diag, "N", 1, 1))
    {
        info = 3;
    }
    else if (*n < 0)
    {
        info = 4;
    }
    else if (*incx == 0)
    {
        info = 7;
    }
    if (info != 0)
    {
        xerbla_("DTPSV ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        return 0;
    }
    nounit = lsame_(diag, "N", 1, 1);
    /* Set up the start point in X if the increment is not unity. This */
    /* will be ( N - 1 )*INCX too small for descending loops. */
    if (*incx <= 0)
    {
        kx = 1 - (*n - 1) * *incx;
    }
    else if (*incx != 1)
    {
        kx = 1;
    }
    /* Start the operations. In this version the elements of AP are */
    /* accessed sequentially with one pass through AP. */
    if (lsame_(trans, "N", 1, 1))
    {
        /* Form x := inv( A )*x. */
        if (lsame_(uplo, "U", 1, 1))
        {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1)
            {
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    if (x[j] != 0.)
                    {
                        if (nounit)
                        {
                            x[j] /= ap[kk];
                        }
                        temp = x[j];
                        k = kk - 1;
                        for (i__ = j - 1;
                                i__ >= 1;
                                --i__)
                        {
                            x[i__] -= temp * ap[k];
                            --k;
                            /* L10: */
                        }
                    }
                    kk -= j;
                    /* L20: */
                }
            }
            else
            {
                jx = kx + (*n - 1) * *incx;
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    if (x[jx] != 0.)
                    {
                        if (nounit)
                        {
                            x[jx] /= ap[kk];
                        }
                        temp = x[jx];
                        ix = jx;
                        i__1 = kk - j + 1;
                        for (k = kk - 1;
                                k >= i__1;
                                --k)
                        {
                            ix -= *incx;
                            x[ix] -= temp * ap[k];
                            /* L30: */
                        }
                    }
                    jx -= *incx;
                    kk -= j;
                    /* L40: */
                }
            }
        }
        else
        {
            kk = 1;
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    if (x[j] != 0.)
                    {
                        if (nounit)
                        {
                            x[j] /= ap[kk];
                        }
                        temp = x[j];
                        k = kk + 1;
                        i__2 = *n;
                        for (i__ = j + 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            x[i__] -= temp * ap[k];
                            ++k;
                            /* L50: */
                        }
                    }
                    kk += *n - j + 1;
                    /* L60: */
                }
            }
            else
            {
                jx = kx;
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    if (x[jx] != 0.)
                    {
                        if (nounit)
                        {
                            x[jx] /= ap[kk];
                        }
                        temp = x[jx];
                        ix = jx;
                        i__2 = kk + *n - j;
                        for (k = kk + 1;
                                k <= i__2;
                                ++k)
                        {
                            ix += *incx;
                            x[ix] -= temp * ap[k];
                            /* L70: */
                        }
                    }
                    jx += *incx;
                    kk += *n - j + 1;
                    /* L80: */
                }
            }
        }
    }
    else
    {
        /* Form x := inv( A' )*x. */
        if (lsame_(uplo, "U", 1, 1))
        {
            kk = 1;
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    temp = x[j];
                    k = kk;
                    i__2 = j - 1;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        temp -= ap[k] * x[i__];
                        ++k;
                        /* L90: */
                    }
                    if (nounit)
                    {
                        temp /= ap[kk + j - 1];
                    }
                    x[j] = temp;
                    kk += j;
                    /* L100: */
                }
            }
            else
            {
                jx = kx;
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    temp = x[jx];
                    ix = kx;
                    i__2 = kk + j - 2;
                    for (k = kk;
                            k <= i__2;
                            ++k)
                    {
                        temp -= ap[k] * x[ix];
                        ix += *incx;
                        /* L110: */
                    }
                    if (nounit)
                    {
                        temp /= ap[kk + j - 1];
                    }
                    x[jx] = temp;
                    jx += *incx;
                    kk += j;
                    /* L120: */
                }
            }
        }
        else
        {
            kk = *n * (*n + 1) / 2;
            if (*incx == 1)
            {
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    temp = x[j];
                    k = kk;
                    i__1 = j + 1;
                    for (i__ = *n;
                            i__ >= i__1;
                            --i__)
                    {
                        temp -= ap[k] * x[i__];
                        --k;
                        /* L130: */
                    }
                    if (nounit)
                    {
                        temp /= ap[kk - *n + j];
                    }
                    x[j] = temp;
                    kk -= *n - j + 1;
                    /* L140: */
                }
            }
            else
            {
                kx += (*n - 1) * *incx;
                jx = kx;
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    temp = x[jx];
                    ix = kx;
                    i__1 = kk - (*n - (j + 1));
                    for (k = kk;
                            k >= i__1;
                            --k)
                    {
                        temp -= ap[k] * x[ix];
                        ix -= *incx;
                        /* L150: */
                    }
                    if (nounit)
                    {
                        temp /= ap[kk - *n + j];
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    kk -= *n - j + 1;
                    /* L160: */
                }
            }
        }
    }
    return 0;
    /* End of DTPSV . */
}
/* dtpsv_ */

