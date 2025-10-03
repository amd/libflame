/* chpmv.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int chpmv_(char *uplo, integer *n, scomplex *alpha, scomplex * ap, scomplex *x, integer *incx, scomplex *beta, scomplex *y, integer * incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1;
    scomplex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    integer info;
    scomplex temp1, temp2;
    integer i__, j, k;
    extern logical lsame_(char *, char *, integer, integer);
    integer kk, ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* CHPMV performs the matrix-vector operation */
    /* y := alpha*A*x + beta*y, */
    /* where alpha and beta are scalars, x and y are n element vectors and */
    /* A is an n by n hermitian matrix, supplied in packed form. */
    /* Parameters */
    /* ========== */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the upper or lower */
    /* triangular part of the matrix A is supplied in the packed */
    /* array AP as follows: */
    /* UPLO = 'U' or 'u' The upper triangular part of A is */
    /* supplied in AP. */
    /* UPLO = 'L' or 'l' The lower triangular part of A is */
    /* supplied in AP. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the order of the matrix A. */
    /* N must be at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX . */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* AP - COMPLEX array of DIMENSION at least */
    /* ( ( n*( n + 1 ) )/2 ). */
    /* Before entry with UPLO = 'U' or 'u', the array AP must */
    /* contain the upper triangular part of the hermitian matrix */
    /* packed sequentially, column by column, so that AP( 1 ) */
    /* contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
    /* and a( 2, 2 ) respectively, and so on. */
    /* Before entry with UPLO = 'L' or 'l', the array AP must */
    /* contain the lower triangular part of the hermitian matrix */
    /* packed sequentially, column by column, so that AP( 1 ) */
    /* contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
    /* and a( 3, 1 ) respectively, and so on. */
    /* Note that the imaginary parts of the diagonal elements need */
    /* not be set and are assumed to be zero. */
    /* Unchanged on exit. */
    /* X - COMPLEX array of dimension at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCX ) ). */
    /* Before entry, the incremented array X must contain the n */
    /* element vector x. */
    /* Unchanged on exit. */
    /* INCX - INTEGER. */
    /* On entry, INCX specifies the increment for the elements of */
    /* X. INCX must not be zero. */
    /* Unchanged on exit. */
    /* BETA - COMPLEX . */
    /* On entry, BETA specifies the scalar beta. When BETA is */
    /* supplied as zero then Y need not be set on input. */
    /* Unchanged on exit. */
    /* Y - COMPLEX array of dimension at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCY ) ). */
    /* Before entry, the incremented array Y must contain the n */
    /* element vector y. On exit, Y is overwritten by the updated */
    /* vector y. */
    /* INCY - INTEGER. */
    /* On entry, INCY specifies the increment for the elements of */
    /* Y. INCY must not be zero. */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --y;
    --x;
    --ap;
    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", 1, 1) && ! lsame_(uplo, "L", 1, 1))
    {
        info = 1;
    }
    else if (*n < 0)
    {
        info = 2;
    }
    else if (*incx == 0)
    {
        info = 6;
    }
    else if (*incy == 0)
    {
        info = 9;
    }
    if (info != 0)
    {
        xerbla_("CHPMV ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || alpha->real == 0.f && alpha->imag == 0.f && (beta->real == 1.f && beta->imag == 0.f))
    {
        return 0;
    }
    /* Set up the start points in X and Y. */
    if (*incx > 0)
    {
        kx = 1;
    }
    else
    {
        kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0)
    {
        ky = 1;
    }
    else
    {
        ky = 1 - (*n - 1) * *incy;
    }
    /* Start the operations. In this version the elements of the array AP */
    /* are accessed sequentially with one pass through AP. */
    /* First form y := beta*y. */
    if (beta->real != 1.f || beta->imag != 0.f)
    {
        if (*incy == 1)
        {
            if (beta->real == 0.f && beta->imag == 0.f)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    y[i__2].real = 0.f, y[i__2].imag = 0.f;
                    /* L10: */
                }
            }
            else
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    i__3 = i__;
                    q__1.real = beta->real * y[i__3].real - beta->imag * y[i__3].imag, q__1.imag = beta->real * y[i__3].imag + beta->imag * y[i__3] .real;
                    y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                    /* L20: */
                }
            }
        }
        else
        {
            iy = ky;
            if (beta->real == 0.f && beta->imag == 0.f)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = iy;
                    y[i__2].real = 0.f, y[i__2].imag = 0.f;
                    iy += *incy;
                    /* L30: */
                }
            }
            else
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = iy;
                    i__3 = iy;
                    q__1.real = beta->real * y[i__3].real - beta->imag * y[i__3].imag, q__1.imag = beta->real * y[i__3].imag + beta->imag * y[i__3] .real;
                    y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                    iy += *incy;
                    /* L40: */
                }
            }
        }
    }
    if (alpha->real == 0.f && alpha->imag == 0.f)
    {
        return 0;
    }
    kk = 1;
    if (lsame_(uplo, "U", 1, 1))
    {
        /* Form y when AP contains the upper triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = q__1.real, temp1.imag = q__1.imag;
                temp2.real = 0.f, temp2.imag = 0.f;
                k = kk;
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = k;
                    q__2.real = temp1.real * ap[i__5].real - temp1.imag * ap[i__5].imag, q__2.imag = temp1.real * ap[i__5].imag + temp1.imag * ap[i__5] .real;
                    q__1.real = y[i__4].real + q__2.real, q__1.imag = y[i__4].imag + q__2.imag;
                    y[i__3].real = q__1.real, y[i__3].imag = q__1.imag;
                    r_cnjg(&q__3, &ap[k]);
                    i__3 = i__;
                    q__2.real = q__3.real * x[i__3].real - q__3.imag * x[i__3].imag, q__2.imag = q__3.real * x[i__3].imag + q__3.imag * x[i__3].real;
                    q__1.real = temp2.real + q__2.real, q__1.imag = temp2.imag + q__2.imag;
                    temp2.real = q__1.real, temp2.imag = q__1.imag;
                    ++k;
                    /* L50: */
                }
                i__2 = j;
                i__3 = j;
                i__4 = kk + j - 1;
                r__1 = ap[i__4].real;
                q__3.real = r__1 * temp1.real, q__3.imag = r__1 * temp1.imag;
                q__2.real = y[i__3].real + q__3.real, q__2.imag = y[i__3].imag + q__3.imag;
                q__4.real = alpha->real * temp2.real - alpha->imag * temp2.imag, q__4.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                q__1.real = q__2.real + q__4.real, q__1.imag = q__2.imag + q__4.imag;
                y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                kk += j;
                /* L60: */
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = jx;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = q__1.real, temp1.imag = q__1.imag;
                temp2.real = 0.f, temp2.imag = 0.f;
                ix = kx;
                iy = ky;
                i__2 = kk + j - 2;
                for (k = kk;
                        k <= i__2;
                        ++k)
                {
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = k;
                    q__2.real = temp1.real * ap[i__5].real - temp1.imag * ap[i__5].imag, q__2.imag = temp1.real * ap[i__5].imag + temp1.imag * ap[i__5] .real;
                    q__1.real = y[i__4].real + q__2.real, q__1.imag = y[i__4].imag + q__2.imag;
                    y[i__3].real = q__1.real, y[i__3].imag = q__1.imag;
                    r_cnjg(&q__3, &ap[k]);
                    i__3 = ix;
                    q__2.real = q__3.real * x[i__3].real - q__3.imag * x[i__3].imag, q__2.imag = q__3.real * x[i__3].imag + q__3.imag * x[i__3].real;
                    q__1.real = temp2.real + q__2.real, q__1.imag = temp2.imag + q__2.imag;
                    temp2.real = q__1.real, temp2.imag = q__1.imag;
                    ix += *incx;
                    iy += *incy;
                    /* L70: */
                }
                i__2 = jy;
                i__3 = jy;
                i__4 = kk + j - 1;
                r__1 = ap[i__4].real;
                q__3.real = r__1 * temp1.real, q__3.imag = r__1 * temp1.imag;
                q__2.real = y[i__3].real + q__3.real, q__2.imag = y[i__3].imag + q__3.imag;
                q__4.real = alpha->real * temp2.real - alpha->imag * temp2.imag, q__4.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                q__1.real = q__2.real + q__4.real, q__1.imag = q__2.imag + q__4.imag;
                y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                jx += *incx;
                jy += *incy;
                kk += j;
                /* L80: */
            }
        }
    }
    else
    {
        /* Form y when AP contains the lower triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = q__1.real, temp1.imag = q__1.imag;
                temp2.real = 0.f, temp2.imag = 0.f;
                i__2 = j;
                i__3 = j;
                i__4 = kk;
                r__1 = ap[i__4].real;
                q__2.real = r__1 * temp1.real, q__2.imag = r__1 * temp1.imag;
                q__1.real = y[i__3].real + q__2.real, q__1.imag = y[i__3].imag + q__2.imag;
                y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                k = kk + 1;
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = k;
                    q__2.real = temp1.real * ap[i__5].real - temp1.imag * ap[i__5].imag, q__2.imag = temp1.real * ap[i__5].imag + temp1.imag * ap[i__5] .real;
                    q__1.real = y[i__4].real + q__2.real, q__1.imag = y[i__4].imag + q__2.imag;
                    y[i__3].real = q__1.real, y[i__3].imag = q__1.imag;
                    r_cnjg(&q__3, &ap[k]);
                    i__3 = i__;
                    q__2.real = q__3.real * x[i__3].real - q__3.imag * x[i__3].imag, q__2.imag = q__3.real * x[i__3].imag + q__3.imag * x[i__3].real;
                    q__1.real = temp2.real + q__2.real, q__1.imag = temp2.imag + q__2.imag;
                    temp2.real = q__1.real, temp2.imag = q__1.imag;
                    ++k;
                    /* L90: */
                }
                i__2 = j;
                i__3 = j;
                q__2.real = alpha->real * temp2.real - alpha->imag * temp2.imag, q__2.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                q__1.real = y[i__3].real + q__2.real, q__1.imag = y[i__3].imag + q__2.imag;
                y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                kk += *n - j + 1;
                /* L100: */
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = jx;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = q__1.real, temp1.imag = q__1.imag;
                temp2.real = 0.f, temp2.imag = 0.f;
                i__2 = jy;
                i__3 = jy;
                i__4 = kk;
                r__1 = ap[i__4].real;
                q__2.real = r__1 * temp1.real, q__2.imag = r__1 * temp1.imag;
                q__1.real = y[i__3].real + q__2.real, q__1.imag = y[i__3].imag + q__2.imag;
                y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                ix = jx;
                iy = jy;
                i__2 = kk + *n - j;
                for (k = kk + 1;
                        k <= i__2;
                        ++k)
                {
                    ix += *incx;
                    iy += *incy;
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = k;
                    q__2.real = temp1.real * ap[i__5].real - temp1.imag * ap[i__5].imag, q__2.imag = temp1.real * ap[i__5].imag + temp1.imag * ap[i__5] .real;
                    q__1.real = y[i__4].real + q__2.real, q__1.imag = y[i__4].imag + q__2.imag;
                    y[i__3].real = q__1.real, y[i__3].imag = q__1.imag;
                    r_cnjg(&q__3, &ap[k]);
                    i__3 = ix;
                    q__2.real = q__3.real * x[i__3].real - q__3.imag * x[i__3].imag, q__2.imag = q__3.real * x[i__3].imag + q__3.imag * x[i__3].real;
                    q__1.real = temp2.real + q__2.real, q__1.imag = temp2.imag + q__2.imag;
                    temp2.real = q__1.real, temp2.imag = q__1.imag;
                    /* L110: */
                }
                i__2 = jy;
                i__3 = jy;
                q__2.real = alpha->real * temp2.real - alpha->imag * temp2.imag, q__2.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                q__1.real = y[i__3].real + q__2.real, q__1.imag = y[i__3].imag + q__2.imag;
                y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                jx += *incx;
                jy += *incy;
                kk += *n - j + 1;
                /* L120: */
            }
        }
    }
    return 0;
    /* End of CHPMV . */
}
/* chpmv_ */

