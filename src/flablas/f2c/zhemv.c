/* zhemv.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zhemv_(char *uplo, integer *n, dcomplex *alpha, dcomplex *a, integer *lda, dcomplex *x, integer *incx, dcomplex *beta, dcomplex *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    dcomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    integer info;
    dcomplex temp1, temp2;
    integer i__, j;
    extern logical lsame_(char *, char *, integer, integer);
    integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* ZHEMV performs the matrix-vector operation */
    /* y := alpha*A*x + beta*y, */
    /* where alpha and beta are scalars, x and y are n element vectors and */
    /* A is an n by n hermitian matrix. */
    /* Parameters */
    /* ========== */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the upper or lower */
    /* triangular part of the array A is to be referenced as */
    /* follows: */
    /* UPLO = 'U' or 'u' Only the upper triangular part of A */
    /* is to be referenced. */
    /* UPLO = 'L' or 'l' Only the lower triangular part of A */
    /* is to be referenced. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the order of the matrix A. */
    /* N must be at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX*16 . */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* A - COMPLEX*16 array of DIMENSION ( LDA, n ). */
    /* Before entry with UPLO = 'U' or 'u', the leading n by n */
    /* upper triangular part of the array A must contain the upper */
    /* triangular part of the hermitian matrix and the strictly */
    /* lower triangular part of A is not referenced. */
    /* Before entry with UPLO = 'L' or 'l', the leading n by n */
    /* lower triangular part of the array A must contain the lower */
    /* triangular part of the hermitian matrix and the strictly */
    /* upper triangular part of A is not referenced. */
    /* Note that the imaginary parts of the diagonal elements need */
    /* not be set and are assumed to be zero. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. LDA must be at least */
    /* fla_max( 1, n ). */
    /* Unchanged on exit. */
    /* X - COMPLEX*16 array of dimension at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCX ) ). */
    /* Before entry, the incremented array X must contain the n */
    /* element vector x. */
    /* Unchanged on exit. */
    /* INCX - INTEGER. */
    /* On entry, INCX specifies the increment for the elements of */
    /* X. INCX must not be zero. */
    /* Unchanged on exit. */
    /* BETA - COMPLEX*16 . */
    /* On entry, BETA specifies the scalar beta. When BETA is */
    /* supplied as zero then Y need not be set on input. */
    /* Unchanged on exit. */
    /* Y - COMPLEX*16 array of dimension at least */
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    --y;
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
    else if (*lda < fla_max(1,*n))
    {
        info = 5;
    }
    else if (*incx == 0)
    {
        info = 7;
    }
    else if (*incy == 0)
    {
        info = 10;
    }
    if (info != 0)
    {
        xerbla_("ZHEMV ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || alpha->real == 0. && alpha->imag == 0. && (beta->real == 1. && beta->imag == 0.))
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
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through the triangular part */
    /* of A. */
    /* First form y := beta*y. */
    if (beta->real != 1. || beta->imag != 0.)
    {
        if (*incy == 1)
        {
            if (beta->real == 0. && beta->imag == 0.)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = i__;
                    y[i__2].real = 0., y[i__2].imag = 0.;
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
                    z__1.real = beta->real * y[i__3].real - beta->imag * y[i__3].imag, z__1.imag = beta->real * y[i__3].imag + beta->imag * y[i__3] .real;
                    y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
                    /* L20: */
                }
            }
        }
        else
        {
            iy = ky;
            if (beta->real == 0. && beta->imag == 0.)
            {
                i__1 = *n;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    i__2 = iy;
                    y[i__2].real = 0., y[i__2].imag = 0.;
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
                    z__1.real = beta->real * y[i__3].real - beta->imag * y[i__3].imag, z__1.imag = beta->real * y[i__3].imag + beta->imag * y[i__3] .real;
                    y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
                    iy += *incy;
                    /* L40: */
                }
            }
        }
    }
    if (alpha->real == 0. && alpha->imag == 0.)
    {
        return 0;
    }
    if (lsame_(uplo, "U", 1, 1))
    {
        /* Form y when A is stored in upper triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                z__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, z__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = z__1.real, temp1.imag = z__1.imag;
                temp2.real = 0., temp2.imag = 0.;
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    z__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag, z__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5] .real;
                    z__1.real = y[i__4].real + z__2.real, z__1.imag = y[i__4].imag + z__2.imag;
                    y[i__3].real = z__1.real, y[i__3].imag = z__1.imag;
                    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
                    i__3 = i__;
                    z__2.real = z__3.real * x[i__3].real - z__3.imag * x[i__3].imag, z__2.imag = z__3.real * x[i__3].imag + z__3.imag * x[i__3].real;
                    z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                    temp2.real = z__1.real, temp2.imag = z__1.imag;
                    /* L50: */
                }
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                d__1 = a[i__4].real;
                z__3.real = d__1 * temp1.real, z__3.imag = d__1 * temp1.imag;
                z__2.real = y[i__3].real + z__3.real, z__2.imag = y[i__3].imag + z__3.imag;
                z__4.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__4.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                z__1.real = z__2.real + z__4.real, z__1.imag = z__2.imag + z__4.imag;
                y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
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
                z__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, z__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = z__1.real, temp1.imag = z__1.imag;
                temp2.real = 0., temp2.imag = 0.;
                ix = kx;
                iy = ky;
                i__2 = j - 1;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    z__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag, z__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5] .real;
                    z__1.real = y[i__4].real + z__2.real, z__1.imag = y[i__4].imag + z__2.imag;
                    y[i__3].real = z__1.real, y[i__3].imag = z__1.imag;
                    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
                    i__3 = ix;
                    z__2.real = z__3.real * x[i__3].real - z__3.imag * x[i__3].imag, z__2.imag = z__3.real * x[i__3].imag + z__3.imag * x[i__3].real;
                    z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                    temp2.real = z__1.real, temp2.imag = z__1.imag;
                    ix += *incx;
                    iy += *incy;
                    /* L70: */
                }
                i__2 = jy;
                i__3 = jy;
                i__4 = j + j * a_dim1;
                d__1 = a[i__4].real;
                z__3.real = d__1 * temp1.real, z__3.imag = d__1 * temp1.imag;
                z__2.real = y[i__3].real + z__3.real, z__2.imag = y[i__3].imag + z__3.imag;
                z__4.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__4.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                z__1.real = z__2.real + z__4.real, z__1.imag = z__2.imag + z__4.imag;
                y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
                jx += *incx;
                jy += *incy;
                /* L80: */
            }
        }
    }
    else
    {
        /* Form y when A is stored in lower triangle. */
        if (*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                z__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, z__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = z__1.real, temp1.imag = z__1.imag;
                temp2.real = 0., temp2.imag = 0.;
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                d__1 = a[i__4].real;
                z__2.real = d__1 * temp1.real, z__2.imag = d__1 * temp1.imag;
                z__1.real = y[i__3].real + z__2.real, z__1.imag = y[i__3].imag + z__2.imag;
                y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    z__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag, z__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5] .real;
                    z__1.real = y[i__4].real + z__2.real, z__1.imag = y[i__4].imag + z__2.imag;
                    y[i__3].real = z__1.real, y[i__3].imag = z__1.imag;
                    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
                    i__3 = i__;
                    z__2.real = z__3.real * x[i__3].real - z__3.imag * x[i__3].imag, z__2.imag = z__3.real * x[i__3].imag + z__3.imag * x[i__3].real;
                    z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                    temp2.real = z__1.real, temp2.imag = z__1.imag;
                    /* L90: */
                }
                i__2 = j;
                i__3 = j;
                z__2.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__2.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                z__1.real = y[i__3].real + z__2.real, z__1.imag = y[i__3].imag + z__2.imag;
                y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
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
                z__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, z__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real;
                temp1.real = z__1.real, temp1.imag = z__1.imag;
                temp2.real = 0., temp2.imag = 0.;
                i__2 = jy;
                i__3 = jy;
                i__4 = j + j * a_dim1;
                d__1 = a[i__4].real;
                z__2.real = d__1 * temp1.real, z__2.imag = d__1 * temp1.imag;
                z__1.real = y[i__3].real + z__2.real, z__1.imag = y[i__3].imag + z__2.imag;
                y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
                ix = jx;
                iy = jy;
                i__2 = *n;
                for (i__ = j + 1;
                        i__ <= i__2;
                        ++i__)
                {
                    ix += *incx;
                    iy += *incy;
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    z__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag, z__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5] .real;
                    z__1.real = y[i__4].real + z__2.real, z__1.imag = y[i__4].imag + z__2.imag;
                    y[i__3].real = z__1.real, y[i__3].imag = z__1.imag;
                    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
                    i__3 = ix;
                    z__2.real = z__3.real * x[i__3].real - z__3.imag * x[i__3].imag, z__2.imag = z__3.real * x[i__3].imag + z__3.imag * x[i__3].real;
                    z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                    temp2.real = z__1.real, temp2.imag = z__1.imag;
                    /* L110: */
                }
                i__2 = jy;
                i__3 = jy;
                z__2.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__2.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                z__1.real = y[i__3].real + z__2.real, z__1.imag = y[i__3].imag + z__2.imag;
                y[i__2].real = z__1.real, y[i__2].imag = z__1.imag;
                jx += *incx;
                jy += *incy;
                /* L120: */
            }
        }
    }
    return 0;
    /* End of ZHEMV . */
}
/* zhemv_ */

