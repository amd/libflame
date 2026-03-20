/* cgeru.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int cgeru_(integer *m, integer *n, scomplex *alpha, scomplex * x, integer *incx, scomplex *y, integer *incy, scomplex *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    scomplex q__1, q__2;
    /* Local variables */
    integer info;
    scomplex temp;
    integer i__, j, ix, jy, kx;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* CGERU performs the rank 1 operation */
    /* A := alpha*x*y' + A, */
    /* where alpha is a scalar, x is an m element vector, y is an n element */
    /* vector and A is an m by n matrix. */
    /* Parameters */
    /* ========== */
    /* M - INTEGER. */
    /* On entry, M specifies the number of rows of the matrix A. */
    /* M must be at least zero. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the number of columns of the matrix A. */
    /* N must be at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX . */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* X - COMPLEX array of dimension at least */
    /* ( 1 + ( m - 1 )*f2c_abs( INCX ) ). */
    /* Before entry, the incremented array X must contain the m */
    /* element vector x. */
    /* Unchanged on exit. */
    /* INCX - INTEGER. */
    /* On entry, INCX specifies the increment for the elements of */
    /* X. INCX must not be zero. */
    /* Unchanged on exit. */
    /* Y - COMPLEX array of dimension at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCY ) ). */
    /* Before entry, the incremented array Y must contain the n */
    /* element vector y. */
    /* Unchanged on exit. */
    /* INCY - INTEGER. */
    /* On entry, INCY specifies the increment for the elements of */
    /* Y. INCY must not be zero. */
    /* Unchanged on exit. */
    /* A - COMPLEX array of DIMENSION ( LDA, n ). */
    /* Before entry, the leading m by n part of the array A must */
    /* contain the matrix of coefficients. On exit, A is */
    /* overwritten by the updated matrix. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. LDA must be at least */
    /* fla_max( 1, m ). */
    /* Unchanged on exit. */
    /* Level 2 Blas routine. */
    /* -- Written on 22-October-1986. */
    /* Jack Dongarra, Argonne National Lab. */
    /* Jeremy Du Croz, Nag Central Office. */
    /* Sven Hammarling, Nag Central Office. */
    /* Richard Hanson, Sandia National Labs. */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. External Subroutines .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    /* Function Body */
    info = 0;
    if (*m < 0)
    {
        info = 1;
    }
    else if (*n < 0)
    {
        info = 2;
    }
    else if (*incx == 0)
    {
        info = 5;
    }
    else if (*incy == 0)
    {
        info = 7;
    }
    else if (*lda < fla_max(1,*m))
    {
        info = 9;
    }
    if (info != 0)
    {
        xerbla_("CGERU ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*m == 0 || *n == 0 || alpha->real == 0.f && alpha->imag == 0.f)
    {
        return 0;
    }
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through A. */
    if (*incy > 0)
    {
        jy = 1;
    }
    else
    {
        jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1)
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = jy;
            if (y[i__2].real != 0.f || y[i__2].imag != 0.f)
            {
                i__2 = jy;
                q__1.real = alpha->real * y[i__2].real - alpha->imag * y[i__2].imag, q__1.imag = alpha->real * y[i__2].imag + alpha->imag * y[i__2].real;
                temp.real = q__1.real, temp.imag = q__1.imag;
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__ + j * a_dim1;
                    i__5 = i__;
                    q__2.real = x[i__5].real * temp.real - x[i__5].imag * temp.imag, q__2.imag = x[i__5].real * temp.imag + x[i__5].imag * temp.real;
                    q__1.real = a[i__4].real + q__2.real, q__1.imag = a[i__4].imag + q__2.imag;
                    a[i__3].real = q__1.real, a[i__3].imag = q__1.imag;
                    /* L10: */
                }
            }
            jy += *incy;
            /* L20: */
        }
    }
    else
    {
        if (*incx > 0)
        {
            kx = 1;
        }
        else
        {
            kx = 1 - (*m - 1) * *incx;
        }
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = jy;
            if (y[i__2].real != 0.f || y[i__2].imag != 0.f)
            {
                i__2 = jy;
                q__1.real = alpha->real * y[i__2].real - alpha->imag * y[i__2].imag, q__1.imag = alpha->real * y[i__2].imag + alpha->imag * y[i__2].real;
                temp.real = q__1.real, temp.imag = q__1.imag;
                ix = kx;
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__ + j * a_dim1;
                    i__5 = ix;
                    q__2.real = x[i__5].real * temp.real - x[i__5].imag * temp.imag, q__2.imag = x[i__5].real * temp.imag + x[i__5].imag * temp.real;
                    q__1.real = a[i__4].real + q__2.real, q__1.imag = a[i__4].imag + q__2.imag;
                    a[i__3].real = q__1.real, a[i__3].imag = q__1.imag;
                    ix += *incx;
                    /* L30: */
                }
            }
            jy += *incy;
            /* L40: */
        }
    }
    return 0;
    /* End of CGERU . */
}
/* cgeru_ */

