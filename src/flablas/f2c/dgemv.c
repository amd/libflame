/* dgemv.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int dgemv_(char *trans, integer *m, integer *n, doublereal * alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    integer info;
    doublereal temp;
    integer lenx, leny, i__, j;
    extern logical lsame_(char *, char *, integer, integer);
    integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* DGEMV performs one of the matrix-vector operations */
    /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y, */
    /* where alpha and beta are scalars, x and y are vectors and A is an */
    /* m by n matrix. */
    /* Parameters */
    /* ========== */
    /* TRANS - CHARACTER*1. */
    /* On entry, TRANS specifies the operation to be performed as */
    /* follows: */
    /* TRANS = 'N' or 'n' y := alpha*A*x + beta*y. */
    /* TRANS = 'T' or 't' y := alpha*A'*x + beta*y. */
    /* TRANS = 'C' or 'c' y := alpha*A'*x + beta*y. */
    /* Unchanged on exit. */
    /* M - INTEGER. */
    /* On entry, M specifies the number of rows of the matrix A. */
    /* M must be at least zero. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the number of columns of the matrix A. */
    /* N must be at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - DOUBLE PRECISION. */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* A - DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
    /* Before entry, the leading m by n part of the array A must */
    /* contain the matrix of coefficients. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. LDA must be at least */
    /* fla_max( 1, m ). */
    /* Unchanged on exit. */
    /* X - DOUBLE PRECISION array of DIMENSION at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCX ) ) when TRANS = 'N' or 'n' */
    /* and at least */
    /* ( 1 + ( m - 1 )*f2c_abs( INCX ) ) otherwise. */
    /* Before entry, the incremented array X must contain the */
    /* vector x. */
    /* Unchanged on exit. */
    /* INCX - INTEGER. */
    /* On entry, INCX specifies the increment for the elements of */
    /* X. INCX must not be zero. */
    /* Unchanged on exit. */
    /* BETA - DOUBLE PRECISION. */
    /* On entry, BETA specifies the scalar beta. When BETA is */
    /* supplied as zero then Y need not be set on input. */
    /* Unchanged on exit. */
    /* Y - DOUBLE PRECISION array of DIMENSION at least */
    /* ( 1 + ( m - 1 )*f2c_abs( INCY ) ) when TRANS = 'N' or 'n' */
    /* and at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCY ) ) otherwise. */
    /* Before entry with BETA non-zero, the incremented array Y */
    /* must contain the vector y. On exit, Y is overwritten by the */
    /* updated vector y. */
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
    if (! lsame_(trans, "N", 1, 1) && ! lsame_(trans, "T", 1, 1) && ! lsame_(trans, "C", 1, 1) )
    {
        info = 1;
    }
    else if (*m < 0)
    {
        info = 2;
    }
    else if (*n < 0)
    {
        info = 3;
    }
    else if (*lda < fla_max(1,*m))
    {
        info = 6;
    }
    else if (*incx == 0)
    {
        info = 8;
    }
    else if (*incy == 0)
    {
        info = 11;
    }
    if (info != 0)
    {
        xerbla_("DGEMV ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.)
    {
        return 0;
    }
    /* Set LENX and LENY, the lengths of the vectors x and y, and set */
    /* up the start points in X and Y. */
    if (lsame_(trans, "N", 1, 1))
    {
        lenx = *n;
        leny = *m;
    }
    else
    {
        lenx = *m;
        leny = *n;
    }
    if (*incx > 0)
    {
        kx = 1;
    }
    else
    {
        kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0)
    {
        ky = 1;
    }
    else
    {
        ky = 1 - (leny - 1) * *incy;
    }
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through A. */
    /* First form y := beta*y. */
    if (*beta != 1.)
    {
        if (*incy == 1)
        {
            if (*beta == 0.)
            {
                i__1 = leny;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    y[i__] = 0.;
                    /* L10: */
                }
            }
            else
            {
                i__1 = leny;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    y[i__] = *beta * y[i__];
                    /* L20: */
                }
            }
        }
        else
        {
            iy = ky;
            if (*beta == 0.)
            {
                i__1 = leny;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    y[iy] = 0.;
                    iy += *incy;
                    /* L30: */
                }
            }
            else
            {
                i__1 = leny;
                for (i__ = 1;
                        i__ <= i__1;
                        ++i__)
                {
                    y[iy] = *beta * y[iy];
                    iy += *incy;
                    /* L40: */
                }
            }
        }
    }
    if (*alpha == 0.)
    {
        return 0;
    }
    if (lsame_(trans, "N", 1, 1))
    {
        /* Form y := alpha*A*x + y. */
        jx = kx;
        if (*incy == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (x[jx] != 0.)
                {
                    temp = *alpha * x[jx];
                    i__2 = *m;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        y[i__] += temp * a[i__ + j * a_dim1];
                        /* L50: */
                    }
                }
                jx += *incx;
                /* L60: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (x[jx] != 0.)
                {
                    temp = *alpha * x[jx];
                    iy = ky;
                    i__2 = *m;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        y[iy] += temp * a[i__ + j * a_dim1];
                        iy += *incy;
                        /* L70: */
                    }
                }
                jx += *incx;
                /* L80: */
            }
        }
    }
    else
    {
        /* Form y := alpha*A'*x + y. */
        jy = ky;
        if (*incx == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                temp = 0.;
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    temp += a[i__ + j * a_dim1] * x[i__];
                    /* L90: */
                }
                y[jy] += *alpha * temp;
                jy += *incy;
                /* L100: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                temp = 0.;
                ix = kx;
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    temp += a[i__ + j * a_dim1] * x[ix];
                    ix += *incx;
                    /* L110: */
                }
                y[jy] += *alpha * temp;
                jy += *incy;
                /* L120: */
            }
        }
    }
    return 0;
    /* End of DGEMV . */
}
/* dgemv_ */

