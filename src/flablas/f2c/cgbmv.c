/* cgbmv.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int cgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, scomplex *alpha, scomplex *a, integer *lda, scomplex *x, integer *incx, scomplex *beta, scomplex *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    integer info;
    scomplex temp;
    integer lenx, leny, i__, j, k;
    extern logical lsame_(char *, char *, integer, integer);
    integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical noconj;
    integer kup1;
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* CGBMV performs one of the matrix-vector operations */
    /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y, or */
    /* y := alpha*conjg( A' )*x + beta*y, */
    /* where alpha and beta are scalars, x and y are vectors and A is an */
    /* m by n band matrix, with kl sub-diagonals and ku super-diagonals. */
    /* Parameters */
    /* ========== */
    /* TRANS - CHARACTER*1. */
    /* On entry, TRANS specifies the operation to be performed as */
    /* follows: */
    /* TRANS = 'N' or 'n' y := alpha*A*x + beta*y. */
    /* TRANS = 'T' or 't' y := alpha*A'*x + beta*y. */
    /* TRANS = 'C' or 'c' y := alpha*conjg( A' )*x + beta*y. */
    /* Unchanged on exit. */
    /* M - INTEGER. */
    /* On entry, M specifies the number of rows of the matrix A. */
    /* M must be at least zero. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the number of columns of the matrix A. */
    /* N must be at least zero. */
    /* Unchanged on exit. */
    /* KL - INTEGER. */
    /* On entry, KL specifies the number of sub-diagonals of the */
    /* matrix A. KL must satisfy 0 .le. KL. */
    /* Unchanged on exit. */
    /* KU - INTEGER. */
    /* On entry, KU specifies the number of super-diagonals of the */
    /* matrix A. KU must satisfy 0 .le. KU. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX . */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* A - COMPLEX array of DIMENSION ( LDA, n ). */
    /* Before entry, the leading ( kl + ku + 1 ) by n part of the */
    /* array A must contain the matrix of coefficients, supplied */
    /* column by column, with the leading diagonal of the matrix in */
    /* row ( ku + 1 ) of the array, the first super-diagonal */
    /* starting at position 2 in row ku, the first sub-diagonal */
    /* starting at position 1 in row ( ku + 2 ), and so on. */
    /* Elements in the array A that do not correspond to elements */
    /* in the band matrix (such as the top left ku by ku triangle) */
    /* are not referenced. */
    /* The following program segment will transfer a band matrix */
    /* from conventional full matrix storage to band storage: */
    /* DO 20, J = 1, N */
    /* K = KU + 1 - J */
    /* DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL ) */
    /* A( K + I, J ) = matrix( I, J ) */
    /* 10 CONTINUE */
    /* 20 CONTINUE */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. LDA must be at least */
    /* ( kl + ku + 1 ). */
    /* Unchanged on exit. */
    /* X - COMPLEX array of DIMENSION at least */
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
    /* BETA - COMPLEX . */
    /* On entry, BETA specifies the scalar beta. When BETA is */
    /* supplied as zero then Y need not be set on input. */
    /* Unchanged on exit. */
    /* Y - COMPLEX array of DIMENSION at least */
    /* ( 1 + ( m - 1 )*f2c_abs( INCY ) ) when TRANS = 'N' or 'n' */
    /* and at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCY ) ) otherwise. */
    /* Before entry, the incremented array Y must contain the */
    /* vector y. On exit, Y is overwritten by the updated vector y. */
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
    else if (*kl < 0)
    {
        info = 4;
    }
    else if (*ku < 0)
    {
        info = 5;
    }
    else if (*lda < *kl + *ku + 1)
    {
        info = 8;
    }
    else if (*incx == 0)
    {
        info = 10;
    }
    else if (*incy == 0)
    {
        info = 13;
    }
    if (info != 0)
    {
        xerbla_("CGBMV ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*m == 0 || *n == 0 || alpha->real == 0.f && alpha->imag == 0.f && (beta->real == 1.f && beta->imag == 0.f))
    {
        return 0;
    }
    noconj = lsame_(trans, "T", 1, 1);
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
    /* accessed sequentially with one pass through the band part of A. */
    /* First form y := beta*y. */
    if (beta->real != 1.f || beta->imag != 0.f)
    {
        if (*incy == 1)
        {
            if (beta->real == 0.f && beta->imag == 0.f)
            {
                i__1 = leny;
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
                i__1 = leny;
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
                i__1 = leny;
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
                i__1 = leny;
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
    kup1 = *ku + 1;
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
                i__2 = jx;
                if (x[i__2].real != 0.f || x[i__2].imag != 0.f)
                {
                    i__2 = jx;
                    q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag, q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2] .real;
                    temp.real = q__1.real, temp.imag = q__1.imag;
                    k = kup1 - j;
                    /* Computing MAX */
                    i__2 = 1, i__3 = j - *ku;
                    /* Computing MIN */
                    i__5 = *m, i__6 = j + *kl;
                    i__4 = fla_min(i__5,i__6);
                    for (i__ = fla_max(i__2,i__3);
                            i__ <= i__4;
                            ++i__)
                    {
                        i__2 = i__;
                        i__3 = i__;
                        i__5 = k + i__ + j * a_dim1;
                        q__2.real = temp.real * a[i__5].real - temp.imag * a[i__5].imag, q__2.imag = temp.real * a[i__5].imag + temp.imag * a[i__5] .real;
                        q__1.real = y[i__3].real + q__2.real, q__1.imag = y[i__3].imag + q__2.imag;
                        y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
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
                i__4 = jx;
                if (x[i__4].real != 0.f || x[i__4].imag != 0.f)
                {
                    i__4 = jx;
                    q__1.real = alpha->real * x[i__4].real - alpha->imag * x[i__4].imag, q__1.imag = alpha->real * x[i__4].imag + alpha->imag * x[i__4] .real;
                    temp.real = q__1.real, temp.imag = q__1.imag;
                    iy = ky;
                    k = kup1 - j;
                    /* Computing MAX */
                    i__4 = 1, i__2 = j - *ku;
                    /* Computing MIN */
                    i__5 = *m, i__6 = j + *kl;
                    i__3 = fla_min(i__5,i__6);
                    for (i__ = fla_max(i__4,i__2);
                            i__ <= i__3;
                            ++i__)
                    {
                        i__4 = iy;
                        i__2 = iy;
                        i__5 = k + i__ + j * a_dim1;
                        q__2.real = temp.real * a[i__5].real - temp.imag * a[i__5].imag, q__2.imag = temp.real * a[i__5].imag + temp.imag * a[i__5] .real;
                        q__1.real = y[i__2].real + q__2.real, q__1.imag = y[i__2].imag + q__2.imag;
                        y[i__4].real = q__1.real, y[i__4].imag = q__1.imag;
                        iy += *incy;
                        /* L70: */
                    }
                }
                jx += *incx;
                if (j > *ku)
                {
                    ky += *incy;
                }
                /* L80: */
            }
        }
    }
    else
    {
        /* Form y := alpha*A'*x + y or y := alpha*conjg( A' )*x + y. */
        jy = ky;
        if (*incx == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                temp.real = 0.f, temp.imag = 0.f;
                k = kup1 - j;
                if (noconj)
                {
                    /* Computing MAX */
                    i__3 = 1, i__4 = j - *ku;
                    /* Computing MIN */
                    i__5 = *m, i__6 = j + *kl;
                    i__2 = fla_min(i__5,i__6);
                    for (i__ = fla_max(i__3,i__4);
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = k + i__ + j * a_dim1;
                        i__4 = i__;
                        q__2.real = a[i__3].real * x[i__4].real - a[i__3].imag * x[i__4] .imag, q__2.imag = a[i__3].real * x[i__4].imag + a[i__3] .imag * x[i__4].real;
                        q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                        temp.real = q__1.real, temp.imag = q__1.imag;
                        /* L90: */
                    }
                }
                else
                {
                    /* Computing MAX */
                    i__2 = 1, i__3 = j - *ku;
                    /* Computing MIN */
                    i__5 = *m, i__6 = j + *kl;
                    i__4 = fla_min(i__5,i__6);
                    for (i__ = fla_max(i__2,i__3);
                            i__ <= i__4;
                            ++i__)
                    {
                        r_cnjg(&q__3, &a[k + i__ + j * a_dim1]);
                        i__2 = i__;
                        q__2.real = q__3.real * x[i__2].real - q__3.imag * x[i__2].imag, q__2.imag = q__3.real * x[i__2].imag + q__3.imag * x[i__2] .real;
                        q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                        temp.real = q__1.real, temp.imag = q__1.imag;
                        /* L100: */
                    }
                }
                i__4 = jy;
                i__2 = jy;
                q__2.real = alpha->real * temp.real - alpha->imag * temp.imag, q__2.imag = alpha->real * temp.imag + alpha->imag * temp.real;
                q__1.real = y[i__2].real + q__2.real, q__1.imag = y[i__2].imag + q__2.imag;
                y[i__4].real = q__1.real, y[i__4].imag = q__1.imag;
                jy += *incy;
                /* L110: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                temp.real = 0.f, temp.imag = 0.f;
                ix = kx;
                k = kup1 - j;
                if (noconj)
                {
                    /* Computing MAX */
                    i__4 = 1, i__2 = j - *ku;
                    /* Computing MIN */
                    i__5 = *m, i__6 = j + *kl;
                    i__3 = fla_min(i__5,i__6);
                    for (i__ = fla_max(i__4,i__2);
                            i__ <= i__3;
                            ++i__)
                    {
                        i__4 = k + i__ + j * a_dim1;
                        i__2 = ix;
                        q__2.real = a[i__4].real * x[i__2].real - a[i__4].imag * x[i__2] .imag, q__2.imag = a[i__4].real * x[i__2].imag + a[i__4] .imag * x[i__2].real;
                        q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                        temp.real = q__1.real, temp.imag = q__1.imag;
                        ix += *incx;
                        /* L120: */
                    }
                }
                else
                {
                    /* Computing MAX */
                    i__3 = 1, i__4 = j - *ku;
                    /* Computing MIN */
                    i__5 = *m, i__6 = j + *kl;
                    i__2 = fla_min(i__5,i__6);
                    for (i__ = fla_max(i__3,i__4);
                            i__ <= i__2;
                            ++i__)
                    {
                        r_cnjg(&q__3, &a[k + i__ + j * a_dim1]);
                        i__3 = ix;
                        q__2.real = q__3.real * x[i__3].real - q__3.imag * x[i__3].imag, q__2.imag = q__3.real * x[i__3].imag + q__3.imag * x[i__3] .real;
                        q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                        temp.real = q__1.real, temp.imag = q__1.imag;
                        ix += *incx;
                        /* L130: */
                    }
                }
                i__2 = jy;
                i__3 = jy;
                q__2.real = alpha->real * temp.real - alpha->imag * temp.imag, q__2.imag = alpha->real * temp.imag + alpha->imag * temp.real;
                q__1.real = y[i__3].real + q__2.real, q__1.imag = y[i__3].imag + q__2.imag;
                y[i__2].real = q__1.real, y[i__2].imag = q__1.imag;
                jy += *incy;
                if (j > *ku)
                {
                    kx += *incx;
                }
                /* L140: */
            }
        }
    }
    return 0;
    /* End of CGBMV . */
}
/* cgbmv_ */

