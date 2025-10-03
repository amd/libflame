/* ctrsv.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int ctrsv_(char *uplo, char *trans, char *diag, integer *n, scomplex *a, integer *lda, scomplex *x, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void c_div(scomplex *, scomplex *, scomplex *), r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    integer info;
    scomplex temp;
    integer i__, j;
    extern logical lsame_(char *, char *, integer, integer);
    integer ix, jx, kx;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical noconj, nounit;
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* CTRSV solves one of the systems of equations */
    /* A*x = b, or A'*x = b, or conjg( A' )*x = b, */
    /* where b and x are n element vectors and A is an n by n unit, or */
    /* non-unit, upper or lower triangular matrix. */
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
    /* TRANS = 'C' or 'c' conjg( A' )*x = b. */
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
    /* A - COMPLEX array of DIMENSION ( LDA, n ). */
    /* Before entry with UPLO = 'U' or 'u', the leading n by n */
    /* upper triangular part of the array A must contain the upper */
    /* triangular matrix and the strictly lower triangular part of */
    /* A is not referenced. */
    /* Before entry with UPLO = 'L' or 'l', the leading n by n */
    /* lower triangular part of the array A must contain the lower */
    /* triangular matrix and the strictly upper triangular part of */
    /* A is not referenced. */
    /* Note that when DIAG = 'U' or 'u', the diagonal elements of */
    /* A are not referenced either, but are assumed to be unity. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. LDA must be at least */
    /* fla_max( 1, n ). */
    /* Unchanged on exit. */
    /* X - COMPLEX array of dimension at least */
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
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
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
    else if (*lda < fla_max(1,*n))
    {
        info = 6;
    }
    else if (*incx == 0)
    {
        info = 8;
    }
    if (info != 0)
    {
        xerbla_("CTRSV ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        return 0;
    }
    noconj = lsame_(trans, "T", 1, 1);
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
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through A. */
    if (lsame_(trans, "N", 1, 1))
    {
        /* Form x := inv( A )*x. */
        if (lsame_(uplo, "U", 1, 1))
        {
            if (*incx == 1)
            {
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    i__1 = j;
                    if (x[i__1].real != 0.f || x[i__1].imag != 0.f)
                    {
                        if (nounit)
                        {
                            i__1 = j;
                            c_div(&q__1, &x[j], &a[j + j * a_dim1]);
                            x[i__1].real = q__1.real, x[i__1].imag = q__1.imag;
                        }
                        i__1 = j;
                        temp.real = x[i__1].real, temp.imag = x[i__1].imag;
                        for (i__ = j - 1;
                                i__ >= 1;
                                --i__)
                        {
                            i__1 = i__;
                            i__2 = i__;
                            i__3 = i__ + j * a_dim1;
                            q__2.real = temp.real * a[i__3].real - temp.imag * a[i__3].imag, q__2.imag = temp.real * a[i__3].imag + temp.imag * a[ i__3].real;
                            q__1.real = x[i__2].real - q__2.real, q__1.imag = x[i__2].imag - q__2.imag;
                            x[i__1].real = q__1.real, x[i__1].imag = q__1.imag;
                            /* L10: */
                        }
                    }
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
                    i__1 = jx;
                    if (x[i__1].real != 0.f || x[i__1].imag != 0.f)
                    {
                        if (nounit)
                        {
                            i__1 = jx;
                            c_div(&q__1, &x[jx], &a[j + j * a_dim1]);
                            x[i__1].real = q__1.real, x[i__1].imag = q__1.imag;
                        }
                        i__1 = jx;
                        temp.real = x[i__1].real, temp.imag = x[i__1].imag;
                        ix = jx;
                        for (i__ = j - 1;
                                i__ >= 1;
                                --i__)
                        {
                            ix -= *incx;
                            i__1 = ix;
                            i__2 = ix;
                            i__3 = i__ + j * a_dim1;
                            q__2.real = temp.real * a[i__3].real - temp.imag * a[i__3].imag, q__2.imag = temp.real * a[i__3].imag + temp.imag * a[ i__3].real;
                            q__1.real = x[i__2].real - q__2.real, q__1.imag = x[i__2].imag - q__2.imag;
                            x[i__1].real = q__1.real, x[i__1].imag = q__1.imag;
                            /* L30: */
                        }
                    }
                    jx -= *incx;
                    /* L40: */
                }
            }
        }
        else
        {
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j;
                    if (x[i__2].real != 0.f || x[i__2].imag != 0.f)
                    {
                        if (nounit)
                        {
                            i__2 = j;
                            c_div(&q__1, &x[j], &a[j + j * a_dim1]);
                            x[i__2].real = q__1.real, x[i__2].imag = q__1.imag;
                        }
                        i__2 = j;
                        temp.real = x[i__2].real, temp.imag = x[i__2].imag;
                        i__2 = *n;
                        for (i__ = j + 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__;
                            i__4 = i__;
                            i__5 = i__ + j * a_dim1;
                            q__2.real = temp.real * a[i__5].real - temp.imag * a[i__5].imag, q__2.imag = temp.real * a[i__5].imag + temp.imag * a[ i__5].real;
                            q__1.real = x[i__4].real - q__2.real, q__1.imag = x[i__4].imag - q__2.imag;
                            x[i__3].real = q__1.real, x[i__3].imag = q__1.imag;
                            /* L50: */
                        }
                    }
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
                    i__2 = jx;
                    if (x[i__2].real != 0.f || x[i__2].imag != 0.f)
                    {
                        if (nounit)
                        {
                            i__2 = jx;
                            c_div(&q__1, &x[jx], &a[j + j * a_dim1]);
                            x[i__2].real = q__1.real, x[i__2].imag = q__1.imag;
                        }
                        i__2 = jx;
                        temp.real = x[i__2].real, temp.imag = x[i__2].imag;
                        ix = jx;
                        i__2 = *n;
                        for (i__ = j + 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            ix += *incx;
                            i__3 = ix;
                            i__4 = ix;
                            i__5 = i__ + j * a_dim1;
                            q__2.real = temp.real * a[i__5].real - temp.imag * a[i__5].imag, q__2.imag = temp.real * a[i__5].imag + temp.imag * a[ i__5].real;
                            q__1.real = x[i__4].real - q__2.real, q__1.imag = x[i__4].imag - q__2.imag;
                            x[i__3].real = q__1.real, x[i__3].imag = q__1.imag;
                            /* L70: */
                        }
                    }
                    jx += *incx;
                    /* L80: */
                }
            }
        }
    }
    else
    {
        /* Form x := inv( A' )*x or x := inv( conjg( A' ) )*x. */
        if (lsame_(uplo, "U", 1, 1))
        {
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j;
                    temp.real = x[i__2].real, temp.imag = x[i__2].imag;
                    if (noconj)
                    {
                        i__2 = j - 1;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * a_dim1;
                            i__4 = i__;
                            q__2.real = a[i__3].real * x[i__4].real - a[i__3].imag * x[ i__4].imag, q__2.imag = a[i__3].real * x[i__4].imag + a[i__3].imag * x[i__4].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            /* L90: */
                        }
                        if (nounit)
                        {
                            c_div(&q__1, &temp, &a[j + j * a_dim1]);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    else
                    {
                        i__2 = j - 1;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            r_cnjg(&q__3, &a[i__ + j * a_dim1]);
                            i__3 = i__;
                            q__2.real = q__3.real * x[i__3].real - q__3.imag * x[i__3].imag, q__2.imag = q__3.real * x[i__3].imag + q__3.imag * x[ i__3].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            /* L100: */
                        }
                        if (nounit)
                        {
                            r_cnjg(&q__2, &a[j + j * a_dim1]);
                            c_div(&q__1, &temp, &q__2);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    i__2 = j;
                    x[i__2].real = temp.real, x[i__2].imag = temp.imag;
                    /* L110: */
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
                    ix = kx;
                    i__2 = jx;
                    temp.real = x[i__2].real, temp.imag = x[i__2].imag;
                    if (noconj)
                    {
                        i__2 = j - 1;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * a_dim1;
                            i__4 = ix;
                            q__2.real = a[i__3].real * x[i__4].real - a[i__3].imag * x[ i__4].imag, q__2.imag = a[i__3].real * x[i__4].imag + a[i__3].imag * x[i__4].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            ix += *incx;
                            /* L120: */
                        }
                        if (nounit)
                        {
                            c_div(&q__1, &temp, &a[j + j * a_dim1]);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    else
                    {
                        i__2 = j - 1;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            r_cnjg(&q__3, &a[i__ + j * a_dim1]);
                            i__3 = ix;
                            q__2.real = q__3.real * x[i__3].real - q__3.imag * x[i__3].imag, q__2.imag = q__3.real * x[i__3].imag + q__3.imag * x[ i__3].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            ix += *incx;
                            /* L130: */
                        }
                        if (nounit)
                        {
                            r_cnjg(&q__2, &a[j + j * a_dim1]);
                            c_div(&q__1, &temp, &q__2);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    i__2 = jx;
                    x[i__2].real = temp.real, x[i__2].imag = temp.imag;
                    jx += *incx;
                    /* L140: */
                }
            }
        }
        else
        {
            if (*incx == 1)
            {
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    i__1 = j;
                    temp.real = x[i__1].real, temp.imag = x[i__1].imag;
                    if (noconj)
                    {
                        i__1 = j + 1;
                        for (i__ = *n;
                                i__ >= i__1;
                                --i__)
                        {
                            i__2 = i__ + j * a_dim1;
                            i__3 = i__;
                            q__2.real = a[i__2].real * x[i__3].real - a[i__2].imag * x[ i__3].imag, q__2.imag = a[i__2].real * x[i__3].imag + a[i__2].imag * x[i__3].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            /* L150: */
                        }
                        if (nounit)
                        {
                            c_div(&q__1, &temp, &a[j + j * a_dim1]);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    else
                    {
                        i__1 = j + 1;
                        for (i__ = *n;
                                i__ >= i__1;
                                --i__)
                        {
                            r_cnjg(&q__3, &a[i__ + j * a_dim1]);
                            i__2 = i__;
                            q__2.real = q__3.real * x[i__2].real - q__3.imag * x[i__2].imag, q__2.imag = q__3.real * x[i__2].imag + q__3.imag * x[ i__2].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            /* L160: */
                        }
                        if (nounit)
                        {
                            r_cnjg(&q__2, &a[j + j * a_dim1]);
                            c_div(&q__1, &temp, &q__2);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    i__1 = j;
                    x[i__1].real = temp.real, x[i__1].imag = temp.imag;
                    /* L170: */
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
                    ix = kx;
                    i__1 = jx;
                    temp.real = x[i__1].real, temp.imag = x[i__1].imag;
                    if (noconj)
                    {
                        i__1 = j + 1;
                        for (i__ = *n;
                                i__ >= i__1;
                                --i__)
                        {
                            i__2 = i__ + j * a_dim1;
                            i__3 = ix;
                            q__2.real = a[i__2].real * x[i__3].real - a[i__2].imag * x[ i__3].imag, q__2.imag = a[i__2].real * x[i__3].imag + a[i__2].imag * x[i__3].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            ix -= *incx;
                            /* L180: */
                        }
                        if (nounit)
                        {
                            c_div(&q__1, &temp, &a[j + j * a_dim1]);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    else
                    {
                        i__1 = j + 1;
                        for (i__ = *n;
                                i__ >= i__1;
                                --i__)
                        {
                            r_cnjg(&q__3, &a[i__ + j * a_dim1]);
                            i__2 = ix;
                            q__2.real = q__3.real * x[i__2].real - q__3.imag * x[i__2].imag, q__2.imag = q__3.real * x[i__2].imag + q__3.imag * x[ i__2].real;
                            q__1.real = temp.real - q__2.real, q__1.imag = temp.imag - q__2.imag;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            ix -= *incx;
                            /* L190: */
                        }
                        if (nounit)
                        {
                            r_cnjg(&q__2, &a[j + j * a_dim1]);
                            c_div(&q__1, &temp, &q__2);
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    i__1 = jx;
                    x[i__1].real = temp.real, x[i__1].imag = temp.imag;
                    jx -= *incx;
                    /* L200: */
                }
            }
        }
    }
    return 0;
    /* End of CTRSV . */
}
/* ctrsv_ */

