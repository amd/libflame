/* stbmv.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int stbmv_(char *uplo, char *trans, char *diag, integer *n, integer *k, real *a, integer *lda, real *x, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    integer info;
    real temp;
    integer i__, j, l;
    extern logical lsame_(char *, char *, integer, integer);
    integer kplus1, ix, jx, kx;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical nounit;
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* STBMV performs one of the matrix-vector operations */
    /* x := A*x, or x := A'*x, */
    /* where x is an n element vector and A is an n by n unit, or non-unit, */
    /* upper or lower triangular band matrix, with ( k + 1 ) diagonals. */
    /* Parameters */
    /* ========== */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the matrix is an upper or */
    /* lower triangular matrix as follows: */
    /* UPLO = 'U' or 'u' A is an upper triangular matrix. */
    /* UPLO = 'L' or 'l' A is a lower triangular matrix. */
    /* Unchanged on exit. */
    /* TRANS - CHARACTER*1. */
    /* On entry, TRANS specifies the operation to be performed as */
    /* follows: */
    /* TRANS = 'N' or 'n' x := A*x. */
    /* TRANS = 'T' or 't' x := A'*x. */
    /* TRANS = 'C' or 'c' x := A'*x. */
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
    /* K - INTEGER. */
    /* On entry with UPLO = 'U' or 'u', K specifies the number of */
    /* super-diagonals of the matrix A. */
    /* On entry with UPLO = 'L' or 'l', K specifies the number of */
    /* sub-diagonals of the matrix A. */
    /* K must satisfy 0 .le. K. */
    /* Unchanged on exit. */
    /* A - REAL array of DIMENSION ( LDA, n ). */
    /* Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
    /* by n part of the array A must contain the upper triangular */
    /* band part of the matrix of coefficients, supplied column by */
    /* column, with the leading diagonal of the matrix in row */
    /* ( k + 1 ) of the array, the first super-diagonal starting at */
    /* position 2 in row k, and so on. The top left k by k triangle */
    /* of the array A is not referenced. */
    /* The following program segment will transfer an upper */
    /* triangular band matrix from conventional full matrix storage */
    /* to band storage: */
    /* DO 20, J = 1, N */
    /* M = K + 1 - J */
    /* DO 10, I = MAX( 1, J - K ), J */
    /* A( M + I, J ) = matrix( I, J ) */
    /* 10 CONTINUE */
    /* 20 CONTINUE */
    /* Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
    /* by n part of the array A must contain the lower triangular */
    /* band part of the matrix of coefficients, supplied column by */
    /* column, with the leading diagonal of the matrix in row 1 of */
    /* the array, the first sub-diagonal starting at position 1 in */
    /* row 2, and so on. The bottom right k by k triangle of the */
    /* array A is not referenced. */
    /* The following program segment will transfer a lower */
    /* triangular band matrix from conventional full matrix storage */
    /* to band storage: */
    /* DO 20, J = 1, N */
    /* M = 1 - J */
    /* DO 10, I = J, MIN( N, J + K ) */
    /* A( M + I, J ) = matrix( I, J ) */
    /* 10 CONTINUE */
    /* 20 CONTINUE */
    /* Note that when DIAG = 'U' or 'u' the elements of the array A */
    /* corresponding to the diagonal elements of the matrix are not */
    /* referenced, but are assumed to be unity. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. LDA must be at least */
    /* ( k + 1 ). */
    /* Unchanged on exit. */
    /* X - REAL array of dimension at least */
    /* ( 1 + ( n - 1 )*f2c_abs( INCX ) ). */
    /* Before entry, the incremented array X must contain the n */
    /* element vector x. On exit, X is overwritten with the */
    /* tranformed vector x. */
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
    else if (*k < 0)
    {
        info = 5;
    }
    else if (*lda < *k + 1)
    {
        info = 7;
    }
    else if (*incx == 0)
    {
        info = 9;
    }
    if (info != 0)
    {
        xerbla_("STBMV ", &info, (ftnlen)6);
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
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through A. */
    if (lsame_(trans, "N", 1, 1))
    {
        /* Form x := A*x. */
        if (lsame_(uplo, "U", 1, 1))
        {
            kplus1 = *k + 1;
            if (*incx == 1)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    if (x[j] != 0.f)
                    {
                        temp = x[j];
                        l = kplus1 - j;
                        /* Computing MAX */
                        i__2 = 1, i__3 = j - *k;
                        i__4 = j - 1;
                        for (i__ = fla_max(i__2,i__3);
                                i__ <= i__4;
                                ++i__)
                        {
                            x[i__] += temp * a[l + i__ + j * a_dim1];
                            /* L10: */
                        }
                        if (nounit)
                        {
                            x[j] *= a[kplus1 + j * a_dim1];
                        }
                    }
                    /* L20: */
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
                    if (x[jx] != 0.f)
                    {
                        temp = x[jx];
                        ix = kx;
                        l = kplus1 - j;
                        /* Computing MAX */
                        i__4 = 1, i__2 = j - *k;
                        i__3 = j - 1;
                        for (i__ = fla_max(i__4,i__2);
                                i__ <= i__3;
                                ++i__)
                        {
                            x[ix] += temp * a[l + i__ + j * a_dim1];
                            ix += *incx;
                            /* L30: */
                        }
                        if (nounit)
                        {
                            x[jx] *= a[kplus1 + j * a_dim1];
                        }
                    }
                    jx += *incx;
                    if (j > *k)
                    {
                        kx += *incx;
                    }
                    /* L40: */
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
                    if (x[j] != 0.f)
                    {
                        temp = x[j];
                        l = 1 - j;
                        /* Computing MIN */
                        i__1 = *n, i__3 = j + *k;
                        i__4 = j + 1;
                        for (i__ = fla_min(i__1,i__3);
                                i__ >= i__4;
                                --i__)
                        {
                            x[i__] += temp * a[l + i__ + j * a_dim1];
                            /* L50: */
                        }
                        if (nounit)
                        {
                            x[j] *= a[j * a_dim1 + 1];
                        }
                    }
                    /* L60: */
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
                    if (x[jx] != 0.f)
                    {
                        temp = x[jx];
                        ix = kx;
                        l = 1 - j;
                        /* Computing MIN */
                        i__4 = *n, i__1 = j + *k;
                        i__3 = j + 1;
                        for (i__ = fla_min(i__4,i__1);
                                i__ >= i__3;
                                --i__)
                        {
                            x[ix] += temp * a[l + i__ + j * a_dim1];
                            ix -= *incx;
                            /* L70: */
                        }
                        if (nounit)
                        {
                            x[jx] *= a[j * a_dim1 + 1];
                        }
                    }
                    jx -= *incx;
                    if (*n - j >= *k)
                    {
                        kx -= *incx;
                    }
                    /* L80: */
                }
            }
        }
    }
    else
    {
        /* Form x := A'*x. */
        if (lsame_(uplo, "U", 1, 1))
        {
            kplus1 = *k + 1;
            if (*incx == 1)
            {
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    temp = x[j];
                    l = kplus1 - j;
                    if (nounit)
                    {
                        temp *= a[kplus1 + j * a_dim1];
                    }
                    /* Computing MAX */
                    i__4 = 1, i__1 = j - *k;
                    i__3 = fla_max(i__4,i__1);
                    for (i__ = j - 1;
                            i__ >= i__3;
                            --i__)
                    {
                        temp += a[l + i__ + j * a_dim1] * x[i__];
                        /* L90: */
                    }
                    x[j] = temp;
                    /* L100: */
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
                    kx -= *incx;
                    ix = kx;
                    l = kplus1 - j;
                    if (nounit)
                    {
                        temp *= a[kplus1 + j * a_dim1];
                    }
                    /* Computing MAX */
                    i__4 = 1, i__1 = j - *k;
                    i__3 = fla_max(i__4,i__1);
                    for (i__ = j - 1;
                            i__ >= i__3;
                            --i__)
                    {
                        temp += a[l + i__ + j * a_dim1] * x[ix];
                        ix -= *incx;
                        /* L110: */
                    }
                    x[jx] = temp;
                    jx -= *incx;
                    /* L120: */
                }
            }
        }
        else
        {
            if (*incx == 1)
            {
                i__3 = *n;
                for (j = 1;
                        j <= i__3;
                        ++j)
                {
                    temp = x[j];
                    l = 1 - j;
                    if (nounit)
                    {
                        temp *= a[j * a_dim1 + 1];
                    }
                    /* Computing MIN */
                    i__1 = *n, i__2 = j + *k;
                    i__4 = fla_min(i__1,i__2);
                    for (i__ = j + 1;
                            i__ <= i__4;
                            ++i__)
                    {
                        temp += a[l + i__ + j * a_dim1] * x[i__];
                        /* L130: */
                    }
                    x[j] = temp;
                    /* L140: */
                }
            }
            else
            {
                jx = kx;
                i__3 = *n;
                for (j = 1;
                        j <= i__3;
                        ++j)
                {
                    temp = x[jx];
                    kx += *incx;
                    ix = kx;
                    l = 1 - j;
                    if (nounit)
                    {
                        temp *= a[j * a_dim1 + 1];
                    }
                    /* Computing MIN */
                    i__1 = *n, i__2 = j + *k;
                    i__4 = fla_min(i__1,i__2);
                    for (i__ = j + 1;
                            i__ <= i__4;
                            ++i__)
                    {
                        temp += a[l + i__ + j * a_dim1] * x[ix];
                        ix += *incx;
                        /* L150: */
                    }
                    x[jx] = temp;
                    jx += *incx;
                    /* L160: */
                }
            }
        }
    }
    return 0;
    /* End of STBMV . */
}
/* stbmv_ */

