/* ctrmm.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int ctrmm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, scomplex *alpha, scomplex *a, integer *lda, scomplex *b, integer *ldb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    scomplex q__1, q__2, q__3;
    /* Builtin functions */
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    integer info;
    scomplex temp;
    integer i__, j, k;
    extern logical lsame_(char *, char *, integer, integer);
    logical lside;
    integer nrowa;
    logical upper;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical noconj, nounit;
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* CTRMM performs one of the matrix-matrix operations */
    /* B := alpha*op( A )*B, or B := alpha*B*op( A ) */
    /* where alpha is a scalar, B is an m by n matrix, A is a unit, or */
    /* non-unit, upper or lower triangular matrix and op( A ) is one of */
    /* op( A ) = A or op( A ) = A' or op( A ) = conjg( A' ). */
    /* Parameters */
    /* ========== */
    /* SIDE - CHARACTER*1. */
    /* On entry, SIDE specifies whether op( A ) multiplies B from */
    /* the left or right as follows: */
    /* SIDE = 'L' or 'l' B := alpha*op( A )*B. */
    /* SIDE = 'R' or 'r' B := alpha*B*op( A ). */
    /* Unchanged on exit. */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the matrix A is an upper or */
    /* lower triangular matrix as follows: */
    /* UPLO = 'U' or 'u' A is an upper triangular matrix. */
    /* UPLO = 'L' or 'l' A is a lower triangular matrix. */
    /* Unchanged on exit. */
    /* TRANSA - CHARACTER*1. */
    /* On entry, TRANSA specifies the form of op( A ) to be used in */
    /* the matrix multiplication as follows: */
    /* TRANSA = 'N' or 'n' op( A ) = A. */
    /* TRANSA = 'T' or 't' op( A ) = A'. */
    /* TRANSA = 'C' or 'c' op( A ) = conjg( A' ). */
    /* Unchanged on exit. */
    /* DIAG - CHARACTER*1. */
    /* On entry, DIAG specifies whether or not A is unit triangular */
    /* as follows: */
    /* DIAG = 'U' or 'u' A is assumed to be unit triangular. */
    /* DIAG = 'N' or 'n' A is not assumed to be unit */
    /* triangular. */
    /* Unchanged on exit. */
    /* M - INTEGER. */
    /* On entry, M specifies the number of rows of B. M must be at */
    /* least zero. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the number of columns of B. N must be */
    /* at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX . */
    /* On entry, ALPHA specifies the scalar alpha. When alpha is */
    /* zero then A is not referenced and B need not be set before */
    /* entry. */
    /* Unchanged on exit. */
    /* A - COMPLEX array of DIMENSION ( LDA, k ), where k is m */
    /* when SIDE = 'L' or 'l' and is n when SIDE = 'R' or 'r'. */
    /* Before entry with UPLO = 'U' or 'u', the leading k by k */
    /* upper triangular part of the array A must contain the upper */
    /* triangular matrix and the strictly lower triangular part of */
    /* A is not referenced. */
    /* Before entry with UPLO = 'L' or 'l', the leading k by k */
    /* lower triangular part of the array A must contain the lower */
    /* triangular matrix and the strictly upper triangular part of */
    /* A is not referenced. */
    /* Note that when DIAG = 'U' or 'u', the diagonal elements of */
    /* A are not referenced either, but are assumed to be unity. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. When SIDE = 'L' or 'l' then */
    /* LDA must be at least fla_max( 1, m ), when SIDE = 'R' or 'r' */
    /* then LDA must be at least fla_max( 1, n ). */
    /* Unchanged on exit. */
    /* B - COMPLEX array of DIMENSION ( LDB, n ). */
    /* Before entry, the leading m by n part of the array B must */
    /* contain the matrix B, and on exit is overwritten by the */
    /* transformed matrix. */
    /* LDB - INTEGER. */
    /* On entry, LDB specifies the first dimension of B as declared */
    /* in the calling (sub) program. LDB must be at least */
    /* fla_max( 1, m ). */
    /* Unchanged on exit. */
    /* Level 3 Blas routine. */
    /* -- Written on 8-February-1989. */
    /* Jack Dongarra, Argonne National Laboratory. */
    /* Iain Duff, AERE Harwell. */
    /* Jeremy Du Croz, Numerical Algorithms Group Ltd. */
    /* Sven Hammarling, Numerical Algorithms Group Ltd. */
    /* .. External Functions .. */
    /* .. External Subroutines .. */
    /* .. Intrinsic Functions .. */
    /* .. Local Scalars .. */
    /* .. Parameters .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    /* Function Body */
    lside = lsame_(side, "L", 1, 1);
    if (lside)
    {
        nrowa = *m;
    }
    else
    {
        nrowa = *n;
    }
    noconj = lsame_(transa, "T", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);
    upper = lsame_(uplo, "U", 1, 1);
    info = 0;
    if (! lside && ! lsame_(side, "R", 1, 1))
    {
        info = 1;
    }
    else if (! upper && ! lsame_(uplo, "L", 1, 1))
    {
        info = 2;
    }
    else if (! lsame_(transa, "N", 1, 1) && ! lsame_(transa, "T", 1, 1) && ! lsame_(transa, "C", 1, 1))
    {
        info = 3;
    }
    else if (! lsame_(diag, "U", 1, 1) && ! lsame_(diag, "N", 1, 1))
    {
        info = 4;
    }
    else if (*m < 0)
    {
        info = 5;
    }
    else if (*n < 0)
    {
        info = 6;
    }
    else if (*lda < fla_max(1,nrowa))
    {
        info = 9;
    }
    else if (*ldb < fla_max(1,*m))
    {
        info = 11;
    }
    if (info != 0)
    {
        xerbla_("CTRMM ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        return 0;
    }
    /* And when alpha.eq.zero. */
    if (alpha->real == 0.f && alpha->imag == 0.f)
    {
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = *m;
            for (i__ = 1;
                    i__ <= i__2;
                    ++i__)
            {
                i__3 = i__ + j * b_dim1;
                b[i__3].real = 0.f, b[i__3].imag = 0.f;
                /* L10: */
            }
            /* L20: */
        }
        return 0;
    }
    /* Start the operations. */
    if (lside)
    {
        if (lsame_(transa, "N", 1, 1))
        {
            /* Form B := alpha*A*B. */
            if (upper)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *m;
                    for (k = 1;
                            k <= i__2;
                            ++k)
                    {
                        i__3 = k + j * b_dim1;
                        if (b[i__3].real != 0.f || b[i__3].imag != 0.f)
                        {
                            i__3 = k + j * b_dim1;
                            q__1.real = alpha->real * b[i__3].real - alpha->imag * b[i__3] .imag, q__1.imag = alpha->real * b[i__3].imag + alpha->imag * b[i__3].real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            i__3 = k - 1;
                            for (i__ = 1;
                                    i__ <= i__3;
                                    ++i__)
                            {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = i__ + k * a_dim1;
                                q__2.real = temp.real * a[i__6].real - temp.imag * a[i__6] .imag, q__2.imag = temp.real * a[i__6].imag + temp.imag * a[i__6].real;
                                q__1.real = b[i__5].real + q__2.real, q__1.imag = b[i__5] .imag + q__2.imag;
                                b[i__4].real = q__1.real, b[i__4].imag = q__1.imag;
                                /* L30: */
                            }
                            if (nounit)
                            {
                                i__3 = k + k * a_dim1;
                                q__1.real = temp.real * a[i__3].real - temp.imag * a[i__3] .imag, q__1.imag = temp.real * a[i__3].imag + temp.imag * a[i__3].real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            i__3 = k + j * b_dim1;
                            b[i__3].real = temp.real, b[i__3].imag = temp.imag;
                        }
                        /* L40: */
                    }
                    /* L50: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    for (k = *m;
                            k >= 1;
                            --k)
                    {
                        i__2 = k + j * b_dim1;
                        if (b[i__2].real != 0.f || b[i__2].imag != 0.f)
                        {
                            i__2 = k + j * b_dim1;
                            q__1.real = alpha->real * b[i__2].real - alpha->imag * b[i__2] .imag, q__1.imag = alpha->real * b[i__2].imag + alpha->imag * b[i__2].real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            i__2 = k + j * b_dim1;
                            b[i__2].real = temp.real, b[i__2].imag = temp.imag;
                            if (nounit)
                            {
                                i__2 = k + j * b_dim1;
                                i__3 = k + j * b_dim1;
                                i__4 = k + k * a_dim1;
                                q__1.real = b[i__3].real * a[i__4].real - b[i__3].imag * a[i__4].imag, q__1.imag = b[i__3].real * a[ i__4].imag + b[i__3].imag * a[i__4].real;
                                b[i__2].real = q__1.real, b[i__2].imag = q__1.imag;
                            }
                            i__2 = *m;
                            for (i__ = k + 1;
                                    i__ <= i__2;
                                    ++i__)
                            {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + k * a_dim1;
                                q__2.real = temp.real * a[i__5].real - temp.imag * a[i__5] .imag, q__2.imag = temp.real * a[i__5].imag + temp.imag * a[i__5].real;
                                q__1.real = b[i__4].real + q__2.real, q__1.imag = b[i__4] .imag + q__2.imag;
                                b[i__3].real = q__1.real, b[i__3].imag = q__1.imag;
                                /* L60: */
                            }
                        }
                        /* L70: */
                    }
                    /* L80: */
                }
            }
        }
        else
        {
            /* Form B := alpha*A'*B or B := alpha*conjg( A' )*B. */
            if (upper)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    for (i__ = *m;
                            i__ >= 1;
                            --i__)
                    {
                        i__2 = i__ + j * b_dim1;
                        temp.real = b[i__2].real, temp.imag = b[i__2].imag;
                        if (noconj)
                        {
                            if (nounit)
                            {
                                i__2 = i__ + i__ * a_dim1;
                                q__1.real = temp.real * a[i__2].real - temp.imag * a[i__2] .imag, q__1.imag = temp.real * a[i__2].imag + temp.imag * a[i__2].real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            i__2 = i__ - 1;
                            for (k = 1;
                                    k <= i__2;
                                    ++k)
                            {
                                i__3 = k + i__ * a_dim1;
                                i__4 = k + j * b_dim1;
                                q__2.real = a[i__3].real * b[i__4].real - a[i__3].imag * b[i__4].imag, q__2.imag = a[i__3].real * b[ i__4].imag + a[i__3].imag * b[i__4].real;
                                q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                                /* L90: */
                            }
                        }
                        else
                        {
                            if (nounit)
                            {
                                r_cnjg(&q__2, &a[i__ + i__ * a_dim1]);
                                q__1.real = temp.real * q__2.real - temp.imag * q__2.imag, q__1.imag = temp.real * q__2.imag + temp.imag * q__2.real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            i__2 = i__ - 1;
                            for (k = 1;
                                    k <= i__2;
                                    ++k)
                            {
                                r_cnjg(&q__3, &a[k + i__ * a_dim1]);
                                i__3 = k + j * b_dim1;
                                q__2.real = q__3.real * b[i__3].real - q__3.imag * b[i__3] .imag, q__2.imag = q__3.real * b[i__3].imag + q__3.imag * b[i__3].real;
                                q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                                /* L100: */
                            }
                        }
                        i__2 = i__ + j * b_dim1;
                        q__1.real = alpha->real * temp.real - alpha->imag * temp.imag, q__1.imag = alpha->real * temp.imag + alpha->imag * temp.real;
                        b[i__2].real = q__1.real, b[i__2].imag = q__1.imag;
                        /* L110: */
                    }
                    /* L120: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *m;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        temp.real = b[i__3].real, temp.imag = b[i__3].imag;
                        if (noconj)
                        {
                            if (nounit)
                            {
                                i__3 = i__ + i__ * a_dim1;
                                q__1.real = temp.real * a[i__3].real - temp.imag * a[i__3] .imag, q__1.imag = temp.real * a[i__3].imag + temp.imag * a[i__3].real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            i__3 = *m;
                            for (k = i__ + 1;
                                    k <= i__3;
                                    ++k)
                            {
                                i__4 = k + i__ * a_dim1;
                                i__5 = k + j * b_dim1;
                                q__2.real = a[i__4].real * b[i__5].real - a[i__4].imag * b[i__5].imag, q__2.imag = a[i__4].real * b[ i__5].imag + a[i__4].imag * b[i__5].real;
                                q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                                /* L130: */
                            }
                        }
                        else
                        {
                            if (nounit)
                            {
                                r_cnjg(&q__2, &a[i__ + i__ * a_dim1]);
                                q__1.real = temp.real * q__2.real - temp.imag * q__2.imag, q__1.imag = temp.real * q__2.imag + temp.imag * q__2.real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            i__3 = *m;
                            for (k = i__ + 1;
                                    k <= i__3;
                                    ++k)
                            {
                                r_cnjg(&q__3, &a[k + i__ * a_dim1]);
                                i__4 = k + j * b_dim1;
                                q__2.real = q__3.real * b[i__4].real - q__3.imag * b[i__4] .imag, q__2.imag = q__3.real * b[i__4].imag + q__3.imag * b[i__4].real;
                                q__1.real = temp.real + q__2.real, q__1.imag = temp.imag + q__2.imag;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                                /* L140: */
                            }
                        }
                        i__3 = i__ + j * b_dim1;
                        q__1.real = alpha->real * temp.real - alpha->imag * temp.imag, q__1.imag = alpha->real * temp.imag + alpha->imag * temp.real;
                        b[i__3].real = q__1.real, b[i__3].imag = q__1.imag;
                        /* L150: */
                    }
                    /* L160: */
                }
            }
        }
    }
    else
    {
        if (lsame_(transa, "N", 1, 1))
        {
            /* Form B := alpha*B*A. */
            if (upper)
            {
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    temp.real = alpha->real, temp.imag = alpha->imag;
                    if (nounit)
                    {
                        i__1 = j + j * a_dim1;
                        q__1.real = temp.real * a[i__1].real - temp.imag * a[i__1].imag, q__1.imag = temp.real * a[i__1].imag + temp.imag * a[i__1] .real;
                        temp.real = q__1.real, temp.imag = q__1.imag;
                    }
                    i__1 = *m;
                    for (i__ = 1;
                            i__ <= i__1;
                            ++i__)
                    {
                        i__2 = i__ + j * b_dim1;
                        i__3 = i__ + j * b_dim1;
                        q__1.real = temp.real * b[i__3].real - temp.imag * b[i__3].imag, q__1.imag = temp.real * b[i__3].imag + temp.imag * b[i__3] .real;
                        b[i__2].real = q__1.real, b[i__2].imag = q__1.imag;
                        /* L170: */
                    }
                    i__1 = j - 1;
                    for (k = 1;
                            k <= i__1;
                            ++k)
                    {
                        i__2 = k + j * a_dim1;
                        if (a[i__2].real != 0.f || a[i__2].imag != 0.f)
                        {
                            i__2 = k + j * a_dim1;
                            q__1.real = alpha->real * a[i__2].real - alpha->imag * a[i__2] .imag, q__1.imag = alpha->real * a[i__2].imag + alpha->imag * a[i__2].real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            i__2 = *m;
                            for (i__ = 1;
                                    i__ <= i__2;
                                    ++i__)
                            {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + k * b_dim1;
                                q__2.real = temp.real * b[i__5].real - temp.imag * b[i__5] .imag, q__2.imag = temp.real * b[i__5].imag + temp.imag * b[i__5].real;
                                q__1.real = b[i__4].real + q__2.real, q__1.imag = b[i__4] .imag + q__2.imag;
                                b[i__3].real = q__1.real, b[i__3].imag = q__1.imag;
                                /* L180: */
                            }
                        }
                        /* L190: */
                    }
                    /* L200: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    temp.real = alpha->real, temp.imag = alpha->imag;
                    if (nounit)
                    {
                        i__2 = j + j * a_dim1;
                        q__1.real = temp.real * a[i__2].real - temp.imag * a[i__2].imag, q__1.imag = temp.real * a[i__2].imag + temp.imag * a[i__2] .real;
                        temp.real = q__1.real, temp.imag = q__1.imag;
                    }
                    i__2 = *m;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * b_dim1;
                        i__4 = i__ + j * b_dim1;
                        q__1.real = temp.real * b[i__4].real - temp.imag * b[i__4].imag, q__1.imag = temp.real * b[i__4].imag + temp.imag * b[i__4] .real;
                        b[i__3].real = q__1.real, b[i__3].imag = q__1.imag;
                        /* L210: */
                    }
                    i__2 = *n;
                    for (k = j + 1;
                            k <= i__2;
                            ++k)
                    {
                        i__3 = k + j * a_dim1;
                        if (a[i__3].real != 0.f || a[i__3].imag != 0.f)
                        {
                            i__3 = k + j * a_dim1;
                            q__1.real = alpha->real * a[i__3].real - alpha->imag * a[i__3] .imag, q__1.imag = alpha->real * a[i__3].imag + alpha->imag * a[i__3].real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                            i__3 = *m;
                            for (i__ = 1;
                                    i__ <= i__3;
                                    ++i__)
                            {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = i__ + k * b_dim1;
                                q__2.real = temp.real * b[i__6].real - temp.imag * b[i__6] .imag, q__2.imag = temp.real * b[i__6].imag + temp.imag * b[i__6].real;
                                q__1.real = b[i__5].real + q__2.real, q__1.imag = b[i__5] .imag + q__2.imag;
                                b[i__4].real = q__1.real, b[i__4].imag = q__1.imag;
                                /* L220: */
                            }
                        }
                        /* L230: */
                    }
                    /* L240: */
                }
            }
        }
        else
        {
            /* Form B := alpha*B*A' or B := alpha*B*conjg( A' ). */
            if (upper)
            {
                i__1 = *n;
                for (k = 1;
                        k <= i__1;
                        ++k)
                {
                    i__2 = k - 1;
                    for (j = 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = j + k * a_dim1;
                        if (a[i__3].real != 0.f || a[i__3].imag != 0.f)
                        {
                            if (noconj)
                            {
                                i__3 = j + k * a_dim1;
                                q__1.real = alpha->real * a[i__3].real - alpha->imag * a[ i__3].imag, q__1.imag = alpha->real * a[i__3] .imag + alpha->imag * a[i__3].real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            else
                            {
                                r_cnjg(&q__2, &a[j + k * a_dim1]);
                                q__1.real = alpha->real * q__2.real - alpha->imag * q__2.imag, q__1.imag = alpha->real * q__2.imag + alpha->imag * q__2.real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            i__3 = *m;
                            for (i__ = 1;
                                    i__ <= i__3;
                                    ++i__)
                            {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = i__ + k * b_dim1;
                                q__2.real = temp.real * b[i__6].real - temp.imag * b[i__6] .imag, q__2.imag = temp.real * b[i__6].imag + temp.imag * b[i__6].real;
                                q__1.real = b[i__5].real + q__2.real, q__1.imag = b[i__5] .imag + q__2.imag;
                                b[i__4].real = q__1.real, b[i__4].imag = q__1.imag;
                                /* L250: */
                            }
                        }
                        /* L260: */
                    }
                    temp.real = alpha->real, temp.imag = alpha->imag;
                    if (nounit)
                    {
                        if (noconj)
                        {
                            i__2 = k + k * a_dim1;
                            q__1.real = temp.real * a[i__2].real - temp.imag * a[i__2].imag, q__1.imag = temp.real * a[i__2].imag + temp.imag * a[ i__2].real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                        else
                        {
                            r_cnjg(&q__2, &a[k + k * a_dim1]);
                            q__1.real = temp.real * q__2.real - temp.imag * q__2.imag, q__1.imag = temp.real * q__2.imag + temp.imag * q__2.real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    if (temp.real != 1.f || temp.imag != 0.f)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + k * b_dim1;
                            i__4 = i__ + k * b_dim1;
                            q__1.real = temp.real * b[i__4].real - temp.imag * b[i__4].imag, q__1.imag = temp.real * b[i__4].imag + temp.imag * b[ i__4].real;
                            b[i__3].real = q__1.real, b[i__3].imag = q__1.imag;
                            /* L270: */
                        }
                    }
                    /* L280: */
                }
            }
            else
            {
                for (k = *n;
                        k >= 1;
                        --k)
                {
                    i__1 = *n;
                    for (j = k + 1;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = j + k * a_dim1;
                        if (a[i__2].real != 0.f || a[i__2].imag != 0.f)
                        {
                            if (noconj)
                            {
                                i__2 = j + k * a_dim1;
                                q__1.real = alpha->real * a[i__2].real - alpha->imag * a[ i__2].imag, q__1.imag = alpha->real * a[i__2] .imag + alpha->imag * a[i__2].real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            else
                            {
                                r_cnjg(&q__2, &a[j + k * a_dim1]);
                                q__1.real = alpha->real * q__2.real - alpha->imag * q__2.imag, q__1.imag = alpha->real * q__2.imag + alpha->imag * q__2.real;
                                temp.real = q__1.real, temp.imag = q__1.imag;
                            }
                            i__2 = *m;
                            for (i__ = 1;
                                    i__ <= i__2;
                                    ++i__)
                            {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + k * b_dim1;
                                q__2.real = temp.real * b[i__5].real - temp.imag * b[i__5] .imag, q__2.imag = temp.real * b[i__5].imag + temp.imag * b[i__5].real;
                                q__1.real = b[i__4].real + q__2.real, q__1.imag = b[i__4] .imag + q__2.imag;
                                b[i__3].real = q__1.real, b[i__3].imag = q__1.imag;
                                /* L290: */
                            }
                        }
                        /* L300: */
                    }
                    temp.real = alpha->real, temp.imag = alpha->imag;
                    if (nounit)
                    {
                        if (noconj)
                        {
                            i__1 = k + k * a_dim1;
                            q__1.real = temp.real * a[i__1].real - temp.imag * a[i__1].imag, q__1.imag = temp.real * a[i__1].imag + temp.imag * a[ i__1].real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                        else
                        {
                            r_cnjg(&q__2, &a[k + k * a_dim1]);
                            q__1.real = temp.real * q__2.real - temp.imag * q__2.imag, q__1.imag = temp.real * q__2.imag + temp.imag * q__2.real;
                            temp.real = q__1.real, temp.imag = q__1.imag;
                        }
                    }
                    if (temp.real != 1.f || temp.imag != 0.f)
                    {
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + k * b_dim1;
                            i__3 = i__ + k * b_dim1;
                            q__1.real = temp.real * b[i__3].real - temp.imag * b[i__3].imag, q__1.imag = temp.real * b[i__3].imag + temp.imag * b[ i__3].real;
                            b[i__2].real = q__1.real, b[i__2].imag = q__1.imag;
                            /* L310: */
                        }
                    }
                    /* L320: */
                }
            }
        }
    }
    return 0;
    /* End of CTRMM . */
}
/* ctrmm_ */

