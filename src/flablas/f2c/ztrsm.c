/* ztrsm.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Table of constant values */
static dcomplex c_b1 =
{
    {1., 0.}
}
;
/* Subroutine */
int ztrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, dcomplex *alpha, dcomplex *a, integer *lda, dcomplex *b, integer *ldb)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    void z_div(dcomplex *, dcomplex *, dcomplex *), d_cnjg( dcomplex *, dcomplex *);
    /* Local variables */
    integer info;
    dcomplex temp;
    integer i__, j, k;
    logical lside;
    extern logical lsame_(char *, char *, integer, integer);
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
    /* ZTRSM solves one of the matrix equations */
    /* op( A )*X = alpha*B, or X*op( A ) = alpha*B, */
    /* where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
    /* non-unit, upper or lower triangular matrix and op( A ) is one of */
    /* op( A ) = A or op( A ) = A' or op( A ) = conjg( A' ). */
    /* The matrix X is overwritten on B. */
    /* Parameters */
    /* ========== */
    /* SIDE - CHARACTER*1. */
    /* On entry, SIDE specifies whether op( A ) appears on the left */
    /* or right of X as follows: */
    /* SIDE = 'L' or 'l' op( A )*X = alpha*B. */
    /* SIDE = 'R' or 'r' X*op( A ) = alpha*B. */
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
    /* ALPHA - COMPLEX*16 . */
    /* On entry, ALPHA specifies the scalar alpha. When alpha is */
    /* zero then A is not referenced and B need not be set before */
    /* entry. */
    /* Unchanged on exit. */
    /* A - COMPLEX*16 array of DIMENSION ( LDA, k ), where k is m */
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
    /* B - COMPLEX*16 array of DIMENSION ( LDB, n ). */
    /* Before entry, the leading m by n part of the array B must */
    /* contain the right-hand side matrix B, and on exit is */
    /* overwritten by the solution matrix X. */
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
        xerbla_("ZTRSM ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0)
    {
        return 0;
    }
    /* And when alpha.eq.zero. */
    if (alpha->real == 0. && alpha->imag == 0.)
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
                b[i__3].real = 0., b[i__3].imag = 0.;
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
            /* Form B := alpha*inv( A )*B. */
            if (upper)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    if (alpha->real != 1. || alpha->imag != 0.)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * b_dim1;
                            i__4 = i__ + j * b_dim1;
                            z__1.real = alpha->real * b[i__4].real - alpha->imag * b[i__4] .imag, z__1.imag = alpha->real * b[i__4].imag + alpha->imag * b[i__4].real;
                            b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                            /* L30: */
                        }
                    }
                    for (k = *m;
                            k >= 1;
                            --k)
                    {
                        i__2 = k + j * b_dim1;
                        if (b[i__2].real != 0. || b[i__2].imag != 0.)
                        {
                            if (nounit)
                            {
                                i__2 = k + j * b_dim1;
                                z_div(&z__1, &b[k + j * b_dim1], &a[k + k * a_dim1]);
                                b[i__2].real = z__1.real, b[i__2].imag = z__1.imag;
                            }
                            i__2 = k - 1;
                            for (i__ = 1;
                                    i__ <= i__2;
                                    ++i__)
                            {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = k + j * b_dim1;
                                i__6 = i__ + k * a_dim1;
                                z__2.real = b[i__5].real * a[i__6].real - b[i__5].imag * a[i__6].imag, z__2.imag = b[i__5].real * a[ i__6].imag + b[i__5].imag * a[i__6].real;
                                z__1.real = b[i__4].real - z__2.real, z__1.imag = b[i__4] .imag - z__2.imag;
                                b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                                /* L40: */
                            }
                        }
                        /* L50: */
                    }
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
                    if (alpha->real != 1. || alpha->imag != 0.)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * b_dim1;
                            i__4 = i__ + j * b_dim1;
                            z__1.real = alpha->real * b[i__4].real - alpha->imag * b[i__4] .imag, z__1.imag = alpha->real * b[i__4].imag + alpha->imag * b[i__4].real;
                            b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                            /* L70: */
                        }
                    }
                    i__2 = *m;
                    for (k = 1;
                            k <= i__2;
                            ++k)
                    {
                        i__3 = k + j * b_dim1;
                        if (b[i__3].real != 0. || b[i__3].imag != 0.)
                        {
                            if (nounit)
                            {
                                i__3 = k + j * b_dim1;
                                z_div(&z__1, &b[k + j * b_dim1], &a[k + k * a_dim1]);
                                b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                            }
                            i__3 = *m;
                            for (i__ = k + 1;
                                    i__ <= i__3;
                                    ++i__)
                            {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = k + j * b_dim1;
                                i__7 = i__ + k * a_dim1;
                                z__2.real = b[i__6].real * a[i__7].real - b[i__6].imag * a[i__7].imag, z__2.imag = b[i__6].real * a[ i__7].imag + b[i__6].imag * a[i__7].real;
                                z__1.real = b[i__5].real - z__2.real, z__1.imag = b[i__5] .imag - z__2.imag;
                                b[i__4].real = z__1.real, b[i__4].imag = z__1.imag;
                                /* L80: */
                            }
                        }
                        /* L90: */
                    }
                    /* L100: */
                }
            }
        }
        else
        {
            /* Form B := alpha*inv( A' )*B */
            /* or B := alpha*inv( conjg( A' ) )*B. */
            if (upper)
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
                        z__1.real = alpha->real * b[i__3].real - alpha->imag * b[i__3].imag, z__1.imag = alpha->real * b[i__3].imag + alpha->imag * b[ i__3].real;
                        temp.real = z__1.real, temp.imag = z__1.imag;
                        if (noconj)
                        {
                            i__3 = i__ - 1;
                            for (k = 1;
                                    k <= i__3;
                                    ++k)
                            {
                                i__4 = k + i__ * a_dim1;
                                i__5 = k + j * b_dim1;
                                z__2.real = a[i__4].real * b[i__5].real - a[i__4].imag * b[i__5].imag, z__2.imag = a[i__4].real * b[ i__5].imag + a[i__4].imag * b[i__5].real;
                                z__1.real = temp.real - z__2.real, z__1.imag = temp.imag - z__2.imag;
                                temp.real = z__1.real, temp.imag = z__1.imag;
                                /* L110: */
                            }
                            if (nounit)
                            {
                                z_div(&z__1, &temp, &a[i__ + i__ * a_dim1]);
                                temp.real = z__1.real, temp.imag = z__1.imag;
                            }
                        }
                        else
                        {
                            i__3 = i__ - 1;
                            for (k = 1;
                                    k <= i__3;
                                    ++k)
                            {
                                d_cnjg(&z__3, &a[k + i__ * a_dim1]);
                                i__4 = k + j * b_dim1;
                                z__2.real = z__3.real * b[i__4].real - z__3.imag * b[i__4] .imag, z__2.imag = z__3.real * b[i__4].imag + z__3.imag * b[i__4].real;
                                z__1.real = temp.real - z__2.real, z__1.imag = temp.imag - z__2.imag;
                                temp.real = z__1.real, temp.imag = z__1.imag;
                                /* L120: */
                            }
                            if (nounit)
                            {
                                d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
                                z_div(&z__1, &temp, &z__2);
                                temp.real = z__1.real, temp.imag = z__1.imag;
                            }
                        }
                        i__3 = i__ + j * b_dim1;
                        b[i__3].real = temp.real, b[i__3].imag = temp.imag;
                        /* L130: */
                    }
                    /* L140: */
                }
            }
            else
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
                        z__1.real = alpha->real * b[i__2].real - alpha->imag * b[i__2].imag, z__1.imag = alpha->real * b[i__2].imag + alpha->imag * b[ i__2].real;
                        temp.real = z__1.real, temp.imag = z__1.imag;
                        if (noconj)
                        {
                            i__2 = *m;
                            for (k = i__ + 1;
                                    k <= i__2;
                                    ++k)
                            {
                                i__3 = k + i__ * a_dim1;
                                i__4 = k + j * b_dim1;
                                z__2.real = a[i__3].real * b[i__4].real - a[i__3].imag * b[i__4].imag, z__2.imag = a[i__3].real * b[ i__4].imag + a[i__3].imag * b[i__4].real;
                                z__1.real = temp.real - z__2.real, z__1.imag = temp.imag - z__2.imag;
                                temp.real = z__1.real, temp.imag = z__1.imag;
                                /* L150: */
                            }
                            if (nounit)
                            {
                                z_div(&z__1, &temp, &a[i__ + i__ * a_dim1]);
                                temp.real = z__1.real, temp.imag = z__1.imag;
                            }
                        }
                        else
                        {
                            i__2 = *m;
                            for (k = i__ + 1;
                                    k <= i__2;
                                    ++k)
                            {
                                d_cnjg(&z__3, &a[k + i__ * a_dim1]);
                                i__3 = k + j * b_dim1;
                                z__2.real = z__3.real * b[i__3].real - z__3.imag * b[i__3] .imag, z__2.imag = z__3.real * b[i__3].imag + z__3.imag * b[i__3].real;
                                z__1.real = temp.real - z__2.real, z__1.imag = temp.imag - z__2.imag;
                                temp.real = z__1.real, temp.imag = z__1.imag;
                                /* L160: */
                            }
                            if (nounit)
                            {
                                d_cnjg(&z__2, &a[i__ + i__ * a_dim1]);
                                z_div(&z__1, &temp, &z__2);
                                temp.real = z__1.real, temp.imag = z__1.imag;
                            }
                        }
                        i__2 = i__ + j * b_dim1;
                        b[i__2].real = temp.real, b[i__2].imag = temp.imag;
                        /* L170: */
                    }
                    /* L180: */
                }
            }
        }
    }
    else
    {
        if (lsame_(transa, "N", 1, 1))
        {
            /* Form B := alpha*B*inv( A ). */
            if (upper)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    if (alpha->real != 1. || alpha->imag != 0.)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * b_dim1;
                            i__4 = i__ + j * b_dim1;
                            z__1.real = alpha->real * b[i__4].real - alpha->imag * b[i__4] .imag, z__1.imag = alpha->real * b[i__4].imag + alpha->imag * b[i__4].real;
                            b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                            /* L190: */
                        }
                    }
                    i__2 = j - 1;
                    for (k = 1;
                            k <= i__2;
                            ++k)
                    {
                        i__3 = k + j * a_dim1;
                        if (a[i__3].real != 0. || a[i__3].imag != 0.)
                        {
                            i__3 = *m;
                            for (i__ = 1;
                                    i__ <= i__3;
                                    ++i__)
                            {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = k + j * a_dim1;
                                i__7 = i__ + k * b_dim1;
                                z__2.real = a[i__6].real * b[i__7].real - a[i__6].imag * b[i__7].imag, z__2.imag = a[i__6].real * b[ i__7].imag + a[i__6].imag * b[i__7].real;
                                z__1.real = b[i__5].real - z__2.real, z__1.imag = b[i__5] .imag - z__2.imag;
                                b[i__4].real = z__1.real, b[i__4].imag = z__1.imag;
                                /* L200: */
                            }
                        }
                        /* L210: */
                    }
                    if (nounit)
                    {
                        z_div(&z__1, &c_b1, &a[j + j * a_dim1]);
                        temp.real = z__1.real, temp.imag = z__1.imag;
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + j * b_dim1;
                            i__4 = i__ + j * b_dim1;
                            z__1.real = temp.real * b[i__4].real - temp.imag * b[i__4].imag, z__1.imag = temp.real * b[i__4].imag + temp.imag * b[ i__4].real;
                            b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                            /* L220: */
                        }
                    }
                    /* L230: */
                }
            }
            else
            {
                for (j = *n;
                        j >= 1;
                        --j)
                {
                    if (alpha->real != 1. || alpha->imag != 0.)
                    {
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + j * b_dim1;
                            i__3 = i__ + j * b_dim1;
                            z__1.real = alpha->real * b[i__3].real - alpha->imag * b[i__3] .imag, z__1.imag = alpha->real * b[i__3].imag + alpha->imag * b[i__3].real;
                            b[i__2].real = z__1.real, b[i__2].imag = z__1.imag;
                            /* L240: */
                        }
                    }
                    i__1 = *n;
                    for (k = j + 1;
                            k <= i__1;
                            ++k)
                    {
                        i__2 = k + j * a_dim1;
                        if (a[i__2].real != 0. || a[i__2].imag != 0.)
                        {
                            i__2 = *m;
                            for (i__ = 1;
                                    i__ <= i__2;
                                    ++i__)
                            {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = k + j * a_dim1;
                                i__6 = i__ + k * b_dim1;
                                z__2.real = a[i__5].real * b[i__6].real - a[i__5].imag * b[i__6].imag, z__2.imag = a[i__5].real * b[ i__6].imag + a[i__5].imag * b[i__6].real;
                                z__1.real = b[i__4].real - z__2.real, z__1.imag = b[i__4] .imag - z__2.imag;
                                b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                                /* L250: */
                            }
                        }
                        /* L260: */
                    }
                    if (nounit)
                    {
                        z_div(&z__1, &c_b1, &a[j + j * a_dim1]);
                        temp.real = z__1.real, temp.imag = z__1.imag;
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + j * b_dim1;
                            i__3 = i__ + j * b_dim1;
                            z__1.real = temp.real * b[i__3].real - temp.imag * b[i__3].imag, z__1.imag = temp.real * b[i__3].imag + temp.imag * b[ i__3].real;
                            b[i__2].real = z__1.real, b[i__2].imag = z__1.imag;
                            /* L270: */
                        }
                    }
                    /* L280: */
                }
            }
        }
        else
        {
            /* Form B := alpha*B*inv( A' ) */
            /* or B := alpha*B*inv( conjg( A' ) ). */
            if (upper)
            {
                for (k = *n;
                        k >= 1;
                        --k)
                {
                    if (nounit)
                    {
                        if (noconj)
                        {
                            z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
                            temp.real = z__1.real, temp.imag = z__1.imag;
                        }
                        else
                        {
                            d_cnjg(&z__2, &a[k + k * a_dim1]);
                            z_div(&z__1, &c_b1, &z__2);
                            temp.real = z__1.real, temp.imag = z__1.imag;
                        }
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + k * b_dim1;
                            i__3 = i__ + k * b_dim1;
                            z__1.real = temp.real * b[i__3].real - temp.imag * b[i__3].imag, z__1.imag = temp.real * b[i__3].imag + temp.imag * b[ i__3].real;
                            b[i__2].real = z__1.real, b[i__2].imag = z__1.imag;
                            /* L290: */
                        }
                    }
                    i__1 = k - 1;
                    for (j = 1;
                            j <= i__1;
                            ++j)
                    {
                        i__2 = j + k * a_dim1;
                        if (a[i__2].real != 0. || a[i__2].imag != 0.)
                        {
                            if (noconj)
                            {
                                i__2 = j + k * a_dim1;
                                temp.real = a[i__2].real, temp.imag = a[i__2].imag;
                            }
                            else
                            {
                                d_cnjg(&z__1, &a[j + k * a_dim1]);
                                temp.real = z__1.real, temp.imag = z__1.imag;
                            }
                            i__2 = *m;
                            for (i__ = 1;
                                    i__ <= i__2;
                                    ++i__)
                            {
                                i__3 = i__ + j * b_dim1;
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + k * b_dim1;
                                z__2.real = temp.real * b[i__5].real - temp.imag * b[i__5] .imag, z__2.imag = temp.real * b[i__5].imag + temp.imag * b[i__5].real;
                                z__1.real = b[i__4].real - z__2.real, z__1.imag = b[i__4] .imag - z__2.imag;
                                b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                                /* L300: */
                            }
                        }
                        /* L310: */
                    }
                    if (alpha->real != 1. || alpha->imag != 0.)
                    {
                        i__1 = *m;
                        for (i__ = 1;
                                i__ <= i__1;
                                ++i__)
                        {
                            i__2 = i__ + k * b_dim1;
                            i__3 = i__ + k * b_dim1;
                            z__1.real = alpha->real * b[i__3].real - alpha->imag * b[i__3] .imag, z__1.imag = alpha->real * b[i__3].imag + alpha->imag * b[i__3].real;
                            b[i__2].real = z__1.real, b[i__2].imag = z__1.imag;
                            /* L320: */
                        }
                    }
                    /* L330: */
                }
            }
            else
            {
                i__1 = *n;
                for (k = 1;
                        k <= i__1;
                        ++k)
                {
                    if (nounit)
                    {
                        if (noconj)
                        {
                            z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
                            temp.real = z__1.real, temp.imag = z__1.imag;
                        }
                        else
                        {
                            d_cnjg(&z__2, &a[k + k * a_dim1]);
                            z_div(&z__1, &c_b1, &z__2);
                            temp.real = z__1.real, temp.imag = z__1.imag;
                        }
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + k * b_dim1;
                            i__4 = i__ + k * b_dim1;
                            z__1.real = temp.real * b[i__4].real - temp.imag * b[i__4].imag, z__1.imag = temp.real * b[i__4].imag + temp.imag * b[ i__4].real;
                            b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                            /* L340: */
                        }
                    }
                    i__2 = *n;
                    for (j = k + 1;
                            j <= i__2;
                            ++j)
                    {
                        i__3 = j + k * a_dim1;
                        if (a[i__3].real != 0. || a[i__3].imag != 0.)
                        {
                            if (noconj)
                            {
                                i__3 = j + k * a_dim1;
                                temp.real = a[i__3].real, temp.imag = a[i__3].imag;
                            }
                            else
                            {
                                d_cnjg(&z__1, &a[j + k * a_dim1]);
                                temp.real = z__1.real, temp.imag = z__1.imag;
                            }
                            i__3 = *m;
                            for (i__ = 1;
                                    i__ <= i__3;
                                    ++i__)
                            {
                                i__4 = i__ + j * b_dim1;
                                i__5 = i__ + j * b_dim1;
                                i__6 = i__ + k * b_dim1;
                                z__2.real = temp.real * b[i__6].real - temp.imag * b[i__6] .imag, z__2.imag = temp.real * b[i__6].imag + temp.imag * b[i__6].real;
                                z__1.real = b[i__5].real - z__2.real, z__1.imag = b[i__5] .imag - z__2.imag;
                                b[i__4].real = z__1.real, b[i__4].imag = z__1.imag;
                                /* L350: */
                            }
                        }
                        /* L360: */
                    }
                    if (alpha->real != 1. || alpha->imag != 0.)
                    {
                        i__2 = *m;
                        for (i__ = 1;
                                i__ <= i__2;
                                ++i__)
                        {
                            i__3 = i__ + k * b_dim1;
                            i__4 = i__ + k * b_dim1;
                            z__1.real = alpha->real * b[i__4].real - alpha->imag * b[i__4] .imag, z__1.imag = alpha->real * b[i__4].imag + alpha->imag * b[i__4].real;
                            b[i__3].real = z__1.real, b[i__3].imag = z__1.imag;
                            /* L370: */
                        }
                    }
                    /* L380: */
                }
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return 0;
    /* End of ZTRSM . */
}
/* ztrsm_ */

