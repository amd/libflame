/* zsyr2k.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zsyr2k_(char *uplo, char *trans, integer *n, integer *k, dcomplex *alpha, dcomplex *a, integer *lda, dcomplex * b, integer *ldb, dcomplex *beta, dcomplex *c__, integer * ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    dcomplex z__1, z__2, z__3, z__4, z__5;
    /* Local variables */
    integer info;
    dcomplex temp1, temp2;
    integer i__, j, l;
    extern logical lsame_(char *, char *, integer, integer);
    integer nrowa;
    logical upper;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* ZSYR2K performs one of the symmetric rank 2k operations */
    /* C := alpha*A*B' + alpha*B*A' + beta*C, */
    /* or */
    /* C := alpha*A'*B + alpha*B'*A + beta*C, */
    /* where alpha and beta are scalars, C is an n by n symmetric matrix */
    /* and A and B are n by k matrices in the first case and k by n */
    /* matrices in the second case. */
    /* Parameters */
    /* ========== */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the upper or lower */
    /* triangular part of the array C is to be referenced as */
    /* follows: */
    /* UPLO = 'U' or 'u' Only the upper triangular part of C */
    /* is to be referenced. */
    /* UPLO = 'L' or 'l' Only the lower triangular part of C */
    /* is to be referenced. */
    /* Unchanged on exit. */
    /* TRANS - CHARACTER*1. */
    /* On entry, TRANS specifies the operation to be performed as */
    /* follows: */
    /* TRANS = 'N' or 'n' C := alpha*A*B' + alpha*B*A' + */
    /* beta*C. */
    /* TRANS = 'T' or 't' C := alpha*A'*B + alpha*B'*A + */
    /* beta*C. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the order of the matrix C. N must be */
    /* at least zero. */
    /* Unchanged on exit. */
    /* K - INTEGER. */
    /* On entry with TRANS = 'N' or 'n', K specifies the number */
    /* of columns of the matrices A and B, and on entry with */
    /* TRANS = 'T' or 't', K specifies the number of rows of the */
    /* matrices A and B. K must be at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX*16 . */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* A - COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is */
    /* k when TRANS = 'N' or 'n', and is n otherwise. */
    /* Before entry with TRANS = 'N' or 'n', the leading n by k */
    /* part of the array A must contain the matrix A, otherwise */
    /* the leading k by n part of the array A must contain the */
    /* matrix A. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. When TRANS = 'N' or 'n' */
    /* then LDA must be at least fla_max( 1, n ), otherwise LDA must */
    /* be at least fla_max( 1, k ). */
    /* Unchanged on exit. */
    /* B - COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is */
    /* k when TRANS = 'N' or 'n', and is n otherwise. */
    /* Before entry with TRANS = 'N' or 'n', the leading n by k */
    /* part of the array B must contain the matrix B, otherwise */
    /* the leading k by n part of the array B must contain the */
    /* matrix B. */
    /* Unchanged on exit. */
    /* LDB - INTEGER. */
    /* On entry, LDB specifies the first dimension of B as declared */
    /* in the calling (sub) program. When TRANS = 'N' or 'n' */
    /* then LDB must be at least fla_max( 1, n ), otherwise LDB must */
    /* be at least fla_max( 1, k ). */
    /* Unchanged on exit. */
    /* BETA - COMPLEX*16 . */
    /* On entry, BETA specifies the scalar beta. */
    /* Unchanged on exit. */
    /* C - COMPLEX*16 array of DIMENSION ( LDC, n ). */
    /* Before entry with UPLO = 'U' or 'u', the leading n by n */
    /* upper triangular part of the array C must contain the upper */
    /* triangular part of the symmetric matrix and the strictly */
    /* lower triangular part of C is not referenced. On exit, the */
    /* upper triangular part of the array C is overwritten by the */
    /* upper triangular part of the updated matrix. */
    /* Before entry with UPLO = 'L' or 'l', the leading n by n */
    /* lower triangular part of the array C must contain the lower */
    /* triangular part of the symmetric matrix and the strictly */
    /* upper triangular part of C is not referenced. On exit, the */
    /* lower triangular part of the array C is overwritten by the */
    /* lower triangular part of the updated matrix. */
    /* LDC - INTEGER. */
    /* On entry, LDC specifies the first dimension of C as declared */
    /* in the calling (sub) program. LDC must be at least */
    /* fla_max( 1, n ). */
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
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    /* Function Body */
    if (lsame_(trans, "N", 1, 1))
    {
        nrowa = *n;
    }
    else
    {
        nrowa = *k;
    }
    upper = lsame_(uplo, "U", 1, 1);
    info = 0;
    if (! upper && ! lsame_(uplo, "L", 1, 1))
    {
        info = 1;
    }
    else if (! lsame_(trans, "N", 1, 1) && ! lsame_(trans, "T", 1, 1))
    {
        info = 2;
    }
    else if (*n < 0)
    {
        info = 3;
    }
    else if (*k < 0)
    {
        info = 4;
    }
    else if (*lda < fla_max(1,nrowa))
    {
        info = 7;
    }
    else if (*ldb < fla_max(1,nrowa))
    {
        info = 9;
    }
    else if (*ldc < fla_max(1,*n))
    {
        info = 12;
    }
    if (info != 0)
    {
        xerbla_("ZSYR2K", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || (alpha->real == 0. && alpha->imag == 0. || *k == 0) && (beta->real == 1. && beta->imag == 0.))
    {
        return 0;
    }
    /* And when alpha.eq.zero. */
    if (alpha->real == 0. && alpha->imag == 0.)
    {
        if (upper)
        {
            if (beta->real == 0. && beta->imag == 0.)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].real = 0., c__[i__3].imag = 0.;
                        /* L10: */
                    }
                    /* L20: */
                }
            }
            else
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.real = beta->real * c__[i__4].real - beta->imag * c__[i__4] .imag, z__1.imag = beta->real * c__[i__4].imag + beta->imag * c__[i__4].real;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                        /* L30: */
                    }
                    /* L40: */
                }
            }
        }
        else
        {
            if (beta->real == 0. && beta->imag == 0.)
            {
                i__1 = *n;
                for (j = 1;
                        j <= i__1;
                        ++j)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].real = 0., c__[i__3].imag = 0.;
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
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.real = beta->real * c__[i__4].real - beta->imag * c__[i__4] .imag, z__1.imag = beta->real * c__[i__4].imag + beta->imag * c__[i__4].real;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                        /* L70: */
                    }
                    /* L80: */
                }
            }
        }
        return 0;
    }
    /* Start the operations. */
    if (lsame_(trans, "N", 1, 1))
    {
        /* Form C := alpha*A*B' + alpha*B*A' + C. */
        if (upper)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (beta->real == 0. && beta->imag == 0.)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].real = 0., c__[i__3].imag = 0.;
                        /* L90: */
                    }
                }
                else if (beta->real != 1. || beta->imag != 0.)
                {
                    i__2 = j;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.real = beta->real * c__[i__4].real - beta->imag * c__[i__4] .imag, z__1.imag = beta->real * c__[i__4].imag + beta->imag * c__[i__4].real;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                        /* L100: */
                    }
                }
                i__2 = *k;
                for (l = 1;
                        l <= i__2;
                        ++l)
                {
                    i__3 = j + l * a_dim1;
                    i__4 = j + l * b_dim1;
                    if (a[i__3].real != 0. || a[i__3].imag != 0. || (b[i__4].real != 0. || b[i__4].imag != 0.))
                    {
                        i__3 = j + l * b_dim1;
                        z__1.real = alpha->real * b[i__3].real - alpha->imag * b[i__3].imag, z__1.imag = alpha->real * b[i__3].imag + alpha->imag * b[ i__3].real;
                        temp1.real = z__1.real, temp1.imag = z__1.imag;
                        i__3 = j + l * a_dim1;
                        z__1.real = alpha->real * a[i__3].real - alpha->imag * a[i__3].imag, z__1.imag = alpha->real * a[i__3].imag + alpha->imag * a[ i__3].real;
                        temp2.real = z__1.real, temp2.imag = z__1.imag;
                        i__3 = j;
                        for (i__ = 1;
                                i__ <= i__3;
                                ++i__)
                        {
                            i__4 = i__ + j * c_dim1;
                            i__5 = i__ + j * c_dim1;
                            i__6 = i__ + l * a_dim1;
                            z__3.real = a[i__6].real * temp1.real - a[i__6].imag * temp1.imag, z__3.imag = a[i__6].real * temp1.imag + a[ i__6].imag * temp1.real;
                            z__2.real = c__[i__5].real + z__3.real, z__2.imag = c__[i__5] .imag + z__3.imag;
                            i__7 = i__ + l * b_dim1;
                            z__4.real = b[i__7].real * temp2.real - b[i__7].imag * temp2.imag, z__4.imag = b[i__7].real * temp2.imag + b[ i__7].imag * temp2.real;
                            z__1.real = z__2.real + z__4.real, z__1.imag = z__2.imag + z__4.imag;
                            c__[i__4].real = z__1.real, c__[i__4].imag = z__1.imag;
                            /* L110: */
                        }
                    }
                    /* L120: */
                }
                /* L130: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                if (beta->real == 0. && beta->imag == 0.)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        c__[i__3].real = 0., c__[i__3].imag = 0.;
                        /* L140: */
                    }
                }
                else if (beta->real != 1. || beta->imag != 0.)
                {
                    i__2 = *n;
                    for (i__ = j;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__1.real = beta->real * c__[i__4].real - beta->imag * c__[i__4] .imag, z__1.imag = beta->real * c__[i__4].imag + beta->imag * c__[i__4].real;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                        /* L150: */
                    }
                }
                i__2 = *k;
                for (l = 1;
                        l <= i__2;
                        ++l)
                {
                    i__3 = j + l * a_dim1;
                    i__4 = j + l * b_dim1;
                    if (a[i__3].real != 0. || a[i__3].imag != 0. || (b[i__4].real != 0. || b[i__4].imag != 0.))
                    {
                        i__3 = j + l * b_dim1;
                        z__1.real = alpha->real * b[i__3].real - alpha->imag * b[i__3].imag, z__1.imag = alpha->real * b[i__3].imag + alpha->imag * b[ i__3].real;
                        temp1.real = z__1.real, temp1.imag = z__1.imag;
                        i__3 = j + l * a_dim1;
                        z__1.real = alpha->real * a[i__3].real - alpha->imag * a[i__3].imag, z__1.imag = alpha->real * a[i__3].imag + alpha->imag * a[ i__3].real;
                        temp2.real = z__1.real, temp2.imag = z__1.imag;
                        i__3 = *n;
                        for (i__ = j;
                                i__ <= i__3;
                                ++i__)
                        {
                            i__4 = i__ + j * c_dim1;
                            i__5 = i__ + j * c_dim1;
                            i__6 = i__ + l * a_dim1;
                            z__3.real = a[i__6].real * temp1.real - a[i__6].imag * temp1.imag, z__3.imag = a[i__6].real * temp1.imag + a[ i__6].imag * temp1.real;
                            z__2.real = c__[i__5].real + z__3.real, z__2.imag = c__[i__5] .imag + z__3.imag;
                            i__7 = i__ + l * b_dim1;
                            z__4.real = b[i__7].real * temp2.real - b[i__7].imag * temp2.imag, z__4.imag = b[i__7].real * temp2.imag + b[ i__7].imag * temp2.real;
                            z__1.real = z__2.real + z__4.real, z__1.imag = z__2.imag + z__4.imag;
                            c__[i__4].real = z__1.real, c__[i__4].imag = z__1.imag;
                            /* L160: */
                        }
                    }
                    /* L170: */
                }
                /* L180: */
            }
        }
    }
    else
    {
        /* Form C := alpha*A'*B + alpha*B'*A + C. */
        if (upper)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    temp1.real = 0., temp1.imag = 0.;
                    temp2.real = 0., temp2.imag = 0.;
                    i__3 = *k;
                    for (l = 1;
                            l <= i__3;
                            ++l)
                    {
                        i__4 = l + i__ * a_dim1;
                        i__5 = l + j * b_dim1;
                        z__2.real = a[i__4].real * b[i__5].real - a[i__4].imag * b[i__5] .imag, z__2.imag = a[i__4].real * b[i__5].imag + a[i__4] .imag * b[i__5].real;
                        z__1.real = temp1.real + z__2.real, z__1.imag = temp1.imag + z__2.imag;
                        temp1.real = z__1.real, temp1.imag = z__1.imag;
                        i__4 = l + i__ * b_dim1;
                        i__5 = l + j * a_dim1;
                        z__2.real = b[i__4].real * a[i__5].real - b[i__4].imag * a[i__5] .imag, z__2.imag = b[i__4].real * a[i__5].imag + b[i__4] .imag * a[i__5].real;
                        z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                        temp2.real = z__1.real, temp2.imag = z__1.imag;
                        /* L190: */
                    }
                    if (beta->real == 0. && beta->imag == 0.)
                    {
                        i__3 = i__ + j * c_dim1;
                        z__2.real = alpha->real * temp1.real - alpha->imag * temp1.imag, z__2.imag = alpha->real * temp1.imag + alpha->imag * temp1.real;
                        z__3.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__3.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    }
                    else
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__3.real = beta->real * c__[i__4].real - beta->imag * c__[i__4] .imag, z__3.imag = beta->real * c__[i__4].imag + beta->imag * c__[i__4].real;
                        z__4.real = alpha->real * temp1.real - alpha->imag * temp1.imag, z__4.imag = alpha->real * temp1.imag + alpha->imag * temp1.real;
                        z__2.real = z__3.real + z__4.real, z__2.imag = z__3.imag + z__4.imag;
                        z__5.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__5.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__5.real, z__1.imag = z__2.imag + z__5.imag;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    }
                    /* L200: */
                }
                /* L210: */
            }
        }
        else
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = *n;
                for (i__ = j;
                        i__ <= i__2;
                        ++i__)
                {
                    temp1.real = 0., temp1.imag = 0.;
                    temp2.real = 0., temp2.imag = 0.;
                    i__3 = *k;
                    for (l = 1;
                            l <= i__3;
                            ++l)
                    {
                        i__4 = l + i__ * a_dim1;
                        i__5 = l + j * b_dim1;
                        z__2.real = a[i__4].real * b[i__5].real - a[i__4].imag * b[i__5] .imag, z__2.imag = a[i__4].real * b[i__5].imag + a[i__4] .imag * b[i__5].real;
                        z__1.real = temp1.real + z__2.real, z__1.imag = temp1.imag + z__2.imag;
                        temp1.real = z__1.real, temp1.imag = z__1.imag;
                        i__4 = l + i__ * b_dim1;
                        i__5 = l + j * a_dim1;
                        z__2.real = b[i__4].real * a[i__5].real - b[i__4].imag * a[i__5] .imag, z__2.imag = b[i__4].real * a[i__5].imag + b[i__4] .imag * a[i__5].real;
                        z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                        temp2.real = z__1.real, temp2.imag = z__1.imag;
                        /* L220: */
                    }
                    if (beta->real == 0. && beta->imag == 0.)
                    {
                        i__3 = i__ + j * c_dim1;
                        z__2.real = alpha->real * temp1.real - alpha->imag * temp1.imag, z__2.imag = alpha->real * temp1.imag + alpha->imag * temp1.real;
                        z__3.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__3.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    }
                    else
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__3.real = beta->real * c__[i__4].real - beta->imag * c__[i__4] .imag, z__3.imag = beta->real * c__[i__4].imag + beta->imag * c__[i__4].real;
                        z__4.real = alpha->real * temp1.real - alpha->imag * temp1.imag, z__4.imag = alpha->real * temp1.imag + alpha->imag * temp1.real;
                        z__2.real = z__3.real + z__4.real, z__2.imag = z__3.imag + z__4.imag;
                        z__5.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__5.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__5.real, z__1.imag = z__2.imag + z__5.imag;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    }
                    /* L230: */
                }
                /* L240: */
            }
        }
    }
    return 0;
    /* End of ZSYR2K. */
}
/* zsyr2k_ */

