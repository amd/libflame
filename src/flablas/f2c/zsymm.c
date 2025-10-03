/* zsymm.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zsymm_(char *side, char *uplo, integer *m, integer *n, dcomplex *alpha, dcomplex *a, integer *lda, dcomplex * b, integer *ldb, dcomplex *beta, dcomplex *c__, integer * ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    dcomplex z__1, z__2, z__3, z__4, z__5;
    /* Local variables */
    integer info;
    dcomplex temp1, temp2;
    integer i__, j, k;
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
    /* ZSYMM performs one of the matrix-matrix operations */
    /* C := alpha*A*B + beta*C, */
    /* or */
    /* C := alpha*B*A + beta*C, */
    /* where alpha and beta are scalars, A is a symmetric matrix and B and */
    /* C are m by n matrices. */
    /* Parameters */
    /* ========== */
    /* SIDE - CHARACTER*1. */
    /* On entry, SIDE specifies whether the symmetric matrix A */
    /* appears on the left or right in the operation as follows: */
    /* SIDE = 'L' or 'l' C := alpha*A*B + beta*C, */
    /* SIDE = 'R' or 'r' C := alpha*B*A + beta*C, */
    /* Unchanged on exit. */
    /* UPLO - CHARACTER*1. */
    /* On entry, UPLO specifies whether the upper or lower */
    /* triangular part of the symmetric matrix A is to be */
    /* referenced as follows: */
    /* UPLO = 'U' or 'u' Only the upper triangular part of the */
    /* symmetric matrix is to be referenced. */
    /* UPLO = 'L' or 'l' Only the lower triangular part of the */
    /* symmetric matrix is to be referenced. */
    /* Unchanged on exit. */
    /* M - INTEGER. */
    /* On entry, M specifies the number of rows of the matrix C. */
    /* M must be at least zero. */
    /* Unchanged on exit. */
    /* N - INTEGER. */
    /* On entry, N specifies the number of columns of the matrix C. */
    /* N must be at least zero. */
    /* Unchanged on exit. */
    /* ALPHA - COMPLEX*16 . */
    /* On entry, ALPHA specifies the scalar alpha. */
    /* Unchanged on exit. */
    /* A - COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is */
    /* m when SIDE = 'L' or 'l' and is n otherwise. */
    /* Before entry with SIDE = 'L' or 'l', the m by m part of */
    /* the array A must contain the symmetric matrix, such that */
    /* when UPLO = 'U' or 'u', the leading m by m upper triangular */
    /* part of the array A must contain the upper triangular part */
    /* of the symmetric matrix and the strictly lower triangular */
    /* part of A is not referenced, and when UPLO = 'L' or 'l', */
    /* the leading m by m lower triangular part of the array A */
    /* must contain the lower triangular part of the symmetric */
    /* matrix and the strictly upper triangular part of A is not */
    /* referenced. */
    /* Before entry with SIDE = 'R' or 'r', the n by n part of */
    /* the array A must contain the symmetric matrix, such that */
    /* when UPLO = 'U' or 'u', the leading n by n upper triangular */
    /* part of the array A must contain the upper triangular part */
    /* of the symmetric matrix and the strictly lower triangular */
    /* part of A is not referenced, and when UPLO = 'L' or 'l', */
    /* the leading n by n lower triangular part of the array A */
    /* must contain the lower triangular part of the symmetric */
    /* matrix and the strictly upper triangular part of A is not */
    /* referenced. */
    /* Unchanged on exit. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. When SIDE = 'L' or 'l' then */
    /* LDA must be at least fla_max( 1, m ), otherwise LDA must be at */
    /* least fla_max( 1, n ). */
    /* Unchanged on exit. */
    /* B - COMPLEX*16 array of DIMENSION ( LDB, n ). */
    /* Before entry, the leading m by n part of the array B must */
    /* contain the matrix B. */
    /* Unchanged on exit. */
    /* LDB - INTEGER. */
    /* On entry, LDB specifies the first dimension of B as declared */
    /* in the calling (sub) program. LDB must be at least */
    /* fla_max( 1, m ). */
    /* Unchanged on exit. */
    /* BETA - COMPLEX*16 . */
    /* On entry, BETA specifies the scalar beta. When BETA is */
    /* supplied as zero then C need not be set on input. */
    /* Unchanged on exit. */
    /* C - COMPLEX*16 array of DIMENSION ( LDC, n ). */
    /* Before entry, the leading m by n part of the array C must */
    /* contain the matrix C, except when beta is zero, in which */
    /* case C need not be set on entry. */
    /* On exit, the array C is overwritten by the m by n updated */
    /* matrix. */
    /* LDC - INTEGER. */
    /* On entry, LDC specifies the first dimension of C as declared */
    /* in the calling (sub) program. LDC must be at least */
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
    /* Set NROWA as the number of rows of A. */
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
    if (lsame_(side, "L", 1, 1))
    {
        nrowa = *m;
    }
    else
    {
        nrowa = *n;
    }
    upper = lsame_(uplo, "U", 1, 1);
    /* Test the input parameters. */
    info = 0;
    if (! lsame_(side, "L", 1, 1) && ! lsame_(side, "R", 1, 1))
    {
        info = 1;
    }
    else if (! upper && ! lsame_(uplo, "L", 1, 1))
    {
        info = 2;
    }
    else if (*m < 0)
    {
        info = 3;
    }
    else if (*n < 0)
    {
        info = 4;
    }
    else if (*lda < fla_max(1,nrowa))
    {
        info = 7;
    }
    else if (*ldb < fla_max(1,*m))
    {
        info = 9;
    }
    else if (*ldc < fla_max(1,*m))
    {
        info = 12;
    }
    if (info != 0)
    {
        xerbla_("ZSYMM ", &info, (ftnlen)6);
        return 0;
    }
    /* Quick return if possible. */
    if (*m == 0 || *n == 0 || alpha->real == 0. && alpha->imag == 0. && (beta->real == 1. && beta->imag == 0.))
    {
        return 0;
    }
    /* And when alpha.eq.zero. */
    if (alpha->real == 0. && alpha->imag == 0.)
    {
        if (beta->real == 0. && beta->imag == 0.)
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
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * c_dim1;
                    i__4 = i__ + j * c_dim1;
                    z__1.real = beta->real * c__[i__4].real - beta->imag * c__[i__4].imag, z__1.imag = beta->real * c__[i__4].imag + beta->imag * c__[ i__4].real;
                    c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    /* L30: */
                }
                /* L40: */
            }
        }
        return 0;
    }
    /* Start the operations. */
    if (lsame_(side, "L", 1, 1))
    {
        /* Form C := alpha*A*B + beta*C. */
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
                    z__1.real = alpha->real * b[i__3].real - alpha->imag * b[i__3].imag, z__1.imag = alpha->real * b[i__3].imag + alpha->imag * b[i__3] .real;
                    temp1.real = z__1.real, temp1.imag = z__1.imag;
                    temp2.real = 0., temp2.imag = 0.;
                    i__3 = i__ - 1;
                    for (k = 1;
                            k <= i__3;
                            ++k)
                    {
                        i__4 = k + j * c_dim1;
                        i__5 = k + j * c_dim1;
                        i__6 = k + i__ * a_dim1;
                        z__2.real = temp1.real * a[i__6].real - temp1.imag * a[i__6].imag, z__2.imag = temp1.real * a[i__6].imag + temp1.imag * a[ i__6].real;
                        z__1.real = c__[i__5].real + z__2.real, z__1.imag = c__[i__5].imag + z__2.imag;
                        c__[i__4].real = z__1.real, c__[i__4].imag = z__1.imag;
                        i__4 = k + j * b_dim1;
                        i__5 = k + i__ * a_dim1;
                        z__2.real = b[i__4].real * a[i__5].real - b[i__4].imag * a[i__5] .imag, z__2.imag = b[i__4].real * a[i__5].imag + b[i__4] .imag * a[i__5].real;
                        z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                        temp2.real = z__1.real, temp2.imag = z__1.imag;
                        /* L50: */
                    }
                    if (beta->real == 0. && beta->imag == 0.)
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + i__ * a_dim1;
                        z__2.real = temp1.real * a[i__4].real - temp1.imag * a[i__4].imag, z__2.imag = temp1.real * a[i__4].imag + temp1.imag * a[ i__4].real;
                        z__3.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__3.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    }
                    else
                    {
                        i__3 = i__ + j * c_dim1;
                        i__4 = i__ + j * c_dim1;
                        z__3.real = beta->real * c__[i__4].real - beta->imag * c__[i__4] .imag, z__3.imag = beta->real * c__[i__4].imag + beta->imag * c__[i__4].real;
                        i__5 = i__ + i__ * a_dim1;
                        z__4.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag, z__4.imag = temp1.real * a[i__5].imag + temp1.imag * a[ i__5].real;
                        z__2.real = z__3.real + z__4.real, z__2.imag = z__3.imag + z__4.imag;
                        z__5.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__5.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__5.real, z__1.imag = z__2.imag + z__5.imag;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    }
                    /* L60: */
                }
                /* L70: */
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
                    z__1.real = alpha->real * b[i__2].real - alpha->imag * b[i__2].imag, z__1.imag = alpha->real * b[i__2].imag + alpha->imag * b[i__2] .real;
                    temp1.real = z__1.real, temp1.imag = z__1.imag;
                    temp2.real = 0., temp2.imag = 0.;
                    i__2 = *m;
                    for (k = i__ + 1;
                            k <= i__2;
                            ++k)
                    {
                        i__3 = k + j * c_dim1;
                        i__4 = k + j * c_dim1;
                        i__5 = k + i__ * a_dim1;
                        z__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag, z__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[ i__5].real;
                        z__1.real = c__[i__4].real + z__2.real, z__1.imag = c__[i__4].imag + z__2.imag;
                        c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                        i__3 = k + j * b_dim1;
                        i__4 = k + i__ * a_dim1;
                        z__2.real = b[i__3].real * a[i__4].real - b[i__3].imag * a[i__4] .imag, z__2.imag = b[i__3].real * a[i__4].imag + b[i__3] .imag * a[i__4].real;
                        z__1.real = temp2.real + z__2.real, z__1.imag = temp2.imag + z__2.imag;
                        temp2.real = z__1.real, temp2.imag = z__1.imag;
                        /* L80: */
                    }
                    if (beta->real == 0. && beta->imag == 0.)
                    {
                        i__2 = i__ + j * c_dim1;
                        i__3 = i__ + i__ * a_dim1;
                        z__2.real = temp1.real * a[i__3].real - temp1.imag * a[i__3].imag, z__2.imag = temp1.real * a[i__3].imag + temp1.imag * a[ i__3].real;
                        z__3.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__3.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
                        c__[i__2].real = z__1.real, c__[i__2].imag = z__1.imag;
                    }
                    else
                    {
                        i__2 = i__ + j * c_dim1;
                        i__3 = i__ + j * c_dim1;
                        z__3.real = beta->real * c__[i__3].real - beta->imag * c__[i__3] .imag, z__3.imag = beta->real * c__[i__3].imag + beta->imag * c__[i__3].real;
                        i__4 = i__ + i__ * a_dim1;
                        z__4.real = temp1.real * a[i__4].real - temp1.imag * a[i__4].imag, z__4.imag = temp1.real * a[i__4].imag + temp1.imag * a[ i__4].real;
                        z__2.real = z__3.real + z__4.real, z__2.imag = z__3.imag + z__4.imag;
                        z__5.real = alpha->real * temp2.real - alpha->imag * temp2.imag, z__5.imag = alpha->real * temp2.imag + alpha->imag * temp2.real;
                        z__1.real = z__2.real + z__5.real, z__1.imag = z__2.imag + z__5.imag;
                        c__[i__2].real = z__1.real, c__[i__2].imag = z__1.imag;
                    }
                    /* L90: */
                }
                /* L100: */
            }
        }
    }
    else
    {
        /* Form C := alpha*B*A + beta*C. */
        i__1 = *n;
        for (j = 1;
                j <= i__1;
                ++j)
        {
            i__2 = j + j * a_dim1;
            z__1.real = alpha->real * a[i__2].real - alpha->imag * a[i__2].imag, z__1.imag = alpha->real * a[i__2].imag + alpha->imag * a[i__2].real;
            temp1.real = z__1.real, temp1.imag = z__1.imag;
            if (beta->real == 0. && beta->imag == 0.)
            {
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * c_dim1;
                    i__4 = i__ + j * b_dim1;
                    z__1.real = temp1.real * b[i__4].real - temp1.imag * b[i__4].imag, z__1.imag = temp1.real * b[i__4].imag + temp1.imag * b[i__4] .real;
                    c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    /* L110: */
                }
            }
            else
            {
                i__2 = *m;
                for (i__ = 1;
                        i__ <= i__2;
                        ++i__)
                {
                    i__3 = i__ + j * c_dim1;
                    i__4 = i__ + j * c_dim1;
                    z__2.real = beta->real * c__[i__4].real - beta->imag * c__[i__4].imag, z__2.imag = beta->real * c__[i__4].imag + beta->imag * c__[ i__4].real;
                    i__5 = i__ + j * b_dim1;
                    z__3.real = temp1.real * b[i__5].real - temp1.imag * b[i__5].imag, z__3.imag = temp1.real * b[i__5].imag + temp1.imag * b[i__5] .real;
                    z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
                    c__[i__3].real = z__1.real, c__[i__3].imag = z__1.imag;
                    /* L120: */
                }
            }
            i__2 = j - 1;
            for (k = 1;
                    k <= i__2;
                    ++k)
            {
                if (upper)
                {
                    i__3 = k + j * a_dim1;
                    z__1.real = alpha->real * a[i__3].real - alpha->imag * a[i__3].imag, z__1.imag = alpha->real * a[i__3].imag + alpha->imag * a[i__3] .real;
                    temp1.real = z__1.real, temp1.imag = z__1.imag;
                }
                else
                {
                    i__3 = j + k * a_dim1;
                    z__1.real = alpha->real * a[i__3].real - alpha->imag * a[i__3].imag, z__1.imag = alpha->real * a[i__3].imag + alpha->imag * a[i__3] .real;
                    temp1.real = z__1.real, temp1.imag = z__1.imag;
                }
                i__3 = *m;
                for (i__ = 1;
                        i__ <= i__3;
                        ++i__)
                {
                    i__4 = i__ + j * c_dim1;
                    i__5 = i__ + j * c_dim1;
                    i__6 = i__ + k * b_dim1;
                    z__2.real = temp1.real * b[i__6].real - temp1.imag * b[i__6].imag, z__2.imag = temp1.real * b[i__6].imag + temp1.imag * b[i__6] .real;
                    z__1.real = c__[i__5].real + z__2.real, z__1.imag = c__[i__5].imag + z__2.imag;
                    c__[i__4].real = z__1.real, c__[i__4].imag = z__1.imag;
                    /* L130: */
                }
                /* L140: */
            }
            i__2 = *n;
            for (k = j + 1;
                    k <= i__2;
                    ++k)
            {
                if (upper)
                {
                    i__3 = j + k * a_dim1;
                    z__1.real = alpha->real * a[i__3].real - alpha->imag * a[i__3].imag, z__1.imag = alpha->real * a[i__3].imag + alpha->imag * a[i__3] .real;
                    temp1.real = z__1.real, temp1.imag = z__1.imag;
                }
                else
                {
                    i__3 = k + j * a_dim1;
                    z__1.real = alpha->real * a[i__3].real - alpha->imag * a[i__3].imag, z__1.imag = alpha->real * a[i__3].imag + alpha->imag * a[i__3] .real;
                    temp1.real = z__1.real, temp1.imag = z__1.imag;
                }
                i__3 = *m;
                for (i__ = 1;
                        i__ <= i__3;
                        ++i__)
                {
                    i__4 = i__ + j * c_dim1;
                    i__5 = i__ + j * c_dim1;
                    i__6 = i__ + k * b_dim1;
                    z__2.real = temp1.real * b[i__6].real - temp1.imag * b[i__6].imag, z__2.imag = temp1.real * b[i__6].imag + temp1.imag * b[i__6] .real;
                    z__1.real = c__[i__5].real + z__2.real, z__1.imag = c__[i__5].imag + z__2.imag;
                    c__[i__4].real = z__1.real, c__[i__4].imag = z__1.imag;
                    /* L150: */
                }
                /* L160: */
            }
            /* L170: */
        }
    }
    return 0;
    /* End of ZSYMM . */
}
/* zsymm_ */

