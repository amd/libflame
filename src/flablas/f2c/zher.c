/* zher.f -- translated by f2c (version 19991025). You must link the resulting object file with the libraries: -lf2c -lm (in that order) */
#include "FLA_f2c.h"
/* Subroutine */
int zher_(char *uplo, integer *n, doublereal *alpha, dcomplex *x, integer *incx, dcomplex *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    dcomplex z__1, z__2;
    /* Builtin functions */
    void d_cnjg(dcomplex *, dcomplex *);
    /* Local variables */
    integer info;
    dcomplex temp;
    integer i__, j;
    extern logical lsame_(char *, char *, integer, integer);
    integer ix, jx, kx;
    extern /* Subroutine */
    int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    /* .. Scalar Arguments .. */
    /* .. Array Arguments .. */
    /* .. */
    /* Purpose */
    /* ======= */
    /* ZHER performs the hermitian rank 1 operation */
    /* A := alpha*x*conjg( x' ) + A, */
    /* where alpha is a real scalar, x is an n element vector and A is an */
    /* n by n hermitian matrix. */
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
    /* ALPHA - DOUBLE PRECISION. */
    /* On entry, ALPHA specifies the scalar alpha. */
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
    /* A - COMPLEX*16 array of DIMENSION ( LDA, n ). */
    /* Before entry with UPLO = 'U' or 'u', the leading n by n */
    /* upper triangular part of the array A must contain the upper */
    /* triangular part of the hermitian matrix and the strictly */
    /* lower triangular part of A is not referenced. On exit, the */
    /* upper triangular part of the array A is overwritten by the */
    /* upper triangular part of the updated matrix. */
    /* Before entry with UPLO = 'L' or 'l', the leading n by n */
    /* lower triangular part of the array A must contain the lower */
    /* triangular part of the hermitian matrix and the strictly */
    /* upper triangular part of A is not referenced. On exit, the */
    /* lower triangular part of the array A is overwritten by the */
    /* lower triangular part of the updated matrix. */
    /* Note that the imaginary parts of the diagonal elements need */
    /* not be set, they are assumed to be zero, and on exit they */
    /* are set to zero. */
    /* LDA - INTEGER. */
    /* On entry, LDA specifies the first dimension of A as declared */
    /* in the calling (sub) program. LDA must be at least */
    /* fla_max( 1, n ). */
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
    --x;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
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
    else if (*incx == 0)
    {
        info = 5;
    }
    else if (*lda < fla_max(1,*n))
    {
        info = 7;
    }
    if (info != 0)
    {
        xerbla_("ZHER ", &info, (ftnlen)5);
        return 0;
    }
    /* Quick return if possible. */
    if (*n == 0 || *alpha == 0.)
    {
        return 0;
    }
    /* Set the start point in X if the increment is not unity. */
    if (*incx <= 0)
    {
        kx = 1 - (*n - 1) * *incx;
    }
    else if (*incx != 1)
    {
        kx = 1;
    }
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through the triangular part */
    /* of A. */
    if (lsame_(uplo, "U", 1, 1))
    {
        /* Form A when A is stored in upper triangle. */
        if (*incx == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                if (x[i__2].real != 0. || x[i__2].imag != 0.)
                {
                    d_cnjg(&z__2, &x[j]);
                    z__1.real = *alpha * z__2.real, z__1.imag = *alpha * z__2.imag;
                    temp.real = z__1.real, temp.imag = z__1.imag;
                    i__2 = j - 1;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = i__;
                        z__2.real = x[i__5].real * temp.real - x[i__5].imag * temp.imag, z__2.imag = x[i__5].real * temp.imag + x[i__5].imag * temp.real;
                        z__1.real = a[i__4].real + z__2.real, z__1.imag = a[i__4].imag + z__2.imag;
                        a[i__3].real = z__1.real, a[i__3].imag = z__1.imag;
                        /* L10: */
                    }
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = j;
                    z__1.real = x[i__4].real * temp.real - x[i__4].imag * temp.imag, z__1.imag = x[i__4].real * temp.imag + x[i__4].imag * temp.real;
                    d__1 = a[i__3].real + z__1.real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
                }
                else
                {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
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
                i__2 = jx;
                if (x[i__2].real != 0. || x[i__2].imag != 0.)
                {
                    d_cnjg(&z__2, &x[jx]);
                    z__1.real = *alpha * z__2.real, z__1.imag = *alpha * z__2.imag;
                    temp.real = z__1.real, temp.imag = z__1.imag;
                    ix = kx;
                    i__2 = j - 1;
                    for (i__ = 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = ix;
                        z__2.real = x[i__5].real * temp.real - x[i__5].imag * temp.imag, z__2.imag = x[i__5].real * temp.imag + x[i__5].imag * temp.real;
                        z__1.real = a[i__4].real + z__2.real, z__1.imag = a[i__4].imag + z__2.imag;
                        a[i__3].real = z__1.real, a[i__3].imag = z__1.imag;
                        ix += *incx;
                        /* L30: */
                    }
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = jx;
                    z__1.real = x[i__4].real * temp.real - x[i__4].imag * temp.imag, z__1.imag = x[i__4].real * temp.imag + x[i__4].imag * temp.real;
                    d__1 = a[i__3].real + z__1.real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
                }
                else
                {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
                }
                jx += *incx;
                /* L40: */
            }
        }
    }
    else
    {
        /* Form A when A is stored in lower triangle. */
        if (*incx == 1)
        {
            i__1 = *n;
            for (j = 1;
                    j <= i__1;
                    ++j)
            {
                i__2 = j;
                if (x[i__2].real != 0. || x[i__2].imag != 0.)
                {
                    d_cnjg(&z__2, &x[j]);
                    z__1.real = *alpha * z__2.real, z__1.imag = *alpha * z__2.imag;
                    temp.real = z__1.real, temp.imag = z__1.imag;
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = j;
                    z__1.real = temp.real * x[i__4].real - temp.imag * x[i__4].imag, z__1.imag = temp.real * x[i__4].imag + temp.imag * x[i__4].real;
                    d__1 = a[i__3].real + z__1.real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
                    i__2 = *n;
                    for (i__ = j + 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = i__;
                        z__2.real = x[i__5].real * temp.real - x[i__5].imag * temp.imag, z__2.imag = x[i__5].real * temp.imag + x[i__5].imag * temp.real;
                        z__1.real = a[i__4].real + z__2.real, z__1.imag = a[i__4].imag + z__2.imag;
                        a[i__3].real = z__1.real, a[i__3].imag = z__1.imag;
                        /* L50: */
                    }
                }
                else
                {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
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
                if (x[i__2].real != 0. || x[i__2].imag != 0.)
                {
                    d_cnjg(&z__2, &x[jx]);
                    z__1.real = *alpha * z__2.real, z__1.imag = *alpha * z__2.imag;
                    temp.real = z__1.real, temp.imag = z__1.imag;
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    i__4 = jx;
                    z__1.real = temp.real * x[i__4].real - temp.imag * x[i__4].imag, z__1.imag = temp.real * x[i__4].imag + temp.imag * x[i__4].real;
                    d__1 = a[i__3].real + z__1.real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
                    ix = jx;
                    i__2 = *n;
                    for (i__ = j + 1;
                            i__ <= i__2;
                            ++i__)
                    {
                        ix += *incx;
                        i__3 = i__ + j * a_dim1;
                        i__4 = i__ + j * a_dim1;
                        i__5 = ix;
                        z__2.real = x[i__5].real * temp.real - x[i__5].imag * temp.imag, z__2.imag = x[i__5].real * temp.imag + x[i__5].imag * temp.real;
                        z__1.real = a[i__4].real + z__2.real, z__1.imag = a[i__4].imag + z__2.imag;
                        a[i__3].real = z__1.real, a[i__3].imag = z__1.imag;
                        /* L70: */
                    }
                }
                else
                {
                    i__2 = j + j * a_dim1;
                    i__3 = j + j * a_dim1;
                    d__1 = a[i__3].real;
                    a[i__2].real = d__1, a[i__2].imag = 0.;
                }
                jx += *incx;
                /* L80: */
            }
        }
    }
    return 0;
    /* End of ZHER . */
}
/* zher_ */

