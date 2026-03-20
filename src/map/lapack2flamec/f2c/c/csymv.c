/* ../netlib/csymv.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b CSYMV computes a matrix-vector product for a scomplex symmetric matrix. */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSYMV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csymv.f
 * "> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csymv.f
 * "> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csymv.f
 * "> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INCX, INCY, LDA, N */
/* COMPLEX ALPHA, BETA */
/* .. */
/* .. Array Arguments .. */
/* COMPLEX A( LDA, * ), X( * ), Y( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYMV performs the matrix-vector operation */
/* > */
/* > y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > On entry, UPLO specifies whether the upper or lower */
/* > triangular part of the array A is to be referenced as */
/* > follows: */
/* > */
/* > UPLO = 'U' or 'u' Only the upper triangular part of A */
/* > is to be referenced. */
/* > */
/* > UPLO = 'L' or 'l' Only the lower triangular part of A */
/* > is to be referenced. */
/* > */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > On entry, N specifies the order of the matrix A. */
/* > N must be at least zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* > ALPHA is COMPLEX */
/* > On entry, ALPHA specifies the scalar alpha. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension ( LDA, N ) */
/* > Before entry, with UPLO = 'U' or 'u', the leading n by n */
/* > upper triangular part of the array A must contain the upper */
/* > triangular part of the symmetric matrix and the strictly */
/* > lower triangular part of A is not referenced. */
/* > Before entry, with UPLO = 'L' or 'l', the leading n by n */
/* > lower triangular part of the array A must contain the lower */
/* > triangular part of the symmetric matrix and the strictly */
/* > upper triangular part of A is not referenced. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > On entry, LDA specifies the first dimension of A as declared */
/* > in the calling (sub) program. LDA must be at least */
/* > fla_max( 1, N ). */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension at least */
/* > ( 1 + ( N - 1 )*f2c_abs( INCX ) ). */
/* > Before entry, the incremented array X must contain the N- */
/* > element vector x. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* > INCX is INTEGER */
/* > On entry, INCX specifies the increment for the elements of */
/* > X. INCX must not be zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* > BETA is COMPLEX */
/* > On entry, BETA specifies the scalar beta. When BETA is */
/* > supplied as zero then Y need not be set on input. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y is COMPLEX array, dimension at least */
/* > ( 1 + ( N - 1 )*f2c_abs( INCY ) ). */
/* > Before entry, the incremented array Y must contain the n */
/* > element vector y. On exit, Y is overwritten by the updated */
/* > vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* > INCY is INTEGER */
/* > On entry, INCY specifies the increment for the elements of */
/* > Y. INCY must not be zero. */
/* > Unchanged on exit. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date September 2012 */
/* > \ingroup complexSYauxiliary */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void csymv_(char *uplo, aocl_int_t *n, scomplex *alpha, scomplex *a, aocl_int_t *lda, scomplex *x,
            aocl_int_t *incx, scomplex *beta, scomplex *y, aocl_int_t *incy)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_csymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t incx_64 = *incx;
    aocl_int64_t incy_64 = *incy;

    aocl_lapack_csymv(uplo, &n_64, alpha, a, &lda_64, x, &incx_64, beta, y, &incy_64);
#endif
}

void aocl_lapack_csymv(char *uplo, aocl_int64_t *n, scomplex *alpha, scomplex *a, aocl_int64_t *lda,
                       scomplex *x, aocl_int64_t *incx, scomplex *beta, scomplex *y,
                       aocl_int64_t *incy)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "csymv inputs: uplo %c, n %lld, lda %lld, incx %lld, incy %lld", *uplo,
             *n, *lda, *incx, *incy);
#else
    snprintf(buffer, 256, "csymv inputs: uplo %c, n %d, lda %d, incx %d, incy %d", *uplo, *n, *lda,
             *incx, *incy);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    scomplex q__1, q__2, q__3, q__4;
    /* Local variables */
    aocl_int64_t i__, j, ix, iy, jx, jy, kx, ky, info;
    scomplex temp1, temp2;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    /* -- LAPACK auxiliary routine (version 3.4.2) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* September 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. */
    /* .. Local Scalars .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;
    /* Function Body */
    info = 0;
    if(!lsame_(uplo, "U", 1, 1) && !lsame_(uplo, "L", 1, 1))
    {
        info = 1;
    }
    else if(*n < 0)
    {
        info = 2;
    }
    else if(*lda < fla_max(1, *n))
    {
        info = 5;
    }
    else if(*incx == 0)
    {
        info = 7;
    }
    else if(*incy == 0)
    {
        info = 10;
    }
    if(info != 0)
    {
        aocl_blas_xerbla("CSYMV ", &info, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible. */
    if(*n == 0 || alpha->real == 0.f && alpha->imag == 0.f && (beta->real == 1.f && beta->imag == 0.f))
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Set up the start points in X and Y. */
    if(*incx > 0)
    {
        kx = 1;
    }
    else
    {
        kx = 1 - (*n - 1) * *incx;
    }
    if(*incy > 0)
    {
        ky = 1;
    }
    else
    {
        ky = 1 - (*n - 1) * *incy;
    }
    /* Start the operations. In this version the elements of A are */
    /* accessed sequentially with one pass through the triangular part */
    /* of A. */
    /* First form y := beta*y. */
    if(beta->real != 1.f || beta->imag != 0.f)
    {
        if(*incy == 1)
        {
            if(beta->real == 0.f && beta->imag == 0.f)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    y[i__2].real = 0.f;
                    y[i__2].imag = 0.f; // , expr subst
                    /* L10: */
                }
            }
            else
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = i__;
                    i__3 = i__;
                    q__1.real = beta->real * y[i__3].real - beta->imag * y[i__3].imag;
                    q__1.imag = beta->real * y[i__3].imag + beta->imag * y[i__3].real; // , expr subst
                    y[i__2].real = q__1.real;
                    y[i__2].imag = q__1.imag; // , expr subst
                    /* L20: */
                }
            }
        }
        else
        {
            iy = ky;
            if(beta->real == 0.f && beta->imag == 0.f)
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = iy;
                    y[i__2].real = 0.f;
                    y[i__2].imag = 0.f; // , expr subst
                    iy += *incy;
                    /* L30: */
                }
            }
            else
            {
                i__1 = *n;
                for(i__ = 1; i__ <= i__1; ++i__)
                {
                    i__2 = iy;
                    i__3 = iy;
                    q__1.real = beta->real * y[i__3].real - beta->imag * y[i__3].imag;
                    q__1.imag = beta->real * y[i__3].imag + beta->imag * y[i__3].real; // , expr subst
                    y[i__2].real = q__1.real;
                    y[i__2].imag = q__1.imag; // , expr subst
                    iy += *incy;
                    /* L40: */
                }
            }
        }
    }
    if(alpha->real == 0.f && alpha->imag == 0.f)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(lsame_(uplo, "U", 1, 1))
    {
        /* Form y when A is stored in upper triangle. */
        if(*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = j;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag;
                q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real; // , expr subst
                temp1.real = q__1.real;
                temp1.imag = q__1.imag; // , expr subst
                temp2.real = 0.f;
                temp2.imag = 0.f; // , expr subst
                i__2 = j - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    q__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag;
                    q__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5].real; // , expr subst
                    q__1.real = y[i__4].real + q__2.real;
                    q__1.imag = y[i__4].imag + q__2.imag; // , expr subst
                    y[i__3].real = q__1.real;
                    y[i__3].imag = q__1.imag; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__;
                    q__2.real = a[i__3].real * x[i__4].real - a[i__3].imag * x[i__4].imag;
                    q__2.imag = a[i__3].real * x[i__4].imag + a[i__3].imag * x[i__4].real; // , expr subst
                    q__1.real = temp2.real + q__2.real;
                    q__1.imag = temp2.imag + q__2.imag; // , expr subst
                    temp2.real = q__1.real;
                    temp2.imag = q__1.imag; // , expr subst
                    /* L50: */
                }
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                q__3.real = temp1.real * a[i__4].real - temp1.imag * a[i__4].imag;
                q__3.imag = temp1.real * a[i__4].imag + temp1.imag * a[i__4].real; // , expr subst
                q__2.real = y[i__3].real + q__3.real;
                q__2.imag = y[i__3].imag + q__3.imag; // , expr subst
                q__4.real = alpha->real * temp2.real - alpha->imag * temp2.imag;
                q__4.imag = alpha->real * temp2.imag + alpha->imag * temp2.real; // , expr subst
                q__1.real = q__2.real + q__4.real;
                q__1.imag = q__2.imag + q__4.imag; // , expr subst
                y[i__2].real = q__1.real;
                y[i__2].imag = q__1.imag; // , expr subst
                /* L60: */
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = jx;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag;
                q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real; // , expr subst
                temp1.real = q__1.real;
                temp1.imag = q__1.imag; // , expr subst
                temp2.real = 0.f;
                temp2.imag = 0.f; // , expr subst
                ix = kx;
                iy = ky;
                i__2 = j - 1;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    q__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag;
                    q__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5].real; // , expr subst
                    q__1.real = y[i__4].real + q__2.real;
                    q__1.imag = y[i__4].imag + q__2.imag; // , expr subst
                    y[i__3].real = q__1.real;
                    y[i__3].imag = q__1.imag; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = ix;
                    q__2.real = a[i__3].real * x[i__4].real - a[i__3].imag * x[i__4].imag;
                    q__2.imag = a[i__3].real * x[i__4].imag + a[i__3].imag * x[i__4].real; // , expr subst
                    q__1.real = temp2.real + q__2.real;
                    q__1.imag = temp2.imag + q__2.imag; // , expr subst
                    temp2.real = q__1.real;
                    temp2.imag = q__1.imag; // , expr subst
                    ix += *incx;
                    iy += *incy;
                    /* L70: */
                }
                i__2 = jy;
                i__3 = jy;
                i__4 = j + j * a_dim1;
                q__3.real = temp1.real * a[i__4].real - temp1.imag * a[i__4].imag;
                q__3.imag = temp1.real * a[i__4].imag + temp1.imag * a[i__4].real; // , expr subst
                q__2.real = y[i__3].real + q__3.real;
                q__2.imag = y[i__3].imag + q__3.imag; // , expr subst
                q__4.real = alpha->real * temp2.real - alpha->imag * temp2.imag;
                q__4.imag = alpha->real * temp2.imag + alpha->imag * temp2.real; // , expr subst
                q__1.real = q__2.real + q__4.real;
                q__1.imag = q__2.imag + q__4.imag; // , expr subst
                y[i__2].real = q__1.real;
                y[i__2].imag = q__1.imag; // , expr subst
                jx += *incx;
                jy += *incy;
                /* L80: */
            }
        }
    }
    else
    {
        /* Form y when A is stored in lower triangle. */
        if(*incx == 1 && *incy == 1)
        {
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = j;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag;
                q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real; // , expr subst
                temp1.real = q__1.real;
                temp1.imag = q__1.imag; // , expr subst
                temp2.real = 0.f;
                temp2.imag = 0.f; // , expr subst
                i__2 = j;
                i__3 = j;
                i__4 = j + j * a_dim1;
                q__2.real = temp1.real * a[i__4].real - temp1.imag * a[i__4].imag;
                q__2.imag = temp1.real * a[i__4].imag + temp1.imag * a[i__4].real; // , expr subst
                q__1.real = y[i__3].real + q__2.real;
                q__1.imag = y[i__3].imag + q__2.imag; // , expr subst
                y[i__2].real = q__1.real;
                y[i__2].imag = q__1.imag; // , expr subst
                i__2 = *n;
                for(i__ = j + 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__ + j * a_dim1;
                    q__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag;
                    q__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5].real; // , expr subst
                    q__1.real = y[i__4].real + q__2.real;
                    q__1.imag = y[i__4].imag + q__2.imag; // , expr subst
                    y[i__3].real = q__1.real;
                    y[i__3].imag = q__1.imag; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = i__;
                    q__2.real = a[i__3].real * x[i__4].real - a[i__3].imag * x[i__4].imag;
                    q__2.imag = a[i__3].real * x[i__4].imag + a[i__3].imag * x[i__4].real; // , expr subst
                    q__1.real = temp2.real + q__2.real;
                    q__1.imag = temp2.imag + q__2.imag; // , expr subst
                    temp2.real = q__1.real;
                    temp2.imag = q__1.imag; // , expr subst
                    /* L90: */
                }
                i__2 = j;
                i__3 = j;
                q__2.real = alpha->real * temp2.real - alpha->imag * temp2.imag;
                q__2.imag = alpha->real * temp2.imag + alpha->imag * temp2.real; // , expr subst
                q__1.real = y[i__3].real + q__2.real;
                q__1.imag = y[i__3].imag + q__2.imag; // , expr subst
                y[i__2].real = q__1.real;
                y[i__2].imag = q__1.imag; // , expr subst
                /* L100: */
            }
        }
        else
        {
            jx = kx;
            jy = ky;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = jx;
                q__1.real = alpha->real * x[i__2].real - alpha->imag * x[i__2].imag;
                q__1.imag = alpha->real * x[i__2].imag + alpha->imag * x[i__2].real; // , expr subst
                temp1.real = q__1.real;
                temp1.imag = q__1.imag; // , expr subst
                temp2.real = 0.f;
                temp2.imag = 0.f; // , expr subst
                i__2 = jy;
                i__3 = jy;
                i__4 = j + j * a_dim1;
                q__2.real = temp1.real * a[i__4].real - temp1.imag * a[i__4].imag;
                q__2.imag = temp1.real * a[i__4].imag + temp1.imag * a[i__4].real; // , expr subst
                q__1.real = y[i__3].real + q__2.real;
                q__1.imag = y[i__3].imag + q__2.imag; // , expr subst
                y[i__2].real = q__1.real;
                y[i__2].imag = q__1.imag; // , expr subst
                ix = jx;
                iy = jy;
                i__2 = *n;
                for(i__ = j + 1; i__ <= i__2; ++i__)
                {
                    ix += *incx;
                    iy += *incy;
                    i__3 = iy;
                    i__4 = iy;
                    i__5 = i__ + j * a_dim1;
                    q__2.real = temp1.real * a[i__5].real - temp1.imag * a[i__5].imag;
                    q__2.imag = temp1.real * a[i__5].imag + temp1.imag * a[i__5].real; // , expr subst
                    q__1.real = y[i__4].real + q__2.real;
                    q__1.imag = y[i__4].imag + q__2.imag; // , expr subst
                    y[i__3].real = q__1.real;
                    y[i__3].imag = q__1.imag; // , expr subst
                    i__3 = i__ + j * a_dim1;
                    i__4 = ix;
                    q__2.real = a[i__3].real * x[i__4].real - a[i__3].imag * x[i__4].imag;
                    q__2.imag = a[i__3].real * x[i__4].imag + a[i__3].imag * x[i__4].real; // , expr subst
                    q__1.real = temp2.real + q__2.real;
                    q__1.imag = temp2.imag + q__2.imag; // , expr subst
                    temp2.real = q__1.real;
                    temp2.imag = q__1.imag; // , expr subst
                    /* L110: */
                }
                i__2 = jy;
                i__3 = jy;
                q__2.real = alpha->real * temp2.real - alpha->imag * temp2.imag;
                q__2.imag = alpha->real * temp2.imag + alpha->imag * temp2.real; // , expr subst
                q__1.real = y[i__3].real + q__2.real;
                q__1.imag = y[i__3].imag + q__2.imag; // , expr subst
                y[i__2].real = q__1.real;
                y[i__2].imag = q__1.imag; // , expr subst
                jx += *incx;
                jy += *incy;
                /* L120: */
            }
        }
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CSYMV */
}
/* csymv_ */
