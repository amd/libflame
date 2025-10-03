/* ../netlib/chetri.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static scomplex c_b2 = {0.f, 0.f};
static aocl_int64_t c__1 = 1;
/* > \brief \b CHETRI */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHETRI + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chetri.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chetri.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chetri.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHETRI computes the inverse of a scomplex Hermitian indefinite matrix */
/* > A using the factorization A = U*D*U**H or A = L*D*L**H computed by */
/* > CHETRF. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**H;
 */
/* > = 'L': Lower triangular, form is A = L*D*L**H. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is COMPLEX array, dimension (LDA,N) */
/* > On entry, the block diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by CHETRF. */
/* > */
/* > On exit, if INFO = 0, the (Hermitian) inverse of the original */
/* > matrix. If UPLO = 'U', the upper triangular part of the */
/* > inverse is formed and the part of A below the diagonal is not */
/* > referenced;
if UPLO = 'L' the lower triangular part of the */
/* > inverse is formed and the part of A above the diagonal is */
/* > not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by CHETRF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) = 0;
the matrix is singular and its */
/* > inverse could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexHEcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chetri_(char *uplo, aocl_int_t *n, scomplex *a, aocl_int_t *lda, aocl_int_t *ipiv,
             scomplex *work, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chetri(uplo, n, a, lda, ipiv, work, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chetri(uplo, &n_64, a, &lda_64, ipiv, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_chetri(char *uplo, aocl_int64_t *n, scomplex *a, aocl_int64_t *lda,
                        aocl_int_t *ipiv, scomplex *work, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "chetri inputs: uplo %c, n %lld, lda %lld", *uplo, *n, *lda);
#else
    snprintf(buffer, 256, "chetri inputs: uplo %c, n %d, lda %d", *uplo, *n, *lda);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    scomplex q__1, q__2;
    /* Builtin functions */
    double c_abs(scomplex *);
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    real d__;
    aocl_int64_t j, k;
    real t, ak;
    aocl_int64_t kp;
    real akp1;
    scomplex temp, akkp1;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t kstep;
    logical upper;
    /* -- LAPACK computational routine (version 3.4.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* November 2011 */
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
    --ipiv;
    --work;
    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CHETRI", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if(upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for(*info = *n; *info >= 1; --(*info))
        {
            i__1 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__1].real == 0.f && a[i__1].imag == 0.f))
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
            /* L10: */
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        i__1 = *n;
        for(*info = 1; *info <= i__1; ++(*info))
        {
            i__2 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__2].real == 0.f && a[i__2].imag == 0.f))
            {
                AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
                return;
            }
            /* L20: */
        }
    }
    *info = 0;
    if(upper)
    {
        /* Compute inv(A) from the factorization A = U*D*U**H. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = 1;
    L30: /* If K > N, exit from loop. */
        if(k > *n)
        {
            goto L50;
        }
        if(ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            i__1 = k + k * a_dim1;
            i__2 = k + k * a_dim1;
            r__1 = 1.f / a[i__2].real;
            a[i__1].real = r__1;
            a[i__1].imag = 0.f; // , expr subst
            /* Compute column K of the inverse. */
            if(k > 1)
            {
                i__1 = k - 1;
                aocl_blas_ccopy(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_chemv(uplo, &i__1, &q__1, &a[a_offset], lda, &work[1], &c__1, &c_b2,
                                &a[k * a_dim1 + 1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = k - 1;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
                r__1 = q__2.real;
                q__1.real = a[i__2].real - r__1;
                q__1.imag = a[i__2].imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            t = c_abs(&a[k + (k + 1) * a_dim1]);
            i__1 = k + k * a_dim1;
            ak = a[i__1].real / t;
            i__1 = k + 1 + (k + 1) * a_dim1;
            akp1 = a[i__1].real / t;
            i__1 = k + (k + 1) * a_dim1;
            q__1.real = a[i__1].real / t;
            q__1.imag = a[i__1].imag / t; // , expr subst
            akkp1.real = q__1.real;
            akkp1.imag = q__1.imag; // , expr subst
            d__ = t * (ak * akp1 - 1.f);
            i__1 = k + k * a_dim1;
            r__1 = akp1 / d__;
            a[i__1].real = r__1;
            a[i__1].imag = 0.f; // , expr subst
            i__1 = k + 1 + (k + 1) * a_dim1;
            r__1 = ak / d__;
            a[i__1].real = r__1;
            a[i__1].imag = 0.f; // , expr subst
            i__1 = k + (k + 1) * a_dim1;
            q__2.real = -akkp1.real;
            q__2.imag = -akkp1.imag; // , expr subst
            q__1.real = q__2.real / d__;
            q__1.imag = q__2.imag / d__; // , expr subst
            a[i__1].real = q__1.real;
            a[i__1].imag = q__1.imag; // , expr subst
            /* Compute columns K and K+1 of the inverse. */
            if(k > 1)
            {
                i__1 = k - 1;
                aocl_blas_ccopy(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_chemv(uplo, &i__1, &q__1, &a[a_offset], lda, &work[1], &c__1, &c_b2,
                                &a[k * a_dim1 + 1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = k - 1;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
                r__1 = q__2.real;
                q__1.real = a[i__2].real - r__1;
                q__1.imag = a[i__2].imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + (k + 1) * a_dim1;
                i__2 = k + (k + 1) * a_dim1;
                i__3 = k - 1;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1],
                           &c__1);
                q__1.real = a[i__2].real - q__2.real;
                q__1.imag = a[i__2].imag - q__2.imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
                i__1 = k - 1;
                aocl_blas_ccopy(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_chemv(uplo, &i__1, &q__1, &a[a_offset], lda, &work[1], &c__1, &c_b2,
                                &a[(k + 1) * a_dim1 + 1], &c__1);
                i__1 = k + 1 + (k + 1) * a_dim1;
                i__2 = k + 1 + (k + 1) * a_dim1;
                i__3 = k - 1;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
                r__1 = q__2.real;
                q__1.real = a[i__2].real - r__1;
                q__1.imag = a[i__2].imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
            }
            kstep = 2;
        }
        kp = (i__1 = ipiv[k], f2c_abs(i__1));
        if(kp != k)
        {
            /* Interchange rows and columns K and KP in the leading */
            /* submatrix A(1:k+1,1:k+1) */
            i__1 = kp - 1;
            aocl_blas_cswap(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
            i__1 = k - 1;
            for(j = kp + 1; j <= i__1; ++j)
            {
                r_cnjg(&q__1, &a[j + k * a_dim1]);
                temp.real = q__1.real;
                temp.imag = q__1.imag; // , expr subst
                i__2 = j + k * a_dim1;
                r_cnjg(&q__1, &a[kp + j * a_dim1]);
                a[i__2].real = q__1.real;
                a[i__2].imag = q__1.imag; // , expr subst
                i__2 = kp + j * a_dim1;
                a[i__2].real = temp.real;
                a[i__2].imag = temp.imag; // , expr subst
                /* L40: */
            }
            i__1 = kp + k * a_dim1;
            r_cnjg(&q__1, &a[kp + k * a_dim1]);
            a[i__1].real = q__1.real;
            a[i__1].imag = q__1.imag; // , expr subst
            i__1 = k + k * a_dim1;
            temp.real = a[i__1].real;
            temp.imag = a[i__1].imag; // , expr subst
            i__1 = k + k * a_dim1;
            i__2 = kp + kp * a_dim1;
            a[i__1].real = a[i__2].real;
            a[i__1].imag = a[i__2].imag; // , expr subst
            i__1 = kp + kp * a_dim1;
            a[i__1].real = temp.real;
            a[i__1].imag = temp.imag; // , expr subst
            if(kstep == 2)
            {
                i__1 = k + (k + 1) * a_dim1;
                temp.real = a[i__1].real;
                temp.imag = a[i__1].imag; // , expr subst
                i__1 = k + (k + 1) * a_dim1;
                i__2 = kp + (k + 1) * a_dim1;
                a[i__1].real = a[i__2].real;
                a[i__1].imag = a[i__2].imag; // , expr subst
                i__1 = kp + (k + 1) * a_dim1;
                a[i__1].real = temp.real;
                a[i__1].imag = temp.imag; // , expr subst
            }
        }
        k += kstep;
        goto L30;
    L50:;
    }
    else
    {
        /* Compute inv(A) from the factorization A = L*D*L**H. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = *n;
    L60: /* If K < 1, exit from loop. */
        if(k < 1)
        {
            goto L80;
        }
        if(ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            i__1 = k + k * a_dim1;
            i__2 = k + k * a_dim1;
            r__1 = 1.f / a[i__2].real;
            a[i__1].real = r__1;
            a[i__1].imag = 0.f; // , expr subst
            /* Compute column K of the inverse. */
            if(k < *n)
            {
                i__1 = *n - k;
                aocl_blas_ccopy(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_chemv(uplo, &i__1, &q__1, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1],
                                &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = *n - k;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
                r__1 = q__2.real;
                q__1.real = a[i__2].real - r__1;
                q__1.imag = a[i__2].imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            t = c_abs(&a[k + (k - 1) * a_dim1]);
            i__1 = k - 1 + (k - 1) * a_dim1;
            ak = a[i__1].real / t;
            i__1 = k + k * a_dim1;
            akp1 = a[i__1].real / t;
            i__1 = k + (k - 1) * a_dim1;
            q__1.real = a[i__1].real / t;
            q__1.imag = a[i__1].imag / t; // , expr subst
            akkp1.real = q__1.real;
            akkp1.imag = q__1.imag; // , expr subst
            d__ = t * (ak * akp1 - 1.f);
            i__1 = k - 1 + (k - 1) * a_dim1;
            r__1 = akp1 / d__;
            a[i__1].real = r__1;
            a[i__1].imag = 0.f; // , expr subst
            i__1 = k + k * a_dim1;
            r__1 = ak / d__;
            a[i__1].real = r__1;
            a[i__1].imag = 0.f; // , expr subst
            i__1 = k + (k - 1) * a_dim1;
            q__2.real = -akkp1.real;
            q__2.imag = -akkp1.imag; // , expr subst
            q__1.real = q__2.real / d__;
            q__1.imag = q__2.imag / d__; // , expr subst
            a[i__1].real = q__1.real;
            a[i__1].imag = q__1.imag; // , expr subst
            /* Compute columns K-1 and K of the inverse. */
            if(k < *n)
            {
                i__1 = *n - k;
                aocl_blas_ccopy(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_chemv(uplo, &i__1, &q__1, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1],
                                &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = *n - k;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
                r__1 = q__2.real;
                q__1.real = a[i__2].real - r__1;
                q__1.imag = a[i__2].imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
                i__1 = k + (k - 1) * a_dim1;
                i__2 = k + (k - 1) * a_dim1;
                i__3 = *n - k;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &a[k + 1 + k * a_dim1], &c__1,
                           &a[k + 1 + (k - 1) * a_dim1], &c__1);
                q__1.real = a[i__2].real - q__2.real;
                q__1.imag = a[i__2].imag - q__2.imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
                i__1 = *n - k;
                aocl_blas_ccopy(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                q__1.real = -1.f;
                q__1.imag = -0.f; // , expr subst
                aocl_blas_chemv(uplo, &i__1, &q__1, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1],
                                &c__1, &c_b2, &a[k + 1 + (k - 1) * a_dim1], &c__1);
                i__1 = k - 1 + (k - 1) * a_dim1;
                i__2 = k - 1 + (k - 1) * a_dim1;
                i__3 = *n - k;
                aocl_lapack_cdotc_f2c(&q__2, &i__3, &work[1], &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
                r__1 = q__2.real;
                q__1.real = a[i__2].real - r__1;
                q__1.imag = a[i__2].imag; // , expr subst
                a[i__1].real = q__1.real;
                a[i__1].imag = q__1.imag; // , expr subst
            }
            kstep = 2;
        }
        kp = (i__1 = ipiv[k], f2c_abs(i__1));
        if(kp != k)
        {
            /* Interchange rows and columns K and KP in the trailing */
            /* submatrix A(k-1:n,k-1:n) */
            if(kp < *n)
            {
                i__1 = *n - kp;
                aocl_blas_cswap(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1],
                                &c__1);
            }
            i__1 = kp - 1;
            for(j = k + 1; j <= i__1; ++j)
            {
                r_cnjg(&q__1, &a[j + k * a_dim1]);
                temp.real = q__1.real;
                temp.imag = q__1.imag; // , expr subst
                i__2 = j + k * a_dim1;
                r_cnjg(&q__1, &a[kp + j * a_dim1]);
                a[i__2].real = q__1.real;
                a[i__2].imag = q__1.imag; // , expr subst
                i__2 = kp + j * a_dim1;
                a[i__2].real = temp.real;
                a[i__2].imag = temp.imag; // , expr subst
                /* L70: */
            }
            i__1 = kp + k * a_dim1;
            r_cnjg(&q__1, &a[kp + k * a_dim1]);
            a[i__1].real = q__1.real;
            a[i__1].imag = q__1.imag; // , expr subst
            i__1 = k + k * a_dim1;
            temp.real = a[i__1].real;
            temp.imag = a[i__1].imag; // , expr subst
            i__1 = k + k * a_dim1;
            i__2 = kp + kp * a_dim1;
            a[i__1].real = a[i__2].real;
            a[i__1].imag = a[i__2].imag; // , expr subst
            i__1 = kp + kp * a_dim1;
            a[i__1].real = temp.real;
            a[i__1].imag = temp.imag; // , expr subst
            if(kstep == 2)
            {
                i__1 = k + (k - 1) * a_dim1;
                temp.real = a[i__1].real;
                temp.imag = a[i__1].imag; // , expr subst
                i__1 = k + (k - 1) * a_dim1;
                i__2 = kp + (k - 1) * a_dim1;
                a[i__1].real = a[i__2].real;
                a[i__1].imag = a[i__2].imag; // , expr subst
                i__1 = kp + (k - 1) * a_dim1;
                a[i__1].real = temp.real;
                a[i__1].imag = temp.imag; // , expr subst
            }
        }
        k -= kstep;
        goto L60;
    L80:;
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHETRI */
}
/* chetri_ */
