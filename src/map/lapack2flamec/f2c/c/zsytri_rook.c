/* ../netlib/zsytri_rook.f -- translated by f2c (version 20160102). You must link the resulting
 object file with libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix
 systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with
 -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c -lm Source for
 libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {1., 0.};
static dcomplex c_b2 = {0., 0.};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZSYTRI_ROOK */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZSYTRI_ROOK + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri_
 * rook.f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri_
 * rook.f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri_
 * rook.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, LDA, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 A( LDA, * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYTRI_ROOK computes the inverse of a scomplex symmetric */
/* > matrix A using the factorization A = U*D*U**T or A = L*D*L**T */
/* > computed by ZSYTRF_ROOK. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > Specifies whether the details of the factorization are stored */
/* > as an upper or lower triangular matrix. */
/* > = 'U': Upper triangular, form is A = U*D*U**T;
 */
/* > = 'L': Lower triangular, form is A = L*D*L**T. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > On entry, the block diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by ZSYTRF_ROOK. */
/* > */
/* > On exit, if INFO = 0, the (symmetric) inverse of the original */
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
/* > as determined by ZSYTRF_ROOK. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (N) */
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
/* > \date December 2016 */
/* > \ingroup complex16SYcomputational */
/* > \par Contributors: */
/* ================== */
/* > */
/* > \verbatim */
/* > */
/* > December 2016, Igor Kozachenko, */
/* > Computer Science Division, */
/* > University of California, Berkeley */
/* > */
/* > September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas, */
/* > School of Mathematics, */
/* > University of Manchester */
/* > */
/* > \endverbatim */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zsytri_rook_(char *uplo, aocl_int_t *n, dcomplex *a, aocl_int_t *lda, aocl_int_t *ipiv,
                  dcomplex *work, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zsytri_rook(uplo, n, a, lda, ipiv, work, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zsytri_rook(uplo, &n_64, a, &lda_64, ipiv, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zsytri_rook(char *uplo, aocl_int64_t *n, dcomplex *a, aocl_int64_t *lda,
                             aocl_int_t *ipiv, dcomplex *work, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zsytri_rook inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n,
                      *lda);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, i__1, i__2, i__3;
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    dcomplex d__;
    aocl_int64_t k;
    dcomplex t, ak;
    aocl_int64_t kp;
    dcomplex akp1, temp, akkp1;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t kstep;
    logical upper;
    /* -- LAPACK computational routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
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
        aocl_blas_xerbla("ZSYTRI_ROOK", &i__1, (ftnlen)11);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Check that the diagonal matrix D is nonsingular. */
    if(upper)
    {
        /* Upper triangular storage: examine D from bottom to top */
        for(*info = *n; *info >= 1; --(*info))
        {
            i__1 = *info + *info * a_dim1;
            if(ipiv[*info] > 0 && (a[i__1].real == 0. && a[i__1].imag == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
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
            if(ipiv[*info] > 0 && (a[i__2].real == 0. && a[i__2].imag == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            /* L20: */
        }
    }
    *info = 0;
    if(upper)
    {
        /* Compute inv(A) from the factorization A = U*D*U**T. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = 1;
    L30: /* If K > N, exit from loop. */
        if(k > *n)
        {
            goto L40;
        }
        if(ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            i__1 = k + k * a_dim1;
            z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            /* Compute column K of the inverse. */
            if(k > 1)
            {
                i__1 = k - 1;
                aocl_blas_zcopy(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zsymv(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1, &c_b2,
                                  &a[k * a_dim1 + 1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            i__1 = k + (k + 1) * a_dim1;
            t.real = a[i__1].real;
            t.imag = a[i__1].imag; // , expr subst
            z_div(&z__1, &a[k + k * a_dim1], &t);
            ak.real = z__1.real;
            ak.imag = z__1.imag; // , expr subst
            z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &t);
            akp1.real = z__1.real;
            akp1.imag = z__1.imag; // , expr subst
            z_div(&z__1, &a[k + (k + 1) * a_dim1], &t);
            akkp1.real = z__1.real;
            akkp1.imag = z__1.imag; // , expr subst
            z__3.real = ak.real * akp1.real - ak.imag * akp1.imag;
            z__3.imag = ak.real * akp1.imag + ak.imag * akp1.real; // , expr subst
            z__2.real = z__3.real - 1.;
            z__2.imag = z__3.imag - 0.; // , expr subst
            z__1.real = t.real * z__2.real - t.imag * z__2.imag;
            z__1.imag = t.real * z__2.imag + t.imag * z__2.real; // , expr subst
            d__.real = z__1.real;
            d__.imag = z__1.imag; // , expr subst
            i__1 = k + k * a_dim1;
            z_div(&z__1, &akp1, &d__);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            i__1 = k + 1 + (k + 1) * a_dim1;
            z_div(&z__1, &ak, &d__);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            i__1 = k + (k + 1) * a_dim1;
            z__2.real = -akkp1.real;
            z__2.imag = -akkp1.imag; // , expr subst
            z_div(&z__1, &z__2, &d__);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            /* Compute columns K and K+1 of the inverse. */
            if(k > 1)
            {
                i__1 = k - 1;
                aocl_blas_zcopy(&i__1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zsymv(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1, &c_b2,
                                  &a[k * a_dim1 + 1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
                i__1 = k + (k + 1) * a_dim1;
                i__2 = k + (k + 1) * a_dim1;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1],
                           &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
                i__1 = k - 1;
                aocl_blas_zcopy(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zsymv(uplo, &i__1, &z__1, &a[a_offset], lda, &work[1], &c__1, &c_b2,
                                  &a[(k + 1) * a_dim1 + 1], &c__1);
                i__1 = k + 1 + (k + 1) * a_dim1;
                i__2 = k + 1 + (k + 1) * a_dim1;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 2;
        }
        if(kstep == 1)
        {
            /* Interchange rows and columns K and IPIV(K) in the leading */
            /* submatrix A(1:k+1,1:k+1) */
            kp = ipiv[k];
            if(kp != k)
            {
                if(kp > 1)
                {
                    i__1 = kp - 1;
                    aocl_blas_zswap(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                i__1 = k - kp - 1;
                aocl_blas_zswap(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1],
                                lda);
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
            }
        }
        else
        {
            /* Interchange rows and columns K and K+1 with -IPIV(K) and */
            /* -IPIV(K+1)in the leading submatrix A(1:k+1,1:k+1) */
            kp = -ipiv[k];
            if(kp != k)
            {
                if(kp > 1)
                {
                    i__1 = kp - 1;
                    aocl_blas_zswap(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                i__1 = k - kp - 1;
                aocl_blas_zswap(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1],
                                lda);
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
            ++k;
            kp = -ipiv[k];
            if(kp != k)
            {
                if(kp > 1)
                {
                    i__1 = kp - 1;
                    aocl_blas_zswap(&i__1, &a[k * a_dim1 + 1], &c__1, &a[kp * a_dim1 + 1], &c__1);
                }
                i__1 = k - kp - 1;
                aocl_blas_zswap(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + (kp + 1) * a_dim1],
                                lda);
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
            }
        }
        ++k;
        goto L30;
    L40:;
    }
    else
    {
        /* Compute inv(A) from the factorization A = L*D*L**T. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        k = *n;
    L50: /* If K < 1, exit from loop. */
        if(k < 1)
        {
            goto L60;
        }
        if(ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            i__1 = k + k * a_dim1;
            z_div(&z__1, &c_b1, &a[k + k * a_dim1]);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            /* Compute column K of the inverse. */
            if(k < *n)
            {
                i__1 = *n - k;
                aocl_blas_zcopy(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zsymv(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1],
                                  &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            i__1 = k + (k - 1) * a_dim1;
            t.real = a[i__1].real;
            t.imag = a[i__1].imag; // , expr subst
            z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &t);
            ak.real = z__1.real;
            ak.imag = z__1.imag; // , expr subst
            z_div(&z__1, &a[k + k * a_dim1], &t);
            akp1.real = z__1.real;
            akp1.imag = z__1.imag; // , expr subst
            z_div(&z__1, &a[k + (k - 1) * a_dim1], &t);
            akkp1.real = z__1.real;
            akkp1.imag = z__1.imag; // , expr subst
            z__3.real = ak.real * akp1.real - ak.imag * akp1.imag;
            z__3.imag = ak.real * akp1.imag + ak.imag * akp1.real; // , expr subst
            z__2.real = z__3.real - 1.;
            z__2.imag = z__3.imag - 0.; // , expr subst
            z__1.real = t.real * z__2.real - t.imag * z__2.imag;
            z__1.imag = t.real * z__2.imag + t.imag * z__2.real; // , expr subst
            d__.real = z__1.real;
            d__.imag = z__1.imag; // , expr subst
            i__1 = k - 1 + (k - 1) * a_dim1;
            z_div(&z__1, &akp1, &d__);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            i__1 = k + k * a_dim1;
            z_div(&z__1, &ak, &d__);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            i__1 = k + (k - 1) * a_dim1;
            z__2.real = -akkp1.real;
            z__2.imag = -akkp1.imag; // , expr subst
            z_div(&z__1, &z__2, &d__);
            a[i__1].real = z__1.real;
            a[i__1].imag = z__1.imag; // , expr subst
            /* Compute columns K-1 and K of the inverse. */
            if(k < *n)
            {
                i__1 = *n - k;
                aocl_blas_zcopy(&i__1, &a[k + 1 + k * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zsymv(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1],
                                  &c__1, &c_b2, &a[k + 1 + k * a_dim1], &c__1);
                i__1 = k + k * a_dim1;
                i__2 = k + k * a_dim1;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + k * a_dim1], &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
                i__1 = k + (k - 1) * a_dim1;
                i__2 = k + (k - 1) * a_dim1;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &a[k + 1 + k * a_dim1], &c__1,
                           &a[k + 1 + (k - 1) * a_dim1], &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
                i__1 = *n - k;
                aocl_blas_zcopy(&i__1, &a[k + 1 + (k - 1) * a_dim1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zsymv(uplo, &i__1, &z__1, &a[k + 1 + (k + 1) * a_dim1], lda, &work[1],
                                  &c__1, &c_b2, &a[k + 1 + (k - 1) * a_dim1], &c__1);
                i__1 = k - 1 + (k - 1) * a_dim1;
                i__2 = k - 1 + (k - 1) * a_dim1;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &a[k + 1 + (k - 1) * a_dim1], &c__1);
                z__1.real = a[i__2].real - z__2.real;
                z__1.imag = a[i__2].imag - z__2.imag; // , expr subst
                a[i__1].real = z__1.real;
                a[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 2;
        }
        if(kstep == 1)
        {
            /* Interchange rows and columns K and IPIV(K) in the trailing */
            /* submatrix A(k-1:n,k-1:n) */
            kp = ipiv[k];
            if(kp != k)
            {
                if(kp < *n)
                {
                    i__1 = *n - kp;
                    aocl_blas_zswap(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1],
                                    &c__1);
                }
                i__1 = kp - k - 1;
                aocl_blas_zswap(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1],
                                lda);
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
            }
        }
        else
        {
            /* Interchange rows and columns K and K-1 with -IPIV(K) and */
            /* -IPIV(K-1) in the trailing submatrix A(k-1:n,k-1:n) */
            kp = -ipiv[k];
            if(kp != k)
            {
                if(kp < *n)
                {
                    i__1 = *n - kp;
                    aocl_blas_zswap(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1],
                                    &c__1);
                }
                i__1 = kp - k - 1;
                aocl_blas_zswap(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1],
                                lda);
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
            --k;
            kp = -ipiv[k];
            if(kp != k)
            {
                if(kp < *n)
                {
                    i__1 = *n - kp;
                    aocl_blas_zswap(&i__1, &a[kp + 1 + k * a_dim1], &c__1, &a[kp + 1 + kp * a_dim1],
                                    &c__1);
                }
                i__1 = kp - k - 1;
                aocl_blas_zswap(&i__1, &a[k + 1 + k * a_dim1], &c__1, &a[kp + (k + 1) * a_dim1],
                                lda);
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
            }
        }
        --k;
        goto L50;
    L60:;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZSYTRI_ROOK */
}
/* zsytri_rook__ */
