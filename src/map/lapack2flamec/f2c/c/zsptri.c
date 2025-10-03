/* ../netlib/zsptri.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static dcomplex c_b1 = {1., 0.};
static dcomplex c_b2 = {0., 0.};
static aocl_int64_t c__1 = 1;
/* > \brief \b ZSPTRI */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZSPTRI + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsptri.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsptri.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsptri.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZSPTRI( UPLO, N, AP, IPIV, WORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX*16 AP( * ), WORK( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSPTRI computes the inverse of a scomplex symmetric indefinite matrix */
/* > A in packed storage using the factorization A = U*D*U**T or */
/* > A = L*D*L**T computed by ZSPTRF. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* > AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* > On entry, the block diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L as computed by ZSPTRF, */
/* > stored as a packed triangular matrix. */
/* > */
/* > On exit, if INFO = 0, the (symmetric) inverse of the original */
/* > matrix, stored as a packed triangular matrix. The j-th column */
/* > of inv(A) is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', */
/* > AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D */
/* > as determined by ZSPTRF. */
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
/* > \date November 2011 */
/* > \ingroup complex16OTHERcomputational */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void zsptri_(char *uplo, aocl_int_t *n, dcomplex *ap, aocl_int_t *ipiv, dcomplex *work,
             aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zsptri(uplo, n, ap, ipiv, work, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zsptri(uplo, &n_64, ap, ipiv, work, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_zsptri(char *uplo, aocl_int64_t *n, dcomplex *ap, aocl_int_t *ipiv,
                        dcomplex *work, aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zsptri inputs: uplo %c, n %" FLA_IS "", *uplo, *n);

    /* System generated locals */
    aocl_int64_t i__1, i__2, i__3;
    dcomplex z__1, z__2, z__3;
    /* Builtin functions */
    void z_div(dcomplex *, dcomplex *, dcomplex *);
    /* Local variables */
    dcomplex d__;
    aocl_int64_t j, k;
    dcomplex t, ak;
    aocl_int64_t kc, kp, kx, kpc, npp;
    dcomplex akp1, temp, akkp1;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t kstep;
    logical upper;
    aocl_int64_t kcnext;
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
    --work;
    --ipiv;
    --ap;
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
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("ZSPTRI", &i__1, (ftnlen)6);
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
        kp = *n * (*n + 1) / 2;
        for(*info = *n; *info >= 1; --(*info))
        {
            i__1 = kp;
            if(ipiv[*info] > 0 && (ap[i__1].real == 0. && ap[i__1].imag == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            kp -= *info;
            /* L10: */
        }
    }
    else
    {
        /* Lower triangular storage: examine D from top to bottom. */
        kp = 1;
        i__1 = *n;
        for(*info = 1; *info <= i__1; ++(*info))
        {
            i__2 = kp;
            if(ipiv[*info] > 0 && (ap[i__2].real == 0. && ap[i__2].imag == 0.))
            {
                AOCL_DTL_TRACE_LOG_EXIT
                return;
            }
            kp = kp + *n - *info + 1;
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
        kc = 1;
    L30: /* If K > N, exit from loop. */
        if(k > *n)
        {
            goto L50;
        }
        kcnext = kc + k;
        if(ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            i__1 = kc + k - 1;
            z_div(&z__1, &c_b1, &ap[kc + k - 1]);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            /* Compute column K of the inverse. */
            if(k > 1)
            {
                i__1 = k - 1;
                aocl_blas_zcopy(&i__1, &ap[kc], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zspmv(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &ap[kc],
                                  &c__1);
                i__1 = kc + k - 1;
                i__2 = kc + k - 1;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &ap[kc], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            i__1 = kcnext + k - 1;
            t.real = ap[i__1].real;
            t.imag = ap[i__1].imag; // , expr subst
            z_div(&z__1, &ap[kc + k - 1], &t);
            ak.real = z__1.real;
            ak.imag = z__1.imag; // , expr subst
            z_div(&z__1, &ap[kcnext + k], &t);
            akp1.real = z__1.real;
            akp1.imag = z__1.imag; // , expr subst
            z_div(&z__1, &ap[kcnext + k - 1], &t);
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
            i__1 = kc + k - 1;
            z_div(&z__1, &akp1, &d__);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            i__1 = kcnext + k;
            z_div(&z__1, &ak, &d__);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            i__1 = kcnext + k - 1;
            z__2.real = -akkp1.real;
            z__2.imag = -akkp1.imag; // , expr subst
            z_div(&z__1, &z__2, &d__);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            /* Compute columns K and K+1 of the inverse. */
            if(k > 1)
            {
                i__1 = k - 1;
                aocl_blas_zcopy(&i__1, &ap[kc], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zspmv(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &ap[kc],
                                  &c__1);
                i__1 = kc + k - 1;
                i__2 = kc + k - 1;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &ap[kc], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
                i__1 = kcnext + k - 1;
                i__2 = kcnext + k - 1;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &ap[kc], &c__1, &ap[kcnext], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
                i__1 = k - 1;
                aocl_blas_zcopy(&i__1, &ap[kcnext], &c__1, &work[1], &c__1);
                i__1 = k - 1;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zspmv(uplo, &i__1, &z__1, &ap[1], &work[1], &c__1, &c_b2, &ap[kcnext],
                                  &c__1);
                i__1 = kcnext + k;
                i__2 = kcnext + k;
                i__3 = k - 1;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &ap[kcnext], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 2;
            kcnext = kcnext + k + 1;
        }
        kp = (i__1 = ipiv[k], f2c_dabs(i__1));
        if(kp != k)
        {
            /* Interchange rows and columns K and KP in the leading */
            /* submatrix A(1:k+1,1:k+1) */
            kpc = (kp - 1) * kp / 2 + 1;
            i__1 = kp - 1;
            aocl_blas_zswap(&i__1, &ap[kc], &c__1, &ap[kpc], &c__1);
            kx = kpc + kp - 1;
            i__1 = k - 1;
            for(j = kp + 1; j <= i__1; ++j)
            {
                kx = kx + j - 1;
                i__2 = kc + j - 1;
                temp.real = ap[i__2].real;
                temp.imag = ap[i__2].imag; // , expr subst
                i__2 = kc + j - 1;
                i__3 = kx;
                ap[i__2].real = ap[i__3].real;
                ap[i__2].imag = ap[i__3].imag; // , expr subst
                i__2 = kx;
                ap[i__2].real = temp.real;
                ap[i__2].imag = temp.imag; // , expr subst
                /* L40: */
            }
            i__1 = kc + k - 1;
            temp.real = ap[i__1].real;
            temp.imag = ap[i__1].imag; // , expr subst
            i__1 = kc + k - 1;
            i__2 = kpc + kp - 1;
            ap[i__1].real = ap[i__2].real;
            ap[i__1].imag = ap[i__2].imag; // , expr subst
            i__1 = kpc + kp - 1;
            ap[i__1].real = temp.real;
            ap[i__1].imag = temp.imag; // , expr subst
            if(kstep == 2)
            {
                i__1 = kc + k + k - 1;
                temp.real = ap[i__1].real;
                temp.imag = ap[i__1].imag; // , expr subst
                i__1 = kc + k + k - 1;
                i__2 = kc + k + kp - 1;
                ap[i__1].real = ap[i__2].real;
                ap[i__1].imag = ap[i__2].imag; // , expr subst
                i__1 = kc + k + kp - 1;
                ap[i__1].real = temp.real;
                ap[i__1].imag = temp.imag; // , expr subst
            }
        }
        k += kstep;
        kc = kcnext;
        goto L30;
    L50:;
    }
    else
    {
        /* Compute inv(A) from the factorization A = L*D*L**T. */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2, depending on the size of the diagonal blocks. */
        npp = *n * (*n + 1) / 2;
        k = *n;
        kc = npp;
    L60: /* If K < 1, exit from loop. */
        if(k < 1)
        {
            goto L80;
        }
        kcnext = kc - (*n - k + 2);
        if(ipiv[k] > 0)
        {
            /* 1 x 1 diagonal block */
            /* Invert the diagonal block. */
            i__1 = kc;
            z_div(&z__1, &c_b1, &ap[kc]);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            /* Compute column K of the inverse. */
            if(k < *n)
            {
                i__1 = *n - k;
                aocl_blas_zcopy(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zspmv(uplo, &i__1, &z__1, &ap[kc + *n - k + 1], &work[1], &c__1, &c_b2,
                                  &ap[kc + 1], &c__1);
                i__1 = kc;
                i__2 = kc;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &ap[kc + 1], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 1;
        }
        else
        {
            /* 2 x 2 diagonal block */
            /* Invert the diagonal block. */
            i__1 = kcnext + 1;
            t.real = ap[i__1].real;
            t.imag = ap[i__1].imag; // , expr subst
            z_div(&z__1, &ap[kcnext], &t);
            ak.real = z__1.real;
            ak.imag = z__1.imag; // , expr subst
            z_div(&z__1, &ap[kc], &t);
            akp1.real = z__1.real;
            akp1.imag = z__1.imag; // , expr subst
            z_div(&z__1, &ap[kcnext + 1], &t);
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
            i__1 = kcnext;
            z_div(&z__1, &akp1, &d__);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            i__1 = kc;
            z_div(&z__1, &ak, &d__);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            i__1 = kcnext + 1;
            z__2.real = -akkp1.real;
            z__2.imag = -akkp1.imag; // , expr subst
            z_div(&z__1, &z__2, &d__);
            ap[i__1].real = z__1.real;
            ap[i__1].imag = z__1.imag; // , expr subst
            /* Compute columns K-1 and K of the inverse. */
            if(k < *n)
            {
                i__1 = *n - k;
                aocl_blas_zcopy(&i__1, &ap[kc + 1], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zspmv(uplo, &i__1, &z__1, &ap[kc + (*n - k + 1)], &work[1], &c__1,
                                  &c_b2, &ap[kc + 1], &c__1);
                i__1 = kc;
                i__2 = kc;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &ap[kc + 1], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
                i__1 = kcnext + 1;
                i__2 = kcnext + 1;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &ap[kc + 1], &c__1, &ap[kcnext + 2], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
                i__1 = *n - k;
                aocl_blas_zcopy(&i__1, &ap[kcnext + 2], &c__1, &work[1], &c__1);
                i__1 = *n - k;
                z__1.real = -1.;
                z__1.imag = -0.; // , expr subst
                aocl_lapack_zspmv(uplo, &i__1, &z__1, &ap[kc + (*n - k + 1)], &work[1], &c__1,
                                  &c_b2, &ap[kcnext + 2], &c__1);
                i__1 = kcnext;
                i__2 = kcnext;
                i__3 = *n - k;
                aocl_lapack_zdotu_f2c(&z__2, &i__3, &work[1], &c__1, &ap[kcnext + 2], &c__1);
                z__1.real = ap[i__2].real - z__2.real;
                z__1.imag = ap[i__2].imag - z__2.imag; // , expr subst
                ap[i__1].real = z__1.real;
                ap[i__1].imag = z__1.imag; // , expr subst
            }
            kstep = 2;
            kcnext -= *n - k + 3;
        }
        kp = (i__1 = ipiv[k], f2c_dabs(i__1));
        if(kp != k)
        {
            /* Interchange rows and columns K and KP in the trailing */
            /* submatrix A(k-1:n,k-1:n) */
            kpc = npp - (*n - kp + 1) * (*n - kp + 2) / 2 + 1;
            if(kp < *n)
            {
                i__1 = *n - kp;
                aocl_blas_zswap(&i__1, &ap[kc + kp - k + 1], &c__1, &ap[kpc + 1], &c__1);
            }
            kx = kc + kp - k;
            i__1 = kp - 1;
            for(j = k + 1; j <= i__1; ++j)
            {
                kx = kx + *n - j + 1;
                i__2 = kc + j - k;
                temp.real = ap[i__2].real;
                temp.imag = ap[i__2].imag; // , expr subst
                i__2 = kc + j - k;
                i__3 = kx;
                ap[i__2].real = ap[i__3].real;
                ap[i__2].imag = ap[i__3].imag; // , expr subst
                i__2 = kx;
                ap[i__2].real = temp.real;
                ap[i__2].imag = temp.imag; // , expr subst
                /* L70: */
            }
            i__1 = kc;
            temp.real = ap[i__1].real;
            temp.imag = ap[i__1].imag; // , expr subst
            i__1 = kc;
            i__2 = kpc;
            ap[i__1].real = ap[i__2].real;
            ap[i__1].imag = ap[i__2].imag; // , expr subst
            i__1 = kpc;
            ap[i__1].real = temp.real;
            ap[i__1].imag = temp.imag; // , expr subst
            if(kstep == 2)
            {
                i__1 = kc - *n + k - 1;
                temp.real = ap[i__1].real;
                temp.imag = ap[i__1].imag; // , expr subst
                i__1 = kc - *n + k - 1;
                i__2 = kc - *n + kp - 1;
                ap[i__1].real = ap[i__2].real;
                ap[i__1].imag = ap[i__2].imag; // , expr subst
                i__1 = kc - *n + kp - 1;
                ap[i__1].real = temp.real;
                ap[i__1].imag = temp.imag; // , expr subst
            }
        }
        k -= kstep;
        kc = kcnext;
        goto L60;
    L80:;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZSPTRI */
}
/* zsptri_ */
