/* ../netlib/chptrf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
 *     Modifications Copyright (c) 2024 Advanced Micro Devices, Inc.  All rights reserved.
 */

#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief \b CHPTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CHPTRF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chptrf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chptrf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chptrf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CHPTRF( UPLO, N, AP, IPIV, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, N */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* COMPLEX AP( * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPTRF computes the factorization of a scomplex Hermitian packed */
/* > matrix A using the Bunch-Kaufman diagonal pivoting method: */
/* > */
/* > A = U*D*U**H or A = L*D*L**H */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is Hermitian and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': Upper triangle of A is stored;
 */
/* > = 'L': Lower triangle of A is stored. */
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
/* > AP is COMPLEX array, dimension (N*(N+1)/2) */
/* > On entry, the upper or lower triangle of the Hermitian matrix */
/* > A, packed columnwise in a linear array. The j-th column of A */
/* > is stored in the array AP as follows: */
/* > if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 */
/* > if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > */
/* > On exit, the block diagonal matrix D and the multipliers used */
/* > to obtain the factor U or L, stored as a packed triangular */
/* > matrix overwriting A (see below for further details). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > Details of the interchanges and the block structure of D. */
/* > If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* > interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* > If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and */
/* > columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* > is a 2-by-2 diagonal block. If UPLO = 'L' and IPIV(k) = */
/* > IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were */
/* > interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, D(i,i) is exactly zero. The factorization */
/* > has been completed, but the block diagonal matrix D is */
/* > exactly singular, and division by zero will occur if it */
/* > is used to solve a system of equations. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexOTHERcomputational */
/* > \par Further Details: */
/* ===================== */
/* > */
/* > \verbatim */
/* > */
/* > If UPLO = 'U', then A = U*D*U**H, where */
/* > U = P(n)*U(n)* ... *P(k)U(k)* ..., */
/* > i.e., U is a product of terms P(k)*U(k), where k decreases from n to */
/* > 1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
/* > and 2-by-2 diagonal blocks D(k). P(k) is a permutation matrix as */
/* > defined by IPIV(k), and U(k) is a unit upper triangular matrix, such */
/* > that if the diagonal block D(k) is of order s (s = 1 or 2), then */
/* > */
/* > ( I v 0 ) k-s */
/* > U(k) = ( 0 I 0 ) s */
/* > ( 0 0 I ) n-k */
/* > k-s s n-k */
/* > */
/* > If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k). */
/* > If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k), */
/* > and A(k,k), and v overwrites A(1:k-2,k-1:k). */
/* > */
/* > If UPLO = 'L', then A = L*D*L**H, where */
/* > L = P(1)*L(1)* ... *P(k)*L(k)* ..., */
/* > i.e., L is a product of terms P(k)*L(k), where k increases from 1 to */
/* > n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1 */
/* > and 2-by-2 diagonal blocks D(k). P(k) is a permutation matrix as */
/* > defined by IPIV(k), and L(k) is a unit lower triangular matrix, such */
/* > that if the diagonal block D(k) is of order s (s = 1 or 2), then */
/* > */
/* > ( I 0 0 ) k-1 */
/* > L(k) = ( 0 I 0 ) s */
/* > ( 0 v I ) n-k-s+1 */
/* > k-1 s n-k-s+1 */
/* > */
/* > If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k). */
/* > If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k), */
/* > and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1). */
/* > \endverbatim */
/* > \par Contributors: */
/* ================== */
/* > */
/* > J. Lewis, Boeing Computer Services Company */
/* > */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void chptrf_(char *uplo, aocl_int_t *n, scomplex *ap, aocl_int_t *ipiv, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_chptrf(uplo, n, ap, ipiv, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t info_64 = *info;

    aocl_lapack_chptrf(uplo, &n_64, ap, ipiv, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_chptrf(char *uplo, aocl_int64_t *n, scomplex *ap, aocl_int_t *ipiv,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "chptrf inputs: uplo %c, n %lld", *uplo, *n);
#else
    snprintf(buffer, 256, "chptrf inputs: uplo %c, n %d", *uplo, *n);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3, r__4;
    scomplex q__1, q__2, q__3, q__4, q__5, q__6;
    /* Builtin functions */
    double sqrt(doublereal), r_imag(scomplex *);
    void r_cnjg(scomplex *, scomplex *);
    /* Local variables */
    real d__;
    aocl_int64_t i__, j, k;
    scomplex t;
    real r1, d11;
    scomplex d12;
    real d22;
    scomplex d21;
    aocl_int64_t kc, kk, kp;
    scomplex wk;
    aocl_int64_t kx;
    real tt;
    aocl_int64_t knc, kpc, npp;
    scomplex wkm1, wkp1;
    aocl_int64_t imax, jmax;
    real alpha;
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    aocl_int64_t kstep;
    logical upper;
    extern real slapy2_(real *, real *);
    real absakk;
    real colmax, rowmax;
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
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    --ipiv;
    --ap;
    /* Function Body */
    *info = 0;
    imax = 0;
    jmax = 0;
    kpc = 0;
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
        aocl_blas_xerbla("CHPTRF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Initialize ALPHA for use in choosing pivot block size. */
    alpha = (sqrt(17.f) + 1.f) / 8.f;
    if(upper)
    {
        /* Factorize A as U*D*U**H using the upper triangle of A */
        /* K is the main loop index, decreasing from N to 1 in steps of */
        /* 1 or 2 */
        k = *n;
        kc = (*n - 1) * *n / 2 + 1;
    L10:
        kpc = knc = kc;
        /* If K < 1, exit from loop */
        if(k < 1)
        {
            goto L110;
        }
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = kc + k - 1;
        absakk = (r__1 = ap[i__1].real, f2c_abs(r__1));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value */
        if(k > 1)
        {
            i__1 = k - 1;
            imax = aocl_blas_icamax(&i__1, &ap[kc], &c__1);
            i__1 = kc + imax - 1;
            colmax = (r__1 = ap[i__1].real, f2c_abs(r__1))
                     + (r__2 = r_imag(&ap[kc + imax - 1]), f2c_abs(r__2));
        }
        else
        {
            colmax = 0.f;
        }
        if(fla_max(absakk, colmax) == 0.f)
        {
            /* Column K is zero: set INFO and continue */
            if(*info == 0)
            {
                *info = k;
            }
            kp = k;
            i__1 = kc + k - 1;
            i__2 = kc + k - 1;
            r__1 = ap[i__2].real;
            ap[i__1].real = r__1;
            ap[i__1].imag = 0.f; // , expr subst
        }
        else
        {
            if(absakk >= alpha * colmax)
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value */
                rowmax = 0.f;
                /* jmax = imax; */
                kx = imax * (imax + 1) / 2 + imax;
                i__1 = k;
                for(j = imax + 1; j <= i__1; ++j)
                {
                    i__2 = kx;
                    if((r__1 = ap[i__2].real, f2c_abs(r__1)) + (r__2 = r_imag(&ap[kx]), f2c_abs(r__2))
                       > rowmax)
                    {
                        i__2 = kx;
                        rowmax = (r__1 = ap[i__2].real, f2c_abs(r__1))
                                 + (r__2 = r_imag(&ap[kx]), f2c_abs(r__2));
                        /* jmax = j; */
                    }
                    kx += j;
                    /* L20: */
                }
                kpc = (imax - 1) * imax / 2 + 1;
                if(imax > 1)
                {
                    i__1 = imax - 1;
                    jmax = aocl_blas_icamax(&i__1, &ap[kpc], &c__1);
                    /* Computing MAX */
                    i__1 = kpc + jmax - 1;
                    r__3 = rowmax;
                    r__4 = (r__1 = ap[i__1].real, f2c_abs(r__1))
                           + (r__2 = r_imag(&ap[kpc + jmax - 1]), f2c_abs(r__2)); // , expr subst
                    rowmax = fla_max(r__3, r__4);
                }
                if(absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else /* if(complicated condition) */
                {
                    i__1 = kpc + imax - 1;
                    if((r__1 = ap[i__1].real, f2c_abs(r__1)) >= alpha * rowmax)
                    {
                        /* interchange rows and columns K and IMAX, use 1-by-1 */
                        /* pivot block */
                        kp = imax;
                    }
                    else
                    {
                        /* interchange rows and columns K-1 and IMAX, use 2-by-2 */
                        /* pivot block */
                        kp = imax;
                        kstep = 2;
                    }
                }
            }
            kk = k - kstep + 1;
            if(kstep == 2)
            {
                knc = knc - k + 1;
            }
            if(kp != kk)
            {
                /* Interchange rows and columns KK and KP in the leading */
                /* submatrix A(1:k,1:k) */
                i__1 = kp - 1;
                aocl_blas_cswap(&i__1, &ap[knc], &c__1, &ap[kpc], &c__1);
                kx = kpc + kp - 1;
                i__1 = kk - 1;
                for(j = kp + 1; j <= i__1; ++j)
                {
                    kx = kx + j - 1;
                    r_cnjg(&q__1, &ap[knc + j - 1]);
                    t.real = q__1.real;
                    t.imag = q__1.imag; // , expr subst
                    i__2 = knc + j - 1;
                    r_cnjg(&q__1, &ap[kx]);
                    ap[i__2].real = q__1.real;
                    ap[i__2].imag = q__1.imag; // , expr subst
                    i__2 = kx;
                    ap[i__2].real = t.real;
                    ap[i__2].imag = t.imag; // , expr subst
                    /* L30: */
                }
                i__1 = kx + kk - 1;
                r_cnjg(&q__1, &ap[kx + kk - 1]);
                ap[i__1].real = q__1.real;
                ap[i__1].imag = q__1.imag; // , expr subst
                i__1 = knc + kk - 1;
                r1 = ap[i__1].real;
                i__1 = knc + kk - 1;
                i__2 = kpc + kp - 1;
                r__1 = ap[i__2].real;
                ap[i__1].real = r__1;
                ap[i__1].imag = 0.f; // , expr subst
                i__1 = kpc + kp - 1;
                ap[i__1].real = r1;
                ap[i__1].imag = 0.f; // , expr subst
                if(kstep == 2)
                {
                    i__1 = kc + k - 1;
                    i__2 = kc + k - 1;
                    r__1 = ap[i__2].real;
                    ap[i__1].real = r__1;
                    ap[i__1].imag = 0.f; // , expr subst
                    i__1 = kc + k - 2;
                    t.real = ap[i__1].real;
                    t.imag = ap[i__1].imag; // , expr subst
                    i__1 = kc + k - 2;
                    i__2 = kc + kp - 1;
                    ap[i__1].real = ap[i__2].real;
                    ap[i__1].imag = ap[i__2].imag; // , expr subst
                    i__1 = kc + kp - 1;
                    ap[i__1].real = t.real;
                    ap[i__1].imag = t.imag; // , expr subst
                }
            }
            else
            {
                i__1 = kc + k - 1;
                i__2 = kc + k - 1;
                r__1 = ap[i__2].real;
                ap[i__1].real = r__1;
                ap[i__1].imag = 0.f; // , expr subst
                if(kstep == 2)
                {
                    i__1 = kc - 1;
                    i__2 = kc - 1;
                    r__1 = ap[i__2].real;
                    ap[i__1].real = r__1;
                    ap[i__1].imag = 0.f; // , expr subst
                }
            }
            /* Update the leading submatrix */
            if(kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k now holds */
                /* W(k) = U(k)*D(k) */
                /* where U(k) is the k-th column of U */
                /* Perform a rank-1 update of A(1:k-1,1:k-1) as */
                /* A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H */
                i__1 = kc + k - 1;
                r1 = 1.f / ap[i__1].real;
                i__1 = k - 1;
                r__1 = -r1;
                aocl_blas_chpr(uplo, &i__1, &r__1, &ap[kc], &c__1, &ap[1]);
                /* Store U(k) in column k */
                i__1 = k - 1;
                aocl_blas_csscal(&i__1, &r1, &ap[kc], &c__1);
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns k and k-1 now hold */
                /* ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */
                /* where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                /* of U */
                /* Perform a rank-2 update of A(1:k-2,1:k-2) as */
                /* A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H */
                /* = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H */
                if(k > 2)
                {
                    i__1 = k - 1 + (k - 1) * k / 2;
                    r__1 = ap[i__1].real;
                    r__2 = r_imag(&ap[k - 1 + (k - 1) * k / 2]);
                    d__ = slapy2_(&r__1, &r__2);
                    i__1 = k - 1 + (k - 2) * (k - 1) / 2;
                    d22 = ap[i__1].real / d__;
                    i__1 = k + (k - 1) * k / 2;
                    d11 = ap[i__1].real / d__;
                    tt = 1.f / (d11 * d22 - 1.f);
                    i__1 = k - 1 + (k - 1) * k / 2;
                    q__1.real = ap[i__1].real / d__;
                    q__1.imag = ap[i__1].imag / d__; // , expr subst
                    d12.real = q__1.real;
                    d12.imag = q__1.imag; // , expr subst
                    d__ = tt / d__;
                    for(j = k - 2; j >= 1; --j)
                    {
                        i__1 = j + (k - 2) * (k - 1) / 2;
                        q__3.real = d11 * ap[i__1].real;
                        q__3.imag = d11 * ap[i__1].imag; // , expr subst
                        r_cnjg(&q__5, &d12);
                        i__2 = j + (k - 1) * k / 2;
                        q__4.real = q__5.real * ap[i__2].real - q__5.imag * ap[i__2].imag;
                        q__4.imag = q__5.real * ap[i__2].imag + q__5.imag * ap[i__2].real; // , expr subst
                        q__2.real = q__3.real - q__4.real;
                        q__2.imag = q__3.imag - q__4.imag; // , expr subst
                        q__1.real = d__ * q__2.real;
                        q__1.imag = d__ * q__2.imag; // , expr subst
                        wkm1.real = q__1.real;
                        wkm1.imag = q__1.imag; // , expr subst
                        i__1 = j + (k - 1) * k / 2;
                        q__3.real = d22 * ap[i__1].real;
                        q__3.imag = d22 * ap[i__1].imag; // , expr subst
                        i__2 = j + (k - 2) * (k - 1) / 2;
                        q__4.real = d12.real * ap[i__2].real - d12.imag * ap[i__2].imag;
                        q__4.imag = d12.real * ap[i__2].imag + d12.imag * ap[i__2].real; // , expr subst
                        q__2.real = q__3.real - q__4.real;
                        q__2.imag = q__3.imag - q__4.imag; // , expr subst
                        q__1.real = d__ * q__2.real;
                        q__1.imag = d__ * q__2.imag; // , expr subst
                        wk.real = q__1.real;
                        wk.imag = q__1.imag; // , expr subst
                        for(i__ = j; i__ >= 1; --i__)
                        {
                            i__1 = i__ + (j - 1) * j / 2;
                            i__2 = i__ + (j - 1) * j / 2;
                            i__3 = i__ + (k - 1) * k / 2;
                            r_cnjg(&q__4, &wk);
                            q__3.real = ap[i__3].real * q__4.real - ap[i__3].imag * q__4.imag;
                            q__3.imag = ap[i__3].real * q__4.imag + ap[i__3].imag * q__4.real; // , expr subst
                            q__2.real = ap[i__2].real - q__3.real;
                            q__2.imag = ap[i__2].imag - q__3.imag; // , expr subst
                            i__4 = i__ + (k - 2) * (k - 1) / 2;
                            r_cnjg(&q__6, &wkm1);
                            q__5.real = ap[i__4].real * q__6.real - ap[i__4].imag * q__6.imag;
                            q__5.imag = ap[i__4].real * q__6.imag + ap[i__4].imag * q__6.real; // , expr subst
                            q__1.real = q__2.real - q__5.real;
                            q__1.imag = q__2.imag - q__5.imag; // , expr subst
                            ap[i__1].real = q__1.real;
                            ap[i__1].imag = q__1.imag; // , expr subst
                            /* L40: */
                        }
                        i__1 = j + (k - 1) * k / 2;
                        ap[i__1].real = wk.real;
                        ap[i__1].imag = wk.imag; // , expr subst
                        i__1 = j + (k - 2) * (k - 1) / 2;
                        ap[i__1].real = wkm1.real;
                        ap[i__1].imag = wkm1.imag; // , expr subst
                        i__1 = j + (j - 1) * j / 2;
                        i__2 = j + (j - 1) * j / 2;
                        r__1 = ap[i__2].real;
                        q__1.real = r__1;
                        q__1.imag = 0.f; // , expr subst
                        ap[i__1].real = q__1.real;
                        ap[i__1].imag = q__1.imag; // , expr subst
                        /* L50: */
                    }
                }
            }
        }
        /* Store details of the interchanges in IPIV */
        if(kstep == 1)
        {
            ipiv[k] = (aocl_int_t)(kp);
        }
        else
        {
            ipiv[k] = (aocl_int_t)(-kp);
            ipiv[k - 1] = (aocl_int_t)(-kp);
        }
        /* Decrease K and return to the start of the main loop */
        k -= kstep;
        kc = knc - k;
        goto L10;
    }
    else
    {
        /* Factorize A as L*D*L**H using the lower triangle of A */
        /* K is the main loop index, increasing from 1 to N in steps of */
        /* 1 or 2 */
        k = 1;
        kc = 1;
        npp = *n * (*n + 1) / 2;
    L60:
        knc = kc;
        /* If K > N, exit from loop */
        if(k > *n)
        {
            goto L110;
        }
        kstep = 1;
        /* Determine rows and columns to be interchanged and whether */
        /* a 1-by-1 or 2-by-2 pivot block will be used */
        i__1 = kc;
        absakk = (r__1 = ap[i__1].real, f2c_abs(r__1));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value */
        if(k < *n)
        {
            i__1 = *n - k;
            imax = k + aocl_blas_icamax(&i__1, &ap[kc + 1], &c__1);
            i__1 = kc + imax - k;
            colmax = (r__1 = ap[i__1].real, f2c_abs(r__1))
                     + (r__2 = r_imag(&ap[kc + imax - k]), f2c_abs(r__2));
        }
        else
        {
            colmax = 0.f;
        }
        if(fla_max(absakk, colmax) == 0.f)
        {
            /* Column K is zero: set INFO and continue */
            if(*info == 0)
            {
                *info = k;
            }
            kp = k;
            i__1 = kc;
            i__2 = kc;
            r__1 = ap[i__2].real;
            ap[i__1].real = r__1;
            ap[i__1].imag = 0.f; // , expr subst
        }
        else
        {
            if(absakk >= alpha * colmax)
            {
                /* no interchange, use 1-by-1 pivot block */
                kp = k;
            }
            else
            {
                /* JMAX is the column-index of the largest off-diagonal */
                /* element in row IMAX, and ROWMAX is its absolute value */
                rowmax = 0.f;
                kx = kc + imax - k;
                i__1 = imax - 1;
                for(j = k; j <= i__1; ++j)
                {
                    i__2 = kx;
                    if((r__1 = ap[i__2].real, f2c_abs(r__1)) + (r__2 = r_imag(&ap[kx]), f2c_abs(r__2))
                       > rowmax)
                    {
                        i__2 = kx;
                        rowmax = (r__1 = ap[i__2].real, f2c_abs(r__1))
                                 + (r__2 = r_imag(&ap[kx]), f2c_abs(r__2));
                        /* jmax = j; */
                    }
                    kx = kx + *n - j;
                    /* L70: */
                }
                kpc = npp - (*n - imax + 1) * (*n - imax + 2) / 2 + 1;
                if(imax < *n)
                {
                    i__1 = *n - imax;
                    jmax = imax + aocl_blas_icamax(&i__1, &ap[kpc + 1], &c__1);
                    /* Computing MAX */
                    i__1 = kpc + jmax - imax;
                    r__3 = rowmax;
                    r__4 = (r__1 = ap[i__1].real, f2c_abs(r__1))
                           + (r__2 = r_imag(&ap[kpc + jmax - imax]), f2c_abs(r__2)); // , expr subst
                    rowmax = fla_max(r__3, r__4);
                }
                if(absakk >= alpha * colmax * (colmax / rowmax))
                {
                    /* no interchange, use 1-by-1 pivot block */
                    kp = k;
                }
                else /* if(complicated condition) */
                {
                    i__1 = kpc;
                    if((r__1 = ap[i__1].real, f2c_abs(r__1)) >= alpha * rowmax)
                    {
                        /* interchange rows and columns K and IMAX, use 1-by-1 */
                        /* pivot block */
                        kp = imax;
                    }
                    else
                    {
                        /* interchange rows and columns K+1 and IMAX, use 2-by-2 */
                        /* pivot block */
                        kp = imax;
                        kstep = 2;
                    }
                }
            }
            kk = k + kstep - 1;
            if(kstep == 2)
            {
                knc = knc + *n - k + 1;
            }
            if(kp != kk)
            {
                /* Interchange rows and columns KK and KP in the trailing */
                /* submatrix A(k:n,k:n) */
                if(kp < *n)
                {
                    i__1 = *n - kp;
                    aocl_blas_cswap(&i__1, &ap[knc + kp - kk + 1], &c__1, &ap[kpc + 1], &c__1);
                }
                kx = knc + kp - kk;
                i__1 = kp - 1;
                for(j = kk + 1; j <= i__1; ++j)
                {
                    kx = kx + *n - j + 1;
                    r_cnjg(&q__1, &ap[knc + j - kk]);
                    t.real = q__1.real;
                    t.imag = q__1.imag; // , expr subst
                    i__2 = knc + j - kk;
                    r_cnjg(&q__1, &ap[kx]);
                    ap[i__2].real = q__1.real;
                    ap[i__2].imag = q__1.imag; // , expr subst
                    i__2 = kx;
                    ap[i__2].real = t.real;
                    ap[i__2].imag = t.imag; // , expr subst
                    /* L80: */
                }
                i__1 = knc + kp - kk;
                r_cnjg(&q__1, &ap[knc + kp - kk]);
                ap[i__1].real = q__1.real;
                ap[i__1].imag = q__1.imag; // , expr subst
                i__1 = knc;
                r1 = ap[i__1].real;
                i__1 = knc;
                i__2 = kpc;
                r__1 = ap[i__2].real;
                ap[i__1].real = r__1;
                ap[i__1].imag = 0.f; // , expr subst
                i__1 = kpc;
                ap[i__1].real = r1;
                ap[i__1].imag = 0.f; // , expr subst
                if(kstep == 2)
                {
                    i__1 = kc;
                    i__2 = kc;
                    r__1 = ap[i__2].real;
                    ap[i__1].real = r__1;
                    ap[i__1].imag = 0.f; // , expr subst
                    i__1 = kc + 1;
                    t.real = ap[i__1].real;
                    t.imag = ap[i__1].imag; // , expr subst
                    i__1 = kc + 1;
                    i__2 = kc + kp - k;
                    ap[i__1].real = ap[i__2].real;
                    ap[i__1].imag = ap[i__2].imag; // , expr subst
                    i__1 = kc + kp - k;
                    ap[i__1].real = t.real;
                    ap[i__1].imag = t.imag; // , expr subst
                }
            }
            else
            {
                i__1 = kc;
                i__2 = kc;
                r__1 = ap[i__2].real;
                ap[i__1].real = r__1;
                ap[i__1].imag = 0.f; // , expr subst
                if(kstep == 2)
                {
                    i__1 = knc;
                    i__2 = knc;
                    r__1 = ap[i__2].real;
                    ap[i__1].real = r__1;
                    ap[i__1].imag = 0.f; // , expr subst
                }
            }
            /* Update the trailing submatrix */
            if(kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k now holds */
                /* W(k) = L(k)*D(k) */
                /* where L(k) is the k-th column of L */
                if(k < *n)
                {
                    /* Perform a rank-1 update of A(k+1:n,k+1:n) as */
                    /* A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H */
                    i__1 = kc;
                    r1 = 1.f / ap[i__1].real;
                    i__1 = *n - k;
                    r__1 = -r1;
                    aocl_blas_chpr(uplo, &i__1, &r__1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1]);
                    /* Store L(k) in column K */
                    i__1 = *n - k;
                    aocl_blas_csscal(&i__1, &r1, &ap[kc + 1], &c__1);
                }
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns K and K+1 now hold */
                /* ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k) */
                /* where L(k) and L(k+1) are the k-th and (k+1)-th columns */
                /* of L */
                if(k < *n - 1)
                {
                    /* Perform a rank-2 update of A(k+2:n,k+2:n) as */
                    /* A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H */
                    /* = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H */
                    /* where L(k) and L(k+1) are the k-th and (k+1)-th */
                    /* columns of L */
                    i__1 = k + 1 + (k - 1) * ((*n << 1) - k) / 2;
                    r__1 = ap[i__1].real;
                    r__2 = r_imag(&ap[k + 1 + (k - 1) * ((*n << 1) - k) / 2]);
                    d__ = slapy2_(&r__1, &r__2);
                    i__1 = k + 1 + k * ((*n << 1) - k - 1) / 2;
                    d11 = ap[i__1].real / d__;
                    i__1 = k + (k - 1) * ((*n << 1) - k) / 2;
                    d22 = ap[i__1].real / d__;
                    tt = 1.f / (d11 * d22 - 1.f);
                    i__1 = k + 1 + (k - 1) * ((*n << 1) - k) / 2;
                    q__1.real = ap[i__1].real / d__;
                    q__1.imag = ap[i__1].imag / d__; // , expr subst
                    d21.real = q__1.real;
                    d21.imag = q__1.imag; // , expr subst
                    d__ = tt / d__;
                    i__1 = *n;
                    for(j = k + 2; j <= i__1; ++j)
                    {
                        i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
                        q__3.real = d11 * ap[i__2].real;
                        q__3.imag = d11 * ap[i__2].imag; // , expr subst
                        i__3 = j + k * ((*n << 1) - k - 1) / 2;
                        q__4.real = d21.real * ap[i__3].real - d21.imag * ap[i__3].imag;
                        q__4.imag = d21.real * ap[i__3].imag + d21.imag * ap[i__3].real; // , expr subst
                        q__2.real = q__3.real - q__4.real;
                        q__2.imag = q__3.imag - q__4.imag; // , expr subst
                        q__1.real = d__ * q__2.real;
                        q__1.imag = d__ * q__2.imag; // , expr subst
                        wk.real = q__1.real;
                        wk.imag = q__1.imag; // , expr subst
                        i__2 = j + k * ((*n << 1) - k - 1) / 2;
                        q__3.real = d22 * ap[i__2].real;
                        q__3.imag = d22 * ap[i__2].imag; // , expr subst
                        r_cnjg(&q__5, &d21);
                        i__3 = j + (k - 1) * ((*n << 1) - k) / 2;
                        q__4.real = q__5.real * ap[i__3].real - q__5.imag * ap[i__3].imag;
                        q__4.imag = q__5.real * ap[i__3].imag + q__5.imag * ap[i__3].real; // , expr subst
                        q__2.real = q__3.real - q__4.real;
                        q__2.imag = q__3.imag - q__4.imag; // , expr subst
                        q__1.real = d__ * q__2.real;
                        q__1.imag = d__ * q__2.imag; // , expr subst
                        wkp1.real = q__1.real;
                        wkp1.imag = q__1.imag; // , expr subst
                        i__2 = *n;
                        for(i__ = j; i__ <= i__2; ++i__)
                        {
                            i__3 = i__ + (j - 1) * ((*n << 1) - j) / 2;
                            i__4 = i__ + (j - 1) * ((*n << 1) - j) / 2;
                            i__5 = i__ + (k - 1) * ((*n << 1) - k) / 2;
                            r_cnjg(&q__4, &wk);
                            q__3.real = ap[i__5].real * q__4.real - ap[i__5].imag * q__4.imag;
                            q__3.imag = ap[i__5].real * q__4.imag + ap[i__5].imag * q__4.real; // , expr subst
                            q__2.real = ap[i__4].real - q__3.real;
                            q__2.imag = ap[i__4].imag - q__3.imag; // , expr subst
                            i__6 = i__ + k * ((*n << 1) - k - 1) / 2;
                            r_cnjg(&q__6, &wkp1);
                            q__5.real = ap[i__6].real * q__6.real - ap[i__6].imag * q__6.imag;
                            q__5.imag = ap[i__6].real * q__6.imag + ap[i__6].imag * q__6.real; // , expr subst
                            q__1.real = q__2.real - q__5.real;
                            q__1.imag = q__2.imag - q__5.imag; // , expr subst
                            ap[i__3].real = q__1.real;
                            ap[i__3].imag = q__1.imag; // , expr subst
                            /* L90: */
                        }
                        i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
                        ap[i__2].real = wk.real;
                        ap[i__2].imag = wk.imag; // , expr subst
                        i__2 = j + k * ((*n << 1) - k - 1) / 2;
                        ap[i__2].real = wkp1.real;
                        ap[i__2].imag = wkp1.imag; // , expr subst
                        i__2 = j + (j - 1) * ((*n << 1) - j) / 2;
                        i__3 = j + (j - 1) * ((*n << 1) - j) / 2;
                        r__1 = ap[i__3].real;
                        q__1.real = r__1;
                        q__1.imag = 0.f; // , expr subst
                        ap[i__2].real = q__1.real;
                        ap[i__2].imag = q__1.imag; // , expr subst
                        /* L100: */
                    }
                }
            }
        }
        /* Store details of the interchanges in IPIV */
        if(kstep == 1)
        {
            ipiv[k] = (aocl_int_t)(kp);
        }
        else
        {
            ipiv[k] = (aocl_int_t)(-kp);
            ipiv[k + 1] = (aocl_int_t)(-kp);
        }
        /* Increase K and return to the start of the main loop */
        k += kstep;
        kc = knc + *n - k + 2;
        goto L60;
    }
L110:
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CHPTRF */
}
/* chptrf_ */
