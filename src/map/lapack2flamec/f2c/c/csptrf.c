/* ../netlib/csptrf.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */

/*
 *     Modifications Copyright (c) 2024 Advanced Micro Devices, Inc.  All rights reserved.
 */

#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
static integer c__1 = 1;
/* > \brief \b CSPTRF */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CSPTRF + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csptrf.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csptrf.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csptrf.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CSPTRF( UPLO, N, AP, IPIV, INFO ) */
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
/* > CSPTRF computes the factorization of a complex symmetric matrix A */
/* > stored in packed format using the Bunch-Kaufman diagonal pivoting */
/* > method: */
/* > */
/* > A = U*D*U**T or A = L*D*L**T */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is symmetric and block diagonal with */
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
/* > On entry, the upper or lower triangle of the symmetric matrix */
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
/* > 5-96 - Based on modifications by J. Lewis, Boeing Computer Services */
/* > Company */
/* > */
/* > If UPLO = 'U', then A = U*D*U**T, where */
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
/* > If UPLO = 'L', then A = L*D*L**T, where */
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
/* > */
/* ===================================================================== */
/* Subroutine */
void csptrf_(char *uplo, integer *n, complex *ap, integer *ipiv, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256, "csptrf inputs: uplo %c, n %lld", *uplo, *n);
#else
    snprintf(buffer, 256, "csptrf inputs: uplo %c, n %d", *uplo, *n);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4;
    /* Builtin functions */
    double sqrt(doublereal), r_imag(complex *);
    void c_div(complex *, complex *, complex *);
    /* Local variables */
    integer i__, j, k;
    complex t, r1, d11, d12, d21, d22;
    integer kc, kk, kp;
    complex wk;
    integer kx, knc, kpc, npp;
    complex wkm1, wkp1;
    integer imax, jmax;
    extern /* Subroutine */
        void
        cspr_(char *, integer *, complex *, complex *, integer *, complex *);
    real alpha;
    extern /* Subroutine */
        void
        cscal_(integer *, complex *, complex *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        cswap_(integer *, complex *, integer *, complex *, integer *);
    integer kstep;
    logical upper;
    real absakk;
    extern integer icamax_(integer *, complex *, integer *);
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
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
    upper = lsame_(uplo, "U", 1, 1);
    imax = 0;
    kpc = 0;
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
        xerbla_("CSPTRF", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Initialize ALPHA for use in choosing pivot block size. */
    alpha = (sqrt(17.f) + 1.f) / 8.f;
    if(upper)
    {
        /* Factorize A as U*D*U**T using the upper triangle of A */
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
        absakk
            = (r__1 = ap[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&ap[kc + k - 1]), f2c_abs(r__2));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value */
        if(k > 1)
        {
            i__1 = k - 1;
            imax = icamax_(&i__1, &ap[kc], &c__1);
            i__1 = kc + imax - 1;
            colmax = (r__1 = ap[i__1].r, f2c_abs(r__1))
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
                rowmax = 0.f;
                /* jmax = imax; */
                kx = imax * (imax + 1) / 2 + imax;
                i__1 = k;
                for(j = imax + 1; j <= i__1; ++j)
                {
                    i__2 = kx;
                    if((r__1 = ap[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&ap[kx]), f2c_abs(r__2))
                       > rowmax)
                    {
                        i__2 = kx;
                        rowmax = (r__1 = ap[i__2].r, f2c_abs(r__1))
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
                    jmax = icamax_(&i__1, &ap[kpc], &c__1);
                    /* Computing MAX */
                    i__1 = kpc + jmax - 1;
                    r__3 = rowmax;
                    r__4 = (r__1 = ap[i__1].r, f2c_abs(r__1))
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
                    if((r__1 = ap[i__1].r, f2c_abs(r__1))
                           + (r__2 = r_imag(&ap[kpc + imax - 1]), f2c_abs(r__2))
                       >= alpha * rowmax)
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
                cswap_(&i__1, &ap[knc], &c__1, &ap[kpc], &c__1);
                kx = kpc + kp - 1;
                i__1 = kk - 1;
                for(j = kp + 1; j <= i__1; ++j)
                {
                    kx = kx + j - 1;
                    i__2 = knc + j - 1;
                    t.r = ap[i__2].r;
                    t.i = ap[i__2].i; // , expr subst
                    i__2 = knc + j - 1;
                    i__3 = kx;
                    ap[i__2].r = ap[i__3].r;
                    ap[i__2].i = ap[i__3].i; // , expr subst
                    i__2 = kx;
                    ap[i__2].r = t.r;
                    ap[i__2].i = t.i; // , expr subst
                    /* L30: */
                }
                i__1 = knc + kk - 1;
                t.r = ap[i__1].r;
                t.i = ap[i__1].i; // , expr subst
                i__1 = knc + kk - 1;
                i__2 = kpc + kp - 1;
                ap[i__1].r = ap[i__2].r;
                ap[i__1].i = ap[i__2].i; // , expr subst
                i__1 = kpc + kp - 1;
                ap[i__1].r = t.r;
                ap[i__1].i = t.i; // , expr subst
                if(kstep == 2)
                {
                    i__1 = kc + k - 2;
                    t.r = ap[i__1].r;
                    t.i = ap[i__1].i; // , expr subst
                    i__1 = kc + k - 2;
                    i__2 = kc + kp - 1;
                    ap[i__1].r = ap[i__2].r;
                    ap[i__1].i = ap[i__2].i; // , expr subst
                    i__1 = kc + kp - 1;
                    ap[i__1].r = t.r;
                    ap[i__1].i = t.i; // , expr subst
                }
            }
            /* Update the leading submatrix */
            if(kstep == 1)
            {
                /* 1-by-1 pivot block D(k): column k now holds */
                /* W(k) = U(k)*D(k) */
                /* where U(k) is the k-th column of U */
                /* Perform a rank-1 update of A(1:k-1,1:k-1) as */
                /* A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T */
                c_div(&q__1, &c_b1, &ap[kc + k - 1]);
                r1.r = q__1.r;
                r1.i = q__1.i; // , expr subst
                i__1 = k - 1;
                q__1.r = -r1.r;
                q__1.i = -r1.i; // , expr subst
                cspr_(uplo, &i__1, &q__1, &ap[kc], &c__1, &ap[1]);
                /* Store U(k) in column k */
                i__1 = k - 1;
                cscal_(&i__1, &r1, &ap[kc], &c__1);
            }
            else
            {
                /* 2-by-2 pivot block D(k): columns k and k-1 now hold */
                /* ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k) */
                /* where U(k) and U(k-1) are the k-th and (k-1)-th columns */
                /* of U */
                /* Perform a rank-2 update of A(1:k-2,1:k-2) as */
                /* A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T */
                /* = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T */
                if(k > 2)
                {
                    i__1 = k - 1 + (k - 1) * k / 2;
                    d12.r = ap[i__1].r;
                    d12.i = ap[i__1].i; // , expr subst
                    c_div(&q__1, &ap[k - 1 + (k - 2) * (k - 1) / 2], &d12);
                    d22.r = q__1.r;
                    d22.i = q__1.i; // , expr subst
                    c_div(&q__1, &ap[k + (k - 1) * k / 2], &d12);
                    d11.r = q__1.r;
                    d11.i = q__1.i; // , expr subst
                    q__3.r = d11.r * d22.r - d11.i * d22.i;
                    q__3.i = d11.r * d22.i + d11.i * d22.r; // , expr subst
                    q__2.r = q__3.r - 1.f;
                    q__2.i = q__3.i - 0.f; // , expr subst
                    c_div(&q__1, &c_b1, &q__2);
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                    c_div(&q__1, &t, &d12);
                    d12.r = q__1.r;
                    d12.i = q__1.i; // , expr subst
                    for(j = k - 2; j >= 1; --j)
                    {
                        i__1 = j + (k - 2) * (k - 1) / 2;
                        q__3.r = d11.r * ap[i__1].r - d11.i * ap[i__1].i;
                        q__3.i = d11.r * ap[i__1].i + d11.i * ap[i__1].r; // , expr subst
                        i__2 = j + (k - 1) * k / 2;
                        q__2.r = q__3.r - ap[i__2].r;
                        q__2.i = q__3.i - ap[i__2].i; // , expr subst
                        q__1.r = d12.r * q__2.r - d12.i * q__2.i;
                        q__1.i = d12.r * q__2.i + d12.i * q__2.r; // , expr subst
                        wkm1.r = q__1.r;
                        wkm1.i = q__1.i; // , expr subst
                        i__1 = j + (k - 1) * k / 2;
                        q__3.r = d22.r * ap[i__1].r - d22.i * ap[i__1].i;
                        q__3.i = d22.r * ap[i__1].i + d22.i * ap[i__1].r; // , expr subst
                        i__2 = j + (k - 2) * (k - 1) / 2;
                        q__2.r = q__3.r - ap[i__2].r;
                        q__2.i = q__3.i - ap[i__2].i; // , expr subst
                        q__1.r = d12.r * q__2.r - d12.i * q__2.i;
                        q__1.i = d12.r * q__2.i + d12.i * q__2.r; // , expr subst
                        wk.r = q__1.r;
                        wk.i = q__1.i; // , expr subst
                        for(i__ = j; i__ >= 1; --i__)
                        {
                            i__1 = i__ + (j - 1) * j / 2;
                            i__2 = i__ + (j - 1) * j / 2;
                            i__3 = i__ + (k - 1) * k / 2;
                            q__3.r = ap[i__3].r * wk.r - ap[i__3].i * wk.i;
                            q__3.i = ap[i__3].r * wk.i + ap[i__3].i * wk.r; // , expr subst
                            q__2.r = ap[i__2].r - q__3.r;
                            q__2.i = ap[i__2].i - q__3.i; // , expr subst
                            i__4 = i__ + (k - 2) * (k - 1) / 2;
                            q__4.r = ap[i__4].r * wkm1.r - ap[i__4].i * wkm1.i;
                            q__4.i = ap[i__4].r * wkm1.i + ap[i__4].i * wkm1.r; // , expr subst
                            q__1.r = q__2.r - q__4.r;
                            q__1.i = q__2.i - q__4.i; // , expr subst
                            ap[i__1].r = q__1.r;
                            ap[i__1].i = q__1.i; // , expr subst
                            /* L40: */
                        }
                        i__1 = j + (k - 1) * k / 2;
                        ap[i__1].r = wk.r;
                        ap[i__1].i = wk.i; // , expr subst
                        i__1 = j + (k - 2) * (k - 1) / 2;
                        ap[i__1].r = wkm1.r;
                        ap[i__1].i = wkm1.i; // , expr subst
                        /* L50: */
                    }
                }
            }
        }
        /* Store details of the interchanges in IPIV */
        if(kstep == 1)
        {
            ipiv[k] = kp;
        }
        else
        {
            ipiv[k] = -kp;
            ipiv[k - 1] = -kp;
        }
        /* Decrease K and return to the start of the main loop */
        k -= kstep;
        kc = knc - k;
        goto L10;
    }
    else
    {
        /* Factorize A as L*D*L**T using the lower triangle of A */
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
        absakk = (r__1 = ap[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&ap[kc]), f2c_abs(r__2));
        /* IMAX is the row-index of the largest off-diagonal element in */
        /* column K, and COLMAX is its absolute value */
        if(k < *n)
        {
            i__1 = *n - k;
            imax = k + icamax_(&i__1, &ap[kc + 1], &c__1);
            i__1 = kc + imax - k;
            colmax = (r__1 = ap[i__1].r, f2c_abs(r__1))
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
                    if((r__1 = ap[i__2].r, f2c_abs(r__1)) + (r__2 = r_imag(&ap[kx]), f2c_abs(r__2))
                       > rowmax)
                    {
                        i__2 = kx;
                        rowmax = (r__1 = ap[i__2].r, f2c_abs(r__1))
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
                    jmax = imax + icamax_(&i__1, &ap[kpc + 1], &c__1);
                    /* Computing MAX */
                    i__1 = kpc + jmax - imax;
                    r__3 = rowmax;
                    r__4 = (r__1 = ap[i__1].r, f2c_abs(r__1))
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
                    if((r__1 = ap[i__1].r, f2c_abs(r__1)) + (r__2 = r_imag(&ap[kpc]), f2c_abs(r__2))
                       >= alpha * rowmax)
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
                    cswap_(&i__1, &ap[knc + kp - kk + 1], &c__1, &ap[kpc + 1], &c__1);
                }
                kx = knc + kp - kk;
                i__1 = kp - 1;
                for(j = kk + 1; j <= i__1; ++j)
                {
                    kx = kx + *n - j + 1;
                    i__2 = knc + j - kk;
                    t.r = ap[i__2].r;
                    t.i = ap[i__2].i; // , expr subst
                    i__2 = knc + j - kk;
                    i__3 = kx;
                    ap[i__2].r = ap[i__3].r;
                    ap[i__2].i = ap[i__3].i; // , expr subst
                    i__2 = kx;
                    ap[i__2].r = t.r;
                    ap[i__2].i = t.i; // , expr subst
                    /* L80: */
                }
                i__1 = knc;
                t.r = ap[i__1].r;
                t.i = ap[i__1].i; // , expr subst
                i__1 = knc;
                i__2 = kpc;
                ap[i__1].r = ap[i__2].r;
                ap[i__1].i = ap[i__2].i; // , expr subst
                i__1 = kpc;
                ap[i__1].r = t.r;
                ap[i__1].i = t.i; // , expr subst
                if(kstep == 2)
                {
                    i__1 = kc + 1;
                    t.r = ap[i__1].r;
                    t.i = ap[i__1].i; // , expr subst
                    i__1 = kc + 1;
                    i__2 = kc + kp - k;
                    ap[i__1].r = ap[i__2].r;
                    ap[i__1].i = ap[i__2].i; // , expr subst
                    i__1 = kc + kp - k;
                    ap[i__1].r = t.r;
                    ap[i__1].i = t.i; // , expr subst
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
                    /* A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T */
                    c_div(&q__1, &c_b1, &ap[kc]);
                    r1.r = q__1.r;
                    r1.i = q__1.i; // , expr subst
                    i__1 = *n - k;
                    q__1.r = -r1.r;
                    q__1.i = -r1.i; // , expr subst
                    cspr_(uplo, &i__1, &q__1, &ap[kc + 1], &c__1, &ap[kc + *n - k + 1]);
                    /* Store L(k) in column K */
                    i__1 = *n - k;
                    cscal_(&i__1, &r1, &ap[kc + 1], &c__1);
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
                    /* A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T */
                    /* = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T */
                    /* where L(k) and L(k+1) are the k-th and (k+1)-th */
                    /* columns of L */
                    i__1 = k + 1 + (k - 1) * ((*n << 1) - k) / 2;
                    d21.r = ap[i__1].r;
                    d21.i = ap[i__1].i; // , expr subst
                    c_div(&q__1, &ap[k + 1 + k * ((*n << 1) - k - 1) / 2], &d21);
                    d11.r = q__1.r;
                    d11.i = q__1.i; // , expr subst
                    c_div(&q__1, &ap[k + (k - 1) * ((*n << 1) - k) / 2], &d21);
                    d22.r = q__1.r;
                    d22.i = q__1.i; // , expr subst
                    q__3.r = d11.r * d22.r - d11.i * d22.i;
                    q__3.i = d11.r * d22.i + d11.i * d22.r; // , expr subst
                    q__2.r = q__3.r - 1.f;
                    q__2.i = q__3.i - 0.f; // , expr subst
                    c_div(&q__1, &c_b1, &q__2);
                    t.r = q__1.r;
                    t.i = q__1.i; // , expr subst
                    c_div(&q__1, &t, &d21);
                    d21.r = q__1.r;
                    d21.i = q__1.i; // , expr subst
                    i__1 = *n;
                    for(j = k + 2; j <= i__1; ++j)
                    {
                        i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
                        q__3.r = d11.r * ap[i__2].r - d11.i * ap[i__2].i;
                        q__3.i = d11.r * ap[i__2].i + d11.i * ap[i__2].r; // , expr subst
                        i__3 = j + k * ((*n << 1) - k - 1) / 2;
                        q__2.r = q__3.r - ap[i__3].r;
                        q__2.i = q__3.i - ap[i__3].i; // , expr subst
                        q__1.r = d21.r * q__2.r - d21.i * q__2.i;
                        q__1.i = d21.r * q__2.i + d21.i * q__2.r; // , expr subst
                        wk.r = q__1.r;
                        wk.i = q__1.i; // , expr subst
                        i__2 = j + k * ((*n << 1) - k - 1) / 2;
                        q__3.r = d22.r * ap[i__2].r - d22.i * ap[i__2].i;
                        q__3.i = d22.r * ap[i__2].i + d22.i * ap[i__2].r; // , expr subst
                        i__3 = j + (k - 1) * ((*n << 1) - k) / 2;
                        q__2.r = q__3.r - ap[i__3].r;
                        q__2.i = q__3.i - ap[i__3].i; // , expr subst
                        q__1.r = d21.r * q__2.r - d21.i * q__2.i;
                        q__1.i = d21.r * q__2.i + d21.i * q__2.r; // , expr subst
                        wkp1.r = q__1.r;
                        wkp1.i = q__1.i; // , expr subst
                        i__2 = *n;
                        for(i__ = j; i__ <= i__2; ++i__)
                        {
                            i__3 = i__ + (j - 1) * ((*n << 1) - j) / 2;
                            i__4 = i__ + (j - 1) * ((*n << 1) - j) / 2;
                            i__5 = i__ + (k - 1) * ((*n << 1) - k) / 2;
                            q__3.r = ap[i__5].r * wk.r - ap[i__5].i * wk.i;
                            q__3.i = ap[i__5].r * wk.i + ap[i__5].i * wk.r; // , expr subst
                            q__2.r = ap[i__4].r - q__3.r;
                            q__2.i = ap[i__4].i - q__3.i; // , expr subst
                            i__6 = i__ + k * ((*n << 1) - k - 1) / 2;
                            q__4.r = ap[i__6].r * wkp1.r - ap[i__6].i * wkp1.i;
                            q__4.i = ap[i__6].r * wkp1.i + ap[i__6].i * wkp1.r; // , expr subst
                            q__1.r = q__2.r - q__4.r;
                            q__1.i = q__2.i - q__4.i; // , expr subst
                            ap[i__3].r = q__1.r;
                            ap[i__3].i = q__1.i; // , expr subst
                            /* L90: */
                        }
                        i__2 = j + (k - 1) * ((*n << 1) - k) / 2;
                        ap[i__2].r = wk.r;
                        ap[i__2].i = wk.i; // , expr subst
                        i__2 = j + k * ((*n << 1) - k - 1) / 2;
                        ap[i__2].r = wkp1.r;
                        ap[i__2].i = wkp1.i; // , expr subst
                        /* L100: */
                    }
                }
            }
        }
        /* Store details of the interchanges in IPIV */
        if(kstep == 1)
        {
            ipiv[k] = kp;
        }
        else
        {
            ipiv[k] = -kp;
            ipiv[k + 1] = -kp;
        }
        /* Increase K and return to the start of the main loop */
        k += kstep;
        kc = knc + *n - k + 2;
        goto L60;
    }
L110:
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CSPTRF */
}
/* csptrf_ */
