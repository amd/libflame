/* ../netlib/cpbrfs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
static integer c__1 = 1;
/* > \brief \b CPBRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CPBRFS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbrfs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbrfs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbrfs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, */
/* LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER UPLO */
/* INTEGER INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* REAL BERR( * ), FERR( * ), RWORK( * ) */
/* COMPLEX AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/* $ WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian positive definite */
/* > and banded, and provides error bounds and backward error estimates */
/* > for the solution. */
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
/* > \param[in] KD */
/* > \verbatim */
/* > KD is INTEGER */
/* > The number of superdiagonals of the matrix A if UPLO = 'U', */
/* > or the number of subdiagonals if UPLO = 'L'. KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > The upper or lower triangle of the Hermitian band matrix A, */
/* > stored in the first KD+1 rows of the array. The j-th column */
/* > of A is stored in the j-th column of the array AB as follows: */
/* > if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j;
 */
/* > if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=fla_min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* > AFB is COMPLEX array, dimension (LDAFB,N) */
/* > The triangular factor U or L from the Cholesky factorization */
/* > A = U**H*U or A = L*L**H of the band matrix A as computed by */
/* > CPBTRF, in the same storage format as A (see AB). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* > LDAFB is INTEGER */
/* > The leading dimension of the array AFB. LDAFB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > The right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* > X is COMPLEX array, dimension (LDX,NRHS) */
/* > On entry, the solution matrix X, as computed by CPBTRS. */
/* > On exit, the improved solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* > FERR is REAL array, dimension (NRHS) */
/* > The estimated forward error bound for each solution vector */
/* > X(j) (the j-th column of the solution matrix X). */
/* > If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* > is an estimated upper bound for the magnitude of the largest */
/* > element in (X(j) - XTRUE) divided by the magnitude of the */
/* > largest element in X(j). The estimate is as reliable as */
/* > the estimate for RCOND, and is almost always a slight */
/* > overestimate of the true error. */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* > BERR is REAL array, dimension (NRHS) */
/* > The componentwise relative backward error of each solution */
/* > vector X(j) (i.e., the smallest relative change in */
/* > any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */
/* > \par Internal Parameters: */
/* ========================= */
/* > */
/* > \verbatim */
/* > ITMAX is the maximum number of steps of iterative refinement. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date November 2011 */
/* > \ingroup complexOTHERcomputational */
/* ===================================================================== */
/* Subroutine */
void cpbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, complex *ab, integer *ldab,
             complex *afb, integer *ldafb, complex *b, integer *ldb, complex *x, integer *ldx,
             real *ferr, real *berr, complex *work, real *rwork, integer *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cpbrfs inputs: uplo %c, n %lld, kd %lld, nrhs %lld, ldab %lld, ldafb %lld, ldb %lld, "
             "ldx %lld",
             *uplo, *n, *kd, *nrhs, *ldab, *ldafb, *ldb, *ldx);
#else
    snprintf(buffer, 256,
             "cpbrfs inputs: uplo %c, n %d, kd %d, nrhs %d, ldab %d, ldafb %d, ldb %d, ldx %d",
             *uplo, *n, *kd, *nrhs, *ldab, *ldafb, *ldb, *ldx);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, x_dim1, x_offset, i__1,
        i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4;
    complex q__1;
    /* Builtin functions */
    double r_imag(complex *);
    /* Local variables */
    integer i__, j, k, l;
    real s, xk;
    integer nz;
    real eps;
    integer kase;
    real safe1, safe2;
    extern /* Subroutine */
        void
        chbmv_(char *, integer *, integer *, complex *, complex *, integer *, complex *, integer *,
               complex *, complex *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    extern /* Subroutine */
        void
        ccopy_(integer *, complex *, integer *, complex *, integer *),
        caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    integer count;
    logical upper;
    extern /* Subroutine */
        void
        clacn2_(integer *, complex *, complex *, real *, integer *, integer *);
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    extern /* Subroutine */
        void
        cpbtrs_(char *, integer *, integer *, integer *, complex *, integer *, complex *, integer *,
                integer *);
    real lstres;
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
    /* .. Local Arrays .. */
    /* .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Statement Functions .. */
    /* .. */
    /* .. Statement Function definitions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    afb_dim1 = *ldafb;
    afb_offset = 1 + afb_dim1;
    afb -= afb_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --rwork;
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
    else if(*kd < 0)
    {
        *info = -3;
    }
    else if(*nrhs < 0)
    {
        *info = -4;
    }
    else if(*ldab < *kd + 1)
    {
        *info = -6;
    }
    else if(*ldafb < *kd + 1)
    {
        *info = -8;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -10;
    }
    else if(*ldx < fla_max(1, *n))
    {
        *info = -12;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CPBRFS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            ferr[j] = 0.f;
            berr[j] = 0.f;
            /* L10: */
        }
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    /* NZ = maximum number of nonzero elements in each row of A, plus 1 */
    /* Computing MIN */
    i__1 = *n + 1;
    i__2 = (*kd << 1) + 2; // , expr subst
    nz = fla_min(i__1, i__2);
    eps = slamch_("Epsilon");
    safmin = slamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;
    /* Do for each right hand side */
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        count = 1;
        lstres = 3.f;
    L20: /* Loop until stopping criterion is satisfied. */
        /* Compute residual R = B - A * X */
        ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
        q__1.r = -1.f;
        q__1.i = -0.f; // , expr subst
        chbmv_(uplo, n, kd, &q__1, &ab[ab_offset], ldab, &x[j * x_dim1 + 1], &c__1, &c_b1, &work[1],
               &c__1);
        /* Compute componentwise relative backward error from formula */
        /* fla_max(i) ( f2c_abs(R(i)) / ( f2c_abs(A)*f2c_abs(X) + f2c_abs(B) )(i) ) */
        /* where f2c_abs(Z) is the componentwise absolute value of the matrix */
        /* or vector Z. If the i-th component of the denominator is less */
        /* than SAFE2, then SAFE1 is added to the i-th components of the */
        /* numerator and denominator before dividing. */
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * b_dim1;
            rwork[i__] = (r__1 = b[i__3].r, f2c_abs(r__1))
                         + (r__2 = r_imag(&b[i__ + j * b_dim1]), f2c_abs(r__2));
            /* L30: */
        }
        /* Compute f2c_abs(A)*f2c_abs(X) + f2c_abs(B). */
        if(upper)
        {
            i__2 = *n;
            for(k = 1; k <= i__2; ++k)
            {
                s = 0.f;
                i__3 = k + j * x_dim1;
                xk = (r__1 = x[i__3].r, f2c_abs(r__1))
                     + (r__2 = r_imag(&x[k + j * x_dim1]), f2c_abs(r__2));
                l = *kd + 1 - k;
                /* Computing MAX */
                i__3 = 1;
                i__4 = k - *kd; // , expr subst
                i__5 = k - 1;
                for(i__ = fla_max(i__3, i__4); i__ <= i__5; ++i__)
                {
                    i__3 = l + i__ + k * ab_dim1;
                    rwork[i__] += ((r__1 = ab[i__3].r, f2c_abs(r__1))
                                   + (r__2 = r_imag(&ab[l + i__ + k * ab_dim1]), f2c_abs(r__2)))
                                  * xk;
                    i__3 = l + i__ + k * ab_dim1;
                    i__4 = i__ + j * x_dim1;
                    s += ((r__1 = ab[i__3].r, f2c_abs(r__1))
                          + (r__2 = r_imag(&ab[l + i__ + k * ab_dim1]), f2c_abs(r__2)))
                         * ((r__3 = x[i__4].r, f2c_abs(r__3))
                            + (r__4 = r_imag(&x[i__ + j * x_dim1]), f2c_abs(r__4)));
                    /* L40: */
                }
                i__5 = *kd + 1 + k * ab_dim1;
                rwork[k] = rwork[k] + (r__1 = ab[i__5].r, f2c_abs(r__1)) * xk + s;
                /* L50: */
            }
        }
        else
        {
            i__2 = *n;
            for(k = 1; k <= i__2; ++k)
            {
                s = 0.f;
                i__5 = k + j * x_dim1;
                xk = (r__1 = x[i__5].r, f2c_abs(r__1))
                     + (r__2 = r_imag(&x[k + j * x_dim1]), f2c_abs(r__2));
                i__5 = k * ab_dim1 + 1;
                rwork[k] += (r__1 = ab[i__5].r, f2c_abs(r__1)) * xk;
                l = 1 - k;
                /* Computing MIN */
                i__3 = *n;
                i__4 = k + *kd; // , expr subst
                i__5 = fla_min(i__3, i__4);
                for(i__ = k + 1; i__ <= i__5; ++i__)
                {
                    i__3 = l + i__ + k * ab_dim1;
                    rwork[i__] += ((r__1 = ab[i__3].r, f2c_abs(r__1))
                                   + (r__2 = r_imag(&ab[l + i__ + k * ab_dim1]), f2c_abs(r__2)))
                                  * xk;
                    i__3 = l + i__ + k * ab_dim1;
                    i__4 = i__ + j * x_dim1;
                    s += ((r__1 = ab[i__3].r, f2c_abs(r__1))
                          + (r__2 = r_imag(&ab[l + i__ + k * ab_dim1]), f2c_abs(r__2)))
                         * ((r__3 = x[i__4].r, f2c_abs(r__3))
                            + (r__4 = r_imag(&x[i__ + j * x_dim1]), f2c_abs(r__4)));
                    /* L60: */
                }
                rwork[k] += s;
                /* L70: */
            }
        }
        s = 0.f;
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(rwork[i__] > safe2)
            {
                /* Computing MAX */
                i__5 = i__;
                r__3 = s;
                r__4 = ((r__1 = work[i__5].r, f2c_abs(r__1))
                        + (r__2 = r_imag(&work[i__]), f2c_abs(r__2)))
                       / rwork[i__]; // , expr subst
                s = fla_max(r__3, r__4);
            }
            else
            {
                /* Computing MAX */
                i__5 = i__;
                r__3 = s;
                r__4 = ((r__1 = work[i__5].r, f2c_abs(r__1))
                        + (r__2 = r_imag(&work[i__]), f2c_abs(r__2)) + safe1)
                       / (rwork[i__] + safe1); // , expr subst
                s = fla_max(r__3, r__4);
            }
            /* L80: */
        }
        berr[j] = s;
        /* Test stopping criterion. Continue iterating if */
        /* 1) The residual BERR(J) is larger than machine epsilon, and */
        /* 2) BERR(J) decreased by at least a factor of 2 during the */
        /* last iteration, and */
        /* 3) At most ITMAX iterations tried. */
        if(berr[j] > eps && berr[j] * 2.f <= lstres && count <= 5)
        {
            /* Update solution and try again. */
            cpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[1], n, info);
            caxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
            lstres = berr[j];
            ++count;
            goto L20;
        }
        /* Bound error from formula */
        /* norm(X - XTRUE) / norm(X) .le. FERR = */
        /* norm( f2c_abs(inv(A))* */
        /* ( f2c_abs(R) + NZ*EPS*( f2c_abs(A)*f2c_abs(X)+f2c_abs(B) ))) / norm(X) */
        /* where */
        /* norm(Z) is the magnitude of the largest component of Z */
        /* inv(A) is the inverse of A */
        /* f2c_abs(Z) is the componentwise absolute value of the matrix or */
        /* vector Z */
        /* NZ is the maximum number of nonzeros in any row of A, plus 1 */
        /* EPS is machine epsilon */
        /* The i-th component of f2c_abs(R)+NZ*EPS*(f2c_abs(A)*f2c_abs(X)+f2c_abs(B)) */
        /* is incremented by SAFE1 if the i-th component of */
        /* f2c_abs(A)*f2c_abs(X) + f2c_abs(B) is less than SAFE2. */
        /* Use CLACN2 to estimate the infinity-norm of the matrix */
        /* inv(A) * diag(W), */
        /* where W = f2c_abs(R) + NZ*EPS*( f2c_abs(A)*f2c_abs(X)+f2c_abs(B) ))) */
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(rwork[i__] > safe2)
            {
                i__5 = i__;
                rwork[i__] = (r__1 = work[i__5].r, f2c_abs(r__1))
                             + (r__2 = r_imag(&work[i__]), f2c_abs(r__2)) + nz * eps * rwork[i__];
            }
            else
            {
                i__5 = i__;
                rwork[i__] = (r__1 = work[i__5].r, f2c_abs(r__1))
                             + (r__2 = r_imag(&work[i__]), f2c_abs(r__2)) + nz * eps * rwork[i__]
                             + safe1;
            }
            /* L90: */
        }
        kase = 0;
    L100:
        clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
        if(kase != 0)
        {
            if(kase == 1)
            {
                /* Multiply by diag(W)*inv(A**H). */
                cpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[1], n, info);
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__5 = i__;
                    i__3 = i__;
                    i__4 = i__;
                    q__1.r = rwork[i__3] * work[i__4].r;
                    q__1.i = rwork[i__3] * work[i__4].i; // , expr subst
                    work[i__5].r = q__1.r;
                    work[i__5].i = q__1.i; // , expr subst
                    /* L110: */
                }
            }
            else if(kase == 2)
            {
                /* Multiply by inv(A)*diag(W). */
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__5 = i__;
                    i__3 = i__;
                    i__4 = i__;
                    q__1.r = rwork[i__3] * work[i__4].r;
                    q__1.i = rwork[i__3] * work[i__4].i; // , expr subst
                    work[i__5].r = q__1.r;
                    work[i__5].i = q__1.i; // , expr subst
                    /* L120: */
                }
                cpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[1], n, info);
            }
            goto L100;
        }
        /* Normalize error. */
        lstres = 0.f;
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            /* Computing MAX */
            i__5 = i__ + j * x_dim1;
            r__3 = lstres;
            r__4 = (r__1 = x[i__5].r, f2c_abs(r__1))
                   + (r__2 = r_imag(&x[i__ + j * x_dim1]), f2c_abs(r__2)); // , expr subst
            lstres = fla_max(r__3, r__4);
            /* L130: */
        }
        if(lstres != 0.f)
        {
            ferr[j] /= lstres;
        }
        /* L140: */
    }
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CPBRFS */
}
/* cpbrfs_ */
