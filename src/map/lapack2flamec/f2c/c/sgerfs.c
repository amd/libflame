/* ../netlib/sgerfs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
static real c_b15 = -1.f;
static real c_b17 = 1.f;
/* > \brief \b SGERFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SGERFS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgerfs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgerfs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgerfs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/* X, LDX, FERR, BERR, WORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER TRANS */
/* INTEGER INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), IWORK( * ) */
/* REAL A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/* $ BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGERFS improves the computed solution to a system of linear */
/* > equations and provides error bounds and backward error estimates for */
/* > the solution. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The order of the matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The original N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* > AF is REAL array, dimension (LDAF,N) */
/* > The factors L and U from the factorization A = P*L*U */
/* > as computed by SGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices from SGETRF;
for 1<=i<=N, row i of the */
/* > matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,NRHS) */
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
/* > X is REAL array, dimension (LDX,NRHS) */
/* > On entry, the solution matrix X, as computed by SGETRS. */
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
/* > WORK is REAL array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK is INTEGER array, dimension (N) */
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
/* > \ingroup realGEcomputational */
/* ===================================================================== */
/* Subroutine */
void sgerfs_(char *trans, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf,
             integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *ferr, real *berr,
             real *work, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgerfs inputs: trans %c ,n %" FLA_IS ",nrhs %" FLA_IS ",lda %" FLA_IS
                      ",ldaf %" FLA_IS ",ldb %" FLA_IS ",ldx %" FLA_IS "",
                      *trans, *n, *nrhs, *lda, *ldaf, *ldb, *ldx);
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2,
        i__3;
    real r__1, r__2, r__3;
    /* Local variables */
    integer i__, j, k;
    real s, xk;
    integer nz;
    real eps;
    integer kase;
    real safe1, safe2;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    extern /* Subroutine */
        void
        sgemv_(char *, integer *, integer *, real *, real *, integer *, real *, integer *, real *,
               real *, integer *);
    integer count;
    extern /* Subroutine */
        void
        scopy_(integer *, real *, integer *, real *, integer *),
        saxpy_(integer *, real *, real *, integer *, real *, integer *),
        slacn2_(integer *, real *, real *, integer *, real *, integer *, integer *);
    extern real slamch_(char *);
    real safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical notran;
    extern /* Subroutine */
        void
        sgetrs_(char *, integer *, integer *, real *, integer *, integer *, real *, integer *,
                integer *);
    char transt[1];
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
    /* .. Executable Statements .. */
    /* Test the input parameters. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", 1, 1);
    if(!notran && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
    {
        *info = -1;
    }
    else if(*n < 0)
    {
        *info = -2;
    }
    else if(*nrhs < 0)
    {
        *info = -3;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -5;
    }
    else if(*ldaf < fla_max(1, *n))
    {
        *info = -7;
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
        xerbla_("SGERFS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
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
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(notran)
    {
        *(unsigned char *)transt = 'T';
    }
    else
    {
        *(unsigned char *)transt = 'N';
    }
    /* NZ = maximum number of nonzero elements in each row of A, plus 1 */
    nz = *n + 1;
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
        /* Compute residual R = B - op(A) * X, */
        /* where op(A) = A, A**T, or A**H, depending on TRANS. */
        scopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
        sgemv_(trans, n, n, &c_b15, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, &c_b17,
               &work[*n + 1], &c__1);
        /* Compute componentwise relative backward error from formula */
        /* fla_max(i) ( f2c_abs(R(i)) / ( f2c_abs(op(A))*f2c_abs(X) + f2c_abs(B) )(i) ) */
        /* where f2c_abs(Z) is the componentwise absolute value of the matrix */
        /* or vector Z. If the i-th component of the denominator is less */
        /* than SAFE2, then SAFE1 is added to the i-th components of the */
        /* numerator and denominator before dividing. */
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            work[i__] = (r__1 = b[i__ + j * b_dim1], f2c_abs(r__1));
            /* L30: */
        }
        /* Compute f2c_abs(op(A))*f2c_abs(X) + f2c_abs(B). */
        if(notran)
        {
            i__2 = *n;
            for(k = 1; k <= i__2; ++k)
            {
                xk = (r__1 = x[k + j * x_dim1], f2c_abs(r__1));
                i__3 = *n;
                for(i__ = 1; i__ <= i__3; ++i__)
                {
                    work[i__] += (r__1 = a[i__ + k * a_dim1], f2c_abs(r__1)) * xk;
                    /* L40: */
                }
                /* L50: */
            }
        }
        else
        {
            i__2 = *n;
            for(k = 1; k <= i__2; ++k)
            {
                s = 0.f;
                i__3 = *n;
                for(i__ = 1; i__ <= i__3; ++i__)
                {
                    s += (r__1 = a[i__ + k * a_dim1], f2c_abs(r__1))
                         * (r__2 = x[i__ + j * x_dim1], f2c_abs(r__2));
                    /* L60: */
                }
                work[k] += s;
                /* L70: */
            }
        }
        s = 0.f;
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(work[i__] > safe2)
            {
                /* Computing MAX */
                r__2 = s;
                r__3 = (r__1 = work[*n + i__], f2c_abs(r__1)) / work[i__]; // , expr subst
                s = fla_max(r__2, r__3);
            }
            else
            {
                /* Computing MAX */
                r__2 = s;
                r__3 = ((r__1 = work[*n + i__], f2c_abs(r__1)) + safe1)
                       / (work[i__] + safe1); // , expr subst
                s = fla_max(r__2, r__3);
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
            sgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n + 1], n, info);
            saxpy_(n, &c_b17, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
            lstres = berr[j];
            ++count;
            goto L20;
        }
        /* Bound error from formula */
        /* norm(X - XTRUE) / norm(X) .le. FERR = */
        /* norm( f2c_abs(inv(op(A)))* */
        /* ( f2c_abs(R) + NZ*EPS*( f2c_abs(op(A))*f2c_abs(X)+f2c_abs(B) ))) / norm(X) */
        /* where */
        /* norm(Z) is the magnitude of the largest component of Z */
        /* inv(op(A)) is the inverse of op(A) */
        /* f2c_abs(Z) is the componentwise absolute value of the matrix or */
        /* vector Z */
        /* NZ is the maximum number of nonzeros in any row of A, plus 1 */
        /* EPS is machine epsilon */
        /* The i-th component of f2c_abs(R)+NZ*EPS*(f2c_abs(op(A))*f2c_abs(X)+f2c_abs(B)) */
        /* is incremented by SAFE1 if the i-th component of */
        /* f2c_abs(op(A))*f2c_abs(X) + f2c_abs(B) is less than SAFE2. */
        /* Use SLACN2 to estimate the infinity-norm of the matrix */
        /* inv(op(A)) * diag(W), */
        /* where W = f2c_abs(R) + NZ*EPS*( f2c_abs(op(A))*f2c_abs(X)+f2c_abs(B) ))) */
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(work[i__] > safe2)
            {
                work[i__] = (r__1 = work[*n + i__], f2c_abs(r__1)) + nz * eps * work[i__];
            }
            else
            {
                work[i__] = (r__1 = work[*n + i__], f2c_abs(r__1)) + nz * eps * work[i__] + safe1;
            }
            /* L90: */
        }
        kase = 0;
    L100:
        slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &kase, isave);
        if(kase != 0)
        {
            if(kase == 1)
            {
                /* Multiply by diag(W)*inv(op(A)**T). */
                sgetrs_(transt, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n + 1], n, info);
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    work[*n + i__] = work[i__] * work[*n + i__];
                    /* L110: */
                }
            }
            else
            {
                /* Multiply by inv(op(A))*diag(W). */
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    work[*n + i__] = work[i__] * work[*n + i__];
                    /* L120: */
                }
                sgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n + 1], n, info);
            }
            goto L100;
        }
        /* Normalize error. */
        lstres = 0.f;
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            /* Computing MAX */
            r__2 = lstres;
            r__3 = (r__1 = x[i__ + j * x_dim1], f2c_abs(r__1)); // , expr subst
            lstres = fla_max(r__2, r__3);
            /* L130: */
        }
        if(lstres != 0.f)
        {
            ferr[j] /= lstres;
        }
        /* L140: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SGERFS */
}
/* sgerfs_ */
