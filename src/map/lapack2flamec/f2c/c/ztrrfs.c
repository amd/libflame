/* ../netlib/ztrrfs.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c__1 = 1;
/* > \brief \b ZTRRFS */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download ZTRRFS + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrrfs.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrrfs.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrrfs.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, */
/* LDX, FERR, BERR, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER DIAG, TRANS, UPLO */
/* INTEGER INFO, LDA, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* DOUBLE PRECISION BERR( * ), FERR( * ), RWORK( * ) */
/* COMPLEX*16 A( LDA, * ), B( LDB, * ), WORK( * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTRRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by ZTRTRS or some other */
/* > means before entering this routine. ZTRRFS does not do iterative */
/* > refinement because doing so cannot improve the backward error. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] UPLO */
/* > \verbatim */
/* > UPLO is CHARACTER*1 */
/* > = 'U': A is upper triangular;
 */
/* > = 'L': A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations: */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* > DIAG is CHARACTER*1 */
/* > = 'N': A is non-unit triangular;
 */
/* > = 'U': A is unit triangular. */
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
/* > A is COMPLEX*16 array, dimension (LDA,N) */
/* > The triangular matrix A. If UPLO = 'U', the leading N-by-N */
/* > upper triangular part of the array A contains the upper */
/* > triangular matrix, and the strictly lower triangular part of */
/* > A is not referenced. If UPLO = 'L', the leading N-by-N lower */
/* > triangular part of the array A contains the lower triangular */
/* > matrix, and the strictly upper triangular part of A is not */
/* > referenced. If DIAG = 'U', the diagonal elements of A are */
/* > also not referenced and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* > The right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* > X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* > The solution matrix X. */
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
/* > FERR is DOUBLE PRECISION array, dimension (NRHS) */
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
/* > BERR is DOUBLE PRECISION array, dimension (NRHS) */
/* > The componentwise relative backward error of each solution */
/* > vector X(j) (i.e., the smallest relative change in */
/* > any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
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
void ztrrfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplex *a,
             integer *lda, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
             doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork,
             integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ztrrfs inputs: uplo %c, trans %c, diag %c, n %" FLA_IS ", nrhs %" FLA_IS
                      ", lda %" FLA_IS ", ldb %" FLA_IS ", ldx %" FLA_IS "",
                      *uplo, *trans, *diag, *n, *nrhs, *lda, *ldb, *ldx);
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    /* Local variables */
    integer i__, j, k;
    doublereal s, xk;
    integer nz;
    doublereal eps;
    integer kase;
    doublereal safe1, safe2;
    extern logical lsame_(char *, char *, integer, integer);
    integer isave[3];
    logical upper;
    extern /* Subroutine */
        void
        zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, integer *),
        zaxpy_(integer *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, integer *),
        ztrmv_(char *, char *, char *, integer *, doublecomplex *, integer *, doublecomplex *,
               integer *),
        ztrsv_(char *, char *, char *, integer *, doublecomplex *, integer *, doublecomplex *,
               integer *),
        zlacn2_(integer *, doublecomplex *, doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *);
    doublereal safmin;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    logical notran;
    char transn[1], transt[1];
    logical nounit;
    doublereal lstres;
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
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
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
    notran = lsame_(trans, "N", 1, 1);
    nounit = lsame_(diag, "N", 1, 1);
    if(!upper && !lsame_(uplo, "L", 1, 1))
    {
        *info = -1;
    }
    else if(!notran && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
    {
        *info = -2;
    }
    else if(!nounit && !lsame_(diag, "U", 1, 1))
    {
        *info = -3;
    }
    else if(*n < 0)
    {
        *info = -4;
    }
    else if(*nrhs < 0)
    {
        *info = -5;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -9;
    }
    else if(*ldx < fla_max(1, *n))
    {
        *info = -11;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("ZTRRFS", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if possible */
    if(*n == 0 || *nrhs == 0)
    {
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            ferr[j] = 0.;
            berr[j] = 0.;
            /* L10: */
        }
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(notran)
    {
        *(unsigned char *)transn = 'N';
        *(unsigned char *)transt = 'C';
    }
    else
    {
        *(unsigned char *)transn = 'C';
        *(unsigned char *)transt = 'N';
    }
    /* NZ = maximum number of nonzero elements in each row of A, plus 1 */
    nz = *n + 1;
    eps = dlamch_("Epsilon");
    safmin = dlamch_("Safe minimum");
    safe1 = nz * safmin;
    safe2 = safe1 / eps;
    /* Do for each right hand side */
    i__1 = *nrhs;
    for(j = 1; j <= i__1; ++j)
    {
        /* Compute residual R = B - op(A) * X, */
        /* where op(A) = A, A**T, or A**H, depending on TRANS. */
        zcopy_(n, &x[j * x_dim1 + 1], &c__1, &work[1], &c__1);
        ztrmv_(uplo, trans, diag, n, &a[a_offset], lda, &work[1], &c__1);
        z__1.r = -1.;
        z__1.i = -0.; // , expr subst
        zaxpy_(n, &z__1, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
        /* Compute componentwise relative backward error from formula */
        /* fla_max(i) ( f2c_dabs(R(i)) / ( f2c_dabs(op(A))*f2c_dabs(X) + f2c_dabs(B) )(i) ) */
        /* where f2c_dabs(Z) is the componentwise absolute value of the matrix */
        /* or vector Z. If the i-th component of the denominator is less */
        /* than SAFE2, then SAFE1 is added to the i-th components of the */
        /* numerator and denominator before dividing. */
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            i__3 = i__ + j * b_dim1;
            rwork[i__] = (d__1 = b[i__3].r, f2c_dabs(d__1))
                         + (d__2 = d_imag(&b[i__ + j * b_dim1]), f2c_dabs(d__2));
            /* L20: */
        }
        if(notran)
        {
            /* Compute f2c_dabs(A)*f2c_dabs(X) + f2c_dabs(B). */
            if(upper)
            {
                if(nounit)
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        i__3 = k + j * x_dim1;
                        xk = (d__1 = x[i__3].r, f2c_dabs(d__1))
                             + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_dabs(d__2));
                        i__3 = k;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            rwork[i__] += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                           + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                          * xk;
                            /* L30: */
                        }
                        /* L40: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        i__3 = k + j * x_dim1;
                        xk = (d__1 = x[i__3].r, f2c_dabs(d__1))
                             + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_dabs(d__2));
                        i__3 = k - 1;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            rwork[i__] += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                           + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                          * xk;
                            /* L50: */
                        }
                        rwork[k] += xk;
                        /* L60: */
                    }
                }
            }
            else
            {
                if(nounit)
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        i__3 = k + j * x_dim1;
                        xk = (d__1 = x[i__3].r, f2c_dabs(d__1))
                             + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_dabs(d__2));
                        i__3 = *n;
                        for(i__ = k; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            rwork[i__] += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                           + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                          * xk;
                            /* L70: */
                        }
                        /* L80: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        i__3 = k + j * x_dim1;
                        xk = (d__1 = x[i__3].r, f2c_dabs(d__1))
                             + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_dabs(d__2));
                        i__3 = *n;
                        for(i__ = k + 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            rwork[i__] += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                           + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                          * xk;
                            /* L90: */
                        }
                        rwork[k] += xk;
                        /* L100: */
                    }
                }
            }
        }
        else
        {
            /* Compute f2c_dabs(A**H)*f2c_dabs(X) + f2c_dabs(B). */
            if(upper)
            {
                if(nounit)
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        s = 0.;
                        i__3 = k;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            i__5 = i__ + j * x_dim1;
                            s += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                  + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                 * ((d__3 = x[i__5].r, f2c_dabs(d__3))
                                    + (d__4 = d_imag(&x[i__ + j * x_dim1]), f2c_dabs(d__4)));
                            /* L110: */
                        }
                        rwork[k] += s;
                        /* L120: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        i__3 = k + j * x_dim1;
                        s = (d__1 = x[i__3].r, f2c_dabs(d__1))
                            + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_dabs(d__2));
                        i__3 = k - 1;
                        for(i__ = 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            i__5 = i__ + j * x_dim1;
                            s += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                  + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                 * ((d__3 = x[i__5].r, f2c_dabs(d__3))
                                    + (d__4 = d_imag(&x[i__ + j * x_dim1]), f2c_dabs(d__4)));
                            /* L130: */
                        }
                        rwork[k] += s;
                        /* L140: */
                    }
                }
            }
            else
            {
                if(nounit)
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        s = 0.;
                        i__3 = *n;
                        for(i__ = k; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            i__5 = i__ + j * x_dim1;
                            s += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                  + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                 * ((d__3 = x[i__5].r, f2c_dabs(d__3))
                                    + (d__4 = d_imag(&x[i__ + j * x_dim1]), f2c_dabs(d__4)));
                            /* L150: */
                        }
                        rwork[k] += s;
                        /* L160: */
                    }
                }
                else
                {
                    i__2 = *n;
                    for(k = 1; k <= i__2; ++k)
                    {
                        i__3 = k + j * x_dim1;
                        s = (d__1 = x[i__3].r, f2c_dabs(d__1))
                            + (d__2 = d_imag(&x[k + j * x_dim1]), f2c_dabs(d__2));
                        i__3 = *n;
                        for(i__ = k + 1; i__ <= i__3; ++i__)
                        {
                            i__4 = i__ + k * a_dim1;
                            i__5 = i__ + j * x_dim1;
                            s += ((d__1 = a[i__4].r, f2c_dabs(d__1))
                                  + (d__2 = d_imag(&a[i__ + k * a_dim1]), f2c_dabs(d__2)))
                                 * ((d__3 = x[i__5].r, f2c_dabs(d__3))
                                    + (d__4 = d_imag(&x[i__ + j * x_dim1]), f2c_dabs(d__4)));
                            /* L170: */
                        }
                        rwork[k] += s;
                        /* L180: */
                    }
                }
            }
        }
        s = 0.;
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(rwork[i__] > safe2)
            {
                /* Computing MAX */
                i__3 = i__;
                d__3 = s;
                d__4 = ((d__1 = work[i__3].r, f2c_dabs(d__1))
                        + (d__2 = d_imag(&work[i__]), f2c_dabs(d__2)))
                       / rwork[i__]; // , expr subst
                s = fla_max(d__3, d__4);
            }
            else
            {
                /* Computing MAX */
                i__3 = i__;
                d__3 = s;
                d__4 = ((d__1 = work[i__3].r, f2c_dabs(d__1))
                        + (d__2 = d_imag(&work[i__]), f2c_dabs(d__2)) + safe1)
                       / (rwork[i__] + safe1); // , expr subst
                s = fla_max(d__3, d__4);
            }
            /* L190: */
        }
        berr[j] = s;
        /* Bound error from formula */
        /* norm(X - XTRUE) / norm(X) .le. FERR = */
        /* norm( f2c_dabs(inv(op(A)))* */
        /* ( f2c_dabs(R) + NZ*EPS*( f2c_dabs(op(A))*f2c_dabs(X)+f2c_dabs(B) ))) / norm(X) */
        /* where */
        /* norm(Z) is the magnitude of the largest component of Z */
        /* inv(op(A)) is the inverse of op(A) */
        /* f2c_dabs(Z) is the componentwise absolute value of the matrix or */
        /* vector Z */
        /* NZ is the maximum number of nonzeros in any row of A, plus 1 */
        /* EPS is machine epsilon */
        /* The i-th component of f2c_dabs(R)+NZ*EPS*(f2c_dabs(op(A))*f2c_dabs(X)+f2c_dabs(B)) */
        /* is incremented by SAFE1 if the i-th component of */
        /* f2c_dabs(op(A))*f2c_dabs(X) + f2c_dabs(B) is less than SAFE2. */
        /* Use ZLACN2 to estimate the infinity-norm of the matrix */
        /* inv(op(A)) * diag(W), */
        /* where W = f2c_dabs(R) + NZ*EPS*( f2c_dabs(op(A))*f2c_dabs(X)+f2c_dabs(B) ))) */
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            if(rwork[i__] > safe2)
            {
                i__3 = i__;
                rwork[i__] = (d__1 = work[i__3].r, f2c_dabs(d__1))
                             + (d__2 = d_imag(&work[i__]), f2c_dabs(d__2)) + nz * eps * rwork[i__];
            }
            else
            {
                i__3 = i__;
                rwork[i__] = (d__1 = work[i__3].r, f2c_dabs(d__1))
                             + (d__2 = d_imag(&work[i__]), f2c_dabs(d__2)) + nz * eps * rwork[i__]
                             + safe1;
            }
            /* L200: */
        }
        kase = 0;
    L210:
        zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
        if(kase != 0)
        {
            if(kase == 1)
            {
                /* Multiply by diag(W)*inv(op(A)**H). */
                ztrsv_(uplo, transt, diag, n, &a[a_offset], lda, &work[1], &c__1);
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__;
                    z__1.r = rwork[i__4] * work[i__5].r;
                    z__1.i = rwork[i__4] * work[i__5].i; // , expr subst
                    work[i__3].r = z__1.r;
                    work[i__3].i = z__1.i; // , expr subst
                    /* L220: */
                }
            }
            else
            {
                /* Multiply by inv(op(A))*diag(W). */
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__;
                    i__4 = i__;
                    i__5 = i__;
                    z__1.r = rwork[i__4] * work[i__5].r;
                    z__1.i = rwork[i__4] * work[i__5].i; // , expr subst
                    work[i__3].r = z__1.r;
                    work[i__3].i = z__1.i; // , expr subst
                    /* L230: */
                }
                ztrsv_(uplo, transn, diag, n, &a[a_offset], lda, &work[1], &c__1);
            }
            goto L210;
        }
        /* Normalize error. */
        lstres = 0.;
        i__2 = *n;
        for(i__ = 1; i__ <= i__2; ++i__)
        {
            /* Computing MAX */
            i__3 = i__ + j * x_dim1;
            d__3 = lstres;
            d__4 = (d__1 = x[i__3].r, f2c_dabs(d__1))
                   + (d__2 = d_imag(&x[i__ + j * x_dim1]), f2c_dabs(d__2)); // , expr subst
            lstres = fla_max(d__3, d__4);
            /* L240: */
        }
        if(lstres != 0.)
        {
            ferr[j] /= lstres;
        }
        /* L250: */
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of ZTRRFS */
}
/* ztrrfs_ */
