/* ../netlib/dsgesv.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static doublereal c_b10 = -1.;
static doublereal c_b11 = 1.;
static aocl_int64_t c__1 = 1;
/* > \brief <b> DSGESV computes the solution to system of linear equations A * X = B for GE
 * matrices</b> (mixe d precision with iterative refinement) */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download DSGESV + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsgesv.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsgesv.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsgesv.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE DSGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK, */
/* SWORK, ITER, INFO ) */
/* .. Scalar Arguments .. */
/* INTEGER INFO, ITER, LDA, LDB, LDX, N, NRHS */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL SWORK( * ) */
/* DOUBLE PRECISION A( LDA, * ), B( LDB, * ), WORK( N, * ), */
/* $ X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSGESV computes the solution to a real system of linear equations */
/* > A * X = B, */
/* > where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */
/* > */
/* > DSGESV first attempts to factorize the matrix in SINGLE PRECISION */
/* > and use this factorization within an iterative refinement procedure */
/* > to produce a solution with DOUBLE PRECISION normwise backward error */
/* > quality (see below). If the approach fails the method switches to a */
/* > DOUBLE PRECISION factorization and solve. */
/* > */
/* > The iterative refinement is not going to be a winning strategy if */
/* > the ratio SINGLE PRECISION performance over DOUBLE PRECISION */
/* > performance is too small. A reasonable strategy should take the */
/* > number of right-hand sides and the size of the matrix into account. */
/* > This might be done with a call to ILAENV in the future. Up to now, we */
/* > always try iterative refinement. */
/* > */
/* > The iterative refinement process is stopped if */
/* > ITER > ITERMAX */
/* > or for all the RHS we have: */
/* > RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX */
/* > where */
/* > o ITER is the number of the current iteration in the iterative */
/* > refinement process */
/* > o RNRM is the infinity-norm of the residual */
/* > o XNRM is the infinity-norm of the solution */
/* > o ANRM is the infinity-operator-norm of the matrix A */
/* > o EPS is the machine epsilon returned by DLAMCH('Epsilon') */
/* > The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 */
/* > respectively. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrix B. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is DOUBLE PRECISION array, */
/* > dimension (LDA,N) */
/* > On entry, the N-by-N coefficient matrix A. */
/* > On exit, if iterative refinement has been successfully used */
/* > (INFO = 0 and ITER >= 0, see description below), then A is */
/* > unchanged, if double precision factorization has been used */
/* > (INFO = 0 and ITER < 0, see description below), then the */
/* > array A contains the factors L and U from the factorization */
/* > A = P*L*U;
the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > The pivot indices that define the permutation matrix P;
 */
/* > row i of the matrix was interchanged with row IPIV(i). */
/* > Corresponds either to the single precision factorization */
/* > (if INFO = 0 and ITER >= 0) or the double precision */
/* > factorization (if INFO = 0 and ITER < 0). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* > B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* > The N-by-NRHS right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB is INTEGER */
/* > The leading dimension of the array B. LDB >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* > X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* > If INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is DOUBLE PRECISION array, dimension (N,NRHS) */
/* > This array is used to hold the residual vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] SWORK */
/* > \verbatim */
/* > SWORK is REAL array, dimension (N*(N+NRHS)) */
/* > This array is used to use the single precision matrix and the */
/* > right-hand sides or solutions in single precision. */
/* > \endverbatim */
/* > */
/* > \param[out] ITER */
/* > \verbatim */
/* > ITER is INTEGER */
/* > < 0: iterative refinement has failed, double precision */
/* > factorization has been performed */
/* > -1 : the routine fell back to full precision for */
/* > implementation- or machine-specific reasons */
/* > -2 : narrowing the precision induced an overflow, */
/* > the routine fell back to full precision */
/* > -3 : failure of SGETRF */
/* > -31: stop the iterative refinement after the 30th */
/* > iterations */
/* > > 0: iterative refinement has been successfully used. */
/* > Returns the number of iterations */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, U(i,i) computed in DOUBLE PRECISION is */
/* > exactly zero. The factorization has been completed, */
/* > but the factor U is exactly singular, so the solution */
/* > could not be computed. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date June 2016 */
/* > \ingroup doubleGEsolve */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void dsgesv_(aocl_int_t *n, aocl_int_t *nrhs, doublereal *a, aocl_int_t *lda, aocl_int_t *ipiv,
             doublereal *b, aocl_int_t *ldb, doublereal *x, aocl_int_t *ldx, doublereal *work,
             real *swork, aocl_int_t *iter, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dsgesv(n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, iter, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t lda_64 = *lda;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t ldx_64 = *ldx;
    aocl_int64_t iter_64 = *iter;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dsgesv(&n_64, &nrhs_64, a, &lda_64, ipiv, b, &ldb_64, x, &ldx_64, work, swork,
                       &iter_64, &info_64);

    *iter = (aocl_int_t)iter_64;
    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_dsgesv(aocl_int64_t *n, aocl_int64_t *nrhs, doublereal *a, aocl_int64_t *lda,
                        aocl_int_t *ipiv, doublereal *b, aocl_int64_t *ldb, doublereal *x,
                        aocl_int64_t *ldx, doublereal *work, real *swork, aocl_int64_t *iter,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dsgesv inputs: n %" FLA_IS ", nrhs %" FLA_IS ", lda %" FLA_IS
                      ", ldb %" FLA_IS ", ldx %" FLA_IS "",
                      *n, *nrhs, *lda, *ldb, *ldx);
    /* System generated locals */
    aocl_int64_t a_dim1, a_offset, b_dim1, b_offset, work_dim1, work_offset, x_dim1, x_offset, i__1;
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    aocl_int64_t i__;
    doublereal cte, eps, anrm;
    aocl_int64_t ptsa;
    doublereal rnrm, xnrm;
    aocl_int64_t ptsx;
    aocl_int64_t iiter;
    /* -- LAPACK driver routine (version 3.8.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* June 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* .. Parameters .. */
    /* .. Local Scalars .. */
    /* .. External Subroutines .. */
    /* .. */
    /* .. External Functions .. */
    /* .. */
    /* .. Intrinsic Functions .. */
    /* .. */
    /* .. Executable Statements .. */
    /* Parameter adjustments */
    work_dim1 = *n;
    work_offset = 1 + work_dim1;
    work -= work_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --swork;
    /* Function Body */
    *info = 0;
    *iter = 0;
    /* Test the input parameters. */
    if(*n < 0)
    {
        *info = -1;
    }
    else if(*nrhs < 0)
    {
        *info = -2;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -4;
    }
    else if(*ldb < fla_max(1, *n))
    {
        *info = -7;
    }
    else if(*ldx < fla_max(1, *n))
    {
        *info = -9;
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("DSGESV", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Quick return if (N.EQ.0). */
    if(*n == 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* Skip single precision iterative refinement if a priori slower */
    /* than double precision factorization. */
    if(FALSE_)
    {
        *iter = -1;
        goto L40;
    }
    /* Compute some constants. */
    anrm = aocl_lapack_dlange("I", n, n, &a[a_offset], lda, &work[work_offset]);
    eps = dlamch_("Epsilon");
    cte = anrm * eps * sqrt((doublereal)(*n)) * 1.;
    /* Set the indices PTSA, PTSX for referencing SA and SX in SWORK. */
    ptsa = 1;
    ptsx = ptsa + *n * *n;
    /* Convert B from double precision to single precision and store the */
    /* result in SX. */
    aocl_lapack_dlag2s(n, nrhs, &b[b_offset], ldb, &swork[ptsx], n, info);
    if(*info != 0)
    {
        *iter = -2;
        goto L40;
    }
    /* Convert A from double precision to single precision and store the */
    /* result in SA. */
    aocl_lapack_dlag2s(n, n, &a[a_offset], lda, &swork[ptsa], n, info);
    if(*info != 0)
    {
        *iter = -2;
        goto L40;
    }
    /* Compute the LU factorization of SA. */
    aocl_lapack_sgetrf(n, n, &swork[ptsa], n, &ipiv[1], info);
    if(*info != 0)
    {
        *iter = -3;
        goto L40;
    }
    /* Solve the system SA*SX = SB. */
    aocl_lapack_sgetrs("No transpose", n, nrhs, &swork[ptsa], n, &ipiv[1], &swork[ptsx], n, info);
    /* Convert SX back to double precision */
    aocl_lapack_slag2d(n, nrhs, &swork[ptsx], n, &x[x_offset], ldx, info);
    /* Compute R = B - AX (R is WORK). */
    aocl_lapack_dlacpy("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n);
    aocl_blas_dgemm("No Transpose", "No Transpose", n, nrhs, n, &c_b10, &a[a_offset], lda,
                    &x[x_offset], ldx, &c_b11, &work[work_offset], n);
    /* Check whether the NRHS normwise backward errors satisfy the */
    /* stopping criterion. If yes, set ITER=0 and return. */
    i__1 = *nrhs;
    for(i__ = 1; i__ <= i__1; ++i__)
    {
        xnrm = (d__1 = x[aocl_blas_idamax(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1],
                f2c_abs(d__1));
        rnrm = (d__1
                = work[aocl_blas_idamax(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * work_dim1],
                f2c_abs(d__1));
        if(rnrm > xnrm * cte)
        {
            goto L10;
        }
    }
    /* If we are here, the NRHS normwise backward errors satisfy the */
    /* stopping criterion. We are good to exit. */
    *iter = 0;
    AOCL_DTL_TRACE_LOG_EXIT
    return;
L10:
    for(iiter = 1; iiter <= 30; ++iiter)
    {
        /* Convert R (in WORK) from double precision to single precision */
        /* and store the result in SX. */
        aocl_lapack_dlag2s(n, nrhs, &work[work_offset], n, &swork[ptsx], n, info);
        if(*info != 0)
        {
            *iter = -2;
            goto L40;
        }
        /* Solve the system SA*SX = SR. */
        aocl_lapack_sgetrs("No transpose", n, nrhs, &swork[ptsa], n, &ipiv[1], &swork[ptsx], n,
                           info);
        /* Convert SX back to double precision and update the current */
        /* iterate. */
        aocl_lapack_slag2d(n, nrhs, &swork[ptsx], n, &work[work_offset], n, info);
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            aocl_blas_daxpy(n, &c_b11, &work[i__ * work_dim1 + 1], &c__1, &x[i__ * x_dim1 + 1],
                            &c__1);
        }
        /* Compute R = B - AX (R is WORK). */
        aocl_lapack_dlacpy("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n);
        aocl_blas_dgemm("No Transpose", "No Transpose", n, nrhs, n, &c_b10, &a[a_offset], lda,
                        &x[x_offset], ldx, &c_b11, &work[work_offset], n);
        /* Check whether the NRHS normwise backward errors satisfy the */
        /* stopping criterion. If yes, set ITER=IITER>0 and return. */
        i__1 = *nrhs;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            xnrm = (d__1 = x[aocl_blas_idamax(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1],
                    f2c_abs(d__1));
            rnrm
                = (d__1
                   = work[aocl_blas_idamax(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * work_dim1],
                   f2c_abs(d__1));
            if(rnrm > xnrm * cte)
            {
                goto L20;
            }
        }
        /* If we are here, the NRHS normwise backward errors satisfy the */
        /* stopping criterion, we are good to exit. */
        *iter = iiter;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    L20: /* L30: */
        ;
    }
    /* If we are at this place of the code, this is because we have */
    /* performed ITER=ITERMAX iterations and never satisfied the */
    /* stopping criterion, set up the ITER flag accordingly and follow up */
    /* on double precision routine. */
    *iter = -31;
L40: /* Single-precision iterative refinement failed to converge to a */
    /* satisfactory solution, so we resort to double precision. */
    aocl_lapack_dgetrf(n, n, &a[a_offset], lda, &ipiv[1], info);
    if(*info != 0)
    {
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    aocl_lapack_dlacpy("All", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    aocl_lapack_dgetrs("No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &x[x_offset], ldx,
                       info);
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of DSGESV. */
}
/* dsgesv_ */
