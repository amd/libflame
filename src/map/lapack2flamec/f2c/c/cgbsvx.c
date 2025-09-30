/* ../netlib/cgbsvx.f -- translated by f2c (version 20100827). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static aocl_int64_t c__1 = 1;
/* > \brief <b> CGBSVX computes the solution to system of linear equations A * X = B for GB
 * matrices</b> */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download CGBSVX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbsvx.
 * f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbsvx.
 * f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbsvx.
 * f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE CGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, */
/* LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, */
/* RCOND, FERR, BERR, WORK, RWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, FACT, TRANS */
/* INTEGER INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS */
/* REAL RCOND */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ) */
/* REAL BERR( * ), C( * ), FERR( * ), R( * ), */
/* $ RWORK( * ) */
/* COMPLEX AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/* $ WORK( * ), X( LDX, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBSVX uses the LU factorization to compute the solution to a scomplex */
/* > system of linear equations A * X = B, A**T * X = B, or A**H * X = B, */
/* > where A is a band matrix of order N with KL subdiagonals and KU */
/* > superdiagonals, and X and B are N-by-NRHS matrices. */
/* > */
/* > Error bounds on the solution and a condition estimate are also */
/* > provided. */
/* > \endverbatim */
/* > \par Description: */
/* ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed by this subroutine: */
/* > */
/* > 1. If FACT = 'E', real scaling factors are computed to equilibrate */
/* > the system: */
/* > TRANS = 'N': diag(R)*A*diag(C) *inv(diag(C))*X = diag(R)*B */
/* > TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B */
/* > TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B */
/* > Whether or not the system will be equilibrated depends on the */
/* > scaling of the matrix A, but if equilibration is used, A is */
/* > overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N') */
/* > or diag(C)*B (if TRANS = 'T' or 'C'). */
/* > */
/* > 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the */
/* > matrix A (after equilibration if FACT = 'E') as */
/* > A = L * U, */
/* > where L is a product of permutation and unit lower triangular */
/* > matrices with KL subdiagonals, and U is upper triangular with */
/* > KL+KU superdiagonals. */
/* > */
/* > 3. If some U(i,i)=0, so that U is exactly singular, then the routine */
/* > returns with INFO = i. Otherwise, the factored form of A is used */
/* > to estimate the condition number of the matrix A. If the */
/* > reciprocal of the condition number is less than machine precision, */
/* > INFO = N+1 is returned as a warning, but the routine still goes on */
/* > to solve for X and compute error bounds as described below. */
/* > */
/* > 4. The system of equations is solved for X using the factored form */
/* > of A. */
/* > */
/* > 5. Iterative refinement is applied to improve the computed solution */
/* > matrix and calculate error bounds and backward error estimates */
/* > for it. */
/* > */
/* > 6. If equilibration was used, the matrix X is premultiplied by */
/* > diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so */
/* > that it solves the original system before equilibration. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \param[in] FACT */
/* > \verbatim */
/* > FACT is CHARACTER*1 */
/* > Specifies whether or not the factored form of the matrix A is */
/* > supplied on entry, and if not, whether the matrix A should be */
/* > equilibrated before it is factored. */
/* > = 'F': On entry, AFB and IPIV contain the factored form of */
/* > A. If EQUED is not 'N', the matrix A has been */
/* > equilibrated with scaling factors given by R and C. */
/* > AB, AFB, and IPIV are not modified. */
/* > = 'N': The matrix A will be copied to AFB and factored. */
/* > = 'E': The matrix A will be equilibrated if necessary, then */
/* > copied to AFB and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* > TRANS is CHARACTER*1 */
/* > Specifies the form of the system of equations. */
/* > = 'N': A * X = B (No transpose) */
/* > = 'T': A**T * X = B (Transpose) */
/* > = 'C': A**H * X = B (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* > N is INTEGER */
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KL */
/* > \verbatim */
/* > KL is INTEGER */
/* > The number of subdiagonals within the band of A. KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* > KU is INTEGER */
/* > The number of superdiagonals within the band of A. KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* > AB is COMPLEX array, dimension (LDAB,N) */
/* > On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* > The j-th column of A is stored in the j-th column of the */
/* > array AB as follows: */
/* > AB(KU+1+i-j,j) = A(i,j) for fla_max(1,j-KU)<=i<=fla_min(N,j+kl) */
/* > */
/* > If FACT = 'F' and EQUED is not 'N', then A must have been */
/* > equilibrated by the scaling factors in R and/or C. AB is not */
/* > modified if FACT = 'F' or 'N', or if FACT = 'E' and */
/* > EQUED = 'N' on exit. */
/* > */
/* > On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* > EQUED = 'R': A := diag(R) * A */
/* > EQUED = 'C': A := A * diag(C) */
/* > EQUED = 'B': A := diag(R) * A * diag(C). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* > LDAB is INTEGER */
/* > The leading dimension of the array AB. LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFB */
/* > \verbatim */
/* > AFB is COMPLEX array, dimension (LDAFB,N) */
/* > If FACT = 'F', then AFB is an input argument and on entry */
/* > contains details of the LU factorization of the band matrix */
/* > A, as computed by CGBTRF. U is stored as an upper triangular */
/* > band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* > and the multipliers used during the factorization are stored */
/* > in rows KL+KU+2 to 2*KL+KU+1. If EQUED .ne. 'N', then AFB is */
/* > the factored form of the equilibrated matrix A. */
/* > */
/* > If FACT = 'N', then AFB is an output argument and on exit */
/* > returns details of the LU factorization of A. */
/* > */
/* > If FACT = 'E', then AFB is an output argument and on exit */
/* > returns details of the LU factorization of the equilibrated */
/* > matrix A (see the description of AB for the form of the */
/* > equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* > LDAFB is INTEGER */
/* > The leading dimension of the array AFB. LDAFB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > If FACT = 'F', then IPIV is an input argument and on entry */
/* > contains the pivot indices from the factorization A = L*U */
/* > as computed by CGBTRF;
row i of the matrix was interchanged */
/* > with row IPIV(i). */
/* > */
/* > If FACT = 'N', then IPIV is an output argument and on exit */
/* > contains the pivot indices from the factorization A = L*U */
/* > of the original matrix A. */
/* > */
/* > If FACT = 'E', then IPIV is an output argument and on exit */
/* > contains the pivot indices from the factorization A = L*U */
/* > of the equilibrated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* > EQUED is CHARACTER*1 */
/* > Specifies the form of equilibration that was done. */
/* > = 'N': No equilibration (always true if FACT = 'N'). */
/* > = 'R': Row equilibration, i.e., A has been premultiplied by */
/* > diag(R). */
/* > = 'C': Column equilibration, i.e., A has been postmultiplied */
/* > by diag(C). */
/* > = 'B': Both row and column equilibration, i.e., A has been */
/* > replaced by diag(R) * A * diag(C). */
/* > EQUED is an input argument if FACT = 'F';
otherwise, it is an */
/* > output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] R */
/* > \verbatim */
/* > R is REAL array, dimension (N) */
/* > The row scale factors for A. If EQUED = 'R' or 'B', A is */
/* > multiplied on the left by diag(R);
if EQUED = 'N' or 'C', R */
/* > is not accessed. R is an input argument if FACT = 'F';
 */
/* > otherwise, R is an output argument. If FACT = 'F' and */
/* > EQUED = 'R' or 'B', each element of R must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* > C is REAL array, dimension (N) */
/* > The column scale factors for A. If EQUED = 'C' or 'B', A is */
/* > multiplied on the right by diag(C);
if EQUED = 'N' or 'R', C */
/* > is not accessed. C is an input argument if FACT = 'F';
 */
/* > otherwise, C is an output argument. If FACT = 'F' and */
/* > EQUED = 'C' or 'B', each element of C must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is COMPLEX array, dimension (LDB,NRHS) */
/* > On entry, the right hand side matrix B. */
/* > On exit, */
/* > if EQUED = 'N', B is not modified;
 */
/* > if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by */
/* > diag(R)*B;
 */
/* > if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is */
/* > overwritten by diag(C)*B. */
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
/* > X is COMPLEX array, dimension (LDX,NRHS) */
/* > If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X */
/* > to the original system of equations. Note that A and B are */
/* > modified on exit if EQUED .ne. 'N', and the solution to the */
/* > equilibrated system is inv(diag(C))*X if TRANS = 'N' and */
/* > EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C' */
/* > and EQUED = 'R' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX is INTEGER */
/* > The leading dimension of the array X. LDX >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* > RCOND is REAL */
/* > The estimate of the reciprocal condition number of the matrix */
/* > A after equilibration (if done). If RCOND is less than the */
/* > machine precision (in particular, if RCOND = 0), the matrix */
/* > is singular to working precision. This condition is */
/* > indicated by a return code of INFO > 0. */
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
/* > On exit, RWORK(1) contains the reciprocal pivot growth */
/* > factor norm(A)/norm(U). The "max absolute element" norm is */
/* > used. If RWORK(1) is much less than 1, then the stability */
/* > of the LU factorization of the (equilibrated) matrix A */
/* > could be poor. This also means that the solution X, condition */
/* > estimator RCOND, and forward error bound FERR could be */
/* > unreliable. If factorization fails with 0<INFO<=N, then */
/* > RWORK(1) contains the reciprocal pivot growth factor for the */
/* > leading INFO columns of A. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO is INTEGER */
/* > = 0: successful exit */
/* > < 0: if INFO = -i, the i-th argument had an illegal value */
/* > > 0: if INFO = i, and i is */
/* > <= N: U(i,i) is exactly zero. The factorization */
/* > has been completed, but the factor U is exactly */
/* > singular, so the solution and error bounds */
/* > could not be computed. RCOND = 0 is returned. */
/* > = N+1: U is nonsingular, but RCOND is less than machine */
/* > precision, meaning that the matrix is singular */
/* > to working precision. Nevertheless, the */
/* > solution and error bounds are computed because */
/* > there are a number of situations where the */
/* > computed solution can be more accurate than the */
/* > value of RCOND would suggest. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date April 2012 */
/* > \ingroup complexGBsolve */
/* ===================================================================== */
/* Subroutine */
/** Generated wrapper function */
void cgbsvx_(char *fact, char *trans, aocl_int_t *n, aocl_int_t *kl, aocl_int_t *ku,
             aocl_int_t *nrhs, scomplex *ab, aocl_int_t *ldab, scomplex *afb, aocl_int_t *ldafb,
             aocl_int_t *ipiv, char *equed, real *r__, real *c__, scomplex *b, aocl_int_t *ldb,
             scomplex *x, aocl_int_t *ldx, real *rcond, real *ferr, real *berr, scomplex *work,
             real *rwork, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r__, c__, b,
                       ldb, x, ldx, rcond, ferr, berr, work, rwork, info);
#else
    aocl_int64_t n_64 = *n;
    aocl_int64_t kl_64 = *kl;
    aocl_int64_t ku_64 = *ku;
    aocl_int64_t nrhs_64 = *nrhs;
    aocl_int64_t ldab_64 = *ldab;
    aocl_int64_t ldafb_64 = *ldafb;
    aocl_int64_t ldb_64 = *ldb;
    aocl_int64_t ldx_64 = *ldx;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgbsvx(fact, trans, &n_64, &kl_64, &ku_64, &nrhs_64, ab, &ldab_64, afb, &ldafb_64,
                       ipiv, equed, r__, c__, b, &ldb_64, x, &ldx_64, rcond, ferr, berr, work,
                       rwork, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

void aocl_lapack_cgbsvx(char *fact, char *trans, aocl_int64_t *n, aocl_int64_t *kl,
                        aocl_int64_t *ku, aocl_int64_t *nrhs, scomplex *ab, aocl_int64_t *ldab,
                        scomplex *afb, aocl_int64_t *ldafb, aocl_int_t *ipiv, char *equed, real *r__,
                        real *c__, scomplex *b, aocl_int64_t *ldb, scomplex *x, aocl_int64_t *ldx,
                        real *rcond, real *ferr, real *berr, scomplex *work, real *rwork,
                        aocl_int64_t *info)
{
    AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);
#if LF_AOCL_DTL_LOG_ENABLE
    char buffer[256];
#if FLA_ENABLE_ILP64
    snprintf(buffer, 256,
             "cgbsvx inputs: fact %c, trans %c, n %lld, kl %lld, ku %lld, nrhs %lld, ldab %lld, "
             "ldafb %lld, equed %c, ldb %lld, ldx %lld",
             *fact, *trans, *n, *kl, *ku, *nrhs, *ldab, *ldafb, *equed, *ldb, *ldx);
#else
    snprintf(buffer, 256,
             "cgbsvx inputs: fact %c, trans %c, n %d, kl %d, ku %d, nrhs %d, ldab %d, ldafb %d, "
             "equed %c, ldb %d, ldx %d",
             *fact, *trans, *n, *kl, *ku, *nrhs, *ldab, *ldafb, *equed, *ldb, *ldx);
#endif
    AOCL_DTL_LOG(AOCL_DTL_LEVEL_TRACE_5, buffer);
#endif
    /* System generated locals */
    aocl_int64_t ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, x_dim1, x_offset, i__1,
        i__2, i__3, i__4, i__5;
    real r__1, r__2;
    scomplex q__1;
    /* Builtin functions */
    double c_abs(scomplex *);
    /* Local variables */
    aocl_int64_t i__, j, j1, j2;
    real amax;
    char norm[1];
    extern logical lsame_(char *, char *, aocl_int64_t, aocl_int64_t);
    real rcmin, rcmax, anorm;
    logical equil;
    real colcnd;
    extern real slamch_(char *);
    logical nofact;
    real bignum;
    aocl_int64_t infequ;
    logical colequ;
    real rowcnd;
    logical notran;
    real smlnum;
    logical rowequ;
    real rpvgrw;
    /* -- LAPACK driver routine (version 3.4.1) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* April 2012 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ===================================================================== */
    /* Moved setting of INFO = N+1 so INFO does not subsequently get */
    /* overwritten. Sven, 17 Mar 05. */
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
    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    afb_dim1 = *ldafb;
    afb_offset = 1 + afb_dim1;
    afb -= afb_offset;
    --ipiv;
    --r__;
    --c__;
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
    nofact = lsame_(fact, "N", 1, 1);
    equil = lsame_(fact, "E", 1, 1);
    notran = lsame_(trans, "N", 1, 1);
    smlnum = 0.f;
    bignum = 0.f;
    if(nofact || equil)
    {
        *(unsigned char *)equed = 'N';
        rowequ = FALSE_;
        colequ = FALSE_;
    }
    else
    {
        rowequ = lsame_(equed, "R", 1, 1) || lsame_(equed, "B", 1, 1);
        colequ = lsame_(equed, "C", 1, 1) || lsame_(equed, "B", 1, 1);
        smlnum = slamch_("Safe minimum");
        bignum = 1.f / smlnum;
    }
    /* Test the input parameters. */
    if(!nofact && !equil && !lsame_(fact, "F", 1, 1))
    {
        *info = -1;
    }
    else if(!notran && !lsame_(trans, "T", 1, 1) && !lsame_(trans, "C", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*kl < 0)
    {
        *info = -4;
    }
    else if(*ku < 0)
    {
        *info = -5;
    }
    else if(*nrhs < 0)
    {
        *info = -6;
    }
    else if(*ldab < *kl + *ku + 1)
    {
        *info = -8;
    }
    else if(*ldafb < (*kl << 1) + *ku + 1)
    {
        *info = -10;
    }
    else if(lsame_(fact, "F", 1, 1) && !(rowequ || colequ || lsame_(equed, "N", 1, 1)))
    {
        *info = -12;
    }
    else
    {
        if(rowequ)
        {
            rcmin = bignum;
            rcmax = 0.f;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MIN */
                r__1 = rcmin;
                r__2 = r__[j]; // , expr subst
                rcmin = fla_min(r__1, r__2);
                /* Computing MAX */
                r__1 = rcmax;
                r__2 = r__[j]; // , expr subst
                rcmax = fla_max(r__1, r__2);
                /* L10: */
            }
            if(rcmin <= 0.f)
            {
                *info = -13;
            }
            else if(*n > 0)
            {
                rowcnd = fla_max(rcmin, smlnum) / fla_min(rcmax, bignum);
            }
            else
            {
                rowcnd = 1.f;
            }
        }
        if(colequ && *info == 0)
        {
            rcmin = bignum;
            rcmax = 0.f;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MIN */
                r__1 = rcmin;
                r__2 = c__[j]; // , expr subst
                rcmin = fla_min(r__1, r__2);
                /* Computing MAX */
                r__1 = rcmax;
                r__2 = c__[j]; // , expr subst
                rcmax = fla_max(r__1, r__2);
                /* L20: */
            }
            if(rcmin <= 0.f)
            {
                *info = -14;
            }
            else if(*n > 0)
            {
                colcnd = fla_max(rcmin, smlnum) / fla_min(rcmax, bignum);
            }
            else
            {
                colcnd = 1.f;
            }
        }
        if(*info == 0)
        {
            if(*ldb < fla_max(1, *n))
            {
                *info = -16;
            }
            else if(*ldx < fla_max(1, *n))
            {
                *info = -18;
            }
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        aocl_blas_xerbla("CGBSVX", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
        return;
    }
    if(equil)
    {
        /* Compute row and column scalings to equilibrate the matrix A. */
        aocl_lapack_cgbequ(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd, &colcnd,
                           &amax, &infequ);
        if(infequ == 0)
        {
            /* Equilibrate the matrix. */
            aocl_lapack_claqgb(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd,
                               &colcnd, &amax, equed);
            rowequ = lsame_(equed, "R", 1, 1) || lsame_(equed, "B", 1, 1);
            colequ = lsame_(equed, "C", 1, 1) || lsame_(equed, "B", 1, 1);
        }
    }
    /* Scale the right hand side. */
    if(notran)
    {
        if(rowequ)
        {
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                i__2 = *n;
                for(i__ = 1; i__ <= i__2; ++i__)
                {
                    i__3 = i__ + j * b_dim1;
                    i__4 = i__;
                    i__5 = i__ + j * b_dim1;
                    q__1.r = r__[i__4] * b[i__5].r;
                    q__1.i = r__[i__4] * b[i__5].i; // , expr subst
                    b[i__3].r = q__1.r;
                    b[i__3].i = q__1.i; // , expr subst
                    /* L30: */
                }
                /* L40: */
            }
        }
    }
    else if(colequ)
    {
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *n;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * b_dim1;
                i__4 = i__;
                i__5 = i__ + j * b_dim1;
                q__1.r = c__[i__4] * b[i__5].r;
                q__1.i = c__[i__4] * b[i__5].i; // , expr subst
                b[i__3].r = q__1.r;
                b[i__3].i = q__1.i; // , expr subst
                /* L50: */
            }
            /* L60: */
        }
    }
    if(nofact || equil)
    {
        /* Compute the LU factorization of the band matrix A. */
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            /* Computing MAX */
            i__2 = j - *ku;
            j1 = fla_max(i__2, 1);
            /* Computing MIN */
            i__2 = j + *kl;
            j2 = fla_min(i__2, *n);
            i__2 = j2 - j1 + 1;
            aocl_blas_ccopy(&i__2, &ab[*ku + 1 - j + j1 + j * ab_dim1], &c__1,
                            &afb[*kl + *ku + 1 - j + j1 + j * afb_dim1], &c__1);
            /* L70: */
        }
        aocl_lapack_cgbtrf(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);
        /* Return if INFO is non-zero. */
        if(*info > 0)
        {
            /* Compute the reciprocal pivot growth factor of the */
            /* leading rank-deficient INFO columns of A. */
            anorm = 0.f;
            i__1 = *info;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MAX */
                i__2 = *ku + 2 - j;
                /* Computing MIN */
                i__4 = *n + *ku + 1 - j;
                i__5 = *kl + *ku + 1; // , expr subst
                i__3 = fla_min(i__4, i__5);
                for(i__ = fla_max(i__2, 1); i__ <= i__3; ++i__)
                {
                    /* Computing MAX */
                    r__1 = anorm;
                    r__2 = c_abs(&ab[i__ + j * ab_dim1]); // , expr subst
                    anorm = fla_max(r__1, r__2);
                    /* L80: */
                }
                /* L90: */
            }
            /* Computing MIN */
            i__3 = *info - 1;
            i__2 = *kl + *ku; // , expr subst
            i__1 = fla_min(i__3, i__2);
            /* Computing MAX */
            i__4 = 1;
            i__5 = *kl + *ku + 2 - *info; // , expr subst
            rpvgrw = aocl_lapack_clantb("M", "U", "N", info, &i__1,
                                        &afb[fla_max(i__4, i__5) + afb_dim1], ldafb, &rwork[1]);
            if(rpvgrw == 0.f)
            {
                rpvgrw = 1.f;
            }
            else
            {
                rpvgrw = anorm / rpvgrw;
            }
            rwork[1] = rpvgrw;
            *rcond = 0.f;
            AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
            return;
        }
    }
    /* Compute the norm of the matrix A and the */
    /* reciprocal pivot growth factor RPVGRW. */
    if(notran)
    {
        *(unsigned char *)norm = '1';
    }
    else
    {
        *(unsigned char *)norm = 'I';
    }
    anorm = aocl_lapack_clangb(norm, n, kl, ku, &ab[ab_offset], ldab, &rwork[1]);
    i__1 = *kl + *ku;
    rpvgrw = aocl_lapack_clantb("M", "U", "N", n, &i__1, &afb[afb_offset], ldafb, &rwork[1]);
    if(rpvgrw == 0.f)
    {
        rpvgrw = 1.f;
    }
    else
    {
        rpvgrw = aocl_lapack_clangb("M", n, kl, ku, &ab[ab_offset], ldab, &rwork[1]) / rpvgrw;
    }
    /* Compute the reciprocal of the condition number of A. */
    aocl_lapack_cgbcon(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond, &work[1],
                       &rwork[1], info);
    /* Compute the solution matrix X. */
    aocl_lapack_clacpy("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    aocl_lapack_cgbtrs(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[x_offset], ldx,
                       info);
    /* Use iterative refinement to improve the computed solution and */
    /* compute error bounds and backward error estimates for it. */
    aocl_lapack_cgbrfs(trans, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], ldafb,
                       &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1],
                       &rwork[1], info);
    /* Transform the solution matrix X to a solution of the original */
    /* system. */
    if(notran)
    {
        if(colequ)
        {
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                i__3 = *n;
                for(i__ = 1; i__ <= i__3; ++i__)
                {
                    i__2 = i__ + j * x_dim1;
                    i__4 = i__;
                    i__5 = i__ + j * x_dim1;
                    q__1.r = c__[i__4] * x[i__5].r;
                    q__1.i = c__[i__4] * x[i__5].i; // , expr subst
                    x[i__2].r = q__1.r;
                    x[i__2].i = q__1.i; // , expr subst
                    /* L100: */
                }
                /* L110: */
            }
            i__1 = *nrhs;
            for(j = 1; j <= i__1; ++j)
            {
                ferr[j] /= colcnd;
                /* L120: */
            }
        }
    }
    else if(rowequ)
    {
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            i__3 = *n;
            for(i__ = 1; i__ <= i__3; ++i__)
            {
                i__2 = i__ + j * x_dim1;
                i__4 = i__;
                i__5 = i__ + j * x_dim1;
                q__1.r = r__[i__4] * x[i__5].r;
                q__1.i = r__[i__4] * x[i__5].i; // , expr subst
                x[i__2].r = q__1.r;
                x[i__2].i = q__1.i; // , expr subst
                /* L130: */
            }
            /* L140: */
        }
        i__1 = *nrhs;
        for(j = 1; j <= i__1; ++j)
        {
            ferr[j] /= rowcnd;
            /* L150: */
        }
    }
    /* Set INFO = N+1 if the matrix is singular to working precision. */
    if(*rcond < slamch_("Epsilon"))
    {
        *info = *n + 1;
    }
    rwork[1] = rpvgrw;
    AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);
    return;
    /* End of CGBSVX */
}
/* cgbsvx_ */
