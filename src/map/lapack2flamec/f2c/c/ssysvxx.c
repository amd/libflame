#ifdef FLA_ENABLE_XBLAS
/* ../netlib/ssysvxx.f -- translated by f2c (version 20160102). You must link the resulting object
 file with libf2c: on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm or, if you install libf2c.a in a
 standard place, with -lf2c -lm -- in that order, at the end of the command line, as in cc *.o -lf2c
 -lm Source for libf2c is in /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* > \brief \b SSYSVXX */
/* =========== DOCUMENTATION =========== */
/* Online html documentation available at */
/* http://www.netlib.org/lapack/explore-html/ */
/* > \htmlonly */
/* > Download SSYSVXX + dependencies */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysvxx
 * .f"> */
/* > [TGZ]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysvxx
 * .f"> */
/* > [ZIP]</a> */
/* > <a
 * href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysvxx
 * .f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */
/* Definition: */
/* =========== */
/* SUBROUTINE SSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, */
/* EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, */
/* N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, */
/* NPARAMS, PARAMS, WORK, IWORK, INFO ) */
/* .. Scalar Arguments .. */
/* CHARACTER EQUED, FACT, UPLO */
/* INTEGER INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, */
/* $ N_ERR_BNDS */
/* REAL RCOND, RPVGRW */
/* .. */
/* .. Array Arguments .. */
/* INTEGER IPIV( * ), IWORK( * ) */
/* REAL A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/* $ X( LDX, * ), WORK( * ) */
/* REAL S( * ), PARAMS( * ), BERR( * ), */
/* $ ERR_BNDS_NORM( NRHS, * ), */
/* $ ERR_BNDS_COMP( NRHS, * ) */
/* .. */
/* > \par Purpose: */
/* ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYSVXX uses the diagonal pivoting factorization to compute the */
/* > solution to a real system of linear equations A * X = B, where A */
/* > is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices. */
/* > */
/* > If requested, both normwise and maximum componentwise error bounds */
/* > are returned. SSYSVXX will return a solution with a tiny */
/* > guaranteed error (O(eps) where eps is the working machine */
/* > precision) unless the matrix is very ill-conditioned, in which */
/* > case a warning is returned. Relevant condition numbers also are */
/* > calculated and returned. */
/* > */
/* > SSYSVXX accepts user-provided factorizations and equilibration */
/* > factors;
see the definitions of the FACT and EQUED options. */
/* > Solving with refinement and using a factorization from a previous */
/* > SSYSVXX call will also produce a solution with either O(eps) */
/* > errors or warnings, but we cannot make that claim for general */
/* > user-provided factorizations and equilibration factors if they */
/* > differ from what SSYSVXX would itself produce. */
/* > \endverbatim */
/* > \par Description: */
/* ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed: */
/* > */
/* > 1. If FACT = 'E', real scaling factors are computed to equilibrate */
/* > the system: */
/* > */
/* > diag(S)*A*diag(S) *inv(diag(S))*X = diag(S)*B */
/* > */
/* > Whether or not the system will be equilibrated depends on the */
/* > scaling of the matrix A, but if equilibration is used, A is */
/* > overwritten by diag(S)*A*diag(S) and B by diag(S)*B. */
/* > */
/* > 2. If FACT = 'N' or 'E', the LU decomposition is used to factor */
/* > the matrix A (after equilibration if FACT = 'E') as */
/* > */
/* > A = U * D * U**T, if UPLO = 'U', or */
/* > A = L * D * L**T, if UPLO = 'L', */
/* > */
/* > where U (or L) is a product of permutation and unit upper (lower) */
/* > triangular matrices, and D is symmetric and block diagonal with */
/* > 1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > 3. If some D(i,i)=0, so that D is exactly singular, then the */
/* > routine returns with INFO = i. Otherwise, the factored form of A */
/* > is used to estimate the condition number of the matrix A (see */
/* > argument RCOND). If the reciprocal of the condition number is */
/* > less than machine precision, the routine still goes on to solve */
/* > for X and compute error bounds as described below. */
/* > */
/* > 4. The system of equations is solved for X using the factored form */
/* > of A. */
/* > */
/* > 5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero), */
/* > the routine will use iterative refinement to try to get a small */
/* > error and error bounds. Refinement calculates the residual to at */
/* > least twice the working precision. */
/* > */
/* > 6. If equilibration was used, the matrix X is premultiplied by */
/* > diag(R) so that it solves the original system before */
/* > equilibration. */
/* > \endverbatim */
/* Arguments: */
/* ========== */
/* > \verbatim */
/* > Some optional parameters are bundled in the PARAMS array. These */
/* > settings determine how refinement is performed, but often the */
/* > defaults are acceptable. If the defaults are acceptable, users */
/* > can pass NPARAMS = 0 which prevents the source code from accessing */
/* > the PARAMS argument. */
/* > \endverbatim */
/* > */
/* > \param[in] FACT */
/* > \verbatim */
/* > FACT is CHARACTER*1 */
/* > Specifies whether or not the factored form of the matrix A is */
/* > supplied on entry, and if not, whether the matrix A should be */
/* > equilibrated before it is factored. */
/* > = 'F': On entry, AF and IPIV contain the factored form of A. */
/* > If EQUED is not 'N', the matrix A has been */
/* > equilibrated with scaling factors given by S. */
/* > A, AF, and IPIV are not modified. */
/* > = 'N': The matrix A will be copied to AF and factored. */
/* > = 'E': The matrix A will be equilibrated if necessary, then */
/* > copied to AF and factored. */
/* > \endverbatim */
/* > */
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
/* > The number of linear equations, i.e., the order of the */
/* > matrix A. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* > NRHS is INTEGER */
/* > The number of right hand sides, i.e., the number of columns */
/* > of the matrices B and X. NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* > A is REAL array, dimension (LDA,N) */
/* > The symmetric matrix A. If UPLO = 'U', the leading N-by-N */
/* > upper triangular part of A contains the upper triangular */
/* > part of the matrix A, and the strictly lower triangular */
/* > part of A is not referenced. If UPLO = 'L', the leading */
/* > N-by-N lower triangular part of A contains the lower */
/* > triangular part of the matrix A, and the strictly upper */
/* > triangular part of A is not referenced. */
/* > */
/* > On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by */
/* > diag(S)*A*diag(S). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* > LDA is INTEGER */
/* > The leading dimension of the array A. LDA >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] AF */
/* > \verbatim */
/* > AF is REAL array, dimension (LDAF,N) */
/* > If FACT = 'F', then AF is an input argument and on entry */
/* > contains the block diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L from the factorization A = */
/* > U*D*U**T or A = L*D*L**T as computed by SSYTRF. */
/* > */
/* > If FACT = 'N', then AF is an output argument and on exit */
/* > returns the block diagonal matrix D and the multipliers */
/* > used to obtain the factor U or L from the factorization A = */
/* > U*D*U**T or A = L*D*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* > LDAF is INTEGER */
/* > The leading dimension of the array AF. LDAF >= fla_max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* > IPIV is INTEGER array, dimension (N) */
/* > If FACT = 'F', then IPIV is an input argument and on entry */
/* > contains details of the interchanges and the block */
/* > structure of D, as determined by SSYTRF. If IPIV(k) > 0, */
/* > then rows and columns k and IPIV(k) were interchanged and */
/* > D(k,k) is a 1-by-1 diagonal block. If UPLO = 'U' and */
/* > IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and */
/* > -IPIV(k) were interchanged and D(k-1:k,k-1:k) is a 2-by-2 */
/* > diagonal block. If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0, */
/* > then rows and columns k+1 and -IPIV(k) were interchanged */
/* > and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > */
/* > If FACT = 'N', then IPIV is an output argument and on exit */
/* > contains details of the interchanges and the block */
/* > structure of D, as determined by SSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* > EQUED is CHARACTER*1 */
/* > Specifies the form of equilibration that was done. */
/* > = 'N': No equilibration (always true if FACT = 'N'). */
/* > = 'Y': Both row and column equilibration, i.e., A has been */
/* > replaced by diag(S) * A * diag(S). */
/* > EQUED is an input argument if FACT = 'F';
otherwise, it is an */
/* > output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] S */
/* > \verbatim */
/* > S is REAL array, dimension (N) */
/* > The scale factors for A. If EQUED = 'Y', A is multiplied on */
/* > the left and right by diag(S). S is an input argument if FACT = */
/* > 'F';
otherwise, S is an output argument. If FACT = 'F' and EQUED */
/* > = 'Y', each element of S must be positive. If S is output, each */
/* > element of S is a power of the radix. If S is input, each element */
/* > of S should be a power of the radix to ensure a reliable solution */
/* > and error estimates. Scaling by powers of the radix does not cause */
/* > rounding errors unless the result underflows or overflows. */
/* > Rounding errors during scaling lead to refining with a matrix that */
/* > is not equivalent to the input matrix, producing error estimates */
/* > that may not be reliable. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* > B is REAL array, dimension (LDB,NRHS) */
/* > On entry, the N-by-NRHS right hand side matrix B. */
/* > On exit, */
/* > if EQUED = 'N', B is not modified;
 */
/* > if EQUED = 'Y', B is overwritten by diag(S)*B;
 */
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
/* > X is REAL array, dimension (LDX,NRHS) */
/* > If INFO = 0, the N-by-NRHS solution matrix X to the original */
/* > system of equations. Note that A and B are modified on exit if */
/* > EQUED .ne. 'N', and the solution to the equilibrated system is */
/* > inv(diag(S))*X. */
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
/* > Reciprocal scaled condition number. This is an estimate of the */
/* > reciprocal Skeel condition number of the matrix A after */
/* > equilibration (if done). If this is less than the machine */
/* > precision (in particular, if it is zero), the matrix is singular */
/* > to working precision. Note that the error may still be small even */
/* > if this number is very small and the matrix appears ill- */
/* > conditioned. */
/* > \endverbatim */
/* > */
/* > \param[out] RPVGRW */
/* > \verbatim */
/* > RPVGRW is REAL */
/* > Reciprocal pivot growth. On exit, this contains the reciprocal */
/* > pivot growth factor norm(A)/norm(U). The "max absolute element" */
/* > norm is used. If this is much less than 1, then the stability of */
/* > the LU factorization of the (equilibrated) matrix A could be poor. */
/* > This also means that the solution X, estimated condition numbers, */
/* > and error bounds could be unreliable. If factorization fails with */
/* > 0<INFO<=N, then this contains the reciprocal pivot growth factor */
/* > for the leading INFO columns of A. */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* > BERR is REAL array, dimension (NRHS) */
/* > Componentwise relative backward error. This is the */
/* > componentwise relative backward error of each solution vector X(j) */
/* > (i.e., the smallest relative change in any element of A or B that */
/* > makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[in] N_ERR_BNDS */
/* > \verbatim */
/* > N_ERR_BNDS is INTEGER */
/* > Number of error bounds to return for each right hand side */
/* > and each type (normwise or componentwise). See ERR_BNDS_NORM and */
/* > ERR_BNDS_COMP below. */
/* > \endverbatim */
/* > */
/* > \param[out] ERR_BNDS_NORM */
/* > \verbatim */
/* > ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) */
/* > For each right-hand side, this array contains information about */
/* > various error bounds and condition numbers corresponding to the */
/* > normwise relative error, which is defined as follows: */
/* > */
/* > Normwise relative error in the ith solution vector: */
/* > max_j (f2c_abs(XTRUE(j,i) - X(j,i))) */
/* > ------------------------------ */
/* > max_j f2c_abs(X(j,i)) */
/* > */
/* > The array is indexed by the type of error information as described */
/* > below. There currently are up to three pieces of information */
/* > returned. */
/* > */
/* > The first index in ERR_BNDS_NORM(i,:) corresponds to the ith */
/* > right-hand side. */
/* > */
/* > The second index in ERR_BNDS_NORM(:,err) contains the following */
/* > three fields: */
/* > err = 1 "Trust/don't trust" boolean. Trust the answer if the */
/* > reciprocal condition number is less than the threshold */
/* > sqrt(n) * slamch('Epsilon'). */
/* > */
/* > err = 2 "Guaranteed" error bound: The estimated forward error, */
/* > almost certainly within a factor of 10 of the true error */
/* > so long as the next entry is greater than the threshold */
/* > sqrt(n) * slamch('Epsilon'). This error bound should only */
/* > be trusted if the previous boolean is true. */
/* > */
/* > err = 3 Reciprocal condition number: Estimated normwise */
/* > reciprocal condition number. Compared with the threshold */
/* > sqrt(n) * slamch('Epsilon') to determine if the error */
/* > estimate is "guaranteed". These reciprocal condition */
/* > numbers are 1 / (norm(Z^{
-1}
,inf) * norm(Z,inf)) for some */
/* > appropriately scaled matrix Z. */
/* > Let Z = S*A, where S scales each row by a power of the */
/* > radix so all absolute row sums of Z are approximately 1. */
/* > */
/* > See Lapack Working Note 165 for further details and extra */
/* > cautions. */
/* > \endverbatim */
/* > */
/* > \param[out] ERR_BNDS_COMP */
/* > \verbatim */
/* > ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) */
/* > For each right-hand side, this array contains information about */
/* > various error bounds and condition numbers corresponding to the */
/* > componentwise relative error, which is defined as follows: */
/* > */
/* > Componentwise relative error in the ith solution vector: */
/* > f2c_abs(XTRUE(j,i) - X(j,i)) */
/* > max_j ---------------------- */
/* > f2c_abs(X(j,i)) */
/* > */
/* > The array is indexed by the right-hand side i (on which the */
/* > componentwise relative error depends), and the type of error */
/* > information as described below. There currently are up to three */
/* > pieces of information returned for each right-hand side. If */
/* > componentwise accuracy is not requested (PARAMS(3) = 0.0), then */
/* > ERR_BNDS_COMP is not accessed. If N_ERR_BNDS < 3, then at most */
/* > the first (:,N_ERR_BNDS) entries are returned. */
/* > */
/* > The first index in ERR_BNDS_COMP(i,:) corresponds to the ith */
/* > right-hand side. */
/* > */
/* > The second index in ERR_BNDS_COMP(:,err) contains the following */
/* > three fields: */
/* > err = 1 "Trust/don't trust" boolean. Trust the answer if the */
/* > reciprocal condition number is less than the threshold */
/* > sqrt(n) * slamch('Epsilon'). */
/* > */
/* > err = 2 "Guaranteed" error bound: The estimated forward error, */
/* > almost certainly within a factor of 10 of the true error */
/* > so long as the next entry is greater than the threshold */
/* > sqrt(n) * slamch('Epsilon'). This error bound should only */
/* > be trusted if the previous boolean is true. */
/* > */
/* > err = 3 Reciprocal condition number: Estimated componentwise */
/* > reciprocal condition number. Compared with the threshold */
/* > sqrt(n) * slamch('Epsilon') to determine if the error */
/* > estimate is "guaranteed". These reciprocal condition */
/* > numbers are 1 / (norm(Z^{
-1}
,inf) * norm(Z,inf)) for some */
/* > appropriately scaled matrix Z. */
/* > Let Z = S*(A*diag(x)), where x is the solution for the */
/* > current right-hand side and S scales each row of */
/* > A*diag(x) by a power of the radix so all absolute row */
/* > sums of Z are approximately 1. */
/* > */
/* > See Lapack Working Note 165 for further details and extra */
/* > cautions. */
/* > \endverbatim */
/* > */
/* > \param[in] NPARAMS */
/* > \verbatim */
/* > NPARAMS is INTEGER */
/* > Specifies the number of parameters set in PARAMS. If <= 0, the */
/* > PARAMS array is never referenced and default values are used. */
/* > \endverbatim */
/* > */
/* > \param[in,out] PARAMS */
/* > \verbatim */
/* > PARAMS is REAL array, dimension NPARAMS */
/* > Specifies algorithm parameters. If an entry is < 0.0, then */
/* > that entry will be filled with default value used for that */
/* > parameter. Only positions up to NPARAMS are accessed;
defaults */
/* > are used for higher-numbered parameters. */
/* > */
/* > PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative */
/* > refinement or not. */
/* > Default: 1.0 */
/* > = 0.0: No refinement is performed, and no error bounds are */
/* > computed. */
/* > = 1.0: Use the double-precision refinement algorithm, */
/* > possibly with doubled-single computations if the */
/* > compilation environment does not support DOUBLE */
/* > PRECISION. */
/* > (other values are reserved for future use) */
/* > */
/* > PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual */
/* > computations allowed for refinement. */
/* > Default: 10 */
/* > Aggressive: Set to 100 to permit convergence using approximate */
/* > factorizations or factorizations other than LU. If */
/* > the factorization uses a technique other than */
/* > Gaussian elimination, the guarantees in */
/* > err_bnds_norm and err_bnds_comp may no longer be */
/* > trustworthy. */
/* > */
/* > PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code */
/* > will attempt to find a solution with small componentwise */
/* > relative error in the double-precision algorithm. Positive */
/* > is true, 0.0 is false. */
/* > Default: 1.0 (attempt componentwise convergence) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK is REAL array, dimension (4*N) */
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
/* > = 0: Successful exit. The solution to every right-hand side is */
/* > guaranteed. */
/* > < 0: If INFO = -i, the i-th argument had an illegal value */
/* > > 0 and <= N: U(INFO,INFO) is exactly zero. The factorization */
/* > has been completed, but the factor U is exactly singular, so */
/* > the solution and error bounds could not be computed. RCOND = 0 */
/* > is returned. */
/* > = N+J: The solution corresponding to the Jth right-hand side is */
/* > not guaranteed. The solutions corresponding to other right- */
/* > hand sides K with K > J may not be guaranteed as well, but */
/* > only the first such right-hand side is reported. If a small */
/* > componentwise error is not requested (PARAMS(3) = 0.0) then */
/* > the Jth right-hand side is the first with a normwise error */
/* > bound that is not guaranteed (the smallest J such */
/* > that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0) */
/* > the Jth right-hand side is the first with either a normwise or */
/* > componentwise error bound that is not guaranteed (the smallest */
/* > J such that either ERR_BNDS_NORM(J,1) = 0.0 or */
/* > ERR_BNDS_COMP(J,1) = 0.0). See the definition of */
/* > ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information */
/* > about all of the right-hand sides check ERR_BNDS_NORM or */
/* > ERR_BNDS_COMP. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */
/* > \date December 2016 */
/* > \ingroup realSYsolve */
/* ===================================================================== */
/* Subroutine */
void ssysvxx_(char *fact, char *uplo, integer *n, integer *nrhs, real *a, integer *lda, real *af,
              integer *ldaf, integer *ipiv, char *equed, real *s, real *b, integer *ldb, real *x,
              integer *ldx, real *rcond, real *rpvgrw, real *berr, integer *n_err_bnds__,
              real *err_bnds_norm__, real *err_bnds_comp__, integer *nparams, real *params,
              real *work, integer *iwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("ssysvxx inputs: fact %c ,uplo %c ,n %" FLA_IS ",nrhs %" FLA_IS
                      ",lda %" FLA_IS ",ldaf %" FLA_IS ",equed %c ,ldb %" FLA_IS ",ldx %" FLA_IS
                      ",nparams %" FLA_IS "",
                      *fact, *uplo, *n, *nrhs, *lda, *ldaf, *equed, *ldb, *ldx, *nparams);
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, x_offset,
        err_bnds_norm_dim1, err_bnds_norm_offset, err_bnds_comp_dim1, err_bnds_comp_offset, i__1;
    real r__1, r__2;
    /* Local variables */
    extern /* Subroutine */
        void
        ssyrfsx_(char *, char *, integer *, integer *, real *, integer *, real *, integer *,
                 integer *, real *, real *, integer *, real *, integer *, real *, real *, integer *,
                 real *, real *, integer *, real *, real *, integer *, integer *);
    integer j;
    real amax, smin, smax;
    extern real sla_syrpvgrw_(char *, integer *, integer *, real *, integer *, real *, integer *,
                              integer *, real *);
    extern logical lsame_(char *, char *, integer, integer);
    real scond;
    logical equil, rcequ;
    extern real slamch_(char *);
    logical nofact;
    extern /* Subroutine */
        void
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real bignum;
    integer infequ;
    extern /* Subroutine */
        void
        slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    real smlnum;
    extern /* Subroutine */
        void
        slaqsy_(char *, integer *, real *, integer *, real *, real *, real *, char *),
        ssytrf_(char *, integer *, real *, integer *, integer *, real *, integer *, integer *),
        slascl2_(integer *, integer *, real *, real *, integer *),
        ssytrs_(char *, integer *, integer *, real *, integer *, integer *, real *, integer *,
                integer *),
        ssyequb_(char *, integer *, real *, integer *, real *, real *, real *, real *, integer *);
    /* -- LAPACK driver routine (version 3.7.0) -- */
    /* -- LAPACK is a software package provided by Univ. of Tennessee, -- */
    /* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /* December 2016 */
    /* .. Scalar Arguments .. */
    /* .. */
    /* .. Array Arguments .. */
    /* .. */
    /* ================================================================== */
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
    err_bnds_comp_dim1 = *nrhs;
    err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
    err_bnds_comp__ -= err_bnds_comp_offset;
    err_bnds_norm_dim1 = *nrhs;
    err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
    err_bnds_norm__ -= err_bnds_norm_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    af_dim1 = *ldaf;
    af_offset = 1 + af_dim1;
    af -= af_offset;
    --ipiv;
    --s;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --berr;
    --params;
    --work;
    --iwork;
    /* Function Body */
    *info = 0;
    nofact = lsame_(fact, "N", 1, 1);
    equil = lsame_(fact, "E", 1, 1);
    smlnum = slamch_("Safe minimum");
    bignum = 1.f / smlnum;
    if(nofact || equil)
    {
        *(unsigned char *)equed = 'N';
        rcequ = FALSE_;
    }
    else
    {
        rcequ = lsame_(equed, "Y", 1, 1);
    }
    /* Default is failure. If an input parameter is wrong or */
    /* factorization fails, make everything look horrible. Only the */
    /* pivot growth is set here, the rest is initialized in SSYRFSX. */
    *rpvgrw = 0.f;
    /* Test the input parameters. PARAMS is not tested until SSYRFSX. */
    if(!nofact && !equil && !lsame_(fact, "F", 1, 1))
    {
        *info = -1;
    }
    else if(!lsame_(uplo, "U", 1, 1) && !lsame_(uplo, "L", 1, 1))
    {
        *info = -2;
    }
    else if(*n < 0)
    {
        *info = -3;
    }
    else if(*nrhs < 0)
    {
        *info = -4;
    }
    else if(*lda < fla_max(1, *n))
    {
        *info = -6;
    }
    else if(*ldaf < fla_max(1, *n))
    {
        *info = -8;
    }
    else if(lsame_(fact, "F", 1, 1) && !(rcequ || lsame_(equed, "N", 1, 1)))
    {
        *info = -10;
    }
    else
    {
        if(rcequ)
        {
            smin = bignum;
            smax = 0.f;
            i__1 = *n;
            for(j = 1; j <= i__1; ++j)
            {
                /* Computing MIN */
                r__1 = smin;
                r__2 = s[j]; // , expr subst
                smin = fla_min(r__1, r__2);
                /* Computing MAX */
                r__1 = smax;
                r__2 = s[j]; // , expr subst
                smax = fla_max(r__1, r__2);
                /* L10: */
            }
            if(smin <= 0.f)
            {
                *info = -11;
            }
            else if(*n > 0)
            {
                scond = fla_max(smin, smlnum) / fla_min(smax, bignum);
            }
            else
            {
                scond = 1.f;
            }
        }
        if(*info == 0)
        {
            if(*ldb < fla_max(1, *n))
            {
                *info = -13;
            }
            else if(*ldx < fla_max(1, *n))
            {
                *info = -15;
            }
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SSYSVXX", &i__1, (ftnlen)7);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(equil)
    {
        /* Compute row and column scalings to equilibrate the matrix A. */
        ssyequb_(uplo, n, &a[a_offset], lda, &s[1], &scond, &amax, &work[1], &infequ);
        if(infequ == 0)
        {
            /* Equilibrate the matrix. */
            slaqsy_(uplo, n, &a[a_offset], lda, &s[1], &scond, &amax, equed);
            rcequ = lsame_(equed, "Y", 1, 1);
        }
    }
    /* Scale the right-hand side. */
    if(rcequ)
    {
        slascl2_(n, nrhs, &s[1], &b[b_offset], ldb);
    }
    if(nofact || equil)
    {
        /* Compute the LDL^T or UDU^T factorization of A. */
        slacpy_(uplo, n, n, &a[a_offset], lda, &af[af_offset], ldaf);
        i__1 = fla_max(1, *n) * 5;
        ssytrf_(uplo, n, &af[af_offset], ldaf, &ipiv[1], &work[1], &i__1, info);
        /* Return if INFO is non-zero. */
        if(*info > 0)
        {
            /* Pivot in column INFO is exactly 0 */
            /* Compute the reciprocal pivot growth factor of the */
            /* leading rank-deficient INFO columns of A. */
            if(*n > 0)
            {
                *rpvgrw = sla_syrpvgrw_(uplo, n, info, &a[a_offset], lda, &af[af_offset], ldaf,
                                        &ipiv[1], &work[1]);
            }
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
    }
    /* Compute the reciprocal pivot growth factor RPVGRW. */
    if(*n > 0)
    {
        *rpvgrw = sla_syrpvgrw_(uplo, n, info, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1],
                                &work[1]);
    }
    /* Compute the solution matrix X. */
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx);
    ssytrs_(uplo, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx, info);
    /* Use iterative refinement to improve the computed solution and */
    /* compute error bounds and backward error estimates for it. */
    ssyrfsx_(uplo, equed, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1], &s[1],
             &b[b_offset], ldb, &x[x_offset], ldx, rcond, &berr[1], n_err_bnds__,
             &err_bnds_norm__[err_bnds_norm_offset], &err_bnds_comp__[err_bnds_comp_offset],
             nparams, &params[1], &work[1], &iwork[1], info);
    /* Scale solutions. */
    if(rcequ)
    {
        slascl2_(n, nrhs, &s[1], &x[x_offset], ldx);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* End of SSYSVXX */
}
/* ssysvxx_ */
#endif
