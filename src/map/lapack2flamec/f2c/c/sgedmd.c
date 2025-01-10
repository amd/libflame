/* ./sgedmd.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static real c_b46 = 1.f;
static real c_b56 = 0.f;
static integer c__2 = 2;
static real c_b154 = -1.f;
/* > \brief \b SGEDMD computes the Dynamic Mode Decomposition (DMD) for a pair of data snapshot
 * matrices. */
/* =========== DOCUMENTATION =========== */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEDMD( JOBS, JOBZ, JOBR, JOBF, WHTSVD, & */
/* M, N, X, LDX, Y, LDY, NRNK, TOL, & */
/* K, REIG, IMEIG, Z, LDZ, RES, & */
/* B, LDB, W, LDW, S, LDS, & */
/* WORK, LWORK, IWORK, LIWORK, INFO ) */
/* ..... */
/* USE iso_fortran_env */
/* IMPLICIT NONE */
/* INTEGER, PARAMETER :: WP = real32 */
/* ..... */
/* Scalar arguments */
/* CHARACTER, INTENT(IN) :: JOBS, JOBZ, JOBR, JOBF */
/* INTEGER, INTENT(IN) :: WHTSVD, M, N, LDX, LDY, & */
/* NRNK, LDZ, LDB, LDW, LDS, & */
/* LWORK, LIWORK */
/* INTEGER, INTENT(OUT) :: K, INFO */
/* REAL, INTENT(IN) :: TOL */
/* Array arguments */
/* REAL, INTENT(INOUT) :: X(LDX,*), Y(LDY,*) */
/* REAL, INTENT(OUT) :: Z(LDZ,*), B(LDB,*), & */
/* W(LDW,*), S(LDS,*) */
/* REAL, INTENT(OUT) :: REIG(*), IMEIG(*), & */
/* RES(*) */
/* REAL, INTENT(OUT) :: WORK(*) */
/* INTEGER, INTENT(OUT) :: IWORK(*) */
/* ............................................................ */
/* > \par Purpose: */
/* ============= */
/* > \verbatim */
/* > SGEDMD computes the Dynamic Mode Decomposition (DMD) for */
/* > a pair of data snapshot matrices. For the input matrices */
/* > X and Y such that Y = A*X with an unaccessible matrix */
/* > A, SGEDMD computes a certain number of Ritz pairs of A using */
/* > the standard Rayleigh-Ritz extraction from a subspace of */
/* > range(X) that is determined using the leading left singular */
/* > vectors of X. Optionally, SGEDMD returns the residuals */
/* > of the computed Ritz pairs, the information needed for */
/* > a refinement of the Ritz vectors, or the eigenvectors of */
/* > the Exact DMD. */
/* > For further details see the references listed */
/* > below. For more details of the implementation see [3]. */
/* > \endverbatim */
/* ............................................................ */
/* > \par References: */
/* ================ */
/* > \verbatim */
/* > [1] P. Schmid: Dynamic mode decomposition of numerical */
/* > and experimental data, */
/* > Journal of Fluid Mechanics 656, 5-28, 2010. */
/* > [2] Z. Drmac, I. Mezic, R. Mohr: Data driven modal */
/* > decompositions: analysis and enhancements, */
/* > SIAM J. on Sci. Comp. 40 (4), A2253-A2285, 2018. */
/* > [3] Z. Drmac: A LAPACK implementation of the Dynamic */
/* > Mode Decomposition I. Technical report. AIMDyn Inc. */
/* > and LAPACK Working Note 298. */
/* > [4] J. Tu, C. W. Rowley, D. M. Luchtenburg, S. L. */
/* > Brunton, N. Kutz: On Dynamic Mode Decomposition: */
/* > Theory and Applications, Journal of Computational */
/* > Dynamics 1(2), 391 -421, 2014. */
/* > \endverbatim */
/* ...................................................................... */
/* > \par Developed and supported by: */
/* ================================ */
/* > \verbatim */
/* > Developed and coded by Zlatko Drmac, Faculty of Science, */
/* > University of Zagreb;
drmac@math.hr */
/* > In cooperation with */
/* > AIMdyn Inc., Santa Barbara, CA. */
/* > and supported by */
/* > - DARPA SBIR project "Koopman Operator-Based Forecasting */
/* > for Nonstationary Processes from Near-Term, Limited */
/* > Observational Data" Contract No: W31P4Q-21-C-0007 */
/* > - DARPA PAI project "Physics-Informed Machine Learning */
/* > Methodologies" Contract No: HR0011-18-9-0033 */
/* > - DARPA MoDyL project "A Data-Driven, Operator-Theoretic */
/* > Framework for Space-Time Analysis of Process Dynamics" */
/* > Contract No: HR0011-16-C-0116 */
/* > Any opinions, findings and conclusions or recommendations */
/* > expressed in this material are those of the author and */
/* > do not necessarily reflect the views of the DARPA SBIR */
/* > Program Office */
/* > \endverbatim */
/* ...................................................................... */
/* > \par Distribution Statement A: */
/* ============================== */
/* > \verbatim */
/* > Distribution Statement A: */
/* > Approved for Public Release, Distribution Unlimited. */
/* > Cleared by DARPA on September 29, 2022 */
/* > \endverbatim */
/* ============================================================ */
/* Arguments */
/* ========= */
/* > \param[in] JOBS */
/* > \verbatim */
/* > JOBS (input) CHARACTER*1 */
/* > Determines whether the initial data snapshots are scaled */
/* > by a diagonal matrix. */
/* > 'S' :: The data snapshots matrices X and Y are multiplied */
/* > with a diagonal matrix D so that X*D has unit */
/* > nonzero columns (in the Euclidean 2-norm) */
/* > 'C' :: The snapshots are scaled as with the 'S' option. */
/* > If it is found that an i-th column of X is zero */
/* > vector and the corresponding i-th column of Y is */
/* > non-zero, then the i-th column of Y is set to */
/* > zero and a warning flag is raised. */
/* > 'Y' :: The data snapshots matrices X and Y are multiplied */
/* > by a diagonal matrix D so that Y*D has unit */
/* > nonzero columns (in the Euclidean 2-norm) */
/* > 'N' :: No data scaling. */
/* > \endverbatim */
/* ..... */
/* > \param[in] JOBZ */
/* > \verbatim */
/* > JOBZ (input) CHARACTER*1 */
/* > Determines whether the eigenvectors (Koopman modes) will */
/* > be computed. */
/* > 'V' :: The eigenvectors (Koopman modes) will be computed */
/* > and returned in the matrix Z. */
/* > See the description of Z. */
/* > 'F' :: The eigenvectors (Koopman modes) will be returned */
/* > in factored form as the product X(:,1:K)*W, where X */
/* > contains a POD basis (leading left singular vectors */
/* > of the data matrix X) and W contains the eigenvectors */
/* > of the corresponding Rayleigh quotient. */
/* > See the descriptions of K, X, W, Z. */
/* > 'N' :: The eigenvectors are not computed. */
/* > \endverbatim */
/* ..... */
/* > \param[in] JOBR */
/* > \verbatim */
/* > JOBR (input) CHARACTER*1 */
/* > Determines whether to compute the residuals. */
/* > 'R' :: The residuals for the computed eigenpairs will be */
/* > computed and stored in the array RES. */
/* > See the description of RES. */
/* > For this option to be legal, JOBZ must be 'V'. */
/* > 'N' :: The residuals are not computed. */
/* > \endverbatim */
/* ..... */
/* > \param[in] JOBF */
/* > \verbatim */
/* > JOBF (input) CHARACTER*1 */
/* > Specifies whether to store information needed for post- */
/* > processing (e.g. computing refined Ritz vectors) */
/* > 'R' :: The matrix needed for the refinement of the Ritz */
/* > vectors is computed and stored in the array B. */
/* > See the description of B. */
/* > 'E' :: The unscaled eigenvectors of the Exact DMD are */
/* > computed and returned in the array B. See the */
/* > description of B. */
/* > 'N' :: No eigenvector refinement data is computed. */
/* > \endverbatim */
/* ..... */
/* > \param[in] WHTSVD */
/* > \verbatim */
/* > WHTSVD (input) INTEGER, WHSTVD in {
1, 2, 3, 4 }
*/
/* > Allows for a selection of the SVD algorithm from the */
/* > LAPACK library. */
/* > 1 :: SGESVD (the QR SVD algorithm) */
/* > 2 :: SGESDD (the Divide and Conquer algorithm;
if enough */
/* > workspace available, this is the fastest option) */
/* > 3 :: SGESVDQ (the preconditioned QR SVD ;
this and 4 */
/* > are the most accurate options) */
/* > 4 :: SGEJSV (the preconditioned Jacobi SVD;
this and 3 */
/* > are the most accurate options) */
/* > For the four methods above, a significant difference in */
/* > the accuracy of small singular values is possible if */
/* > the snapshots vary in norm so that X is severely */
/* > ill-conditioned. If small (smaller than EPS*||X||) */
/* > singular values are of interest and JOBS=='N', then */
/* > the options (3, 4) give the most accurate results, where */
/* > the option 4 is slightly better and with stronger */
/* > theoretical background. */
/* > If JOBS=='S', i.e. the columns of X will be normalized, */
/* > then all methods give nearly equally accurate results. */
/* > \endverbatim */
/* ..... */
/* > \param[in] M */
/* > \verbatim */
/* > M (input) INTEGER, M>= 0 */
/* > The state space dimension (the row dimension of X, Y). */
/* > \endverbatim */
/* ..... */
/* > \param[in] N */
/* > \verbatim */
/* > N (input) INTEGER, 0 <= N <= M */
/* > The number of data snapshot pairs */
/* > (the number of columns of X and Y). */
/* > \endverbatim */
/* ..... */
/* > \param[in,out] X */
/* > \verbatim */
/* > X (input/output) REAL M-by-N array */
/* > > On entry, X contains the data snapshot matrix X. It is */
/* > assumed that the column norms of X are in the range of */
/* > the normalized floating point numbers. */
/* > < On exit, the leading K columns of X contain a POD basis, */
/* > i.e. the leading K left singular vectors of the input */
/* > data matrix X, U(:,1:K). All N columns of X contain all */
/* > left singular vectors of the input matrix X. */
/* > See the descriptions of K, Z and W. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX (input) INTEGER, LDX >= M */
/* > The leading dimension of the array X. */
/* > \endverbatim */
/* ..... */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y (input/workspace/output) REAL M-by-N array */
/* > > On entry, Y contains the data snapshot matrix Y */
/* > < On exit, */
/* > If JOBR == 'R', the leading K columns of Y contain */
/* > the residual vectors for the computed Ritz pairs. */
/* > See the description of RES. */
/* > If JOBR == 'N', Y contains the original input data, */
/* > scaled according to the value of JOBS. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDY */
/* > \verbatim */
/* > LDY (input) INTEGER , LDY >= M */
/* > The leading dimension of the array Y. */
/* > \endverbatim */
/* ..... */
/* > \param[in] NRNK */
/* > \verbatim */
/* > NRNK (input) INTEGER */
/* > Determines the mode how to compute the numerical rank, */
/* > i.e. how to truncate small singular values of the input */
/* > matrix X. On input, if */
/* > NRNK = -1 :: i-th singular value sigma(i) is truncated */
/* > if sigma(i) <= TOL*sigma(1) */
/* > This option is recommended. */
/* > NRNK = -2 :: i-th singular value sigma(i) is truncated */
/* > if sigma(i) <= TOL*sigma(i-1) */
/* > This option is included for R&D purposes. */
/* > It requires highly accurate SVD, which */
/* > may not be feasible. */
/* > The numerical rank can be enforced by using positive */
/* > value of NRNK as follows: */
/* > 0 < NRNK <= N :: at most NRNK largest singular values */
/* > will be used. If the number of the computed nonzero */
/* > singular values is less than NRNK, then only those */
/* > nonzero values will be used and the actually used */
/* > dimension is less than NRNK. The actual number of */
/* > the nonzero singular values is returned in the variable */
/* > K. See the descriptions of TOL and K. */
/* > \endverbatim */
/* ..... */
/* > \param[in] TOL */
/* > \verbatim */
/* > TOL (input) REAL, 0 <= TOL < 1 */
/* > The tolerance for truncating small singular values. */
/* > See the description of NRNK. */
/* > \endverbatim */
/* ..... */
/* > \param[out] K */
/* > \verbatim */
/* > K (output) INTEGER, 0 <= K <= N */
/* > The dimension of the POD basis for the data snapshot */
/* > matrix X and the number of the computed Ritz pairs. */
/* > The value of K is determined according to the rule set */
/* > by the parameters NRNK and TOL. */
/* > See the descriptions of NRNK and TOL. */
/* > \endverbatim */
/* ..... */
/* > \param[out] REIG */
/* > \verbatim */
/* > REIG (output) REAL N-by-1 array */
/* > The leading K (K<=N) entries of REIG contain */
/* > the real parts of the computed eigenvalues */
/* > REIG(1:K) + sqrt(-1)*IMEIG(1:K). */
/* > See the descriptions of K, IMEIG, and Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] IMEIG */
/* > \verbatim */
/* > IMEIG (output) REAL N-by-1 array */
/* > The leading K (K<=N) entries of IMEIG contain */
/* > the imaginary parts of the computed eigenvalues */
/* > REIG(1:K) + sqrt(-1)*IMEIG(1:K). */
/* > The eigenvalues are determined as follows: */
/* > If IMEIG(i) == 0, then the corresponding eigenvalue is */
/* > real, LAMBDA(i) = REIG(i). */
/* > If IMEIG(i)>0, then the corresponding complex */
/* > conjugate pair of eigenvalues reads */
/* > LAMBDA(i) = REIG(i) + sqrt(-1)*IMAG(i) */
/* > LAMBDA(i+1) = REIG(i) - sqrt(-1)*IMAG(i) */
/* > That is, complex conjugate pairs have consecutive */
/* > indices (i,i+1), with the positive imaginary part */
/* > listed first. */
/* > See the descriptions of K, REIG, and Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] Z */
/* > \verbatim */
/* > Z (workspace/output) REAL M-by-N array */
/* > If JOBZ =='V' then */
/* > Z contains real Ritz vectors as follows: */
/* > If IMEIG(i)=0, then Z(:,i) is an eigenvector of */
/* > the i-th Ritz value;
||Z(:,i)||_2=1. */
/* > If IMEIG(i) > 0 (and IMEIG(i+1) < 0) then */
/* > [Z(:,i) Z(:,i+1)] span an invariant subspace and */
/* > the Ritz values extracted from this subspace are */
/* > REIG(i) + sqrt(-1)*IMEIG(i) and */
/* > REIG(i) - sqrt(-1)*IMEIG(i). */
/* > The corresponding eigenvectors are */
/* > Z(:,i) + sqrt(-1)*Z(:,i+1) and */
/* > Z(:,i) - sqrt(-1)*Z(:,i+1), respectively. */
/* > || Z(:,i:i+1)||_F = 1. */
/* > If JOBZ == 'F', then the above descriptions hold for */
/* > the columns of X(:,1:K)*W(1:K,1:K), where the columns */
/* > of W(1:k,1:K) are the computed eigenvectors of the */
/* > K-by-K Rayleigh quotient. The columns of W(1:K,1:K) */
/* > are similarly structured: If IMEIG(i) == 0 then */
/* > X(:,1:K)*W(:,i) is an eigenvector, and if IMEIG(i)>0 */
/* > then X(:,1:K)*W(:,i)+sqrt(-1)*X(:,1:K)*W(:,i+1) and */
/* > X(:,1:K)*W(:,i)-sqrt(-1)*X(:,1:K)*W(:,i+1) */
/* > are the eigenvectors of LAMBDA(i), LAMBDA(i+1). */
/* > See the descriptions of REIG, IMEIG, X and W. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDZ */
/* > \verbatim */
/* > LDZ (input) INTEGER , LDZ >= M */
/* > The leading dimension of the array Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] RES */
/* > \verbatim */
/* > RES (output) REAL N-by-1 array */
/* > RES(1:K) contains the residuals for the K computed */
/* > Ritz pairs. */
/* > If LAMBDA(i) is real, then */
/* > RES(i) = || A * Z(:,i) - LAMBDA(i)*Z(:,i))||_2. */
/* > If [LAMBDA(i), LAMBDA(i+1)] is a complex conjugate pair */
/* > then */
/* > RES(i)=RES(i+1) = || A * Z(:,i:i+1) - Z(:,i:i+1) *B||_F */
/* > where B = [ real(LAMBDA(i)) imag(LAMBDA(i)) ] */
/* > [-imag(LAMBDA(i)) real(LAMBDA(i)) ]. */
/* > It holds that */
/* > RES(i) = || A*ZC(:,i) - LAMBDA(i) *ZC(:,i) ||_2 */
/* > RES(i+1) = || A*ZC(:,i+1) - LAMBDA(i+1)*ZC(:,i+1) ||_2 */
/* > where ZC(:,i) = Z(:,i) + sqrt(-1)*Z(:,i+1) */
/* > ZC(:,i+1) = Z(:,i) - sqrt(-1)*Z(:,i+1) */
/* > See the description of REIG, IMEIG and Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] B */
/* > \verbatim */
/* > B (output) REAL M-by-N array. */
/* > IF JOBF =='R', B(1:M,1:K) contains A*U(:,1:K), and can */
/* > be used for computing the refined vectors;
see further */
/* > details in the provided references. */
/* > If JOBF == 'E', B(1:M,1;
K) contains */
/* > A*U(:,1:K)*W(1:K,1:K), which are the vectors from the */
/* > Exact DMD, up to scaling by the inverse eigenvalues. */
/* > If JOBF =='N', then B is not referenced. */
/* > See the descriptions of X, W, K. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB (input) INTEGER, LDB >= M */
/* > The leading dimension of the array B. */
/* > \endverbatim */
/* ..... */
/* > \param[out] W */
/* > \verbatim */
/* > W (workspace/output) REAL N-by-N array */
/* > On exit, W(1:K,1:K) contains the K computed */
/* > eigenvectors of the matrix Rayleigh quotient (real and */
/* > imaginary parts for each complex conjugate pair of the */
/* > eigenvalues). The Ritz vectors (returned in Z) are the */
/* > product of X (containing a POD basis for the input */
/* > matrix X) and W. See the descriptions of K, S, X and Z. */
/* > W is also used as a workspace to temporarily store the */
/* > left singular vectors of X. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDW */
/* > \verbatim */
/* > LDW (input) INTEGER, LDW >= N */
/* > The leading dimension of the array W. */
/* > \endverbatim */
/* ..... */
/* > \param[out] S */
/* > \verbatim */
/* > S (workspace/output) REAL N-by-N array */
/* > The array S(1:K,1:K) is used for the matrix Rayleigh */
/* > quotient. This content is overwritten during */
/* > the eigenvalue decomposition by SGEEV. */
/* > See the description of K. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDS */
/* > \verbatim */
/* > LDS (input) INTEGER, LDS >= N */
/* > The leading dimension of the array S. */
/* > \endverbatim */
/* ..... */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK (workspace/output) REAL LWORK-by-1 array */
/* > On exit, WORK(1:N) contains the singular values of */
/* > X (for JOBS=='N') or column scaled X (JOBS=='S', 'C'). */
/* > If WHTSVD==4, then WORK(N+1) and WORK(N+2) contain */
/* > scaling factor WORK(N+2)/WORK(N+1) used to scale X */
/* > and Y to avoid overflow in the SVD of X. */
/* > This may be of interest if the scaling option is off */
/* > and as many as possible smallest eigenvalues are */
/* > desired to the highest feasible accuracy. */
/* > If the call to SGEDMD is only workspace query, then */
/* > WORK(1) contains the minimal workspace length and */
/* > WORK(2) is the optimal workspace length. Hence, the */
/* > length of work is at least 2. */
/* > See the description of LWORK. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LWORK */
/* > \verbatim */
/* > LWORK (input) INTEGER */
/* > The minimal length of the workspace vector WORK. */
/* > LWORK is calculated as follows: */
/* > If WHTSVD == 1 :: */
/* > If JOBZ == 'V', then */
/* > LWORK >= MAX(2, N + LWORK_SVD, N+MAX(1,4*N)). */
/* > If JOBZ == 'N' then */
/* > LWORK >= MAX(2, N + LWORK_SVD, N+MAX(1,3*N)). */
/* > Here LWORK_SVD = MAX(1,3*N+M,5*N) is the minimal */
/* > workspace length of SGESVD. */
/* > If WHTSVD == 2 :: */
/* > If JOBZ == 'V', then */
/* > LWORK >= MAX(2, N + LWORK_SVD, N+MAX(1,4*N)) */
/* > If JOBZ == 'N', then */
/* > LWORK >= MAX(2, N + LWORK_SVD, N+MAX(1,3*N)) */
/* > Here LWORK_SVD = MAX(M, 5*N*N+4*N)+3*N*N is the */
/* > minimal workspace length of SGESDD. */
/* > If WHTSVD == 3 :: */
/* > If JOBZ == 'V', then */
/* > LWORK >= MAX(2, N+LWORK_SVD,N+MAX(1,4*N)) */
/* > If JOBZ == 'N', then */
/* > LWORK >= MAX(2, N+LWORK_SVD,N+MAX(1,3*N)) */
/* > Here LWORK_SVD = N+M+MAX(3*N+1, */
/* > MAX(1,3*N+M,5*N),MAX(1,N)) */
/* > is the minimal workspace length of SGESVDQ. */
/* > If WHTSVD == 4 :: */
/* > If JOBZ == 'V', then */
/* > LWORK >= MAX(2, N+LWORK_SVD,N+MAX(1,4*N)) */
/* > If JOBZ == 'N', then */
/* > LWORK >= MAX(2, N+LWORK_SVD,N+MAX(1,3*N)) */
/* > Here LWORK_SVD = MAX(7,2*M+N,6*N+2*N*N) is the */
/* > minimal workspace length of SGEJSV. */
/* > The above expressions are not simplified in order to */
/* > make the usage of WORK more transparent, and for */
/* > easier checking. In any case, LWORK >= 2. */
/* > If on entry LWORK = -1, then a workspace query is */
/* > assumed and the procedure only computes the minimal */
/* > and the optimal workspace lengths for both WORK and */
/* > IWORK. See the descriptions of WORK and IWORK. */
/* > \endverbatim */
/* ..... */
/* > \param[out] IWORK */
/* > \verbatim */
/* > IWORK (workspace/output) INTEGER LIWORK-by-1 array */
/* > Workspace that is required only if WHTSVD equals */
/* > 2 , 3 or 4. (See the description of WHTSVD). */
/* > If on entry LWORK =-1 or LIWORK=-1, then the */
/* > minimal length of IWORK is computed and returned in */
/* > IWORK(1). See the description of LIWORK. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LIWORK */
/* > \verbatim */
/* > LIWORK (input) INTEGER */
/* > The minimal length of the workspace vector IWORK. */
/* > If WHTSVD == 1, then only IWORK(1) is used;
LIWORK >=1 */
/* > If WHTSVD == 2, then LIWORK >= MAX(1,8*MIN(M,N)) */
/* > If WHTSVD == 3, then LIWORK >= MAX(1,M+N-1) */
/* > If WHTSVD == 4, then LIWORK >= MAX(3,M+3*N) */
/* > If on entry LIWORK = -1, then a workspace query is */
/* > assumed and the procedure only computes the minimal */
/* > and the optimal workspace lengths for both WORK and */
/* > IWORK. See the descriptions of WORK and IWORK. */
/* > \endverbatim */
/* ..... */
/* > \param[out] INFO */
/* > \verbatim */
/* > INFO (output) INTEGER */
/* > -i < 0 :: On entry, the i-th argument had an */
/* > illegal value */
/* > = 0 :: Successful return. */
/* > = 1 :: Void input. Quick exit (M=0 or N=0). */
/* > = 2 :: The SVD computation of X did not converge. */
/* > Suggestion: Check the input data and/or */
/* > repeat with different WHTSVD. */
/* > = 3 :: The computation of the eigenvalues did not */
/* > converge. */
/* > = 4 :: If data scaling was requested on input and */
/* > the procedure found inconsistency in the data */
/* > such that for some column index i, */
/* > X(:,i) = 0 but Y(:,i) /= 0, then Y(:,i) is set */
/* > to zero if JOBS=='C'. The computation proceeds */
/* > with original or modified data and warning */
/* > flag is set with INFO=4. */
/* > \endverbatim */
/* Authors: */
/* ======== */
/* > \author Zlatko Drmac */
/* > \ingroup gedmd */
/* ............................................................. */
/* ............................................................. */
/* Subroutine */
void sgedmd_(char *jobs, char *jobz, char *jobr, char *jobf, integer *whtsvd, integer *m,
             integer *n, real *x, integer *ldx, real *y, integer *ldy, integer *nrnk, real *tol,
             integer *k, real *reig, real *imeig, real *z__, integer *ldz, real *res, real *b,
             integer *ldb, real *w, integer *ldw, real *s, integer *lds, real *work, integer *lwork,
             integer *iwork, integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgedmd inputs: jobs %c ,jobz %c ,jobr %c ,jobf %c ,whtsvd %" FLA_IS
                      ",m %" FLA_IS ",n %" FLA_IS ",ldx %" FLA_IS ",ldy %" FLA_IS ",nrnk %" FLA_IS
                      ",ldz %" FLA_IS ", ldb %" FLA_IS ",ldw %" FLA_IS ",lds %" FLA_IS
                      ",lwork %" FLA_IS ",liwork %" FLA_IS "",
                      *jobs, *jobz, *jobr, *jobf, *whtsvd, *m, *n, *ldx, *ldy, *nrnk, *ldz, *ldb,
                      *ldw, *lds, *lwork, *liwork);

    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, z_dim1, z_offset, b_dim1, b_offset, w_dim1,
        w_offset, s_dim1, s_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    integer i__, j;
    real ab[4] /* was [2][2] */
        ,
        ofl, ssum;
    integer info1, info2;
    real xscl1, xscl2;
    extern real snrm2_(integer *, real *, integer *);
    real scale;
    extern logical lsame_(char *, char *, integer, integer);
    extern /* Subroutine */
        void
        sscal_(integer *, real *, real *, integer *);
    logical badxy;
    real small_val;
    extern /* Subroutine */
        void
        sgemm_(char *, char *, integer *, integer *, integer *, real *, real *, integer *, real *,
               integer *, real *, real *, integer *),
        sgeev_(char *, char *, integer *, real *, integer *, real *, real *, real *, integer *,
               real *, integer *, real *, integer *, integer *);
    char jobzl[1];
    extern /* Subroutine */
        void
        saxpy_(integer *, real *, real *, integer *, real *, integer *);
    logical wntex;
    extern real slamch_(char *), slange_(char *, integer *, integer *, real *, integer *, real *);
    extern /* Subroutine */
        void
        sgesdd_(char *, integer *, integer *, real *, integer *, real *, real *, integer *, real *,
                integer *, real *, integer *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    char t_or_n__[1];
    extern /* Subroutine */
        void
        slascl_(char *, integer *, integer *, real *, real *, integer *, integer *, real *,
                integer *, integer *);
    extern integer isamax_(integer *, real *, integer *);
    logical sccolx, sccoly;
    extern logical sisnan_(real *);
    extern /* Subroutine */
        void
        sgesvd_(char *, char *, integer *, integer *, real *, integer *, real *, real *, integer *,
                real *, integer *, real *, integer *, integer *);
    integer lwrsdd, mwrsdd;
    extern /* Subroutine */
        void
        sgejsv_(char *, char *, char *, char *, char *, char *, integer *, integer *, real *,
                integer *, real *, real *, integer *, real *, integer *, real *, integer *,
                integer *, integer *),
        slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *);
    integer iminwr;
    logical wntref, wntvec;
    real rootsc;
    integer lwrkev, mlwork, mwrkev, numrnk, olwork;
    real rdummy[2];
    integer lwrsvd, mwrsvd;
    logical lquery, wntres;
    char jsvopt[1];
    extern /* Subroutine */
        void
        slassq_(integer *, real *, integer *, real *, real *);
    integer mwrsvj, lwrsvq, mwrsvq;
    real rdummy2[2];
    extern /* Subroutine */
        void
        sgesvdq_(char *, char *, char *, char *, char *, integer *, integer *, real *, integer *,
                 real *, real *, integer *, real *, integer *, integer *, integer *, integer *,
                 real *, integer *, real *, integer *, integer *);
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 04:19:00 10/24/24 */
    /* ...Switches: */
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by University of -- */
    /* -- Tennessee, University of California Berkeley, University of -- */
    /* -- Colorado Denver and NAG Ltd.. -- */
    /* ..... */
    /* Scalar arguments */
    /* ~~~~~~~~~~~~~~~~ */
    /* Array arguments */
    /* ~~~~~~~~~~~~~~~ */
    /* Parameters */
    /* ~~~~~~~~~~ */
    /* Local scalars */
    /* ~~~~~~~~~~~~~ */
    /* Local arrays */
    /* ~~~~~~~~~~~~ */
    /* External functions (BLAS and LAPACK) */
    /* ~~~~~~~~~~~~~~~~~ */
    /* ............................................................ */
    /* Test the input arguments */
    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    y_dim1 = *ldy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --reig;
    --imeig;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --res;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    w_dim1 = *ldw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    s_dim1 = *lds;
    s_offset = 1 + s_dim1;
    s -= s_offset;
    --work;
    --iwork;
    /* Function Body */
    wntres = lsame_(jobr, "R", 1, 1);
    sccolx = lsame_(jobs, "S", 1, 1) || lsame_(jobs, "C", 1, 1);
    sccoly = lsame_(jobs, "Y", 1, 1);
    wntvec = lsame_(jobz, "V", 1, 1);
    wntref = lsame_(jobf, "R", 1, 1);
    wntex = lsame_(jobf, "E", 1, 1);
    *info = 0;
    lquery = *lwork == -1 || *liwork == -1;
    if(!(sccolx || sccoly || lsame_(jobs, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(wntvec || lsame_(jobz, "N", 1, 1) || lsame_(jobz, "F", 1, 1)))
    {
        *info = -2;
    }
    else if(!(wntres || lsame_(jobr, "N", 1, 1)) || wntres && !wntvec)
    {
        *info = -3;
    }
    else if(!(wntref || wntex || lsame_(jobf, "N", 1, 1)))
    {
        *info = -4;
    }
    else if(!(*whtsvd == 1 || *whtsvd == 2 || *whtsvd == 3 || *whtsvd == 4))
    {
        *info = -5;
    }
    else if(*m < 0)
    {
        *info = -6;
    }
    else if(*n < 0 || *n > *m)
    {
        *info = -7;
    }
    else if(*ldx < *m)
    {
        *info = -9;
    }
    else if(*ldy < *m)
    {
        *info = -11;
    }
    else if(!(*nrnk == -2 || *nrnk == -1 || *nrnk >= 1 && *nrnk <= *n))
    {
        *info = -12;
    }
    else if(*tol < 0.f || *tol >= 1.f)
    {
        *info = -13;
    }
    else if(*ldz < *m)
    {
        *info = -18;
    }
    else if((wntref || wntex) && *ldb < *m)
    {
        *info = -21;
    }
    else if(*ldw < *n)
    {
        *info = -23;
    }
    else if(*lds < *n)
    {
        *info = -25;
    }
    if(*info == 0)
    {
        /* Compute the minimal and the optimal workspace */
        /* requirements. Simulate running the code and */
        /* determine minimal and optimal sizes of the */
        /* workspace at any moment of the run. */
        if(*n == 0)
        {
            /* Quick return. All output except K is void. */
            /* INFO=1 signals the void input. */
            /* In case of a workspace query, the default */
            /* minimal workspace lengths are returned. */
            if(lquery)
            {
                iwork[1] = 1;
                work[1] = 2.f;
                work[2] = 2.f;
            }
            else
            {
                *k = 0;
            }
            *info = 1;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        mlwork = fla_max(2, *n);
        olwork = fla_max(2, *n);
        iminwr = 1;
        if(*whtsvd == 1)
        {
            /* The following is specified as the minimal */
            /* length of WORK in the definition of SGESVD: */
            /* MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) */
            /* Computing MAX */
            /* Computing MAX */
            i__3 = 1;
            i__4 = fla_min(*m, *n) * 3 + fla_max(*m, *n); // , expr subst
            i__1 = fla_max(i__3, i__4);
            i__2 = fla_min(*m, *n) * 5; // , expr subst
            mwrsvd = fla_max(i__1, i__2);
            /* Computing MAX */
            i__1 = mlwork;
            i__2 = *n + mwrsvd; // , expr subst
            mlwork = fla_max(i__1, i__2);
            if(lquery)
            {
                sgesvd_("O", "S", m, n, &x[x_offset], ldx, &work[1], &b[b_offset], ldb,
                        &w[w_offset], ldw, rdummy, &c_n1, &info1);
                /* Computing MAX */
                i__1 = mwrsvd;
                i__2 = (integer)rdummy[0]; // , expr subst
                lwrsvd = fla_max(i__1, i__2);
                /* Computing MAX */
                i__1 = olwork;
                i__2 = *n + lwrsvd; // , expr subst
                olwork = fla_max(i__1, i__2);
            }
        }
        else if(*whtsvd == 2)
        {
            /* The following is specified as the minimal */
            /* length of WORK in the definition of SGESDD: */
            /* MWRSDD = 3*MIN(M,N)*MIN(M,N) + */
            /* MAX( MAX(M,N),5*MIN(M,N)*MIN(M,N)+4*MIN(M,N) ) */
            /* IMINWR = 8*MIN(M,N) */
            /* Computing MAX */
            i__1 = fla_max(*m, *n);
            i__2 = fla_min(*m, *n) * 5 * fla_min(*m, *n) + (fla_min(*m, *n) << 2); // , expr subst
            mwrsdd = fla_min(*m, *n) * 3 * fla_min(*m, *n) + fla_max(i__1, i__2);
            /* Computing MAX */
            i__1 = mlwork;
            i__2 = *n + mwrsdd; // , expr subst
            mlwork = fla_max(i__1, i__2);
            iminwr = fla_min(*m, *n) << 3;
            if(lquery)
            {
                sgesdd_("O", m, n, &x[x_offset], ldx, &work[1], &b[b_offset], ldb, &w[w_offset],
                        ldw, rdummy, &c_n1, &iwork[1], &info1);
                /* Computing MAX */
                i__1 = mwrsdd;
                i__2 = (integer)rdummy[0]; // , expr subst
                lwrsdd = fla_max(i__1, i__2);
                /* Computing MAX */
                i__1 = olwork;
                i__2 = *n + lwrsdd; // , expr subst
                olwork = fla_max(i__1, i__2);
            }
        }
        else if(*whtsvd == 3)
        {
            /* LWQP3 = 3*N+1 */
            /* LWORQ = MAX(N, 1) */
            /* MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) */
            /* MWRSVQ = N + MAX( LWQP3, MWRSVD, LWORQ )+ MAX(M,2) */
            /* MLWORK = N + MWRSVQ */
            /* IMINWR = M+N-1 */
            sgesvdq_("H", "P", "N", "R", "R", m, n, &x[x_offset], ldx, &work[1], &z__[z_offset],
                     ldz, &w[w_offset], ldw, &numrnk, &iwork[1], &c_n1, rdummy, &c_n1, rdummy2,
                     &c_n1, &info1);
            iminwr = iwork[1];
            mwrsvq = (integer)rdummy[1];
            /* Computing MAX */
            i__1 = mlwork;
            i__2 = *n + mwrsvq + (integer)rdummy2[0]; // , expr subst
            mlwork = fla_max(i__1, i__2);
            if(lquery)
            {
                lwrsvq = (integer)rdummy[0];
                /* Computing MAX */
                i__1 = olwork;
                i__2 = *n + lwrsvq + (integer)rdummy2[0]; // , expr subst
                olwork = fla_max(i__1, i__2);
            }
        }
        else if(*whtsvd == 4)
        {
            *(unsigned char *)jsvopt = 'J';
            /* MWRSVJ = MAX( 7, 2*M+N, 6*N+2*N*N )! for JSVOPT='V' */
            /* Computing MAX */
            /* Computing MAX */
            /* Computing MAX */
            i__5 = 7;
            i__6 = (*m << 1) + *n; // , expr subst
            i__3 = fla_max(i__5, i__6);
            i__4 = (*n << 2) + *n * *n; // , expr subst
            i__1 = fla_max(i__3, i__4);
            i__2 = (*n << 1) + *n * *n + 6; // , expr subst
            mwrsvj = fla_max(i__1, i__2);
            /* Computing MAX */
            i__1 = mlwork;
            i__2 = *n + mwrsvj; // , expr subst
            mlwork = fla_max(i__1, i__2);
            /* Computing MAX */
            i__1 = 3;
            i__2 = *m + *n * 3; // , expr subst
            iminwr = fla_max(i__1, i__2);
            if(lquery)
            {
                /* Computing MAX */
                i__1 = olwork;
                i__2 = *n + mwrsvj; // , expr subst
                olwork = fla_max(i__1, i__2);
            }
        }
        if(wntvec || wntex || lsame_(jobz, "F", 1, 1))
        {
            *(unsigned char *)jobzl = 'V';
        }
        else
        {
            *(unsigned char *)jobzl = 'N';
        }
        /* Workspace calculation to the SGEEV call */
        if(lsame_(jobzl, "V", 1, 1))
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = *n << 2; // , expr subst
            mwrkev = fla_max(i__1, i__2);
        }
        else
        {
            /* Computing MAX */
            i__1 = 1;
            i__2 = *n * 3; // , expr subst
            mwrkev = fla_max(i__1, i__2);
        }
        /* Computing MAX */
        i__1 = mlwork;
        i__2 = *n + mwrkev; // , expr subst
        mlwork = fla_max(i__1, i__2);
        if(lquery)
        {
            sgeev_("N", jobzl, n, &s[s_offset], lds, &reig[1], &imeig[1], &w[w_offset], ldw,
                   &w[w_offset], ldw, rdummy, &c_n1, &info1);
            /* Computing MAX */
            i__1 = mwrkev;
            i__2 = (integer)rdummy[0]; // , expr subst
            lwrkev = fla_max(i__1, i__2);
            /* Computing MAX */
            i__1 = olwork;
            i__2 = *n + lwrkev; // , expr subst
            olwork = fla_max(i__1, i__2);
        }
        if(*liwork < iminwr && !lquery)
        {
            *info = -29;
        }
        if(*lwork < mlwork && !lquery)
        {
            *info = -27;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGEDMD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        /* Return minimal and optimal workspace sizes */
        iwork[1] = iminwr;
        work[1] = (real)mlwork;
        work[2] = (real)olwork;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ............................................................ */
    ofl = slamch_("O");
    small_val = slamch_("S");
    badxy = FALSE_;
    /* <1> Optional scaling of the snapshots (columns of X, Y) */
    /* ========================================================== */
    if(sccolx)
    {
        /* The columns of X will be normalized. */
        /* To prevent overflows, the column norms of X are */
        /* carefully computed using SLASSQ. */
        *k = 0;
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* WORK(i) = DNRM2( M, X(1,i), 1 ) */
            scale = 0.f;
            slassq_(m, &x[i__ * x_dim1 + 1], &c__1, &scale, &ssum);
            if(sisnan_(&scale) || sisnan_(&ssum))
            {
                *k = 0;
                *info = -8;
                i__2 = -(*info);
                xerbla_("SGEDMD", &i__2, (ftnlen)6);
            }
            if(scale != 0.f && ssum != 0.f)
            {
                rootsc = sqrt(ssum);
                if(scale >= ofl / rootsc)
                {
                    /* Norm of X(:,i) overflows. First, X(:,i) */
                    /* is scaled by */
                    /* ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2. */
                    /* Next, the norm of X(:,i) is stored without */
                    /* overflow as WORK(i) = - SCALE * (ROOTSC/M), */
                    /* the minus sign indicating the 1/M factor. */
                    /* Scaling is performed without overflow, and */
                    /* underflow may occur in the smallest entries */
                    /* of X(:,i). The relative backward and forward */
                    /* errors are small in the ell_2 norm. */
                    r__1 = 1.f / rootsc;
                    slascl_("G", &c__0, &c__0, &scale, &r__1, m, &c__1, &x[i__ * x_dim1 + 1], m,
                            &info2);
                    work[i__] = -scale * (rootsc / (real)(*m));
                }
                else
                {
                    /* X(:,i) will be scaled to unit 2-norm */
                    work[i__] = scale * rootsc;
                    slascl_("G", &c__0, &c__0, &work[i__], &c_b46, m, &c__1, &x[i__ * x_dim1 + 1],
                            m, &info2);
                    /* X(1:M,i) = (ONE/WORK(i)) * X(1:M,i) ! INTRINSIC */
                    /* LAPACK CALL */
                }
            }
            else
            {
                work[i__] = 0.f;
                ++(*k);
            }
        }
        if(*k == *n)
        {
            /* All columns of X are zero. Return error code -8. */
            /* (the 8th input variable had an illegal value) */
            *k = 0;
            *info = -8;
            i__1 = -(*info);
            xerbla_("SGEDMD", &i__1, (ftnlen)6);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Now, apply the same scaling to the columns of Y. */
            if(work[i__] > 0.f)
            {
                /* ! BLAS CALL */
                r__1 = 1.f / work[i__];
                sscal_(m, &r__1, &y[i__ * y_dim1 + 1], &c__1);
                /* Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i) ! INTRINSIC */
            }
            else if(work[i__] < 0.f)
            {
                r__1 = -work[i__];
                r__2 = 1.f / (real)(*m);
                slascl_("G", &c__0, &c__0, &r__1, &r__2, m, &c__1, &y[i__ * y_dim1 + 1], m, &info2);
                /* LAPACK CALL */
            }
            else if(y[isamax_(m, &y[i__ * y_dim1 + 1], &c__1) + i__ * y_dim1] != 0.f)
            {
                /* X(:,i) is zero vector. For consistency, */
                /* Y(:,i) should also be zero. If Y(:,i) is not */
                /* zero, then the data might be inconsistent or */
                /* corrupted. If JOBS == 'C', Y(:,i) is set to */
                /* zero and a warning flag is raised. */
                /* The computation continues but the */
                /* situation will be reported in the output. */
                badxy = TRUE_;
                /* ! BLAS CALL */
                if(lsame_(jobs, "C", 1, 1))
                {
                    sscal_(m, &c_b56, &y[i__ * y_dim1 + 1], &c__1);
                }
            }
        }
    }
    if(sccoly)
    {
        /* The columns of Y will be normalized. */
        /* To prevent overflows, the column norms of Y are */
        /* carefully computed using SLASSQ. */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* WORK(i) = DNRM2( M, Y(1,i), 1 ) */
            scale = 0.f;
            slassq_(m, &y[i__ * y_dim1 + 1], &c__1, &scale, &ssum);
            if(sisnan_(&scale) || sisnan_(&ssum))
            {
                *k = 0;
                *info = -10;
                i__2 = -(*info);
                xerbla_("SGEDMD", &i__2, (ftnlen)6);
            }
            if(scale != 0.f && ssum != 0.f)
            {
                rootsc = sqrt(ssum);
                if(scale >= ofl / rootsc)
                {
                    /* Norm of Y(:,i) overflows. First, Y(:,i) */
                    /* is scaled by */
                    /* ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2. */
                    /* Next, the norm of Y(:,i) is stored without */
                    /* overflow as WORK(i) = - SCALE * (ROOTSC/M), */
                    /* the minus sign indicating the 1/M factor. */
                    /* Scaling is performed without overflow, and */
                    /* underflow may occur in the smallest entries */
                    /* of Y(:,i). The relative backward and forward */
                    /* errors are small in the ell_2 norm. */
                    r__1 = 1.f / rootsc;
                    slascl_("G", &c__0, &c__0, &scale, &r__1, m, &c__1, &y[i__ * y_dim1 + 1], m,
                            &info2);
                    work[i__] = -scale * (rootsc / (real)(*m));
                }
                else
                {
                    /* X(:,i) will be scaled to unit 2-norm */
                    work[i__] = scale * rootsc;
                    slascl_("G", &c__0, &c__0, &work[i__], &c_b46, m, &c__1, &y[i__ * y_dim1 + 1],
                            m, &info2);
                    /* Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i) ! INTRINSIC */
                    /* LAPACK CALL */
                }
            }
            else
            {
                work[i__] = 0.f;
            }
        }
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Now, apply the same scaling to the columns of X. */
            if(work[i__] > 0.f)
            {
                /* ! BLAS CALL */
                r__1 = 1.f / work[i__];
                sscal_(m, &r__1, &x[i__ * x_dim1 + 1], &c__1);
                /* X(1:M,i) = (ONE/WORK(i)) * X(1:M,i) ! INTRINSIC */
            }
            else if(work[i__] < 0.f)
            {
                r__1 = -work[i__];
                r__2 = 1.f / (real)(*m);
                slascl_("G", &c__0, &c__0, &r__1, &r__2, m, &c__1, &x[i__ * x_dim1 + 1], m, &info2);
                /* LAPACK CALL */
            }
            else if(x[isamax_(m, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1] != 0.f)
            {
                /* Y(:,i) is zero vector. If X(:,i) is not */
                /* zero, then a warning flag is raised. */
                /* The computation continues but the */
                /* situation will be reported in the output. */
                badxy = TRUE_;
            }
        }
    }
    /* <2> SVD of the data snapshot matrix X. */
    /* ===================================== */
    /* The left singular vectors are stored in the array X. */
    /* The right singular vectors are in the array W. */
    /* The array W will later on contain the eigenvectors */
    /* of a Rayleigh quotient. */
    numrnk = *n;
    if(*whtsvd == 1)
    {
        i__1 = *lwork - *n;
        sgesvd_("O", "S", m, n, &x[x_offset], ldx, &work[1], &b[b_offset], ldb, &w[w_offset], ldw,
                &work[*n + 1], &i__1, &info1);
        /* LAPACK CALL */
        *(unsigned char *)t_or_n__ = 'T';
    }
    else if(*whtsvd == 2)
    {
        i__1 = *lwork - *n;
        sgesdd_("O", m, n, &x[x_offset], ldx, &work[1], &b[b_offset], ldb, &w[w_offset], ldw,
                &work[*n + 1], &i__1, &iwork[1], &info1);
        /* LAPACK CALL */
        *(unsigned char *)t_or_n__ = 'T';
    }
    else if(*whtsvd == 3)
    {
        /* ! LAPACK CALL */
        i__1 = *lwork - *n - fla_max(2, *m);
        i__2 = fla_max(2, *m);
        sgesvdq_("H", "P", "N", "R", "R", m, n, &x[x_offset], ldx, &work[1], &z__[z_offset], ldz,
                 &w[w_offset], ldw, &numrnk, &iwork[1], liwork, &work[*n + fla_max(2, *m) + 1],
                 &i__1, &work[*n + 1], &i__2, &info1);
        /* ! LAPACK CALL */
        slacpy_("A", m, &numrnk, &z__[z_offset], ldz, &x[x_offset], ldx);
        *(unsigned char *)t_or_n__ = 'T';
    }
    else if(*whtsvd == 4)
    {
        /* ! LAPACK CALL */
        i__1 = *lwork - *n;
        sgejsv_("F", "U", jsvopt, "N", "N", "P", m, n, &x[x_offset], ldx, &work[1], &z__[z_offset],
                ldz, &w[w_offset], ldw, &work[*n + 1], &i__1, &iwork[1], &info1);
        slacpy_("A", m, n, &z__[z_offset], ldz, &x[x_offset], ldx);
        /* LAPACK CALL */
        *(unsigned char *)t_or_n__ = 'N';
        xscl1 = work[*n + 1];
        xscl2 = work[*n + 2];
        /* This is an exceptional situation. If the */
        /* data matrices are not scaled and the */
        /* largest singular value of X overflows. */
        /* In that case SGEJSV can return the SVD */
        /* in scaled form. The scaling factor can be used */
        /* to rescale the data (X and Y). */
        if(xscl1 != xscl2)
        {
            slascl_("G", &c__0, &c__0, &xscl1, &xscl2, m, n, &y[y_offset], ldy, &info2);
        }
    }
    if(info1 > 0)
    {
        /* The SVD selected subroutine did not converge. */
        /* Return with an error code. */
        *info = 2;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    if(work[1] == 0.f)
    {
        /* The largest computed singular value of (scaled) */
        /* X is zero. Return error code -8 */
        /* (the 8th input variable had an illegal value). */
        *k = 0;
        *info = -8;
        i__1 = -(*info);
        xerbla_("SGEDMD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* <3> Determine the numerical rank of the data */
    /* snapshots matrix X. This depends on the */
    /* parameters NRNK and TOL. */
    if(*nrnk == -1)
    {
        *k = 1;
        i__1 = numrnk;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            if(work[i__] <= work[1] * *tol || work[i__] <= small_val)
            {
                goto L2;
            }
            ++(*k);
        }
    L2:;
    }
    else if(*nrnk == -2)
    {
        *k = 1;
        i__1 = numrnk - 1;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            if(work[i__ + 1] <= work[i__] * *tol || work[i__] <= small_val)
            {
                goto L3;
            }
            ++(*k);
        }
    L3:;
    }
    else
    {
        *k = 1;
        i__1 = *nrnk;
        for(i__ = 2; i__ <= i__1; ++i__)
        {
            if(work[i__] <= small_val)
            {
                goto L4;
            }
            ++(*k);
        }
    L4:;
    }
    /* Now, U = X(1:M,1:K) is the SVD/POD basis for the */
    /* snapshot data in the input matrix X. */
    /* <4> Compute the Rayleigh quotient S = U^T * A * U. */
    /* Depending on the requested outputs, the computation */
    /* is organized to compute additional auxiliary */
    /* matrices (for the residuals and refinements). */
    /* In all formulas below, we need V_k*Sigma_k^(-1) */
    /* where either V_k is in W(1:N,1:K), or V_k^T is in */
    /* W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)). */
    if(lsame_(t_or_n__, "N", 1, 1))
    {
        i__1 = *k;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* ! BLAS CALL */
            r__1 = 1.f / work[i__];
            sscal_(n, &r__1, &w[i__ * w_dim1 + 1], &c__1);
            /* W(1:N,i) = (ONE/WORK(i)) * W(1:N,i) ! INTRINSIC */
        }
    }
    else
    {
        /* This non-unit stride access is due to the fact */
        /* that SGESVD, SGESVDQ and SGESDD return the */
        /* transposed matrix of the right singular vectors. */
        /* DO i = 1, K */
        /* CALL SSCAL( N, ONE/WORK(i), W(i,1), LDW ) ! BLAS CALL */
        /* ! W(i,1:N) = (ONE/WORK(i)) * W(i,1:N) ! INTRINSIC */
        /* END DO */
        i__1 = *k;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            work[*n + i__] = 1.f / work[i__];
        }
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *k;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                w[i__ + j * w_dim1] = work[*n + i__] * w[i__ + j * w_dim1];
            }
        }
    }
    if(wntref)
    {
        /* Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K))) */
        /* for computing the refined Ritz vectors */
        /* (optionally, outside SGEDMD). */
        sgemm_("N", t_or_n__, m, k, n, &c_b46, &y[y_offset], ldy, &w[w_offset], ldw, &c_b56,
               &z__[z_offset], ldz);
        /* Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(W(1:K,1:N))) ! INTRI */
        /* Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K)) ! INTRI */
        /* At this point Z contains */
        /* A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and */
        /* this is needed for computing the residuals. */
        /* This matrix is returned in the array B and */
        /* it can be used to compute refined Ritz vectors. */
        /* BLAS CALL */
        slacpy_("A", m, k, &z__[z_offset], ldz, &b[b_offset], ldb);
        /* B(1:M,1:K) = Z(1:M,1:K) ! INTRINSIC */
        /* BLAS CALL */
        sgemm_("T", "N", k, k, m, &c_b46, &x[x_offset], ldx, &z__[z_offset], ldz, &c_b56,
               &s[s_offset], lds);
        /* S(1:K,1:K) = MATMUL(TANSPOSE(X(1:M,1:K)),Z(1:M,1:K)) ! INTRI */
        /* At this point S = U^T * A * U is the Rayleigh quotient. */
        /* BLAS CALL */
    }
    else
    {
        /* A * U(:,1:K) is not explicitly needed and the */
        /* computation is organized differently. The Rayleigh */
        /* quotient is computed more efficiently. */
        sgemm_("T", "N", k, n, m, &c_b46, &x[x_offset], ldx, &y[y_offset], ldy, &c_b56,
               &z__[z_offset], ldz);
        /* Z(1:K,1:N) = MATMUL( TRANSPOSE(X(1:M,1:K)), Y(1:M,1:N) ) ! IN */
        /* In the two SGEMM calls here, can use K for LDZ */
        /* BLAS CALL */
        sgemm_("N", t_or_n__, k, k, n, &c_b46, &z__[z_offset], ldz, &w[w_offset], ldw, &c_b56,
               &s[s_offset], lds);
        /* S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(W(1:K,1:N))) ! INTRIN */
        /* S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K))) ! INTRIN */
        /* At this point S = U^T * A * U is the Rayleigh quotient. */
        /* If the residuals are requested, save scaled V_k into Z. */
        /* Recall that V_k or V_k^T is stored in W. */
        /* BLAS CALL */
        if(wntres || wntex)
        {
            if(lsame_(t_or_n__, "N", 1, 1))
            {
                slacpy_("A", n, k, &w[w_offset], ldw, &z__[z_offset], ldz);
            }
            else
            {
                slacpy_("A", k, n, &w[w_offset], ldw, &z__[z_offset], ldz);
            }
        }
    }
    /* <5> Compute the Ritz values and (if requested) the */
    /* right eigenvectors of the Rayleigh quotient. */
    i__1 = *lwork - *n;
    sgeev_("N", jobzl, k, &s[s_offset], lds, &reig[1], &imeig[1], &w[w_offset], ldw, &w[w_offset],
           ldw, &work[*n + 1], &i__1, &info1);
    /* W(1:K,1:K) contains the eigenvectors of the Rayleigh */
    /* quotient. Even in the case of complex spectrum, all */
    /* computation is done in real arithmetic. REIG and */
    /* IMEIG are the real and the imaginary parts of the */
    /* eigenvalues, so that the spectrum is given as */
    /* REIG(:) + sqrt(-1)*IMEIG(:). Complex conjugate pairs */
    /* are listed at consecutive positions. For such a */
    /* complex conjugate pair of the eigenvalues, the */
    /* corresponding eigenvectors are also a complex */
    /* conjugate pair with the real and imaginary parts */
    /* stored column-wise in W at the corresponding */
    /* consecutive column indices. See the description of Z. */
    /* Also, see the description of SGEEV. */
    /* LAPACK CALL */
    if(info1 > 0)
    {
        /* SGEEV failed to compute the eigenvalues and */
        /* eigenvectors of the Rayleigh quotient. */
        *info = 3;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* <6> Compute the eigenvectors (if requested) and, */
    /* the residuals (if requested). */
    if(wntvec || wntex)
    {
        if(wntres)
        {
            if(wntref)
            {
                /* Here, if the refinement is requested, we have */
                /* A*U(:,1:K) already computed and stored in Z. */
                /* For the residuals, need Y = A * U(:,1;
               K) * W. */
                sgemm_("N", "N", m, k, k, &c_b46, &z__[z_offset], ldz, &w[w_offset], ldw, &c_b56,
                       &y[y_offset], ldy);
                /* Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K) ! INTRINSIC */
                /* This frees Z;
                Y contains A * U(:,1:K) * W. */
                /* BLAS CALL */
            }
            else
            {
                /* Compute S = V_k * Sigma_k^(-1) * W, where */
                /* V_k * Sigma_k^(-1) is stored in Z */
                sgemm_(t_or_n__, "N", n, k, k, &c_b46, &z__[z_offset], ldz, &w[w_offset], ldw,
                       &c_b56, &s[s_offset], lds);
                /* Then, compute Z = Y * S = */
                /* = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) = */
                /* = A * U(:,1:K) * W(1:K,1:K) */
                sgemm_("N", "N", m, k, n, &c_b46, &y[y_offset], ldy, &s[s_offset], lds, &c_b56,
                       &z__[z_offset], ldz);
                /* Save a copy of Z into Y and free Z for holding */
                /* the Ritz vectors. */
                slacpy_("A", m, k, &z__[z_offset], ldz, &y[y_offset], ldy);
                if(wntex)
                {
                    slacpy_("A", m, k, &z__[z_offset], ldz, &b[b_offset], ldb);
                }
            }
        }
        else if(wntex)
        {
            /* Compute S = V_k * Sigma_k^(-1) * W, where */
            /* V_k * Sigma_k^(-1) is stored in Z */
            sgemm_(t_or_n__, "N", n, k, k, &c_b46, &z__[z_offset], ldz, &w[w_offset], ldw, &c_b56,
                   &s[s_offset], lds);
            /* Then, compute Z = Y * S = */
            /* = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) = */
            /* = A * U(:,1:K) * W(1:K,1:K) */
            sgemm_("N", "N", m, k, n, &c_b46, &y[y_offset], ldy, &s[s_offset], lds, &c_b56,
                   &b[b_offset], ldb);
            /* The above call replaces the following two calls */
            /* that were used in the developing-testing phase. */
            /* CALL SGEMM( 'N', 'N', M, K, N, ONE, Y, LDY, S, & */
            /* LDS, ZERO, Z, LDZ) */
            /* Save a copy of Z into B and free Z for holding */
            /* the Ritz vectors. */
            /* CALL SLACPY( 'A', M, K, Z, LDZ, B, LDB ) */
        }
        /* Compute the real form of the Ritz vectors */
        if(wntvec)
        {
            sgemm_("N", "N", m, k, k, &c_b46, &x[x_offset], ldx, &w[w_offset], ldw, &c_b56,
                   &z__[z_offset], ldz);
        }
        /* Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K)) ! INTRINSIC */
        /* BLAS CALL */
        if(wntres)
        {
            i__ = 1;
        L5:
            if(i__ <= *k)
            {
                if(imeig[i__] == 0.f)
                {
                    /* have a real eigenvalue with real eigenvector */
                    /* ! BLAS CALL */
                    r__1 = -reig[i__];
                    saxpy_(m, &r__1, &z__[i__ * z_dim1 + 1], &c__1, &y[i__ * y_dim1 + 1], &c__1);
                    /* Y(1:M,i) = Y(1:M,i) - REIG(i) * Z(1:M,i) ! */
                    res[i__] = snrm2_(m, &y[i__ * y_dim1 + 1], &c__1);
                    /* BLAS CALL */
                    ++i__;
                }
                else
                {
                    /* Have a complex conjugate pair */
                    /* REIG(i) +- sqrt(-1)*IMEIG(i). */
                    /* Since all computation is done in real */
                    /* arithmetic, the formula for the residual */
                    /* is recast for real representation of the */
                    /* complex conjugate eigenpair. See the */
                    /* description of RES. */
                    ab[0] = reig[i__];
                    ab[1] = -imeig[i__];
                    ab[2] = imeig[i__];
                    ab[3] = reig[i__];
                    sgemm_("N", "N", m, &c__2, &c__2, &c_b154, &z__[i__ * z_dim1 + 1], ldz, ab,
                           &c__2, &c_b46, &y[i__ * y_dim1 + 1], ldy);
                    /* Y(1:M,i:i+1) = Y(1:M,i:i+1) - Z(1:M,i:i+1) * AB ! INT */
                    /* ! LAPACK CALL */
                    /* BLAS CALL */
                    res[i__] = slange_("F", m, &c__2, &y[i__ * y_dim1 + 1], ldy, &work[*n + 1]);
                    res[i__ + 1] = res[i__];
                    i__ += 2;
                }
                goto L5;
            }
        }
    }
    if(*whtsvd == 4)
    {
        work[*n + 1] = xscl1;
        work[*n + 2] = xscl2;
    }
    /* Successful exit. */
    if(!badxy)
    {
        *info = 0;
    }
    else
    {
        /* A warning on possible data inconsistency. */
        /* This should be a rare event. */
        *info = 4;
    }
    /* ............................................................ */
    AOCL_DTL_TRACE_LOG_EXIT
    return;
    /* ...... */
}
/* sgedmd_ */
