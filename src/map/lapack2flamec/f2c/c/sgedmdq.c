/* ./sgedmdq.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static integer c_n1 = -1;
static real c_b30 = 0.f;
/* > \brief \b SGEDMDQ computes the Dynamic Mode Decomposition (DMD) for a pair of data snapshot
 * matrices. */
/* =========== DOCUMENTATION =========== */
/* Definition: */
/* =========== */
/* SUBROUTINE SGEDMDQ( JOBS, JOBZ, JOBR, JOBQ, JOBT, JOBF, & */
/* WHTSVD, M, N, F, LDF, X, LDX, Y, & */
/* LDY, NRNK, TOL, K, REIG, IMEIG, & */
/* Z, LDZ, RES, B, LDB, V, LDV, & */
/* S, LDS, WORK, LWORK, IWORK, LIWORK, INFO ) */
/* ..... */
/* USE iso_fortran_env */
/* IMPLICIT NONE */
/* INTEGER, PARAMETER :: WP = real32 */
/* ..... */
/* Scalar arguments */
/* CHARACTER, INTENT(IN) :: JOBS, JOBZ, JOBR, JOBQ, & */
/* JOBT, JOBF */
/* INTEGER, INTENT(IN) :: WHTSVD, M, N, LDF, LDX, & */
/* LDY, NRNK, LDZ, LDB, LDV, & */
/* LDS, LWORK, LIWORK */
/* INTEGER, INTENT(OUT) :: INFO, K */
/* REAL, INTENT(IN) :: TOL */
/* Array arguments */
/* REAL, INTENT(INOUT) :: F(LDF,*) */
/* REAL, INTENT(OUT) :: X(LDX,*), Y(LDY,*), & */
/* Z(LDZ,*), B(LDB,*), & */
/* V(LDV,*), S(LDS,*) */
/* REAL, INTENT(OUT) :: REIG(*), IMEIG(*), & */
/* RES(*) */
/* REAL, INTENT(OUT) :: WORK(*) */
/* INTEGER, INTENT(OUT) :: IWORK(*) */
/* ............................................................ */
/* > \par Purpose: */
/* ============= */
/* > \verbatim */
/* > SGEDMDQ computes the Dynamic Mode Decomposition (DMD) for */
/* > a pair of data snapshot matrices, using a QR factorization */
/* > based compression of the data. For the input matrices */
/* > X and Y such that Y = A*X with an unaccessible matrix */
/* > A, SGEDMDQ computes a certain number of Ritz pairs of A using */
/* > the standard Rayleigh-Ritz extraction from a subspace of */
/* > range(X) that is determined using the leading left singular */
/* > vectors of X. Optionally, SGEDMDQ returns the residuals */
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
/* > Program Office. */
/* > \endverbatim */
/* ...................................................................... */
/* > \par Distribution Statement A: */
/* ============================== */
/* > \verbatim */
/* > Approved for Public Release, Distribution Unlimited. */
/* > Cleared by DARPA on September 29, 2022 */
/* > \endverbatim */
/* ...................................................................... */
/* Arguments */
/* ========= */
/* > \param[in] JOBS */
/* > \verbatim */
/* > JOBS (input) CHARACTER*1 */
/* > Determines whether the initial data snapshots are scaled */
/* > by a diagonal matrix. The data snapshots are the columns */
/* > of F. The leading N-1 columns of F are denoted X and the */
/* > trailing N-1 columns are denoted Y. */
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
/* > in factored form as the product Z*V, where Z */
/* > is orthonormal and V contains the eigenvectors */
/* > of the corresponding Rayleigh quotient. */
/* > See the descriptions of F, V, Z. */
/* > 'Q' :: The eigenvectors (Koopman modes) will be returned */
/* > in factored form as the product Q*Z, where Z */
/* > contains the eigenvectors of the compression of the */
/* > underlying discretized operator onto the span of */
/* > the data snapshots. See the descriptions of F, V, Z. */
/* > Q is from the initial QR factorization. */
/* > 'N' :: The eigenvectors are not computed. */
/* > \endverbatim */
/* ..... */
/* > \param[in] JOBR */
/* > \verbatim */
/* > JOBR (input) CHARACTER*1 */
/* > Determines whether to compute the residuals. */
/* > 'R' :: The residuals for the computed eigenpairs will */
/* > be computed and stored in the array RES. */
/* > See the description of RES. */
/* > For this option to be legal, JOBZ must be 'V'. */
/* > 'N' :: The residuals are not computed. */
/* > \endverbatim */
/* ..... */
/* > \param[in] JOBQ */
/* > \verbatim */
/* > JOBQ (input) CHARACTER*1 */
/* > Specifies whether to explicitly compute and return the */
/* > orthogonal matrix from the QR factorization. */
/* > 'Q' :: The matrix Q of the QR factorization of the data */
/* > snapshot matrix is computed and stored in the */
/* > array F. See the description of F. */
/* > 'N' :: The matrix Q is not explicitly computed. */
/* > \endverbatim */
/* ..... */
/* > \param[in] JOBT */
/* > \verbatim */
/* > JOBT (input) CHARACTER*1 */
/* > Specifies whether to return the upper triangular factor */
/* > from the QR factorization. */
/* > 'R' :: The matrix R of the QR factorization of the data */
/* > snapshot matrix F is returned in the array Y. */
/* > See the description of Y and Further details. */
/* > 'N' :: The matrix R is not returned. */
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
/* > To be useful on exit, this option needs JOBQ='Q'. */
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
/* > M (input) INTEGER, M >= 0 */
/* > The state space dimension (the number of rows of F) */
/* > \endverbatim */
/* ..... */
/* > \param[in] N */
/* > \verbatim */
/* > N (input) INTEGER, 0 <= N <= M */
/* > The number of data snapshots from a single trajectory, */
/* > taken at equidistant discrete times. This is the */
/* > number of columns of F. */
/* > \endverbatim */
/* ..... */
/* > \param[in,out] F */
/* > \verbatim */
/* > F (input/output) REAL M-by-N array */
/* > > On entry, */
/* > the columns of F are the sequence of data snapshots */
/* > from a single trajectory, taken at equidistant discrete */
/* > times. It is assumed that the column norms of F are */
/* > in the range of the normalized floating point numbers. */
/* > < On exit, */
/* > If JOBQ == 'Q', the array F contains the orthogonal */
/* > matrix/factor of the QR factorization of the initial */
/* > data snapshots matrix F. See the description of JOBQ. */
/* > If JOBQ == 'N', the entries in F strictly below the main */
/* > diagonal contain, column-wise, the information on the */
/* > Householder vectors, as returned by SGEQRF. The */
/* > remaining information to restore the orthogonal matrix */
/* > of the initial QR factorization is stored in WORK(1:N). */
/* > See the description of WORK. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDF */
/* > \verbatim */
/* > LDF (input) INTEGER, LDF >= M */
/* > The leading dimension of the array F. */
/* > \endverbatim */
/* ..... */
/* > \param[in,out] X */
/* > \verbatim */
/* > X (workspace/output) REAL MIN(M,N)-by-(N-1) array */
/* > X is used as workspace to hold representations of the */
/* > leading N-1 snapshots in the orthonormal basis computed */
/* > in the QR factorization of F. */
/* > On exit, the leading K columns of X contain the leading */
/* > K left singular vectors of the above described content */
/* > of X. To lift them to the space of the left singular */
/* > vectors U(:,1:K)of the input data, pre-multiply with the */
/* > Q factor from the initial QR factorization. */
/* > See the descriptions of F, K, V and Z. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDX */
/* > \verbatim */
/* > LDX (input) INTEGER, LDX >= N */
/* > The leading dimension of the array X */
/* > \endverbatim */
/* ..... */
/* > \param[in,out] Y */
/* > \verbatim */
/* > Y (workspace/output) REAL MIN(M,N)-by-(N-1) array */
/* > Y is used as workspace to hold representations of the */
/* > trailing N-1 snapshots in the orthonormal basis computed */
/* > in the QR factorization of F. */
/* > On exit, */
/* > If JOBT == 'R', Y contains the MIN(M,N)-by-N upper */
/* > triangular factor from the QR factorization of the data */
/* > snapshot matrix F. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDY */
/* > \verbatim */
/* > LDY (input) INTEGER , LDY >= N */
/* > The leading dimension of the array Y */
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
/* > 0 < NRNK <= N-1 :: at most NRNK largest singular values */
/* > will be used. If the number of the computed nonzero */
/* > singular values is less than NRNK, then only those */
/* > nonzero values will be used and the actually used */
/* > dimension is less than NRNK. The actual number of */
/* > the nonzero singular values is returned in the variable */
/* > K. See the description of K. */
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
/* > The dimension of the SVD/POD basis for the leading N-1 */
/* > data snapshots (columns of F) and the number of the */
/* > computed Ritz pairs. The value of K is determined */
/* > according to the rule set by the parameters NRNK and */
/* > TOL. See the descriptions of NRNK and TOL. */
/* > \endverbatim */
/* ..... */
/* > \param[out] REIG */
/* > \verbatim */
/* > REIG (output) REAL (N-1)-by-1 array */
/* > The leading K (K<=N) entries of REIG contain */
/* > the real parts of the computed eigenvalues */
/* > REIG(1:K) + sqrt(-1)*IMEIG(1:K). */
/* > See the descriptions of K, IMEIG, Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] IMEIG */
/* > \verbatim */
/* > IMEIG (output) REAL (N-1)-by-1 array */
/* > The leading K (K<N) entries of REIG contain */
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
/* > See the descriptions of K, REIG, Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] Z */
/* > \verbatim */
/* > Z (workspace/output) REAL M-by-(N-1) array */
/* > If JOBZ =='V' then */
/* > Z contains real Ritz vectors as follows: */
/* > If IMEIG(i)=0, then Z(:,i) is an eigenvector of */
/* > the i-th Ritz value. */
/* > If IMEIG(i) > 0 (and IMEIG(i+1) < 0) then */
/* > [Z(:,i) Z(:,i+1)] span an invariant subspace and */
/* > the Ritz values extracted from this subspace are */
/* > REIG(i) + sqrt(-1)*IMEIG(i) and */
/* > REIG(i) - sqrt(-1)*IMEIG(i). */
/* > The corresponding eigenvectors are */
/* > Z(:,i) + sqrt(-1)*Z(:,i+1) and */
/* > Z(:,i) - sqrt(-1)*Z(:,i+1), respectively. */
/* > If JOBZ == 'F', then the above descriptions hold for */
/* > the columns of Z*V, where the columns of V are the */
/* > eigenvectors of the K-by-K Rayleigh quotient, and Z is */
/* > orthonormal. The columns of V are similarly structured: */
/* > If IMEIG(i) == 0 then Z*V(:,i) is an eigenvector, and if */
/* > IMEIG(i) > 0 then Z*V(:,i)+sqrt(-1)*Z*V(:,i+1) and */
/* > Z*V(:,i)-sqrt(-1)*Z*V(:,i+1) */
/* > are the eigenvectors of LAMBDA(i), LAMBDA(i+1). */
/* > See the descriptions of REIG, IMEIG, X and V. */
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
/* > RES (output) REAL (N-1)-by-1 array */
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
/* > See the description of Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] B */
/* > \verbatim */
/* > B (output) REAL MIN(M,N)-by-(N-1) array. */
/* > IF JOBF =='R', B(1:N,1:K) contains A*U(:,1:K), and can */
/* > be used for computing the refined vectors;
see further */
/* > details in the provided references. */
/* > If JOBF == 'E', B(1:N,1;
K) contains */
/* > A*U(:,1:K)*W(1:K,1:K), which are the vectors from the */
/* > Exact DMD, up to scaling by the inverse eigenvalues. */
/* > In both cases, the content of B can be lifted to the */
/* > original dimension of the input data by pre-multiplying */
/* > with the Q factor from the initial QR factorization. */
/* > Here A denotes a compression of the underlying operator. */
/* > See the descriptions of F and X. */
/* > If JOBF =='N', then B is not referenced. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDB */
/* > \verbatim */
/* > LDB (input) INTEGER, LDB >= MIN(M,N) */
/* > The leading dimension of the array B. */
/* > \endverbatim */
/* ..... */
/* > \param[out] V */
/* > \verbatim */
/* > V (workspace/output) REAL (N-1)-by-(N-1) array */
/* > On exit, V(1:K,1:K) contains the K eigenvectors of */
/* > the Rayleigh quotient. The eigenvectors of a complex */
/* > conjugate pair of eigenvalues are returned in real form */
/* > as explained in the description of Z. The Ritz vectors */
/* > (returned in Z) are the product of X and V;
see */
/* > the descriptions of X and Z. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDV */
/* > \verbatim */
/* > LDV (input) INTEGER, LDV >= N-1 */
/* > The leading dimension of the array V. */
/* > \endverbatim */
/* ..... */
/* > \param[out] S */
/* > \verbatim */
/* > S (output) REAL (N-1)-by-(N-1) array */
/* > The array S(1:K,1:K) is used for the matrix Rayleigh */
/* > quotient. This content is overwritten during */
/* > the eigenvalue decomposition by SGEEV. */
/* > See the description of K. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDS */
/* > \verbatim */
/* > LDS (input) INTEGER, LDS >= N-1 */
/* > The leading dimension of the array S. */
/* > \endverbatim */
/* ..... */
/* > \param[out] WORK */
/* > \verbatim */
/* > WORK (workspace/output) REAL LWORK-by-1 array */
/* > On exit, */
/* > WORK(1:MIN(M,N)) contains the scalar factors of the */
/* > elementary reflectors as returned by SGEQRF of the */
/* > M-by-N input matrix F. */
/* > WORK(MIN(M,N)+1:MIN(M,N)+N-1) contains the singular values of */
/* > the input submatrix F(1:M,1:N-1). */
/* > If the call to SGEDMDQ is only workspace query, then */
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
/* > Let MLWQR = N (minimal workspace for SGEQRF[M,N]) */
/* > MLWDMD = minimal workspace for SGEDMD (see the */
/* > description of LWORK in SGEDMD) for */
/* > snapshots of dimensions MIN(M,N)-by-(N-1) */
/* > MLWMQR = N (minimal workspace for */
/* > SORMQR['L','N',M,N,N]) */
/* > MLWGQR = N (minimal workspace for SORGQR[M,N,N]) */
/* > Then */
/* > LWORK = MAX(N+MLWQR, N+MLWDMD) */
/* > is updated as follows: */
/* > if JOBZ == 'V' or JOBZ == 'F' THEN */
/* > LWORK = MAX( LWORK,MIN(M,N)+N-1 +MLWMQR ) */
/* > if JOBQ == 'Q' THEN */
/* > LWORK = MAX( LWORK,MIN(M,N)+N-1+MLWGQR) */
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
/* > Let M1=MIN(M,N), N1=N-1. Then */
/* > If WHTSVD == 2, then LIWORK >= MAX(1,8*MIN(M1,N1)) */
/* > If WHTSVD == 3, then LIWORK >= MAX(1,M1+N1-1) */
/* > If WHTSVD == 4, then LIWORK >= MAX(3,M1+3*N1) */
/* > If on entry LIWORK = -1, then a worskpace query is */
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
void sgedmdq_(char *jobs, char *jobz, char *jobr, char *jobq, char *jobt, char *jobf,
              integer *whtsvd, integer *m, integer *n, real *f, integer *ldf, real *x, integer *ldx,
              real *y, integer *ldy, integer *nrnk, real *tol, integer *k, real *reig, real *imeig,
              real *z__, integer *ldz, real *res, real *b, integer *ldb, real *v, integer *ldv,
              real *s, integer *lds, real *work, integer *lwork, integer *iwork, integer *liwork,
              integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgedmdq inputs: jobs %\c ,jobz %\c ,jobr %\c ,jobq %\c ,jobt %\c ,jobf %\c "
                      ",whtsvd %" FLA_IS ",m %" FLA_IS ",n %" FLA_IS ",ldf %" FLA_IS ",ldx %" FLA_IS
                      ",ldy %" FLA_IS ",nrnk %" FLA_IS ",ldz %" FLA_IS ",ldb %" FLA_IS
                      ",ldv %" FLA_IS ",lds %" FLA_IS ",lwork %" FLA_IS ",liwork %" FLA_IS "",
                      *jobs, *jobz, *jobr, *jobq, *jobt, *jobf, *whtsvd, *m, *n, *ldf, *ldx, *ldy,
                      *nrnk, *ldz, *ldb, *ldv, *lds, *lwork, *liwork);
    /* System generated locals */
    integer f_dim1, f_offset, x_dim1, x_offset, y_dim1, y_offset, z_dim1, z_offset, b_dim1,
        b_offset, v_dim1, v_offset, s_dim1, s_offset, i__1, i__2;
    /* Local variables */
    integer info1;
    extern logical lsame_(char *, char *, integer, integer);
    char jobvl[1];
    integer minmn;
    logical wantq;
    integer mlwqr, olwqr;
    logical wntex;
    extern /* Subroutine */
        void
        sgedmd_(char *, char *, char *, char *, integer *, integer *, integer *, real *, integer *,
                real *, integer *, integer *, real *, integer *, real *, real *, real *, integer *,
                real *, real *, integer *, real *, integer *, real *, integer *, real *, integer *,
                integer *, integer *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    integer mlwdmd, olwdmd;
    extern /* Subroutine */
        void
        sgeqrf_(integer *, integer *, real *, integer *, real *, real *, integer *, integer *);
    logical sccolx, sccoly;
    extern /* Subroutine */
        void
        slacpy_(char *, integer *, integer *, real *, integer *, real *, integer *),
        slaset_(char *, integer *, integer *, real *, real *, real *, integer *);
    integer iminwr;
    logical wntvec, wntvcf;
    integer mlwgqr;
    logical wntref;
    integer mlwork, olwgqr, olwork;
    real rdummy[2];
    integer mlwmqr, olwmqr;
    logical lquery, wntres, wnttrf, wntvcq;
    extern /* Subroutine */
        void
        sorgqr_(integer *, integer *, integer *, real *, integer *, real *, real *, integer *,
                integer *),
        sormqr_(char *, char *, integer *, integer *, integer *, real *, integer *, real *, real *,
                integer *, real *, integer *, integer *);
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 04:20:24 10/24/24 */
    /* ...Switches: */
    /* -- LAPACK driver routine -- */
    /* -- LAPACK is a software package provided by University of -- */
    /* -- Tennessee, University of California Berkeley, University of -- */
    /* -- Colorado Denver and NAG Ltd.. -- */
    /* ..... */
    /* Scalar arguments */
    /* ~~~~~~~~~~~~~~~~ */
    /* Array arguments */
    /* ~~~~~~~~~~~~~~~~ */
    /* Parameters */
    /* ~~~~~~~~~~ */
    /* Local scalars */
    /* ~~~~~~~~~~~~~ */
    /* Local array */
    /* ~~~~~~~~~~~ */
    /* External functions (BLAS and LAPACK) */
    /* ~~~~~~~~~~~~~~~~~ */
    /* .......................................................... */
    /* Test the input arguments */
    /* Parameter adjustments */
    f_dim1 = *ldf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
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
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
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
    wntvcf = lsame_(jobz, "F", 1, 1);
    wntvcq = lsame_(jobz, "Q", 1, 1);
    wntref = lsame_(jobf, "R", 1, 1);
    wntex = lsame_(jobf, "E", 1, 1);
    wantq = lsame_(jobq, "Q", 1, 1);
    wnttrf = lsame_(jobt, "R", 1, 1);
    minmn = fla_min(*m, *n);
    *info = 0;
    lquery = *lwork == -1 || *liwork == -1;
    if(!(sccolx || sccoly || lsame_(jobs, "N", 1, 1)))
    {
        *info = -1;
    }
    else if(!(wntvec || wntvcf || wntvcq || lsame_(jobz, "N", 1, 1)))
    {
        *info = -2;
    }
    else if(!(wntres || lsame_(jobr, "N", 1, 1)) || wntres && lsame_(jobz, "N", 1, 1))
    {
        *info = -3;
    }
    else if(!(wantq || lsame_(jobq, "N", 1, 1)))
    {
        *info = -4;
    }
    else if(!(wnttrf || lsame_(jobt, "N", 1, 1)))
    {
        *info = -5;
    }
    else if(!(wntref || wntex || lsame_(jobf, "N", 1, 1)))
    {
        *info = -6;
    }
    else if(!(*whtsvd == 1 || *whtsvd == 2 || *whtsvd == 3 || *whtsvd == 4))
    {
        *info = -7;
    }
    else if(*m < 0)
    {
        *info = -8;
    }
    else if(*n < 0 || *n > *m + 1)
    {
        *info = -9;
    }
    else if(*ldf < *m)
    {
        *info = -11;
    }
    else if(*ldx < minmn)
    {
        *info = -13;
    }
    else if(*ldy < minmn)
    {
        *info = -15;
    }
    else if(!(*nrnk == -2 || *nrnk == -1 || *nrnk >= 1 && *nrnk <= *n))
    {
        *info = -16;
    }
    else if(*tol < 0.f || *tol >= 1.f)
    {
        *info = -17;
    }
    else if(*ldz < *m)
    {
        *info = -22;
    }
    else if((wntref || wntex) && *ldb < minmn)
    {
        *info = -25;
    }
    else if(*ldv < *n - 1)
    {
        *info = -27;
    }
    else if(*lds < *n - 1)
    {
        *info = -29;
    }
    if(wntvec || wntvcf)
    {
        *(unsigned char *)jobvl = 'V';
    }
    else
    {
        *(unsigned char *)jobvl = 'N';
    }
    if(*info == 0)
    {
        /* Compute the minimal and the optimal workspace */
        /* requirements. Simulate running the code and */
        /* determine minimal and optimal sizes of the */
        /* workspace at any moment of the run. */
        if(*n == 0 || *n == 1)
        {
            /* All output except K is void. INFO=1 signals */
            /* the void input. In case of a workspace query, */
            /* the minimal workspace lengths are returned. */
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
        mlwqr = fla_max(1, *n);
        /* Minimal workspace length for SGEQRF. */
        mlwork = fla_min(*m, *n) + mlwqr;
        if(lquery)
        {
            sgeqrf_(m, n, &f[f_offset], ldf, &work[1], rdummy, &c_n1, &info1);
            olwqr = (integer)rdummy[0];
            olwork = fla_min(*m, *n) + olwqr;
        }
        i__1 = *n - 1;
        sgedmd_(jobs, jobvl, jobr, jobf, whtsvd, &minmn, &i__1, &x[x_offset], ldx, &y[y_offset],
                ldy, nrnk, tol, k, &reig[1], &imeig[1], &z__[z_offset], ldz, &res[1], &b[b_offset],
                ldb, &v[v_offset], ldv, &s[s_offset], lds, &work[1], &c_n1, &iwork[1], liwork,
                &info1);
        mlwdmd = (integer)work[1];
        /* Computing MAX */
        i__1 = mlwork;
        i__2 = minmn + mlwdmd; // , expr subst
        mlwork = fla_max(i__1, i__2);
        iminwr = iwork[1];
        if(lquery)
        {
            olwdmd = (integer)work[2];
            /* Computing MAX */
            i__1 = olwork;
            i__2 = minmn + olwdmd; // , expr subst
            olwork = fla_max(i__1, i__2);
        }
        if(wntvec || wntvcf)
        {
            mlwmqr = fla_max(1, *n);
            /* Computing MAX */
            i__1 = mlwork;
            i__2 = minmn + *n - 1 + mlwmqr; // , expr subst
            mlwork = fla_max(i__1, i__2);
            if(lquery)
            {
                sormqr_("L", "N", m, n, &minmn, &f[f_offset], ldf, &work[1], &z__[z_offset], ldz,
                        &work[1], &c_n1, &info1);
                olwmqr = (integer)work[1];
                /* Computing MAX */
                i__1 = olwork;
                i__2 = minmn + *n - 1 + olwmqr; // , expr subst
                olwork = fla_max(i__1, i__2);
            }
        }
        if(wantq)
        {
            mlwgqr = *n;
            /* Computing MAX */
            i__1 = mlwork;
            i__2 = minmn + *n - 1 + mlwgqr; // , expr subst
            mlwork = fla_max(i__1, i__2);
            if(lquery)
            {
                sorgqr_(m, &minmn, &minmn, &f[f_offset], ldf, &work[1], &work[1], &c_n1, &info1);
                olwgqr = (integer)work[1];
                /* Computing MAX */
                i__1 = olwork;
                i__2 = minmn + *n - 1 + olwgqr; // , expr subst
                olwork = fla_max(i__1, i__2);
            }
        }
        iminwr = fla_max(1, iminwr);
        mlwork = fla_max(2, mlwork);
        if(*lwork < mlwork && !lquery)
        {
            *info = -31;
        }
        if(*liwork < iminwr && !lquery)
        {
            *info = -33;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("SGEDMDQ", &i__1, (ftnlen)7);
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
    /* ..... */
    /* Initial QR factorization that is used to represent the */
    /* snapshots as elements of lower dimensional subspace. */
    /* For large scale computation with M >>N , at this place */
    /* one can use an out of core QRF. */
    i__1 = *lwork - minmn;
    sgeqrf_(m, n, &f[f_offset], ldf, &work[1], &work[minmn + 1], &i__1, &info1);
    /* Define X and Y as the snapshots representations in the */
    /* orthogonal basis computed in the QR factorization. */
    /* X corresponds to the leading N-1 and Y to the trailing */
    /* N-1 snapshots. */
    i__1 = *n - 1;
    slaset_("L", &minmn, &i__1, &c_b30, &c_b30, &x[x_offset], ldx);
    i__1 = *n - 1;
    slacpy_("U", &minmn, &i__1, &f[f_offset], ldf, &x[x_offset], ldx);
    i__1 = *n - 1;
    slacpy_("A", &minmn, &i__1, &f[(f_dim1 << 1) + 1], ldf, &y[y_offset], ldy);
    if(*m >= 3)
    {
        i__1 = minmn - 2;
        i__2 = *n - 2;
        slaset_("L", &i__1, &i__2, &c_b30, &c_b30, &y[y_dim1 + 3], ldy);
    }
    /* Compute the DMD of the projected snapshot pairs (X,Y) */
    i__1 = *n - 1;
    i__2 = *lwork - minmn;
    sgedmd_(jobs, jobvl, jobr, jobf, whtsvd, &minmn, &i__1, &x[x_offset], ldx, &y[y_offset], ldy,
            nrnk, tol, k, &reig[1], &imeig[1], &z__[z_offset], ldz, &res[1], &b[b_offset], ldb,
            &v[v_offset], ldv, &s[s_offset], lds, &work[minmn + 1], &i__2, &iwork[1], liwork,
            &info1);
    if(info1 == 2 || info1 == 3)
    {
        /* Return with error code. */
        *info = info1;
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else
    {
        *info = info1;
    }
    /* The Ritz vectors (Koopman modes) can be explicitly */
    /* formed or returned in factored form. */
    if(wntvec)
    {
        /* Compute the eigenvectors explicitly. */
        if(*m > minmn)
        {
            i__1 = *m - minmn;
            slaset_("A", &i__1, k, &c_b30, &c_b30, &z__[minmn + 1 + z_dim1], ldz);
        }
        i__1 = *lwork - (minmn + *n - 1);
        sormqr_("L", "N", m, k, &minmn, &f[f_offset], ldf, &work[1], &z__[z_offset], ldz,
                &work[minmn + *n], &i__1, &info1);
    }
    else if(wntvcf)
    {
        /* Return the Ritz vectors (eigenvectors) in factored */
        /* form Z*V, where Z contains orthonormal matrix (the */
        /* product of Q from the initial QR factorization and */
        /* the SVD/POD_basis returned by SGEDMD in X) and the */
        /* second factor (the eigenvectors of the Rayleigh */
        /* quotient) is in the array V, as returned by SGEDMD. */
        slacpy_("A", n, k, &x[x_offset], ldx, &z__[z_offset], ldz);
        if(*m > *n)
        {
            i__1 = *m - *n;
            slaset_("A", &i__1, k, &c_b30, &c_b30, &z__[*n + 1 + z_dim1], ldz);
        }
        i__1 = *lwork - (minmn + *n - 1);
        sormqr_("L", "N", m, k, &minmn, &f[f_offset], ldf, &work[1], &z__[z_offset], ldz,
                &work[minmn + *n], &i__1, &info1);
    }
    /* Some optional output variables: */
    /* The upper triangular factor in the initial QR */
    /* factorization is optionally returned in the array Y. */
    /* This is useful if this call to SGEDMDQ is to be */
    /* followed by a streaming DMD that is implemented in a */
    /* QR compressed form. */
    if(wnttrf)
    {
        /* Return the upper triangular R in Y */
        slaset_("A", &minmn, n, &c_b30, &c_b30, &y[y_offset], ldy);
        slacpy_("U", &minmn, n, &f[f_offset], ldf, &y[y_offset], ldy);
    }
    if(wantq)
    {
        i__1 = *lwork - (minmn + *n - 1);
        sorgqr_(m, &minmn, &minmn, &f[f_offset], ldf, &work[1], &work[minmn + *n], &i__1, &info1);
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
/* sgedmdq_ */
