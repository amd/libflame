/* ./cgedmd.f -- translated by f2c (version 20190311). You must link the resulting object file with
 libf2c: on Microsoft Windows system, link with libf2c.lib; on Linux or Unix systems, link with
 .../path/to/libf2c.a -lm or, if you install libf2c.a in a standard place, with -lf2c -lm -- in that
 order, at the end of the command line, as in cc *.o -lf2c -lm Source for libf2c is in
 /netlib/f2c/libf2c.zip, e.g., http://www.netlib.org/f2c/libf2c.zip */
#include "FLA_f2c.h" /* Table of constant values */
static complex c_b1 = {1.f, 0.f};
static complex c_b2 = {0.f, 0.f};
static integer c_n1 = -1;
static integer c__1 = 1;
static integer c__0 = 0;
static real c_b55 = 1.f;
static real c_b65 = 0.f;
/* > \brief \b CGEDMD computes the Dynamic Mode Decomposition (DMD) for a pair of data snapshot
 * matrices. */
/* =========== DOCUMENTATION =========== */
/* Definition: */
/* =========== */
/* SUBROUTINE CGEDMD( JOBS, JOBZ, JOBR, JOBF, WHTSVD, & */
/* M, N, X, LDX, Y, LDY, NRNK, TOL, & */
/* K, EIGS, Z, LDZ, RES, B, LDB, & */
/* W, LDW, S, LDS, ZWORK, LZWORK, & */
/* RWORK, LRWORK, IWORK, LIWORK, INFO ) */
/* ..... */
/* USE iso_fortran_env */
/* IMPLICIT NONE */
/* INTEGER, PARAMETER :: WP = real32 */
/* ..... */
/* Scalar arguments */
/* CHARACTER, INTENT(IN) :: JOBS, JOBZ, JOBR, JOBF */
/* INTEGER, INTENT(IN) :: WHTSVD, M, N, LDX, LDY, & */
/* NRNK, LDZ, LDB, LDW, LDS, & */
/* LIWORK, LRWORK, LZWORK */
/* INTEGER, INTENT(OUT) :: K, INFO */
/* REAL, INTENT(IN) :: TOL */
/* Array arguments */
/* COMPLEX, INTENT(INOUT) :: X(LDX,*), Y(LDY,*) */
/* COMPLEX, INTENT(OUT) :: Z(LDZ,*), B(LDB,*), & */
/* W(LDW,*), S(LDS,*) */
/* COMPLEX, INTENT(OUT) :: EIGS(*) */
/* COMPLEX, INTENT(OUT) :: ZWORK(*) */
/* REAL, INTENT(OUT) :: RES(*) */
/* REAL, INTENT(OUT) :: RWORK(*) */
/* INTEGER, INTENT(OUT) :: IWORK(*) */
/* ............................................................ */
/* > \par Purpose: */
/* ============= */
/* > \verbatim */
/* > CGEDMD computes the Dynamic Mode Decomposition (DMD) for */
/* > a pair of data snapshot matrices. For the input matrices */
/* > X and Y such that Y = A*X with an unaccessible matrix */
/* > A, CGEDMD computes a certain number of Ritz pairs of A using */
/* > the standard Rayleigh-Ritz extraction from a subspace of */
/* > range(X) that is determined using the leading left singular */
/* > vectors of X. Optionally, CGEDMD returns the residuals */
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
/* > 1 :: CGESVD (the QR SVD algorithm) */
/* > 2 :: CGESDD (the Divide and Conquer algorithm;
if enough */
/* > workspace available, this is the fastest option) */
/* > 3 :: CGESVDQ (the preconditioned QR SVD ;
this and 4 */
/* > are the most accurate options) */
/* > 4 :: CGEJSV (the preconditioned Jacobi SVD;
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
/* > X (input/output) COMPLEX M-by-N array */
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
/* > Y (input/workspace/output) COMPLEX M-by-N array */
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
/* > \param[out] EIGS */
/* > \verbatim */
/* > EIGS (output) COMPLEX N-by-1 array */
/* > The leading K (K<=N) entries of EIGS contain */
/* > the computed eigenvalues (Ritz values). */
/* > See the descriptions of K, and Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] Z */
/* > \verbatim */
/* > Z (workspace/output) COMPLEX M-by-N array */
/* > If JOBZ =='V' then Z contains the Ritz vectors. Z(:,i) */
/* > is an eigenvector of the i-th Ritz value;
||Z(:,i)||_2=1. */
/* > If JOBZ == 'F', then the Z(:,i)'s are given implicitly as */
/* > the columns of X(:,1:K)*W(1:K,1:K), i.e. X(:,1:K)*W(:,i) */
/* > is an eigenvector corresponding to EIGS(i). The columns */
/* > of W(1:k,1:K) are the computed eigenvectors of the */
/* > K-by-K Rayleigh quotient. */
/* > See the descriptions of EIGS, X and W. */
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
/* > Ritz pairs, */
/* > RES(i) = || A * Z(:,i) - EIGS(i)*Z(:,i))||_2. */
/* > See the description of EIGS and Z. */
/* > \endverbatim */
/* ..... */
/* > \param[out] B */
/* > \verbatim */
/* > B (output) COMPLEX M-by-N array. */
/* > IF JOBF =='R', B(1:M,1:K) contains A*U(:,1:K), and can */
/* > be used for computing the refined vectors;
see further */
/* > details in the provided references. */
/* > If JOBF == 'E', B(1:M,1:K) contains */
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
/* > W (workspace/output) COMPLEX N-by-N array */
/* > On exit, W(1:K,1:K) contains the K computed */
/* > eigenvectors of the matrix Rayleigh quotient. */
/* > The Ritz vectors (returned in Z) are the */
/* > product of X (containing a POD basis for the input */
/* > matrix X) and W. See the descriptions of K, S, X and Z. */
/* > W is also used as a workspace to temporarily store the */
/* > right singular vectors of X. */
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
/* > S (workspace/output) COMPLEX N-by-N array */
/* > The array S(1:K,1:K) is used for the matrix Rayleigh */
/* > quotient. This content is overwritten during */
/* > the eigenvalue decomposition by CGEEV. */
/* > See the description of K. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LDS */
/* > \verbatim */
/* > LDS (input) INTEGER, LDS >= N */
/* > The leading dimension of the array S. */
/* > \endverbatim */
/* ..... */
/* > \param[out] ZWORK */
/* > \verbatim */
/* > ZWORK (workspace/output) COMPLEX LZWORK-by-1 array */
/* > ZWORK is used as complex workspace in the complex SVD, as */
/* > specified by WHTSVD (1,2, 3 or 4) and for CGEEV for computing */
/* > the eigenvalues of a Rayleigh quotient. */
/* > If the call to CGEDMD is only workspace query, then */
/* > ZWORK(1) contains the minimal complex workspace length and */
/* > ZWORK(2) is the optimal complex workspace length. */
/* > Hence, the length of work is at least 2. */
/* > See the description of LZWORK. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LZWORK */
/* > \verbatim */
/* > LZWORK (input) INTEGER */
/* > The minimal length of the workspace vector ZWORK. */
/* > LZWORK is calculated as MAX(LZWORK_SVD, LZWORK_CGEEV), */
/* > where LZWORK_CGEEV = MAX( 1, 2*N ) and the minimal */
/* > LZWORK_SVD is calculated as follows */
/* > If WHTSVD == 1 :: CGESVD :: */
/* > LZWORK_SVD = MAX(1,2*MIN(M,N)+MAX(M,N)) */
/* > If WHTSVD == 2 :: CGESDD :: */
/* > LZWORK_SVD = 2*MIN(M,N)*MIN(M,N)+2*MIN(M,N)+MAX(M,N) */
/* > If WHTSVD == 3 :: CGESVDQ :: */
/* > LZWORK_SVD = obtainable by a query */
/* > If WHTSVD == 4 :: CGEJSV :: */
/* > LZWORK_SVD = obtainable by a query */
/* > If on entry LZWORK = -1, then a workspace query is */
/* > assumed and the procedure only computes the minimal */
/* > and the optimal workspace lengths and returns them in */
/* > LZWORK(1) and LZWORK(2), respectively. */
/* > \endverbatim */
/* ..... */
/* > \param[out] RWORK */
/* > \verbatim */
/* > RWORK (workspace/output) REAL LRWORK-by-1 array */
/* > On exit, RWORK(1:N) contains the singular values of */
/* > X (for JOBS=='N') or column scaled X (JOBS=='S', 'C'). */
/* > If WHTSVD==4, then RWORK(N+1) and RWORK(N+2) contain */
/* > scaling factor RWORK(N+2)/RWORK(N+1) used to scale X */
/* > and Y to avoid overflow in the SVD of X. */
/* > This may be of interest if the scaling option is off */
/* > and as many as possible smallest eigenvalues are */
/* > desired to the highest feasible accuracy. */
/* > If the call to CGEDMD is only workspace query, then */
/* > RWORK(1) contains the minimal workspace length. */
/* > See the description of LRWORK. */
/* > \endverbatim */
/* ..... */
/* > \param[in] LRWORK */
/* > \verbatim */
/* > LRWORK (input) INTEGER */
/* > The minimal length of the workspace vector RWORK. */
/* > LRWORK is calculated as follows: */
/* > LRWORK = MAX(1, N+LRWORK_SVD,N+LRWORK_CGEEV), where */
/* > LRWORK_CGEEV = MAX(1,2*N) and RWORK_SVD is the real workspace */
/* > for the SVD subroutine determined by the input parameter */
/* > WHTSVD. */
/* > If WHTSVD == 1 :: CGESVD :: */
/* > LRWORK_SVD = 5*MIN(M,N) */
/* > If WHTSVD == 2 :: CGESDD :: */
/* > LRWORK_SVD = MAX(5*MIN(M,N)*MIN(M,N)+7*MIN(M,N), */
/* > 2*MAX(M,N)*MIN(M,N)+2*MIN(M,N)*MIN(M,N)+MIN(M,N) ) ) */
/* > If WHTSVD == 3 :: CGESVDQ :: */
/* > LRWORK_SVD = obtainable by a query */
/* > If WHTSVD == 4 :: CGEJSV :: */
/* > LRWORK_SVD = obtainable by a query */
/* > If on entry LRWORK = -1, then a workspace query is */
/* > assumed and the procedure only computes the minimal */
/* > real workspace length and returns it in RWORK(1). */
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
/* > and the optimal workspace lengths for ZWORK, RWORK and */
/* > IWORK. See the descriptions of ZWORK, RWORK and IWORK. */
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
void cgedmd_(char *jobs, char *jobz, char *jobr, char *jobf, integer *whtsvd, integer *m,
             integer *n, complex *x, integer *ldx, complex *y, integer *ldy, integer *nrnk,
             real *tol, integer *k, complex *eigs, complex *z__, integer *ldz, real *res,
             complex *b, integer *ldb, complex *w, integer *ldw, complex *s, integer *lds,
             complex *zwork, integer *lzwork, real *rwork, integer *lrwork, integer *iwork,
             integer *liwork, integer *info)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgedmd inputs: jobs %c ,jobz %c ,jobr %c ,jobf %c ,whtsvd %" FLA_IS
                      ",m %" FLA_IS ",n %" FLA_IS ",ldx %" FLA_IS ",ldy %" FLA_IS ",nrnk %" FLA_IS
                      ",ldz %" FLA_IS ", ldb %" FLA_IS ",ldw %" FLA_IS ",lds %" FLA_IS
                      ", lzwork %" FLA_IS ",lrwork %" FLA_IS ",liwork %" FLA_IS "",
                      *jobs, *jobz, *jobr, *jobf, *whtsvd, *m, *n, *ldx, *ldy, *nrnk, *ldz, *ldb,
                      *ldw, *lds, *lzwork, *lrwork, *liwork);
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, z_dim1, z_offset, b_dim1, b_offset, w_dim1,
        w_offset, s_dim1, s_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1, r__2;
    complex q__1, q__2;
    /* Builtin functions */
    double sqrt(doublereal), c_abs(complex *);
    /* Local variables */
    integer i__, j;
    real ofl, ssum;
    integer info1, info2;
    real xscl1, xscl2, scale;
    extern /* Subroutine */
        void
        cgemm_(char *, char *, integer *, integer *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, complex *, integer *),
        cgeev_(char *, char *, integer *, complex *, integer *, complex *, complex *, integer *,
               complex *, integer *, complex *, integer *, real *, integer *);
    extern logical lsame_(char *, char *, integer, integer);
    logical badxy;
    real small_val;
    char jobzl[1];
    extern /* Subroutine */
        void
        caxpy_(integer *, complex *, complex *, integer *, complex *, integer *);
    logical wntex;
    extern real scnrm2_(integer *, complex *, integer *);
    extern /* Subroutine */
        void
        cgesdd_(char *, integer *, integer *, complex *, integer *, real *, complex *, integer *,
                complex *, integer *, complex *, integer *, real *, integer *, integer *),
        clascl_(char *, integer *, integer *, real *, real *, integer *, integer *, complex *,
                integer *, integer *);
    extern integer icamax_(integer *, complex *, integer *);
    extern real slamch_(char *);
    extern /* Subroutine */
        void
        csscal_(integer *, real *, complex *, integer *),
        cgesvd_(char *, char *, integer *, integer *, complex *, integer *, real *, complex *,
                integer *, complex *, integer *, complex *, integer *, real *, integer *),
        clacpy_(char *, integer *, integer *, complex *, integer *, complex *, integer *),
        xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    char t_or_n__[1];
    extern /* Subroutine */
        void
        cgejsv_(char *, char *, char *, char *, char *, char *, integer *, integer *, complex *,
                integer *, real *, complex *, integer *, complex *, integer *, complex *, integer *,
                real *, integer *, integer *, integer *),
        classq_(integer *, complex *, integer *, real *, real *);
    logical sccolx, sccoly;
    extern logical sisnan_(real *);
    integer lwrsdd, mwrsdd, iminwr;
    logical wntref, wntvec;
    real rootsc;
    integer lwrkev, mlwork, mwrkev, numrnk, olwork, lwrsvd, mwrsvd, mlrwrk;
    logical lquery, wntres;
    char jsvopt[1];
    integer lwrsvj, mwrsvj;
    real rdummy[2];
    integer lwrsvq, mwrsvq;
    extern /* Subroutine */
        void
        cgesvdq_(char *, char *, char *, char *, char *, integer *, integer *, complex *, integer *,
                 real *, complex *, integer *, complex *, integer *, integer *, integer *,
                 integer *, complex *, integer *, real *, integer *, integer *);
    /* ...Translated by Pacific-Sierra Research vf90 Personal 3.4N3 02:36:14 10/25/24 */
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
    --eigs;
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
    --zwork;
    --rwork;
    --iwork;
    /* Function Body */
    wntres = lsame_(jobr, "R", 1, 1);
    sccolx = lsame_(jobs, "S", 1, 1) || lsame_(jobs, "C", 1, 1);
    sccoly = lsame_(jobs, "Y", 1, 1);
    wntvec = lsame_(jobz, "V", 1, 1);
    wntref = lsame_(jobf, "R", 1, 1);
    wntex = lsame_(jobf, "E", 1, 1);
    *info = 0;
    lquery = *lzwork == -1 || *liwork == -1 || *lrwork == -1;
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
        *info = -17;
    }
    else if((wntref || wntex) && *ldb < *m)
    {
        *info = -20;
    }
    else if(*ldw < *n)
    {
        *info = -22;
    }
    else if(*lds < *n)
    {
        *info = -24;
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
                rwork[1] = 1.f;
                zwork[1].r = 2.f;
                zwork[1].i = 0.f; // , expr subst
                zwork[2].r = 2.f;
                zwork[2].i = 0.f; // , expr subst
            }
            else
            {
                *k = 0;
            }
            *info = 1;
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        iminwr = 1;
        mlrwrk = fla_max(1, *n);
        mlwork = 2;
        olwork = 2;
        if(*whtsvd == 1)
        {
            /* The following is specified as the minimal */
            /* length of WORK in the definition of CGESVD: */
            /* MWRSVD = MAX(1,2*MIN(M,N)+MAX(M,N)) */
            /* Computing MAX */
            i__1 = 1;
            i__2 = (fla_min(*m, *n) << 1) + fla_max(*m, *n); // , expr subst
            mwrsvd = fla_max(i__1, i__2);
            mlwork = fla_max(2, mwrsvd);
            /* Computing MAX */
            i__1 = mlrwrk;
            i__2 = *n + fla_min(*m, *n) * 5; // , expr subst
            mlrwrk = fla_max(i__1, i__2);
            if(lquery)
            {
                cgesvd_("O", "S", m, n, &x[x_offset], ldx, &rwork[1], &b[b_offset], ldb,
                        &w[w_offset], ldw, &zwork[1], &c_n1, rdummy, &info1);
                lwrsvd = (integer)zwork[1].r;
                olwork = fla_max(2, lwrsvd);
            }
        }
        else if(*whtsvd == 2)
        {
            /* The following is specified as the minimal */
            /* length of WORK in the definition of CGESDD: */
            /* MWRSDD = 2*fla_min(M,N)*fla_min(M,N)+2*fla_min(M,N)+max(M,N). */
            /* RWORK length: 5*MIN(M,N)*MIN(M,N)+7*MIN(M,N) */
            /* In LAPACK 3.10.1 RWORK is defined differently. */
            /* Below we take max over the two versions. */
            /* IMINWR = 8*MIN(M,N) */
            mwrsdd = (fla_min(*m, *n) << 1) * fla_min(*m, *n) + (fla_min(*m, *n) << 1)
                     + fla_max(*m, *n);
            mlwork = fla_max(2, mwrsdd);
            iminwr = fla_min(*m, *n) << 3;
            /* Computing MAX */
            /* Computing MAX */
            /* Computing MAX */
            i__5 = fla_min(*m, *n) * 5 * fla_min(*m, *n) + fla_min(*m, *n) * 7;
            i__6 = fla_min(*m, *n) * 5 * fla_min(*m, *n) + fla_min(*m, *n) * 5; // , expr subst
            i__3 = fla_max(i__5, i__6);
            i__4 = (fla_max(*m, *n) << 1) * fla_min(*m, *n)
                   + (fla_min(*m, *n) << 1) * fla_min(*m, *n) + fla_min(*m, *n); // , expr subst
            i__1 = mlrwrk;
            i__2 = *n + fla_max(i__3, i__4); // , expr subst
            mlrwrk = fla_max(i__1, i__2);
            if(lquery)
            {
                cgesdd_("O", m, n, &x[x_offset], ldx, &rwork[1], &b[b_offset], ldb, &w[w_offset],
                        ldw, &zwork[1], &c_n1, rdummy, &iwork[1], &info1);
                /* Computing MAX */
                i__1 = mwrsdd;
                i__2 = (integer)zwork[1].r; // , expr subst
                lwrsdd = fla_max(i__1, i__2);
                olwork = fla_max(2, lwrsdd);
            }
        }
        else if(*whtsvd == 3)
        {
            cgesvdq_("H", "P", "N", "R", "R", m, n, &x[x_offset], ldx, &rwork[1], &z__[z_offset],
                     ldz, &w[w_offset], ldw, &numrnk, &iwork[1], &c_n1, &zwork[1], &c_n1, rdummy,
                     &c_n1, &info1);
            iminwr = iwork[1];
            mwrsvq = (integer)zwork[2].r;
            mlwork = fla_max(2, mwrsvq);
            /* Computing MAX */
            i__1 = mlrwrk;
            i__2 = *n + (integer)rdummy[0]; // , expr subst
            mlrwrk = fla_max(i__1, i__2);
            if(lquery)
            {
                lwrsvq = (integer)zwork[1].r;
                olwork = fla_max(2, lwrsvq);
            }
        }
        else if(*whtsvd == 4)
        {
            *(unsigned char *)jsvopt = 'J';
            cgejsv_("F", "U", jsvopt, "N", "N", "P", m, n, &x[x_offset], ldx, &rwork[1],
                    &z__[z_offset], ldz, &w[w_offset], ldw, &zwork[1], &c_n1, rdummy, &c_n1,
                    &iwork[1], &info1);
            iminwr = iwork[1];
            mwrsvj = (integer)zwork[2].r;
            mlwork = fla_max(2, mwrsvj);
            /* Computing MAX */
            /* Computing MAX */
            i__3 = 7;
            i__4 = (integer)rdummy[0]; // , expr subst
            i__1 = mlrwrk;
            i__2 = *n + fla_max(i__3, i__4); // , expr subst
            mlrwrk = fla_max(i__1, i__2);
            if(lquery)
            {
                lwrsvj = (integer)zwork[1].r;
                olwork = fla_max(2, lwrsvj);
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
        /* Workspace calculation to the CGEEV call */
        /* Computing MAX */
        i__1 = 1;
        i__2 = *n << 1; // , expr subst
        mwrkev = fla_max(i__1, i__2);
        mlwork = fla_max(mlwork, mwrkev);
        /* Computing MAX */
        i__1 = mlrwrk;
        i__2 = *n + (*n << 1); // , expr subst
        mlrwrk = fla_max(i__1, i__2);
        if(lquery)
        {
            cgeev_("N", jobzl, n, &s[s_offset], lds, &eigs[1], &w[w_offset], ldw, &w[w_offset], ldw,
                   &zwork[1], &c_n1, &rwork[1], &info1);
            /* LAPACK CALL */
            lwrkev = (integer)zwork[1].r;
            olwork = fla_max(olwork, lwrkev);
            olwork = fla_max(2, olwork);
        }
        if(*liwork < iminwr && !lquery)
        {
            *info = -30;
        }
        if(*lrwork < mlrwrk && !lquery)
        {
            *info = -28;
        }
        if(*lzwork < mlwork && !lquery)
        {
            *info = -26;
        }
    }
    if(*info != 0)
    {
        i__1 = -(*info);
        xerbla_("CGEDMD", &i__1, (ftnlen)6);
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    else if(lquery)
    {
        /* Return minimal and optimal workspace sizes */
        iwork[1] = iminwr;
        rwork[1] = (real)mlrwrk;
        zwork[1].r = (real)mlwork;
        zwork[1].i = 0.f; // , expr subst
        zwork[2].r = (real)olwork;
        zwork[2].i = 0.f; // , expr subst
        AOCL_DTL_TRACE_LOG_EXIT
        return;
    }
    /* ............................................................ */
    ofl = slamch_("O") * slamch_("P");
    small_val = slamch_("S");
    badxy = FALSE_;
    /* <1> Optional scaling of the snapshots (columns of X, Y) */
    /* ========================================================== */
    if(sccolx)
    {
        /* The columns of X will be normalized. */
        /* To prevent overflows, the column norms of X are */
        /* carefully computed using CLASSQ. */
        *k = 0;
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* WORK(i) = SCNRM2( M, X(1,i), 1 ) */
            scale = 0.f;
            classq_(m, &x[i__ * x_dim1 + 1], &c__1, &scale, &ssum);
            if(sisnan_(&scale) || sisnan_(&ssum))
            {
                *k = 0;
                *info = -8;
                i__2 = -(*info);
                xerbla_("CGEDMD", &i__2, (ftnlen)6);
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
                    clascl_("G", &c__0, &c__0, &scale, &r__1, m, &c__1, &x[i__ * x_dim1 + 1], ldx,
                            &info2);
                    rwork[i__] = -scale * (rootsc / (real)(*m));
                }
                else
                {
                    /* X(:,i) will be scaled to unit 2-norm */
                    rwork[i__] = scale * rootsc;
                    clascl_("G", &c__0, &c__0, &rwork[i__], &c_b55, m, &c__1, &x[i__ * x_dim1 + 1],
                            ldx, &info2);
                    /* X(1:M,i) = (ONE/RWORK(i)) * X(1:M,i) ! INTRINSIC */
                    /* LAPACK CALL */
                }
            }
            else
            {
                rwork[i__] = 0.f;
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
            xerbla_("CGEDMD", &i__1, (ftnlen)6);
            AOCL_DTL_TRACE_LOG_EXIT
            return;
        }
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Now, apply the same scaling to the columns of Y. */
            if(rwork[i__] > 0.f)
            {
                /* ! BLAS CALL */
                r__1 = 1.f / rwork[i__];
                csscal_(m, &r__1, &y[i__ * y_dim1 + 1], &c__1);
                /* Y(1:M,i) = (ONE/RWORK(i)) * Y(1:M,i) ! INTRINSIC */
            }
            else if(rwork[i__] < 0.f)
            {
                r__1 = -rwork[i__];
                r__2 = 1.f / (real)(*m);
                clascl_("G", &c__0, &c__0, &r__1, &r__2, m, &c__1, &y[i__ * y_dim1 + 1], ldy,
                        &info2);
                /* LAPACK CALL */
            }
            else if(c_abs(&y[icamax_(m, &y[i__ * y_dim1 + 1], &c__1) + i__ * y_dim1]) != 0.f)
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
                    csscal_(m, &c_b65, &y[i__ * y_dim1 + 1], &c__1);
                }
            }
        }
    }
    if(sccoly)
    {
        /* The columns of Y will be normalized. */
        /* To prevent overflows, the column norms of Y are */
        /* carefully computed using CLASSQ. */
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* RWORK(i) = SCNRM2( M, Y(1,i), 1 ) */
            scale = 0.f;
            classq_(m, &y[i__ * y_dim1 + 1], &c__1, &scale, &ssum);
            if(sisnan_(&scale) || sisnan_(&ssum))
            {
                *k = 0;
                *info = -10;
                i__2 = -(*info);
                xerbla_("CGEDMD", &i__2, (ftnlen)6);
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
                    /* overflow as RWORK(i) = - SCALE * (ROOTSC/M), */
                    /* the minus sign indicating the 1/M factor. */
                    /* Scaling is performed without overflow, and */
                    /* underflow may occur in the smallest entries */
                    /* of Y(:,i). The relative backward and forward */
                    /* errors are small in the ell_2 norm. */
                    r__1 = 1.f / rootsc;
                    clascl_("G", &c__0, &c__0, &scale, &r__1, m, &c__1, &y[i__ * y_dim1 + 1], ldy,
                            &info2);
                    rwork[i__] = -scale * (rootsc / (real)(*m));
                }
                else
                {
                    /* Y(:,i) will be scaled to unit 2-norm */
                    rwork[i__] = scale * rootsc;
                    clascl_("G", &c__0, &c__0, &rwork[i__], &c_b55, m, &c__1, &y[i__ * y_dim1 + 1],
                            ldy, &info2);
                    /* Y(1:M,i) = (ONE/RWORK(i)) * Y(1:M,i) ! INTRINSIC */
                    /* LAPACK CALL */
                }
            }
            else
            {
                rwork[i__] = 0.f;
            }
        }
        i__1 = *n;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* Now, apply the same scaling to the columns of X. */
            if(rwork[i__] > 0.f)
            {
                /* ! BLAS CALL */
                r__1 = 1.f / rwork[i__];
                csscal_(m, &r__1, &x[i__ * x_dim1 + 1], &c__1);
                /* X(1:M,i) = (ONE/RWORK(i)) * X(1:M,i) ! INTRINSIC */
            }
            else if(rwork[i__] < 0.f)
            {
                r__1 = -rwork[i__];
                r__2 = 1.f / (real)(*m);
                clascl_("G", &c__0, &c__0, &r__1, &r__2, m, &c__1, &x[i__ * x_dim1 + 1], ldx,
                        &info2);
                /* LAPACK CALL */
            }
            else if(c_abs(&x[icamax_(m, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1]) != 0.f)
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
        cgesvd_("O", "S", m, n, &x[x_offset], ldx, &rwork[1], &b[b_offset], ldb, &w[w_offset], ldw,
                &zwork[1], lzwork, &rwork[*n + 1], &info1);
        /* LAPACK CALL */
        *(unsigned char *)t_or_n__ = 'C';
    }
    else if(*whtsvd == 2)
    {
        cgesdd_("O", m, n, &x[x_offset], ldx, &rwork[1], &b[b_offset], ldb, &w[w_offset], ldw,
                &zwork[1], lzwork, &rwork[*n + 1], &iwork[1], &info1);
        /* LAPACK CALL */
        *(unsigned char *)t_or_n__ = 'C';
    }
    else if(*whtsvd == 3)
    {
        i__1 = *lrwork - *n;
        cgesvdq_("H", "P", "N", "R", "R", m, n, &x[x_offset], ldx, &rwork[1], &z__[z_offset], ldz,
                 &w[w_offset], ldw, &numrnk, &iwork[1], liwork, &zwork[1], lzwork, &rwork[*n + 1],
                 &i__1, &info1);
        /* ! LAPACK CALL */
        /* LAPACK CALL */
        clacpy_("A", m, &numrnk, &z__[z_offset], ldz, &x[x_offset], ldx);
        *(unsigned char *)t_or_n__ = 'C';
    }
    else if(*whtsvd == 4)
    {
        i__1 = *lrwork - *n;
        cgejsv_("F", "U", jsvopt, "N", "N", "P", m, n, &x[x_offset], ldx, &rwork[1], &z__[z_offset],
                ldz, &w[w_offset], ldw, &zwork[1], lzwork, &rwork[*n + 1], &i__1, &iwork[1],
                &info1);
        /* LAPACK CALL */
        clacpy_("A", m, n, &z__[z_offset], ldz, &x[x_offset], ldx);
        /* LAPACK CALL */
        *(unsigned char *)t_or_n__ = 'N';
        xscl1 = rwork[*n + 1];
        xscl2 = rwork[*n + 2];
        /* This is an exceptional situation. If the */
        /* data matrices are not scaled and the */
        /* largest singular value of X overflows. */
        /* In that case CGEJSV can return the SVD */
        /* in scaled form. The scaling factor can be used */
        /* to rescale the data (X and Y). */
        if(xscl1 != xscl2)
        {
            clascl_("G", &c__0, &c__0, &xscl1, &xscl2, m, n, &y[y_offset], ldy, &info2);
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
    if(rwork[1] == 0.f)
    {
        /* The largest computed singular value of (scaled) */
        /* X is zero. Return error code -8 */
        /* (the 8th input variable had an illegal value). */
        *k = 0;
        *info = -8;
        i__1 = -(*info);
        xerbla_("CGEDMD", &i__1, (ftnlen)6);
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
            if(rwork[i__] <= rwork[1] * *tol || rwork[i__] <= small_val)
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
            if(rwork[i__ + 1] <= rwork[i__] * *tol || rwork[i__] <= small_val)
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
            if(rwork[i__] <= small_val)
            {
                goto L4;
            }
            ++(*k);
        }
    L4:;
    }
    /* Now, U = X(1:M,1:K) is the SVD/POD basis for the */
    /* snapshot data in the input matrix X. */
    /* <4> Compute the Rayleigh quotient S = U^H * A * U. */
    /* Depending on the requested outputs, the computation */
    /* is organized to compute additional auxiliary */
    /* matrices (for the residuals and refinements). */
    /* In all formulas below, we need V_k*Sigma_k^(-1) */
    /* where either V_k is in W(1:N,1:K), or V_k^H is in */
    /* W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)). */
    if(lsame_(t_or_n__, "N", 1, 1))
    {
        i__1 = *k;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            /* ! BLAS CALL */
            r__1 = 1.f / rwork[i__];
            csscal_(n, &r__1, &w[i__ * w_dim1 + 1], &c__1);
            /* W(1:N,i) = (ONE/RWORK(i)) * W(1:N,i) ! INTRINSIC */
        }
    }
    else
    {
        /* This non-unit stride access is due to the fact */
        /* that CGESVD, CGESVDQ and CGESDD return the */
        /* adjoint matrix of the right singular vectors. */
        /* DO i = 1, K */
        /* CALL DSCAL( N, ONE/RWORK(i), W(i,1), LDW ) ! BLAS CALL */
        /* ! W(i,1:N) = (ONE/RWORK(i)) * W(i,1:N) ! INTRINSIC */
        /* END DO */
        i__1 = *k;
        for(i__ = 1; i__ <= i__1; ++i__)
        {
            rwork[*n + i__] = 1.f / rwork[i__];
        }
        i__1 = *n;
        for(j = 1; j <= i__1; ++j)
        {
            i__2 = *k;
            for(i__ = 1; i__ <= i__2; ++i__)
            {
                i__3 = i__ + j * w_dim1;
                i__4 = *n + i__;
                q__2.r = rwork[i__4];
                q__2.i = 0.f; // , expr subst
                i__5 = i__ + j * w_dim1;
                q__1.r = q__2.r * w[i__5].r - q__2.i * w[i__5].i;
                q__1.i = q__2.r * w[i__5].i + q__2.i * w[i__5].r; // , expr subst
                w[i__3].r = q__1.r;
                w[i__3].i = q__1.i; // , expr subst
            }
        }
    }
    if(wntref)
    {
        /* Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K))) */
        /* for computing the refined Ritz vectors */
        /* (optionally, outside CGEDMD). */
        cgemm_("N", t_or_n__, m, k, n, &c_b1, &y[y_offset], ldy, &w[w_offset], ldw, &c_b2,
               &z__[z_offset], ldz);
        /* Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(W(1:K,1:N))) ! INTRI */
        /* Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K)) ! INTRI */
        /* At this point Z contains */
        /* A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and */
        /* this is needed for computing the residuals. */
        /* This matrix is returned in the array B and */
        /* it can be used to compute refined Ritz vectors. */
        /* BLAS CALL */
        clacpy_("A", m, k, &z__[z_offset], ldz, &b[b_offset], ldb);
        /* B(1:M,1:K) = Z(1:M,1:K) ! INTRINSIC */
        /* BLAS CALL */
        cgemm_("C", "N", k, k, m, &c_b1, &x[x_offset], ldx, &z__[z_offset], ldz, &c_b2,
               &s[s_offset], lds);
        /* S(1:K,1:K) = MATMUL(TANSPOSE(X(1:M,1:K)),Z(1:M,1:K)) ! INTRI */
        /* At this point S = U^H * A * U is the Rayleigh quotient. */
        /* BLAS CALL */
    }
    else
    {
        /* A * U(:,1:K) is not explicitly needed and the */
        /* computation is organized differently. The Rayleigh */
        /* quotient is computed more efficiently. */
        cgemm_("C", "N", k, n, m, &c_b1, &x[x_offset], ldx, &y[y_offset], ldy, &c_b2,
               &z__[z_offset], ldz);
        /* Z(1:K,1:N) = MATMUL( TRANSPOSE(X(1:M,1:K)), Y(1:M,1:N) ) ! IN */
        /* BLAS CALL */
        cgemm_("N", t_or_n__, k, k, n, &c_b1, &z__[z_offset], ldz, &w[w_offset], ldw, &c_b2,
               &s[s_offset], lds);
        /* S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(W(1:K,1:N))) ! INTRIN */
        /* S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K))) ! INTRIN */
        /* At this point S = U^H * A * U is the Rayleigh quotient. */
        /* If the residuals are requested, save scaled V_k into Z. */
        /* Recall that V_k or V_k^H is stored in W. */
        /* BLAS CALL */
        if(wntres || wntex)
        {
            if(lsame_(t_or_n__, "N", 1, 1))
            {
                clacpy_("A", n, k, &w[w_offset], ldw, &z__[z_offset], ldz);
            }
            else
            {
                clacpy_("A", k, n, &w[w_offset], ldw, &z__[z_offset], ldz);
            }
        }
    }
    /* <5> Compute the Ritz values and (if requested) the */
    /* right eigenvectors of the Rayleigh quotient. */
    cgeev_("N", jobzl, k, &s[s_offset], lds, &eigs[1], &w[w_offset], ldw, &w[w_offset], ldw,
           &zwork[1], lzwork, &rwork[*n + 1], &info1);
    /* W(1:K,1:K) contains the eigenvectors of the Rayleigh */
    /* quotient. See the description of Z. */
    /* Also, see the description of CGEEV. */
    /* LAPACK CALL */
    if(info1 > 0)
    {
        /* CGEEV failed to compute the eigenvalues and */
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
                cgemm_("N", "N", m, k, k, &c_b1, &z__[z_offset], ldz, &w[w_offset], ldw, &c_b2,
                       &y[y_offset], ldy);
                /* Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K) ! INTRINSIC */
                /* This frees Z;
                Y contains A * U(:,1:K) * W. */
                /* BLAS CALL */
            }
            else
            {
                /* Compute S = V_k * Sigma_k^(-1) * W, where */
                /* V_k * Sigma_k^(-1) (or its adjoint) is stored in Z */
                cgemm_(t_or_n__, "N", n, k, k, &c_b1, &z__[z_offset], ldz, &w[w_offset], ldw, &c_b2,
                       &s[s_offset], lds);
                /* Then, compute Z = Y * S = */
                /* = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) = */
                /* = A * U(:,1:K) * W(1:K,1:K) */
                cgemm_("N", "N", m, k, n, &c_b1, &y[y_offset], ldy, &s[s_offset], lds, &c_b2,
                       &z__[z_offset], ldz);
                /* Save a copy of Z into Y and free Z for holding */
                /* the Ritz vectors. */
                clacpy_("A", m, k, &z__[z_offset], ldz, &y[y_offset], ldy);
                if(wntex)
                {
                    clacpy_("A", m, k, &z__[z_offset], ldz, &b[b_offset], ldb);
                }
            }
        }
        else if(wntex)
        {
            /* Compute S = V_k * Sigma_k^(-1) * W, where */
            /* V_k * Sigma_k^(-1) is stored in Z */
            cgemm_(t_or_n__, "N", n, k, k, &c_b1, &z__[z_offset], ldz, &w[w_offset], ldw, &c_b2,
                   &s[s_offset], lds);
            /* Then, compute Z = Y * S = */
            /* = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) = */
            /* = A * U(:,1:K) * W(1:K,1:K) */
            cgemm_("N", "N", m, k, n, &c_b1, &y[y_offset], ldy, &s[s_offset], lds, &c_b2,
                   &b[b_offset], ldb);
            /* The above call replaces the following two calls */
            /* that were used in the developing-testing phase. */
            /* CALL CGEMM( 'N', 'N', M, K, N, ZONE, Y, LDY, S, & */
            /* LDS, ZZERO, Z, LDZ) */
            /* Save a copy of Z into Y and free Z for holding */
            /* the Ritz vectors. */
            /* CALL CLACPY( 'A', M, K, Z, LDZ, B, LDB ) */
        }
        /* Compute the Ritz vectors */
        if(wntvec)
        {
            cgemm_("N", "N", m, k, k, &c_b1, &x[x_offset], ldx, &w[w_offset], ldw, &c_b2,
                   &z__[z_offset], ldz);
        }
        /* Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K)) ! INTRIN */
        /* BLAS CALL */
        if(wntres)
        {
            i__1 = *k;
            for(i__ = 1; i__ <= i__1; ++i__)
            {
                /* ! BLAS CALL */
                i__2 = i__;
                q__1.r = -eigs[i__2].r;
                q__1.i = -eigs[i__2].i; // , expr subst
                caxpy_(m, &q__1, &z__[i__ * z_dim1 + 1], &c__1, &y[i__ * y_dim1 + 1], &c__1);
                /* Y(1:M,i) = Y(1:M,i) - EIGS(i) * Z(1:M,i) ! */
                res[i__] = scnrm2_(m, &y[i__ * y_dim1 + 1], &c__1);
                /* BLAS CALL */
            }
        }
    }
    if(*whtsvd == 4)
    {
        rwork[*n + 1] = xscl1;
        rwork[*n + 2] = xscl2;
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
/* cgedmd_ */
