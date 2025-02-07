/******************************************************************************
 * Copyright (C) 2021-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface_lin.hh
 *  libflame_interface.hh defines all the Linear solve routines for Libflame CPP templated public
 *  interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_LIN_HH
#define LIBFLAME_INTERFACE_LIN_HH

#include "libflame.hh"

namespace libflame
{
     /** @defgroup LinearSolve Linear solve, AX = B
     * @ingroup LAPACK
     * @{
     */

    /** @defgroup Cholesky Cholesky
     * @ingroup LinearSolve
     * @{
     */

    /** @defgroup potrf potrf
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief Cholesky factorization of a real symmetric positive definite matrix a
        *
        * @details
        * \b Purpose:
        * \verbatim
            Cholesky factorization of a real symmetric positive definite matrix a.
            The factorization has the form
                A = U**T * U,  if uplo = 'U', or
                A = L * L**T,  if uplo = 'L',
            where U is an upper triangular matrix and L is lower triangular.

                This is the block version of the algorithm, calling Level 3 BLAS.
            \endverbatim

            * @param[in] uplo
                      uplo is char*. \n
                      uplo specifies output format \n
                      = 'U': Output is upper triangular factorization of A \n
                      = 'L': Output is lower triangular factorization of A \n
            * @param[in] n
                      n is int*. \n
                      The order of the matrix a. n >= 0 \n
            * @param[in,out] a
                      a is REAL/DOUBLE PRECISION/COMPLEX/COMPLEX*16 array, dimension (lda,n). \n
                      On entry, the symmetric matrix a.  If uplo = 'U', the leading
                      n-by-n upper triangular part of a contains the upper
                      triangular part of the matrix a, and the strictly lower
                      triangular part of a is not referenced.  If uplo = 'L', the
                      leading n-by-n lower triangular part of a contains the lower
                      triangular part of the matrix a, and the strictly upper
                      triangular part of a is not referenced. \n
                      On exit, if info = 0, the factor U or L from the Cholesky
                      factorization A = U**T*U or A = L*L**T. \n
            * @param[in] lda
                      lda is int*. \n
                      The leading dimension of the matrix a, lda >= fla_max(1,n) \n
            * @param[out]	INFO
                      INFO is INTEGER \n
                      = 0:  successful exit \n
                      < 0:  if INFO = -i, the i-th argument had an illegal value \n
                      > 0:  if INFO = i, the leading minor of order i is not
                            positive definite, and the factorization could not be
                            completed. \n

            *     *  */
    template <typename T>
    void potrf(char *uplo, integer *n, T *a, integer *lda, integer *info)
    {
        potrf(uplo, n, a, lda, info);
    }
    /** @}*/ // end of potrf

    /** @defgroup potrf2 potrf2
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief POTRF2 computes the Cholesky factorization of a real symmetric \n
     positive definite matrix A using the recursive algorithm.

 * @details
 * \b Purpose:
    \verbatim
     POTRF2 computes the Cholesky factorization of a real symmetric
     positive definite matrix A using the recursive algorithm.

     The factorization has the form
        A = U**T * U,  if UPLO = 'U', or
        A = L  * L**T,  if UPLO = 'L',
     where U is an upper triangular matrix and L is lower triangular.

     This is the recursive version of the algorithm. It divides
     the matrix into four submatrices:

            [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
        A = [ -----|----- ]  with n1 = n/2
            [  A21 | A22  ]       n2 = n-n1

     The subroutine calls itself to factor A11. Update and scale A21
     or A12, update A22 then call itself to factor A22.
    \endverbatim

  * @param[in] UPLO
           UPLO is CHARACTER*1 \n
           = 'U':  Upper triangle of A is stored; \n
           = 'L':  Lower triangle of A is stored. \n
  * @param[in] N
           N is INTEGER \n
           The order of the matrix A.  N >= 0. \n
  * @param[in,out] A
           A is REAL array, dimension (LDA,N) \n
           On entry, the symmetric matrix A.  If UPLO = 'U', the leading
           N-by-N upper triangular part of A contains the upper
           triangular part of the matrix A, and the strictly lower
           triangular part of A is not referenced.  If UPLO = 'L', the
           leading N-by-N lower triangular part of A contains the lower
           triangular part of the matrix A, and the strictly upper
           triangular part of A is not referenced. \n
  \n
           On exit, if INFO = 0, the factor U or L from the Cholesky
           factorization A = U**T*U or A = L*L**T. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,N). \n
  * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the leading minor of order i is not
                positive definite, and the factorization could not be
                completed. \n

 *  * */
    template <typename T>
    void potrf2(char *uplo, integer *n, T *a, integer *lda, integer *info)
    {
        potrf2(uplo, n, a, lda, info);
    }
    /** @}*/ // end of potrf2

    /** @defgroup potf2 potf2
     * @{
     * @ingroup Cholesky Cholesky
     */

    /*! @brief Cholesky factorization of a real symmetric
    *         positive definite matrix a
    *
    * @details
    * \b Purpose:
    * \verbatim
        Cholesky factorization of a real symmetric positive definite matrix a
        The factorization has the form
            A = U**T * U,  if uplo = 'U', or
            A = L * L**T,  if uplo = 'L',
        where U is an upper triangular matrix and L is lower triangular.

        This is the unblocked version of the algorithm, calling Level 2 BLAS.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              uplo specifies output format \n
              = 'U': Output is upper triangular factorization of A \n
              = 'L': Output is lower triangular factorization of A \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a. n >= 0 \n
    * @param[in,out] a
              a is REAL/DOUBLE PRECISION/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the symmetric matrix a.  If uplo = 'U', the leading
              n-by-n upper triangular part of a contains the upper
              triangular part of the matrix a, and the strictly lower
              triangular part of a is not referenced.  If uplo = 'L', the
              leading n-by-n lower triangular part of a contains the lower
              triangular part of the matrix a, and the strictly upper
              triangular part of a is not referenced. \n
              On exit, if info = 0, the factor U or L from the Cholesky
              factorization A = U**T *U  or A = L*L**T. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,n) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -k, the k-th argument had an illegal value \n
              > 0: if INFO = k, the leading minor of order k is not
                   positive definite, and the factorization could not be
                   completed. \n

    *     *  */
    template <typename T>
    void potf2(char *uplo, integer *n, T *a, integer *lda, integer *info)
    {
        potf2(uplo, n, a, lda, info);
    }
    /** @}*/ // end of potf2

    /** @defgroup pstrf pstrf
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PSTRF computes the Cholesky factorization with complete pivoting of a real symmetric
 positive semidefinite matrix.

 * @details
 * \b Purpose:
    \verbatim
     PSTRF computes the Cholesky factorization with complete
     pivoting of a real symmetric positive semidefinite matrix A.

     The factorization has the form
        P**T * A * P = U**T * U ,  if UPLO = 'U',
        P**T * A * P = L  * L**T,  if UPLO = 'L',
     where U is an upper triangular matrix and L is lower triangular, and
     P is stored as vector PIV.

     This algorithm does not attempt to check that A is positive
     semidefinite. This version of the algorithm calls level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A. If UPLO = 'U', the leading
          n by n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n by n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, if INFO = 0, the factor U or L from the Cholesky
          factorization as above. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] PIV
          PIV is INTEGER array, dimension (N) \n
          PIV is such that the nonzero entries are P( PIV(K), K) = 1. \n
 * @param[out] RANK
          RANK is INTEGER \n
          The rank of A given by the number of steps the algorithm
          completed. \n
 * @param[in] TOL
          TOL is REAL \n
          User defined tolerance. If TOL < 0, then N*U*MAX( A(K,K))
          will be used. The algorithm terminates at the (K-1)st step
          if the pivot <= TOL. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (2*N) \n
          Work space. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          < 0: If INFO = -K, the K-th argument had an illegal value, \n
          = 0: algorithm completed successfully, and \n
          > 0: the matrix A is either rank deficient with computed rank
               as returned in RANK, or is not positive semidefinite. See
               Section 7 of LAPACK Working Note #161 for further
               information. \n

 *  * */
    template <typename T>
    void pstrf(char *uplo, integer *n, T *a, integer *lda, integer *piv, integer *rank, T *tol,
               T *work, integer *info)
    {
        pstrf(uplo, n, a, lda, piv, rank, tol, work, info);
    }
    template <typename T, typename Ta>
    void pstrf(char *uplo, integer *n, T *a, integer *lda, integer *piv, integer *rank, Ta *tol,
               Ta *work, integer *info)
    {
        pstrf(uplo, n, a, lda, piv, rank, tol, work, info);
    }
    /** @}*/ // end of pstrf

    /** @defgroup pstf2 pstf2
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PSTF2 computes the Cholesky factorization with complete pivoting \n
     of a real symmetric positive semidefinite matrix
 * @details
 * \b Purpose:
    \verbatim
    PSTF2 computes the Cholesky factorization with complete
    pivoting of a real symmetric positive semidefinite matrix A.

    The factorization has the form
       P**T * A * P = U**T * U , if UPLO = 'U',
       P**T * A * P = L  * L**T, if UPLO = 'L',
    where U is an upper triangular matrix and L is lower triangular, and
    P is stored as vector PIV.

    This algorithm does not attempt to check that A is positive
    semidefinite. This version of the algorithm calls level 2 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          n by n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n by n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, if INFO = 0, the factor U or L from the Cholesky
          factorization as above. \n
 * @param[out] PIV
          PIV is INTEGER array, dimension (N) \n
          PIV is such that the nonzero entries are P(PIV(K), K) = 1. \n
 * @param[out] RANK
          RANK is INTEGER \n
          The rank of A given by the number of steps the algorithm
          completed. \n
 * @param[in] TOL
          TOL is REAL \n
          User defined tolerance. If TOL < 0, then N*U*MAX(A(K,K))
          will be used. The algorithm terminates at the (K-1)st step
          if the pivot <= TOL. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] WORK
          WORK is REAL array, dimension (2*N) \n
          Work space. \n
 * @param[out] INFO
          INFO is INTEGER \n
          < 0: If INFO = -K, the K-th argument had an illegal value, \n
          = 0: algorithm completed successfully, and \n
          > 0: the matrix A is either rank deficient with computed rank
               as   returned in RANK, or is not positive semidefinite. See
               Section 7 of LAPACK Working Note #161 for further
               information. \n

 *  * */
    template <typename T>
    void pstf2(char *uplo, integer *n, T *a, integer *lda, integer *piv, integer *rank, T *tol,
               T *work, integer *info)
    {
        pstf2(uplo, n, a, lda, piv, rank, tol, work, info);
    }
    template <typename T, typename Ta>
    void pstf2(char *uplo, integer *n, T *a, integer *lda, integer *piv, integer *rank, Ta *tol,
               Ta *work, integer *info)
    {
        pstf2(uplo, n, a, lda, piv, rank, tol, work, info);
    }
    /** @}*/ // end of pstf2

    /** @defgroup potrs potrs
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief POTRS solves a system of linear equations A*X = B

 * @details
 * \b Purpose:
    \verbatim
     POTRS solves a system of linear equations A*X = B with a symmetric
     positive definite matrix A using the Cholesky factorization
     A = U**T*U or A = L*L**T computed by SPOTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The triangular factor U or L from the Cholesky factorization
          A = U**T*U or A = L*L**T, as computed by SPOTRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void potrs(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *b, integer *ldb,
               integer *info)
    {
        potrs(uplo, n, nrhs, a, lda, b, ldb, info);
    }

    /** @}*/ // end of potrs

    /** @defgroup potri potri
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief Inverse of a real symmetric positive definite matrix.
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of the inverse of a real symmetric positive definite matrix a using the
        Cholesky factorization A = U**T*U or A = L*L**T computed by SPOTRF.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  Upper triangle of a is stored; \n
              = 'L':  Lower triangle of a is stored. \n
    * @param[in] n
              n is integer*  \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the triangular factor U or L from the Cholesky
              factorization A = U**T*U or A = L*L**T, as computed by
              SPOTRF. \n
              On exit, the upper or lower triangle of the (symmetric)
              inverse of A, overwriting the input factor U or L. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,n) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  if INFO = i, the (i,i) element of the factor U or L is
                    zero, and the inverse could not be computed. \n

    *     *  */
    template <typename T>
    void potri(char *uplo, integer *n, T *buff_A, integer *ldim_A, integer *info)
    {
        potri(uplo, n, buff_A, ldim_A, info);
    }
    /** @}*/ // end of potri

    /** @defgroup porfs porfs
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PORFS improves the computed solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
     PORFS improves the computed solution to a system of linear
     equations when the coefficient matrix is symmetric positive definite,
     and provides error bounds and backward error estimates for the
     solution.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of A contains the upper triangular part
          of the matrix A, and the strictly lower triangular part of A
          is not referenced.  If UPLO = 'L', the leading N-by-N lower
          triangular part of A contains the lower triangular part of
          the matrix A, and the strictly upper triangular part of A is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The triangular factor U or L from the Cholesky factorization
          A = U**T*U or A = L*L**T, as computed by SPOTRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by SPOTRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X). \n
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void porfs(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af, integer *ldaf,
               T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr, T *work, integer *iwork,
               integer *info)
    {
        porfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void porfs(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af, integer *ldaf,
               T *b, integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork,
               integer *info)
    {
        porfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }

    /** @}*/ // end of porfs

    /** @defgroup porfsx porfsx
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PORFSX improves the computed solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
     PORFSX improves the computed solution to a system of linear
     equations when the coefficient matrix is symmetric positive
     definite, and provides error bounds and backward error estimates
     for the solution.  In addition to normwise error bound, the code
     provides maximum componentwise error bound if possible.  See
     comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the
     error bounds.

     The original system of linear equations may have been equilibrated
     before calling this routine, as described by arguments EQUED and S
     below. In this case, the solution and error bounds returned are
     for the original unequilibrated system.
    \endverbatim

  * @param[in] UPLO
           UPLO is CHARACTER*1 \n
           = 'U':  Upper triangle of A is stored; \n
           = 'L':  Lower triangle of A is stored. \n
  * @param[in] EQUED
           EQUED is CHARACTER*1 \n
           Specifies the form of equilibration that was done to A
           before calling this routine. This is needed to compute
           the solution and error bounds correctly. \n
             = 'N':  No equilibration \n
             = 'Y':  Both row and column equilibration, i.e., A has been
                     replaced by diag(S) * A * diag(S).
                     The right hand side B has been changed accordingly. \n
  * @param[in] N
           N is INTEGER \n
           The order of the matrix A.  N >= 0. \n
  * @param[in] NRHS
           NRHS is INTEGER \n
           The number of right hand sides, i.e., the number of columns
           of the matrices B and X.  NRHS >= 0. \n
  * @param[in] A
           A is REAL array, dimension (LDA,N) \n
           The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
           upper triangular part of A contains the upper triangular part
           of the matrix A, and the strictly lower triangular part of A
           is not referenced.  If UPLO = 'L', the leading N-by-N lower
           triangular part of A contains the lower triangular part of
           the matrix A, and the strictly upper triangular part of A is
           not referenced. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,N). \n
  * @param[in] AF
           AF is REAL array, dimension (LDAF,N) \n
           The triangular factor U or L from the Cholesky factorization
           A = U**T*U or A = L*L**T, as computed by SPOTRF. \n
  * @param[in] LDAF
           LDAF is INTEGER \n
           The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
  * @param[in,out] S
           S is REAL array, dimension (N) \n
           The scale factors for A.  If EQUED = 'Y', A is multiplied on
           the left and right by diag(S).  S is an input argument if FACT =
           'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED
           = 'Y', each element of S must be positive.  If S is output, each
           element of S is a power of the radix. If S is input, each element
           of S should be a power of the radix to ensure a reliable solution
           and error estimates. Scaling by powers of the radix does not cause
           rounding errors unless the result underflows or overflows.
           Rounding errors during scaling lead to refining with a matrix that
           is not equivalent to the input matrix, producing error estimates
           that may not be reliable. \n
  * @param[in] B
           B is REAL array, dimension (LDB,NRHS) \n
           The right hand side matrix B. \n
  * @param[in] LDB
           LDB is INTEGER \n
           The leading dimension of the array B.  LDB >= fla_max(1,N). \n
  * @param[in,out] X
           X is REAL array, dimension (LDX,NRHS) \n
           On entry, the solution matrix X, as computed by SGETRS.
           On exit, the improved solution matrix X. \n
  * @param[in] LDX
           LDX is INTEGER \n
           The leading dimension of the array X.  LDX >= fla_max(1,N). \n
  * @param[out] RCOND
           RCOND is REAL \n
           Reciprocal scaled condition number. This is an estimate of the
           reciprocal Skeel condition number of the matrix A after
           equilibration (if done).  If this is less than the machine
           precision (in particular, if it is zero), the matrix is singular
           to working precision.  Note that the error may still be small even
           if this number is very small and the matrix appears ill-
           conditioned. \n
  * @param[out] BERR
           BERR is REAL array, dimension (NRHS) \n
           Componentwise relative backward error. This is the
           componentwise relative backward error of each solution vector X(j)
           (i.e., the smallest relative change in any element of A or B that
           makes X(j) an exact solution). \n
  * @param[in] N_ERR_BNDS
           N_ERR_BNDS is INTEGER \n
           Number of error bounds to return for each right hand side
           and each type (normwise or componentwise).  See ERR_BNDS_NORM and
           ERR_BNDS_COMP below. \n
  * @param[out] ERR_BNDS_NORM
           ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) \n
           For each right-hand side, this array contains information about
           various error bounds and condition numbers corresponding to the
           normwise relative error, which is defined as follows: \n
 \n
           Normwise relative error in the ith solution vector: \n

           \f[ {\frac{{max\_j}\ {(abs(XTRUE(j,i)} - X(j,i)))}{{max\_j}\ {abs(X(j,i))}}}\f] \n

           The array is indexed by the type of error information as described
           below. There currently are up to three pieces of information
           returned. \n
 \n
           The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
           right-hand side. \n
 \n
           The second index in ERR_BNDS_NORM(:,err) contains the following
           three fields: \n
           err = 1 "Trust/don't trust" boolean. Trust the answer if the
                    reciprocal condition number is less than the threshold
                    sqrt(n) * slamch('Epsilon'). \n
 \n
           err = 2 "Guaranteed" error bound: The estimated forward error,
                    almost certainly within a factor of 10 of the true error
                    so long as the next entry is greater than the threshold
                    sqrt(n) * slamch('Epsilon'). This error bound should only
                    be trusted if the previous boolean is true. \n
 \n
           err = 3  Reciprocal condition number: Estimated normwise
                    reciprocal condition number.  Compared with the threshold
                    sqrt(n) * slamch('Epsilon') to determine if the error
                    estimate is "guaranteed". These reciprocal condition
                    numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                    appropriately scaled matrix Z.
                    Let Z = S*A, where S scales each row by a power of the
                    radix so all absolute row sums of Z are approximately 1. \n
 \n
           See Lapack Working Note 165 for further details and extra
           cautions. \n
  * @param[out] ERR_BNDS_COMP
           ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
           For each right-hand side, this array contains information about
           various error bounds and condition numbers corresponding to the
           componentwise relative error, which is defined as follows: \n
 \n
           Componentwise relative error in the ith solution vector: \n

           \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))}\f] \n

 \n
           The array is indexed by the right-hand side i (on which the
           componentwise relative error depends), and the type of error
           information as described below. There currently are up to three
           pieces of information returned for each right-hand side. If
           componentwise accuracy is not requested (PARAMS(3) = 0.0), then
           ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
           the first (:,N_ERR_BNDS) entries are returned. \n
 \n
           The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
           right-hand side. \n
 \n
           The second index in ERR_BNDS_COMP(:,err) contains the following
           three fields: \n
           err = 1 "Trust/don't trust" boolean. Trust the answer if the
                    reciprocal condition number is less than the threshold
                    sqrt(n) * slamch('Epsilon'). \n
 \n
           err = 2 "Guaranteed" error bound: The estimated forward error,
                    almost certainly within a factor of 10 of the true error
                    so long as the next entry is greater than the threshold
                    sqrt(n) * slamch('Epsilon'). This error bound should only
                    be trusted if the previous boolean is true. \n
 \n
           err = 3  Reciprocal condition number: Estimated componentwise
                    reciprocal condition number.  Compared with the threshold
                    sqrt(n) * slamch('Epsilon') to determine if the error
                    estimate is "guaranteed". These reciprocal condition
                    numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                    appropriately scaled matrix Z. \n
                    Let Z = S*(A*diag(x)), where x is the solution for the
                    current right-hand side and S scales each row of
                    A*diag(x) by a power of the radix so all absolute row
                    sums of Z are approximately 1. \n
 \n
           See Lapack Working Note 165 for further details and extra
           cautions. \n
  * @param[in] NPARAMS
           NPARAMS is INTEGER \n
           Specifies the number of parameters set in PARAMS.  If <= 0, the
           PARAMS array is never referenced and default values are used. \n
  * @param[in,out] PARAMS
           PARAMS is REAL array, dimension NPARAMS \n
           Specifies algorithm parameters.  If an entry is < 0.0, then
           that entry will be filled with default value used for that
           parameter.  Only positions up to NPARAMS are accessed; defaults
           are used for higher-numbered parameters. \n
 \n
              PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
                   refinement or not. \n
                Default: 1.0 \n
                   = 0.0:  No refinement is performed, and no error bounds are
                           computed. \n
                   = 1.0:  Use the double-precision refinement algorithm,
                           possibly with doubled-single computations if the
                           compilation environment does not support DOUBLE
                           PRECISION. \n
                     (other values are reserved for future use) \n
 \n
              PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
                   computations allowed for refinement. \n
                Default: 10 \n
                Aggressive: Set to 100 to permit convergence using approximate
                            factorizations or factorizations other than LU. If
                            the factorization uses a technique other than
                            Gaussian elimination, the guarantees in
                            err_bnds_norm and err_bnds_comp may no longer be
                            trustworthy. \n
 \n
              PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
                   will attempt to find a solution with small componentwise
                   relative error in the double-precision algorithm.  Positive
                   is true, 0.0 is false. \n
                Default: 1.0 (attempt componentwise convergence) \n
  * @param[out]	WORK
          WORK is REAL array, dimension (4*N) \n
  * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
  * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  Successful exit. The solution to every right-hand side is
           guaranteed. \n
          < 0:  If INFO = -i, the i-th argument had an illegal value \n
          > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
           has been completed, but the factor U is exactly singular, so
           the solution and error bounds could not be computed. RCOND = 0
           is returned. \n
          = N+J: The solution corresponding to the Jth right-hand side is
           not guaranteed. The solutions corresponding to other right-
           hand sides K with K > J may not be guaranteed as well, but
           only the first such right-hand side is reported. If a small
           componentwise error is not requested (PARAMS(3) = 0.0) then
           the Jth right-hand side is the first with a normwise error
           bound that is not guaranteed (the smallest J such
           that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
           the Jth right-hand side is the first with either a normwise or
           componentwise error bound that is not guaranteed (the smallest
           J such that either ERR_BNDS_NORM(J,1) = 0.0 or
           ERR_BNDS_COMP(J,1) = 0.0). See the definition of
           ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
           about all of the right-hand sides check ERR_BNDS_NORM or
           ERR_BNDS_COMP. \n

 *  * */
    template <typename T>
    void porfsx(char *uplo, char *equed, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, T *s, T *b, integer *ldb, T *x, integer *ldx, T *rcond, T *berr,
                integer *n_err_bnds, T *err_bnds_norm, T *err_bnds_comp, integer *nparams,
                T *params, T *work, integer *iwork, integer *info)
    {
        porfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds,
               err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
    }
    template <typename T, typename Ta>
    void porfsx(char *uplo, char *equed, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, Ta *s, T *b, integer *ldb, T *x, integer *ldx, Ta *rcond, Ta *berr,
                integer *n_err_bnds, Ta *err_bnds_norm, Ta *err_bnds_comp, integer *nparams,
                Ta *params, T *work, Ta *rwork, integer *info)
    {
        porfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds,
               err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    }
    /** @}*/ // end of porfsx

    /** @defgroup poequ poequ
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief POEQU computes row and column scalings

 * @details
 * \b Purpose:
    \verbatim
     POEQU computes row and column scalings intended to equilibrate a
     symmetric positive definite matrix A and reduce its condition number
     (with respect to the two-norm).  S contains the scale factors,
     S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with
     elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This
     choice of S puts the condition number of B within a factor N of the
     smallest possible condition number over all possible diagonal
     scalings.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The N-by-N symmetric positive definite matrix whose scaling
          factors are to be computed.  Only the diagonal elements of A
          are referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] S
          S is REAL array, dimension (N) \n
          If INFO = 0, S contains the scale factors for A. \n
 * @param[out] SCOND
          SCOND is REAL \n
          If INFO = 0, S contains the ratio of the smallest S(i) to
          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
          large nor too small, it is not worth scaling by S. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the i-th diagonal element is nonpositive. \n

 *  * */
    template <typename T>
    void poequ(integer *n, T *a, integer *lda, T *s, T *scond, T *amax, integer *info)
    {
        poequ(n, a, lda, s, scond, amax, info);
    }
    template <typename T, typename Ta>
    void poequ(integer *n, T *a, integer *lda, Ta *s, Ta *scond, Ta *amax, integer *info)
    {
        poequ(n, a, lda, s, scond, amax, info);
    }
    /** @}*/ // end of poequ

    /** @defgroup  poequb poequb
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief POEQUB computes row and column scalings

 * @details
 * \b Purpose:
    \verbatim
     POEQUB computes row and column scalings intended to equilibrate a
     symmetric positive definite matrix A and reduce its condition number
     (with respect to the two-norm).  S contains the scale factors,
     S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with
     elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This
     choice of S puts the condition number of B within a factor N of the
     smallest possible condition number over all possible diagonal
     scalings.

     This routine differs from SPOEQU by restricting the scaling factors
     to a power of the radix.  Barring over- and underflow, scaling by
     these factors introduces no additional rounding errors.  However, the
     scaled diagonal entries are no longer approximately 1 but lie
     between sqrt(radix) and 1/sqrt(radix).
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The N-by-N symmetric positive definite matrix whose scaling
          factors are to be computed.  Only the diagonal elements of A
          are referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] S
          S is REAL array, dimension (N) \n
          If INFO = 0, S contains the scale factors for A. \n
 * @param[out] SCOND
          SCOND is REAL \n
          If INFO = 0, S contains the ratio of the smallest S(i) to
          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
          large nor too small, it is not worth scaling by S. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the i-th diagonal element is nonpositive. \n

 *  * */
    template <typename T>
    void poequb(integer *n, T *a, integer *lda, T *s, T *scond, T *amax, integer *info)
    {
        poequb(n, a, lda, s, scond, amax, info);
    }
    template <typename T, typename Ta>
    void poequb(integer *n, T *a, integer *lda, Ta *s, Ta *scond, Ta *amax, integer *info)
    {
        poequb(n, a, lda, s, scond, amax, info);
    }

    /** @}*/ // end of poequb
    /** @defgroup laqhe laqhe
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief LAQHE scales a Hermitian matrix

 * @details
 * \b Purpose:
    \verbatim
    LAQHE equilibrates a Hermitian matrix A using the scaling factors
    in the vector S.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          Hermitian matrix A is stored. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
          n by n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n by n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, if EQUED = 'Y', the equilibrated matrix:
          diag(S) * A * diag(S). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(N,1). \n
 * @param[in] S
          S is REAL array, dimension (N) \n
          The scale factors for A. \n
 * @param[in] SCOND
          SCOND is REAL \n
          Ratio of the smallest S(i) to the largest S(i). \n
 * @param[in] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix entry. \n
 * @param[out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies whether or not equilibration was done. \n
          = 'N':  No equilibration. \n
          = 'Y':  Equilibration was done, i.e., A has been replaced by
                  diag(S) * A * diag(S).   \n

 * */
    template <typename T, typename Ta>
    void laqhe(char *uplo, integer *n, T *a, integer *lda, Ta *s, Ta *scond, Ta *amax, char *equed)
    {
        laqhe(uplo, n, a, lda, s, scond, amax, equed);
    }
    /** @}*/ // end of laqhe

    /** @defgroup la_porcond la_porcond
     * @{
     * @ingroup Cholesky Cholesky
     */

    /*! @brief LA_PORCOND estimates the Skeel condition number for a symmetric positive-definite
     matrix

     * @details
     * \b Purpose:
        \verbatim
        LA_PORCOND Estimates the Skeel condition number of  op(A) * op2(C)
        where op2 is determined by CMODE as follows
        CMODE =  1    op2(C) = C
        CMODE =  0    op2(C) = I
        CMODE = -1    op2(C) = inv(C)
        The Skeel condition number  cond(A) = norminf(|inv(A)||A|)
        is computed by computing scaling factors R such that
        diag(R)*A*op2(C) is row equilibrated and computing the standard
        infinity-norm condition number.
        \endverbatim

     * @param[in] UPLO
              UPLO is CHARACTER*1 \n
              = 'U':  Upper triangle of A is stored; \n
              = 'L':  Lower triangle of A is stored. \n
     * @param[in] N
              N is INTEGER \n
              The number of linear equations, i.e., the order of the
              matrix A.  N >= 0. \n
     * @param[in] A
              A is REAL array, dimension (LDA,N) \n
              On entry, the N-by-N matrix A. \n
     * @param[in] LDA
              LDA is INTEGER \n
              The leading dimension of the array A.  LDA >= fla_max(1,N). \n
     * @param[in] AF
              AF is REAL array, dimension (LDAF,N) \n
              The triangular factor U or L from the Cholesky factorization
              A = U**T*U or A = L*L**T, as computed by SPOTRF. \n
     * @param[in] LDAF
              LDAF is INTEGER \n
              The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
     * @param[in] CMODE
              CMODE is INTEGER \n
              Determines op2(C) in the formula op(A) * op2(C) as follows: \n
              CMODE =  1    op2(C) = C \n
              CMODE =  0    op2(C) = I \n
              CMODE = -1    op2(C) = inv(C) \n
     * @param[in] C
              C is REAL array, dimension (N) \n
              The vector C in the formula op(A) * op2(C). \n
     * @param[out] INFO
              INFO is INTEGER \n
              = 0:  Successful exit. \n
              i > 0:  The ith argument is invalid. \n
     * @param[out] WORK
              WORK is REAL array, dimension (3*N). \n
              Workspace. \n
     * @param[out] IWORK
              IWORK is INTEGER array, dimension (N). \n
              Workspace. \n

     *  * */
    template <typename T>
    T la_porcond(char *uplo, integer *n, T *a, integer *lda, T *af, integer *ldaf, integer *cmode,
                 T *c, integer *info, T *work, integer *iwork)
    {
        return la_porcond(uplo, n, a, lda, af, ldaf, cmode, c, info, work, iwork);
    }
    /** @}*/ // end of la_porcond

    /** @defgroup la_porcondx la_porcondx
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief LA_PORCOND_X computes the infinity norm condition number of op(A)*diag(x) for
 Hermitian positive-definite matrices

 * @details
 * \b Purpose:
    \verbatim
    LA_PORCOND_X Computes the infinity norm condition number of
    op(A) * diag(X) where X is a COMPLEX vector.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is COMPLEX array, dimension (LDAF,N) \n
          The triangular factor U or L from the Cholesky factorization
          A = U**H*U or A = L*L**H, as computed by CPOTRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] X
          X is COMPLEX array, dimension (N) \n
          The vector X in the formula op(A) * diag(X). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          i > 0:  The ith argument is invalid. \n
 * @param[out] WORK
          WORK is COMPLEX array, dimension (2*N). \n
          Workspace. \n
 * @param[out] RWORK
          RWORK is REAL array, dimension (N). \n
          Workspace. \n

 *  * */
    template <typename T, typename Ta>
    Ta la_porcond_x(char *uplo, integer *n, T *a, integer *lda, T *af, integer *ldaf, T *x,
                    integer *info, T *work, Ta *rwork)
    {
        return la_porcond_x(uplo, n, a, lda, af, ldaf, x, info, work, rwork);
    }
    /** @}*/ // end of la_porcondx

    /** @defgroup la_porpvgrw la_porpvgrw
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief LA_PORPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a
 symmetric or Hermitian positive-definite matrix

 * @details
 * \b Purpose:
    \verbatim
    LA_PORPVGRW computes the reciprocal pivot growth factor
    norm(A)/norm(U). The "max absolute element" norm is used. If this is
    much less than 1, the stability of the LU factorization of the
    (equilibrated) matrix A could be poor. This also means that the
    solution X, estimated condition numbers, and error bounds could be
    unreliable.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] NCOLS
          NCOLS is INTEGER \n
          The number of columns of the matrix A. NCOLS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The triangular factor U or L from the Cholesky factorization
          A = U**T*U or A = L*L**T, as computed by SPOTRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[out] WORK
          WORK is REAL array, dimension (2*N) \n

 *  * */
    template <typename T>
    T la_porpvgrw(char *uplo, integer *ncols, T *a, integer *lda, T *af, integer *ldaf, T *work)
    {
        return la_porpvgrw(uplo, ncols, a, lda, af, ldaf, work);
    }
    template <typename T, typename Ta>
    Ta la_porpvgrw(char *uplo, integer *ncols, T *a, integer *lda, T *af, integer *ldaf, Ta *work)
    {
        return la_porpvgrw(uplo, ncols, a, lda, af, ldaf, work);
    }
    /** @}*/ // end of la_porpvgrw

    /** @defgroup ppcon ppcon
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PPCON estimates the reciprocal of the condition number

 * @details
 * \b Purpose:
    \verbatim
     PPCON estimates the reciprocal of the condition number (in the
     1-norm) of a real symmetric positive definite packed matrix using
     the Cholesky factorization A = U**T*U or A = L*L**T computed by
     PPTRF.

     An estimate is obtained for norm(inv(A)), and the reciprocal of the
     condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The triangular factor U or L from the Cholesky factorization
          A = U**T*U or A = L*L**T, packed columnwise in a linear
          array.  The j-th column of U or L is stored in the array AP
          as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm (or infinity-norm) of the symmetric matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void ppcon(char *uplo, integer *n, T *ap, T *anorm, T *rcond, T *work, integer *iwork,
               integer *info)
    {
        ppcon(uplo, n, ap, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void ppcon(char *uplo, integer *n, T *ap, Ta *anorm, Ta *rcond, T *work, Ta *rwork,
               integer *info)
    {
        ppcon(uplo, n, ap, anorm, rcond, work, rwork, info);
    }
    /** @}*/ // end of ppcon ppcon

    /** @defgroup pptrf pptrf
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PPTRF computes the Cholesky factorization of a real symmetric matrix

 * @details
 * \b Purpose:
    \verbatim
     PPTRF computes the Cholesky factorization of a real symmetric
     positive definite matrix A stored in packed format.

     The factorization has the form
        A = U**T * U,  if UPLO = 'U', or
        A = L  * L**T,  if UPLO = 'L',
     where U is an upper triangular matrix and L is lower triangular.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
          See below for further details. \n
 \n
          On exit, if INFO = 0, the triangular factor U or L from the
          Cholesky factorization A = U**T*U or A = L*L**T, in the same
          storage format as A. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the leading minor of order i is not
                positive definite, and the factorization could not be
                completed. \n

 *  * */
    template <typename T>
    void pptrf(char *uplo, integer *n, T *ap, integer *info)
    {
        pptrf(uplo, n, ap, info);
    }
    /** @}*/ // end of pptrf

    /** @defgroup pptrs pptrs
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PPTRS solves a system of linear equations A*X = B with a symmetric matrix

 * @details
 * \b Purpose:
    \verbatim
     PPTRS solves a system of linear equations A*X = B with a symmetric
     positive definite matrix A in packed storage using the Cholesky
     factorization A = U**T*U or A = L*L**T computed by PPTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The triangular factor U or L from the Cholesky factorization
          A = U**T*U or A = L*L**T, packed columnwise in a linear
          array.  The j-th column of U or L is stored in the array AP
          as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void pptrs(char *uplo, integer *n, integer *nrhs, T *ap, T *b, integer *ldb, integer *info)
    {
        pptrs(uplo, n, nrhs, ap, b, ldb, info);
    }
    /** @}*/ // end of pptrs

    /** @defgroup pptri pptri
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PPTRI computes the inverse of a real symmetric matrix

 * @details
 * \b Purpose:
    \verbatim
     PPTRI computes the inverse of a real symmetric positive definite
     matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
     computed by PPTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangular factor is stored in AP; \n
          = 'L':  Lower triangular factor is stored in AP. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the triangular factor U or L from the Cholesky
          factorization A = U**T*U or A = L*L**T, packed columnwise as
          a linear array. The j-th column of U or L is stored in the
          array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n. \n
 \n
          On exit, the upper or lower triangle of the (symmetric)
          inverse of A, overwriting the input factor U or L. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the (i,i) element of the factor U or L is
                zero, and the inverse could not be computed. \n

 *  * */
    template <typename T>
    void pptri(char *uplo, integer *n, T *ap, integer *info)
    {
        pptri(uplo, n, ap, info);
    }
    /** @}*/ // end of pptri

    /** @defgroup pprfs pprfs
   * @ingroup Cholesky Cholesky
   * @{
   * /
* @details
* \b Purpose:
  \verbatim
   PPRFS improves the computed solution to a system of linear
   equations when the coefficient matrix is symmetric positive definite
   and packed, and provides error bounds and backward error estimates
   for the solution.
  \endverbatim

* @param[in] UPLO
        UPLO is CHARACTER*1 \n
        = 'U':  Upper triangle of A is stored; \n
        = 'L':  Lower triangle of A is stored. \n
* @param[in] N
        N is INTEGER \n
        The order of the matrix A.  N >= 0. \n
* @param[in] NRHS
        NRHS is INTEGER \n
        The number of right hand sides, i.e., the number of columns
        of the matrices B and X.  NRHS >= 0. \n
* @param[in] AP
        AP is REAL array, dimension (N*(N+1)/2) \n
        The upper or lower triangle of the symmetric matrix A, packed
        columnwise in a linear array.  The j-th column of A is stored
        in the array AP as follows: \n
        if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
        if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
* @param[in] AFP
        AFP is REAL array, dimension (N*(N+1)/2) \n
        The triangular factor U or L from the Cholesky factorization
        A = U**T*U or A = L*L**T, as computed by SPPTRF/CPPTRF,
        packed columnwise in a linear array in the same format as A
        (see AP). \n
* @param[in] B
        B is REAL array, dimension (LDB,NRHS) \n
        The right hand side matrix B. \n
* @param[in] LDB
        LDB is INTEGER \n
        The leading dimension of the array B.  LDB >= fla_max(1,N). \n
* @param[in,out] X
        X is REAL array, dimension (LDX,NRHS) \n
        On entry, the solution matrix X, as computed by SPPTRS.
        On exit, the improved solution matrix X. \n
* @param[in] LDX
        LDX is INTEGER \n
        The leading dimension of the array X.  LDX >= fla_max(1,N). \n
* @param[out] FERR
        FERR is REAL array, dimension (NRHS) \n
        The estimated forward error bound for each solution vector
        X(j) (the j-th column of the solution matrix X). \n
        If XTRUE is the true solution corresponding to X(j), FERR(j)
        is an estimated upper bound for the magnitude of the largest
        element in (X(j) - XTRUE) divided by the magnitude of the
        largest element in X(j).  The estimate is as reliable as
        the estimate for RCOND, and is almost always a slight
        overestimate of the true error. \n
* @param[out] BERR
        BERR is REAL array, dimension (NRHS) \n
        The componentwise relative backward error of each solution
        vector X(j) (i.e., the smallest relative change in
        any element of A or B that makes X(j) an exact solution). \n
* @param[out]	WORK
        WORK is REAL array, dimension (3*N) \n
* @param[out]	IWORK
        IWORK is INTEGER array, dimension (N) \n
* @param[out]	INFO
        INFO is INTEGER \n
        = 0:  successful exit \n
        < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void pprfs(char *uplo, integer *n, integer *nrhs, T *ap, T *afp, T *b, integer *ldb, T *x,
               integer *ldx, T *ferr, T *berr, T *work, integer *iwork, integer *info)
    {
        pprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void pprfs(char *uplo, integer *n, integer *nrhs, T *ap, T *afp, T *b, integer *ldb, T *x,
               integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork, integer *info)
    {
        pprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }
    /** @}*/ // end of pprfs

    /** @defgroup ppequ ppequ
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PPEQU computes row and column scalings

 * @details
 * \b Purpose:
    \verbatim
     PPEQU computes row and column scalings intended to equilibrate a
     symmetric positive definite matrix A in packed storage and reduce
     its condition number (with respect to the two-norm).  S contains the
     scale factors, S(i)=1/sqrt(A(i,i)), chosen so that the scaled matrix
     B with elements B(i,j)=S(i)*A(i,j)*S(j) has ones on the diagonal.
     This choice of S puts the condition number of B within a factor N of
     the smallest possible condition number over all possible diagonal
     scalings.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangle of the symmetric matrix A, packed
          columnwise in a linear array.  The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[out] S
          S is REAL array, dimension (N) \n
          If INFO = 0, S contains the scale factors for A. \n
 * @param[out] SCOND
          SCOND is REAL \n
          If INFO = 0, S contains the ratio of the smallest S(i) to
          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
          large nor too small, it is not worth scaling by S. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the i-th diagonal element is nonpositive. \n

 *  * */
    template <typename T>
    void ppequ(char *uplo, integer *n, T *ap, T *s, T *scond, T *amax, integer *info)
    {
        ppequ(uplo, n, ap, s, scond, amax, info);
    }
    template <typename T, typename Ta>
    void ppequ(char *uplo, integer *n, T *ap, Ta *s, Ta *scond, Ta *amax, integer *info)
    {
        ppequ(uplo, n, ap, s, scond, amax, info);
    }
    /** @}*/ // end of ppequ

    /** @defgroup laqhp laqhp
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief LAQHP scales a Hermitian matrix stored in packed form

 * @details
 * \b Purpose:
    \verbatim
    LAQHP equilibrates a Hermitian matrix A using the scaling factors
    in the vector S.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          Hermitian matrix A is stored. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] AP
          AP is COMPLEX array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the Hermitian matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 \n
          On exit, the equilibrated matrix:  diag(S) * A * diag(S), in
          the same storage format as A. \n
 * @param[in] S
          S is REAL array, dimension (N) \n
          The scale factors for A. \n
 * @param[in] SCOND
          SCOND is REAL \n
          Ratio of the smallest S(i) to the largest S(i). \n
 * @param[in] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix entry. \n
 * @param[out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies whether or not equilibration was done. \n
          = 'N':  No equilibration. \n
          = 'Y':  Equilibration was done, i.e., A has been replaced by
                  diag(S) * A * diag(S).  \n

 * */
    template <typename T, typename Ta>
    void laqhp(char *uplo, integer *n, T *ap, Ta *s, Ta *scond, Ta *amax, char *equed)
    {
        laqhp(uplo, n, ap, s, scond, amax, equed);
    }

    /** @}*/ // end of laqhp

    /** @defgroup petrf petrf
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PFTRF computes the Cholesky factorization

 * @details
 * \b Purpose:
    \verbatim
     PFTRF computes the Cholesky factorization of a real symmetric
     positive definite matrix A.

     The factorization has the form
        A = U**T * U,  if UPLO = 'U', or
        A = L  * L**T,  if UPLO = 'L',
     where U is an upper triangular matrix and L is lower triangular.

     This is the block version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  The Normal TRANSR of RFP A is stored; \n
          = 'T':  The Transpose TRANSR of RFP A is stored. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of RFP A is stored; \n
          = 'L':  Lower triangle of RFP A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension ( N*(N+1)/2); \n
          On entry, the symmetric matrix A in RFP format. RFP format is
          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'
          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is
          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'T' then RFP is
          the transpose of RFP A as defined when
          TRANSR = 'N'. The contents of RFP A are defined by UPLO as
          follows: If UPLO = 'U' the RFP A contains the NT elements of
          upper packed A. If UPLO = 'L' the RFP A contains the elements
          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR =
          'T'. When TRANSR is 'N' the LDA is N+1 when N is even and N
          is odd. See the Note below for more details. \n
 \n
          On exit, if INFO = 0, the factor U or L from the Cholesky
          factorization RFP A = U**T*U or RFP A = L*L**T. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the leading minor of order i is not
                positive definite, and the factorization could not be
                completed. \n

 *  * */
    template <typename T>
    void pftrf(char *transr, char *uplo, integer *n, T *a, integer *info)
    {
        pftrf(transr, uplo, n, a, info);
    }
    /** @}*/ // end of pftrf

    /** @defgroup pfrs pftrs
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PFTRS solves a system of linear equations A*X = B with a symmetric positive definite
 matrix

 * @details
 * \b Purpose:
    \verbatim
     PFTRS solves a system of linear equations A*X = B with a symmetric
     positive definite matrix A using the Cholesky factorization
     A = U**T*U or A = L*L**T computed by PFTRF.
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  The Normal TRANSR of RFP A is stored; \n
          = 'T':  The Transpose TRANSR of RFP A is stored. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of RFP A is stored; \n
          = 'L':  Lower triangle of RFP A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension ( N*(N+1)/2) \n
          The triangular factor U or L from the Cholesky factorization
          of RFP A = U**H*U or RFP A = L*L**T, as computed by SPFTRF.
          See note below for more details about RFP A. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void pftrs(char *transr, char *uplo, integer *n, integer *nrhs, T *a, T *b, integer *ldb,
               integer *info)
    {
        pftrs(transr, uplo, n, nrhs, a, b, ldb, info);
    }
    /** @}*/ // end of pftrs

    /** @defgroup pftri pftri
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PFTRI computes the inverse of a real (symmetric) positive definite matrix

 * @details
 * \b Purpose:
    \verbatim
     PFTRI computes the inverse of a real (symmetric) positive definite
     matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
     computed by PFTRF.
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  The Normal TRANSR of RFP A is stored; \n
          = 'T':  The Transpose TRANSR of RFP A is stored. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension ( N*(N+1)/2) \n
          On entry, the symmetric matrix A in RFP format. RFP format is
          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'
          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is
          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'T' then RFP is
          the transpose of RFP A as defined when
          TRANSR = 'N'. The contents of RFP A are defined by UPLO as
          follows: If UPLO = 'U' the RFP A contains the nt elements of
          upper packed A. If UPLO = 'L' the RFP A contains the elements
          of lower packed A. The LDA of RFP A is (N+1)/2 when TRANSR =
          'T'. When TRANSR is 'N' the LDA is N+1 when N is even and N
          is odd. See the Note below for more details. \n
 \n
          On exit, the symmetric inverse of the original matrix, in the
          same storage format. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the (i,i) element of the factor U or L is
                zero, and the inverse could not be computed. \n

 *  * */
    template <typename T>
    void pftri(char *transr, char *uplo, integer *n, T *a, integer *info)
    {
        pftri(transr, uplo, n, a, info);
    }
    /** @}*/ // end of pftri

    /** @defgroup pbcon pbcon
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PBCON estimates the reciprocal of the condition number

 * @details
 * \b Purpose:
    \verbatim
     PBCON estimates the reciprocal of the condition number (in the
     1-norm) of a real symmetric positive definite band matrix using the
     Cholesky factorization A = U**T*U or A = L*L**T computed by SPBTRF.

     An estimate is obtained for norm(inv(A)), and the reciprocal of the
     condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangular factor stored in AB; \n
          = 'L':  Lower triangular factor stored in AB. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The triangular factor U or L from the Cholesky factorization
          A = U**T*U or A = L*L**T of the band matrix A, stored in the
          first KD+1 rows of the array.  The j-th column of U or L is
          stored in the j-th column of the array AB as follows:
          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for fla_max(1,j-kd)<=i<=j;
          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm (or infinity-norm) of the symmetric band matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void pbcon(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, T *anorm, T *rcond,
               T *work, integer *iwork, integer *info)
    {
        pbcon(uplo, n, kd, ab, ldab, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void pbcon(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, Ta *anorm, Ta *rcond,
               T *work, Ta *rwork, integer *info)
    {
        pbcon(uplo, n, kd, ab, ldab, anorm, rcond, work, rwork, info);
    }

    /** @}*/ // end of pbcon

    /** @defgroup  pbtrf pbtrf
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PBTRF computes the Cholesky factorization of a real symmetric  \n
     positive definite band matrix A.

 * @details
 * \b Purpose:
    \verbatim
     PBTRF computes the Cholesky factorization of a real symmetric
     positive definite band matrix A.

     The factorization has the form
        A = U**T * U,  if UPLO = 'U', or
        A = L  * L**T,  if UPLO = 'L',
     where U is an upper triangular matrix and L is lower triangular.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix A, stored in the first KD+1 rows of the array.  The
          j-th column of A is stored in the j-th column of the array AB
          as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
 \n
          On exit, if INFO = 0, the triangular factor U or L from the
          Cholesky factorization A = U**T*U or A = L*L**T of the band
          matrix A, in the same storage format as A. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the leading minor of order i is not
                positive definite, and the factorization could not be
                completed. \n

 *  * */
    template <typename T>
    void pbtrf(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, integer *info)
    {
        pbtrf(uplo, n, kd, ab, ldab, info);
    }

    /** @}*/ // end of pbtrf

    /** @defgroup pbequ pbequ
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PBEQU computes row and column scalings

 * @details
 * \b Purpose:
    \verbatim
     PBEQU computes row and column scalings intended to equilibrate a
     symmetric positive definite band matrix A and reduce its condition
     number (with respect to the two-norm).  S contains the scale factors,
     S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with
     elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This
     choice of S puts the condition number of B within a factor N of the
     smallest possible condition number over all possible diagonal
     scalings.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangular of A is stored; \n
          = 'L':  Lower triangular of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of superdiagonals of the matrix A if UPLO = 'U',
          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The upper or lower triangle of the symmetric band matrix A,
          stored in the first KD+1 rows of the array.  The j-th column
          of A is stored in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array A.  LDAB >= KD+1. \n
 * @param[out] S
          S is REAL array, dimension (N) \n
          If INFO = 0, S contains the scale factors for A. \n
 * @param[out] SCOND
          SCOND is REAL \n
          If INFO = 0, S contains the ratio of the smallest S(i) to
          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too
          large nor too small, it is not worth scaling by S. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, the i-th diagonal element is nonpositive. \n

 *  * */
    template <typename T>
    void pbequ(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, T *s, T *scond, T *amax,
               integer *info)
    {
        pbequ(uplo, n, kd, ab, ldab, s, scond, amax, info);
    }
    template <typename T, typename Ta>
    void pbequ(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, T *s, Ta *scond, Ta *amax,
               integer *info)
    {
        pbequ(uplo, n, kd, ab, ldab, s, scond, amax, info);
    }
    /** @}*/ // end of pbequ

    /** @defgroup laqhb laqhb
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief LAQHB scales a Hermitian band matrix, using scaling factors computed by cpbequ

 * @details
 * \b Purpose:
    \verbatim
    LAQHB equilibrates an Hermitian band matrix A using the scaling
    factors in the vector S.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of super-diagonals of the matrix A if UPLO = 'U',
          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0. \n
 * @param[in,out] AB
          AB is COMPLEX array, dimension (LDAB,N) \n
          On entry, the upper or lower triangle of the symmetric band
          matrix A, stored in the first KD+1 rows of the array.  The
          j-th column of A is stored in the j-th column of the array AB
          as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
 \n
          On exit, if INFO = 0, the triangular factor U or L from the
          Cholesky factorization A = U**H *U or A = L*L**H of the band
          matrix A, in the same storage format as A. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[out] S
          S is REAL array, dimension (N) \n
          The scale factors for A. \n
 * @param[in] SCOND
          SCOND is REAL \n
          Ratio of the smallest S(i) to the largest S(i). \n
 * @param[in] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix entry. \n
 * @param[out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies whether or not equilibration was done. \n
          = 'N':  No equilibration. \n
          = 'Y':  Equilibration was done, i.e., A has been replaced by
                  diag(S) * A * diag(S).  \n

 * */
    template <typename T, typename Ta>
    void laqhb(char *uplo, integer *n, integer *kd, T *ab, integer *ldab, Ta *s, Ta *scond,
               Ta *amax, char *equed)
    {
        laqhb(uplo, n, kd, ab, ldab, s, scond, amax, equed);
    }
    /** @}*/ // end of laqhb

    /** @defgroup ptcon ptcon
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PTCON computes the reciprocal of the condition number

 * @details
 * \b Purpose:
    \verbatim
     PTCON computes the reciprocal of the condition number (in the
     1-norm) of a real symmetric positive definite tridiagonal matrix
     using the factorization A = L*D*L**T or A = U**T*D*U computed by
     SPTTRF.

     Norm(inv(A)) is computed by a direct method, and the reciprocal of
     the condition number is computed as
           RCOND = 1 / (ANORM * norm(inv(A))).
     \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the diagonal matrix D from the
          factorization of A, as computed by SPTTRF. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) off-diagonal elements of the unit bidiagonal factor
          U or L from the factorization of A,  as computed by SPTTRF. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the
          1-norm of inv(A) computed in this routine. \n
 * @param[out] RWORK
          RWORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void ptcon(integer *n, T *d, T *e, T *anorm, T *rcond, T *rwork, integer *info)
    {
        ptcon(n, d, e, anorm, rcond, rwork, info);
    }
    template <typename T, typename Ta>
    void ptcon(integer *n, Ta *d, T *e, Ta *anorm, Ta *rcond, Ta *rwork, integer *info)
    {
        ptcon(n, d, e, anorm, rcond, rwork, info);
    }
    /** @}*/ // end of ptcon

    /** @defgroup pttrf pptrf
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PTTRF computes the L*D*L**T factorization of a real symmetric matrix

 * @details
 * \b Purpose:
    \verbatim
     PTTRF computes the L*D*L**T factorization of a real symmetric
     positive definite tridiagonal matrix A.  The factorization may also
     be regarded as having the form A = U**T*D*U.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the n diagonal elements of the tridiagonal matrix
          A.  On exit, the n diagonal elements of the diagonal matrix
          D from the L*D*L**T factorization of A. \n
 * @param[in,out] E
          E is REAL array, dimension (N-1) \n
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix A.  On exit, the (n-1) subdiagonal elements of the
          unit bidiagonal factor L from the L*D*L**T factorization of A. \n
          E can also be regarded as the superdiagonal of the unit
          bidiagonal factor U from the U**T*D*U factorization of A. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -k, the k-th argument had an illegal value \n
          > 0: if INFO = k, the leading minor of order k is not
               positive definite; if k < N, the factorization could not
               be completed, while if k = N, the factorization was
               completed, but D(N) <= 0. \n

 *  * */
    template <typename T>
    void pttrf(integer *n, T *d, T *e, integer *info)
    {
        pttrf(n, d, e, info);
    }
    template <typename T, typename Ta>
    void pttrf(integer *n, Ta *d, T *e, integer *info)
    {
        pttrf(n, d, e, info);
    }

    /** @}*/ // end of pttrf

    /** @defgroup pttrs pttrs
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PTTRS solves a tridiagonal system of the form A * X = B

 * @details
 * \b Purpose:
    \verbatim
     PTTRS solves a tridiagonal system of the form
        A * X = B
     using the L*D*L**T factorization of A computed by SPTTRF.  D is a
     diagonal matrix specified in the vector D, L is a unit bidiagonal
     matrix whose subdiagonal is specified in the vector E, and X and B
     are N by NRHS matrices.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the tridiagonal matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the diagonal matrix D from the
          L*D*L**T factorization of A. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of the unit bidiagonal factor
          L from the L*D*L**T factorization of A.  E can also be regarded
          as the superdiagonal of the unit bidiagonal factor U from the
          factorization A = U**T*D*U. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side vectors B for the system of
          linear equations. \n
          On exit, the solution vectors, X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N).  \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -k, the k-th argument had an illegal value \n

 *  * */
    template <typename T>
    void pttrs(integer *n, integer *nrhs, T *d, T *e, T *b, integer *ldb, integer *info)
    {
        pttrs(n, nrhs, d, e, b, ldb, info);
    }
    template <typename T, typename Ta>
    void pttrs(char *uplo, integer *n, integer *nrhs, Ta *d, T *e, T *b, integer *ldb,
               integer *info)
    {
        pttrs(uplo, n, nrhs, d, e, b, ldb, info);
    }
    /** @}*/ // end of pttrs

    /** @defgroup ptts2 ptts2
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PTTS2 solves a tridiagonal system of the form AX=B using \n
     the L D LH factorization computed by spttrf
 * @details
 * \b Purpose:
    \verbatim
    PTTS2 solves a tridiagonal system of the form
       A * X = B
    using the L*D*L**T factorization of A computed by SPTTRF.  D is a
    diagonal matrix specified in the vector D, L is a unit bidiagonal
    matrix whose subdiagonal is specified in the vector E, and X and B
    are N by NRHS matrices.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the tridiagonal matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the diagonal matrix D from the
          L*D*L**T factorization of A. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of the unit bidiagonal factor
          L from the L*D*L**T factorization of A.  E can also be regarded
          as the superdiagonal of the unit bidiagonal factor U from the
          factorization A = U**T*D*U. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side vectors B for the system of
          linear equations. \n
          On exit, the solution vectors, X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n

 *  * */
    template <typename T>
    void ptts2(integer *n, integer *nrhs, T *d, T *e, T *b, integer *ldb)
    {
        ptts2(n, nrhs, d, e, b, ldb);
    }
    template <typename T, typename Ta>
    void ptts2(integer *iuplo, integer *n, integer *nrhs, Ta *d, T *e, T *b, integer *ldb)
    {
        ptts2(iuplo, n, nrhs, d, e, b, ldb);
    }
    /** @}*/ // end of ptts2

    /** @defgroup  ptrfs ptrfs
     * @{
     * @ingroup Cholesky Cholesky
     */
    /*! @brief PTRFS improves the computed solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
    PTRFS improves the computed solution to a system of linear
    equations when the coefficient matrix is symmetric positive definite
    and tridiagonal, and provides error bounds and backward error
    estimates for the solution.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the tridiagonal matrix A. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of the tridiagonal matrix A. \n
 * @param[in] DF
          DF is REAL array, dimension (N) \n
          The n diagonal elements of the diagonal matrix D from the
          factorization computed by SPTTRF. \n
 * @param[in] EF
          EF is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of the unit bidiagonal factor
          L from the factorization computed by SPTTRF. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by SPTTRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X). \n
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j). \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void ptrfs(integer *n, integer *nrhs, T *d, T *e, T *df, T *ef, T *b, integer *ldb, T *x,
               integer *ldx, T *ferr, T *berr, T *work, integer *info)
    {
        ptrfs(n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, info);
    }
    template <typename T, typename Ta>
    void ptrfs(char *uplo, integer *n, integer *nrhs, Ta *d, T *e, Ta *df, T *ef, T *b,
               integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork,
               integer *info)
    {
        ptrfs(uplo, n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }
    /** @}*/ // end of ptrfs
    /** @}*/ // end of Cholesky

    /** @defgroup LU LU General matrix, driver
     * @ingroup LinearSolve
     * @{
     */
    /** @defgroup gesv gesv
     * @{
     * @ingroup LU LU
     */
    /*! @brief GESV computes the solution to a real system of linear equations
 * @details
 * \b Purpose:
    \verbatim
     GESV computes the solution to a real system of linear equations
        A * X = B,
     where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

     The LU decomposition with partial pivoting and row interchanges is
     used to factor A as
        A = P * L * U,
     where P is a permutation matrix, L is unit lower triangular, and U is
     upper triangular.  The factored form of A is then used to solve the
     system of equations A * X = B.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @parm[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N coefficient matrix A.
          On exit, the factors L and U from the factorization
          A = P*L*U; the unit diagonal elements of L are not stored. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices that define the permutation matrix P;
          row i of the matrix was interchanged with row IPIV(i). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS matrix of right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
                has been completed, but the factor U is exactly
                singular, so the solution could not be computed. \n

 *  * */
    template <typename T>
    void gesv(integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b, integer *ldb,
              integer *info)
    {
        gesv(n, nrhs, a, lda, ipiv, b, ldb, info);
    }

    /** @}*/ // end of gesv

    /** @defgroup gesvx gesvx
     * @{
     * @ingroup LU LU
     */
    /*! @brief GESVX computes the solution to system of linear equations A * X = B for GE matrices
 * @details
 * \b Purpose:
    \verbatim
     GESVX uses the LU factorization to compute the solution to a real
     system of linear equations
        A * X = B,
     where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

     Error bounds on the solution and a condition estimate are also
     provided.
    \endverbatim

 * @param[in] FACT
          FACT is CHARACTER*1 \n
          Specifies whether or not the factored form of the matrix A is
          supplied on entry, and if not, whether the matrix A should be
          equilibrated before it is factored. \n
          = 'F':  On entry, AF and IPIV contain the factored form of A.
                  If EQUED is not 'N', the matrix A has been
                  equilibrated with scaling factors given by R and C.
                  A, AF, and IPIV are not modified. \n
          = 'N':  The matrix A will be copied to AF and factored. \n
          = 'E':  The matrix A will be equilibrated if necessary, then
                  copied to AF and factored. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
          not 'N', then A must have been equilibrated by the scaling
          factors in R and/or C.  A is not modified if FACT = 'F' or
          'N', or if FACT = 'E' and EQUED = 'N' on exit. \n
 \n
          On exit, if EQUED .ne. 'N', A is scaled as follows: \n
          EQUED = 'R':  A := diag(R) * A \n
          EQUED = 'C':  A := A * diag(C) \n
          EQUED = 'B':  A := diag(R) * A * diag(C). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] AF
          AF is REAL array, dimension (LDAF,N) \n
          If FACT = 'F', then AF is an input argument and on entry
          contains the factors L and U from the factorization
          A = P*L*U as computed by GETRF.  If EQUED .ne. 'N', then
          AF is the factored form of the equilibrated matrix A. \n
 \n
          If FACT = 'N', then AF is an output argument and on exit
          returns the factors L and U from the factorization A = P*L*U
          of the original matrix A. \n
 \n
          If FACT = 'E', then AF is an output argument and on exit
          returns the factors L and U from the factorization A = P*L*U
          of the equilibrated matrix A (see the description of A for
          the form of the equilibrated matrix). \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          If FACT = 'F', then IPIV is an input argument and on entry
          contains the pivot indices from the factorization A = P*L*U
          as computed by GETRF; row i of the matrix was interchanged
          with row IPIV(i). \n
 \n
          If FACT = 'N', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = P*L*U
          of the original matrix A. \n
 \n
          If FACT = 'E', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = P*L*U
          of the equilibrated matrix A. \n
 * @param[in,out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done. \n
          = 'N':  No equilibration (always true if FACT = 'N'). \n
          = 'R':  Row equilibration, i.e., A has been premultiplied by
                  diag(R). \n
          = 'C':  Column equilibration, i.e., A has been postmultiplied
                  by diag(C). \n
          = 'B':  Both row and column equilibration, i.e., A has been
                  replaced by diag(R) * A * diag(C). \n
          EQUED is an input argument if FACT = 'F'; otherwise, it is an
          output argument. \n
 * @param[in,out] R
          R is REAL array, dimension (N) \n
          The row scale factors for A.  If EQUED = 'R' or 'B', A is
          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
          is not accessed.  R is an input argument if FACT = 'F';
          otherwise, R is an output argument.  If FACT = 'F' and
          EQUED = 'R' or 'B', each element of R must be positive. \n
 * @param[in,out] C
          C is REAL array, dimension (N) \n
          The column scale factors for A.  If EQUED = 'C' or 'B', A is
          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
          is not accessed.  C is an input argument if FACT = 'F';
          otherwise, C is an output argument.  If FACT = 'F' and
          EQUED = 'C' or 'B', each element of C must be positive. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B. \n
          On exit, \n
          if EQUED = 'N', B is not modified; \n
          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
          diag(R)*B; \n
          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
          overwritten by diag(C)*B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
          X is REAL array, dimension (LDX,NRHS) \n
          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
          to the original system of equations.  Note that A and B are
          modified on exit if EQUED .ne. 'N', and the solution to the
          equilibrated system is inv(diag(C))*X if TRANS = 'N' and
          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
          and EQUED = 'R' or 'B'. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          The estimate of the reciprocal condition number of the matrix
          A after equilibration (if done).  If RCOND is less than the
          machine precision (in particular, if RCOND = 0), the matrix
          is singular to working precision.  This condition is
          indicated by a return code of INFO > 0. \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*N) \n
          On exit, RWORK(1) contains the reciprocal pivot growth
          factor norm(A)/norm(U). The "max absolute element" norm is
          used. If RWORK(1) is much less than 1, then the stability
          of the LU factorization of the (equilibrated) matrix A
          could be poor. This also means that the solution X, condition
          estimator RCOND, and forward error bound FERR could be
          unreliable. If factorization fails with 0<INFO<=N, then
          RWORK(1) contains the reciprocal pivot growth factor for the
          leading INFO columns of A. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, and i is \n
                <= N:  U(i,i) is exactly zero.  The factorization has
                       been completed, but the factor U is exactly
                       singular, so the solution and error bounds
                       could not be computed. RCOND = 0 is returned. \n
                = N+1: U is nonsingular, but RCOND is less than machine
                       precision, meaning that the matrix is singular
                       to working precision.  Nevertheless, the
                       solution and error bounds are computed because
                       there are a number of situations where the
                       computed solution can be more accurate than the
                       value of RCOND would suggest. \n

 *  * */
    template <typename T>
    void gesvx(char *fact, char *trans, integer *n, integer *nrhs, T *a, integer *lda, T *af,
               integer *ldaf, integer *ipiv, char *equed, T *r, T *c, T *b, integer *ldb, T *x,
               integer *ldx, T *rcond, T *ferr, T *berr, T *work, integer *iwork, integer *info)
    {
        gesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond,
              ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void gesvx(char *fact, char *trans, integer *n, integer *nrhs, T *a, integer *lda, T *af,
               integer *ldaf, integer *ipiv, char *equed, Ta *r, Ta *c, T *b, integer *ldb, T *x,
               integer *ldx, Ta *rcond, Ta *ferr, Ta *berr, T *work, Ta *rwork, integer *info)
    {
        gesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond,
              ferr, berr, work, rwork, info);
    }

    /** @}*/ // end of gesvx

    /** @defgroup gesvxx gesvxx
     * @{
     * @ingroup LU LU
     */
    /*! @brief GESVXX computes the solution to system of linear equations A * X = B for GE matrices

 * @details
 * \b Purpose:
    \verbatim
     GESVXX uses the LU factorization to compute the solution to a
     real system of linear equations  A * X = B,  where A is an
     N-by-N matrix and X and B are N-by-NRHS matrices.

     If requested, both normwise and maximum componentwise error bounds
     are returned. GESVXX will return a solution with a tiny
     guaranteed error (O(eps) where eps is the working machine
     precision) unless the matrix is very ill-conditioned, in which
     case a warning is returned. Relevant condition numbers also are
     calculated and returned.

     GESVXX accepts user-provided factorizations and equilibration
     factors; see the definitions of the FACT and EQUED options.
     Solving with refinement and using a factorization from a previous
     GESVXX call will also produce a solution with either O(eps)
     errors or warnings, but we cannot make that claim for general
     user-provided factorizations and equilibration factors if they
     differ from what GESVXX would itself produce.
    \endverbatim

 * @param[in] FACT
          FACT is CHARACTER*1 \n
          Specifies whether or not the factored form of the matrix A is
          supplied on entry, and if not, whether the matrix A should be
          equilibrated before it is factored. \n
           = 'F':  On entry, AF and IPIV contain the factored form of A.
                   If EQUED is not 'N', the matrix A has been
                   equilibrated with scaling factors given by R and C.
                   A, AF, and IPIV are not modified. \n
           = 'N':  The matrix A will be copied to AF and factored. \n
           = 'E':  The matrix A will be equilibrated if necessary, then
                   copied to AF and factored. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
           = 'N':  A * X = B     (No transpose) \n
           = 'T':  A**T * X = B  (Transpose) \n
           = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
          not 'N', then A must have been equilibrated by the scaling
          factors in R and/or C.  A is not modified if FACT = 'F' or
          'N', or if FACT = 'E' and EQUED = 'N' on exit. \n
 \n
          On exit, if EQUED .ne. 'N', A is scaled as follows: \n
          EQUED = 'R':  A := diag(R) * A \n
          EQUED = 'C':  A := A * diag(C) \n
          EQUED = 'B':  A := diag(R) * A * diag(C). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] AF
          AF is REAL array, dimension (LDAF,N) \n
          If FACT = 'F', then AF is an input argument and on entry
          contains the factors L and U from the factorization
          A = P*L*U as computed by GETRF.  If EQUED .ne. 'N', then
          AF is the factored form of the equilibrated matrix A. \n
 \n
          If FACT = 'N', then AF is an output argument and on exit
          returns the factors L and U from the factorization A = P*L*U
          of the original matrix A. \n
 \n
          If FACT = 'E', then AF is an output argument and on exit
          returns the factors L and U from the factorization A = P*L*U
          of the equilibrated matrix A (see the description of A for
          the form of the equilibrated matrix). \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          If FACT = 'F', then IPIV is an input argument and on entry
          contains the pivot indices from the factorization A = P*L*U
          as computed by GETRF; row i of the matrix was interchanged
          with row IPIV(i). \n
 \n
          If FACT = 'N', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = P*L*U
          of the original matrix A. \n
 \n
          If FACT = 'E', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = P*L*U
          of the equilibrated matrix A. \n
 * @param[in,out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done. \n
           = 'N':  No equilibration (always true if FACT = 'N'). \n
           = 'R':  Row equilibration, i.e., A has been premultiplied by
                   diag(R). \n
           = 'C':  Column equilibration, i.e., A has been postmultiplied
                   by diag(C). \n
           = 'B':  Both row and column equilibration, i.e., A has been
                   replaced by diag(R) * A * diag(C). \n
          EQUED is an input argument if FACT = 'F'; otherwise, it is an
          output argument. \n
 * @param[in,out] R
          R is REAL array, dimension (N) \n
          The row scale factors for A.  If EQUED = 'R' or 'B', A is
          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
          is not accessed.  R is an input argument if FACT = 'F';
          otherwise, R is an output argument.  If FACT = 'F' and
          EQUED = 'R' or 'B', each element of R must be positive. \n
          If R is output, each element of R is a power of the radix.
          If R is input, each element of R should be a power of the radix
          to ensure a reliable solution and error estimates. Scaling by
          powers of the radix does not cause rounding errors unless the
          result underflows or overflows. Rounding errors during scaling
          lead to refining with a matrix that is not equivalent to the
          input matrix, producing error estimates that may not be
          reliable. \n
 * @param[in,out] C
          C is REAL array, dimension (N) \n
          The column scale factors for A. If EQUED = 'C' or 'B', A is
          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
          is not accessed.  C is an input argument if FACT = 'F';
          otherwise, C is an output argument.  If FACT = 'F' and
          EQUED = 'C' or 'B', each element of C must be positive. \n
          If C is output, each element of C is a power of the radix. \n
          If C is input, each element of C should be a power of the radix
          to ensure a reliable solution and error estimates. Scaling by
          powers of the radix does not cause rounding errors unless the
          result underflows or overflows. Rounding errors during scaling
          lead to refining with a matrix that is not equivalent to the
          input matrix, producing error estimates that may not be
          reliable. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B. \n
          On exit, \n
          if EQUED = 'N', B is not modified; \n
          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
            diag(R)*B; \n
          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
            overwritten by diag(C)*B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
          X is REAL array, dimension (LDX,NRHS) \n
          If INFO = 0, the N-by-NRHS solution matrix X to the original
          system of equations. Note that A and B are modified on exit
          if EQUED .ne. 'N', and the solution to the equilibrated system is
          inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or
          inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          Reciprocal scaled condition number. This is an estimate of the
          reciprocal Skeel condition number of the matrix A after
          equilibration (if done). If this is less than the machine
          precision (in particular, if it is zero), the matrix is singular
          to working precision.  Note that the error may still be small even
          if this number is very small and the matrix appears ill-
          conditioned. \n
 * @param[out] RPVGRW
          RPVGRW is REAL \n
          Reciprocal pivot growth.  On exit, this contains the reciprocal
          pivot growth factor norm(A)/norm(U). The "max absolute element"
          norm is used.  If this is much less than 1, then the stability of
          the LU factorization of the (equilibrated) matrix A could be poor. \n
          This also means that the solution X, estimated condition numbers,
          and error bounds could be unreliable. If factorization fails with
          0<INFO<=N, then this contains the reciprocal pivot growth factor
          for the leading INFO columns of A.  In GESVX, this quantity is
          returned in WORK(1). \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS)
          Componentwise relative backward error. This is the
          componentwise relative backward error of each solution vector X(j)
          (i.e., the smallest relative change in any element of A or B that
          makes X(j) an exact solution). \n
 * @param[in] N_ERR_BNDS
          N_ERR_BNDS is INTEGER \n
          Number of error bounds to return for each right hand side
          and each type (normwise or componentwise). See ERR_BNDS_NORM and
          ERR_BNDS_COMP below. \n
 * @param[out] ERR_BNDS_NORM
          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          normwise relative error, which is defined as follows: \n
 \n
          Normwise relative error in the ith solution vector: \n

          \f[\frac{{max\_j}\  (abs(XTRUE(j,i) - X(j,i)))} {{max\_j}\ abs(X(j,i))} \f] \n

 \n
          The array is indexed by the type of error information as described
          below. There currently are up to three pieces of information
          returned. \n
 \n
          The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_NORM(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated normwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*A, where S scales each row by a power of the
                  radix so all absolute row sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[out] ERR_BNDS_COMP
          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          componentwise relative error, which is defined as follows: \n
 \n
          Componentwise relative error in the ith solution vector: \n

          \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

          The array is indexed by the right-hand side i (on which the
          componentwise relative error depends), and the type of error
          information as described below. There currently are up to three
          pieces of information returned for each right-hand side. If
          componentwise accuracy is not requested (PARAMS(3) = 0.0), then
          ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
          the first (:,N_ERR_BNDS) entries are returned. \n
 \n
          The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_COMP(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated componentwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z. \n
                  Let Z = S*(A*diag(x)), where x is the solution for the
                  current right-hand side and S scales each row of
                  A*diag(x) by a power of the radix so all absolute row
                  sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in] NPARAMS
          NPARAMS is INTEGER \n
          Specifies the number of parameters set in PARAMS.  If <= 0, the
          PARAMS array is never referenced and default values are used. \n
 * @param[in,out] PARAMS
          PARAMS is REAL array, dimension NPARAMS \n
          Specifies algorithm parameters.  If an entry is < 0.0, then
          that entry will be filled with default value used for that
          parameter.  Only positions up to NPARAMS are accessed; defaults
          are used for higher-numbered parameters. \n
 \n
           PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
                refinement or not. \n
             Default: 1.0 \n
                = 0.0:  No refinement is performed, and no error bounds are
                        computed. \n
                = 1.0:  Use the double-precision refinement algorithm,
                        possibly with doubled-single computations if the
                        compilation environment does not support DOUBLE
                        PRECISION.
                  (other values are reserved for future use) \n
 \n
           PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
                computations allowed for refinement. \n
             Default: 10 \n
             Aggressive: Set to 100 to permit convergence using approximate
                         factorizations or factorizations other than LU. If
                         the factorization uses a technique other than
                         Gaussian elimination, the guarantees in
                         err_bnds_norm and err_bnds_comp may no longer be
                         trustworthy. \n
 \n
           PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
                will attempt to find a solution with small componentwise
                relative error in the double-precision algorithm.  Positive
                is true, 0.0 is false. \n
             Default: 1.0 (attempt componentwise convergence) \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  Successful exit. The solution to every right-hand side is
           guaranteed. \n
          < 0:  If INFO = -i, the i-th argument had an illegal value \n
          > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
           has been completed, but the factor U is exactly singular, so
           the solution and error bounds could not be computed. RCOND = 0
           is returned. \n
          = N+J: The solution corresponding to the Jth right-hand side is
           not guaranteed. The solutions corresponding to other right-
           hand sides K with K > J may not be guaranteed as well, but
           only the first such right-hand side is reported. If a small
           componentwise error is not requested (PARAMS(3) = 0.0) then
           the Jth right-hand side is the first with a normwise error
           bound that is not guaranteed (the smallest J such
           that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
           the Jth right-hand side is the first with either a normwise or
           componentwise error bound that is not guaranteed (the smallest
           J such that either ERR_BNDS_NORM(J,1) = 0.0 or
           ERR_BNDS_COMP(J,1) = 0.0). See the definition of
           ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
           about all of the right-hand sides check ERR_BNDS_NORM or
           ERR_BNDS_COMP. \n

 *  * */
    template <typename T>
    void gesvxx(char *fact, char *trans, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, char *equed, T *r, T *c, T *b, integer *ldb, T *x,
                integer *ldx, T *rcond, T *rpvgrw, T *berr, integer *n_err_bnds, T *err_bnds_norm,
                T *err_bnds_comp, integer *nparams, T *params, T *work, integer *iwork,
                integer *info)
    {
        gesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond,
               rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork,
               info);
    }
    template <typename T, typename Ta>
    void gesvxx(char *fact, char *trans, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, char *equed, Ta *r, Ta *c, T *b, integer *ldb, T *x,
                integer *ldx, Ta *rcond, Ta *rpvgrw, Ta *berr, integer *n_err_bnds,
                Ta *err_bnds_norm, Ta *err_bnds_comp, integer *nparams, Ta *params, T *work,
                Ta *rwork, integer *info)
    {
        gesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, rcond,
               rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork,
               info);
    }

    /** @}*/ // end of gesvxx

    /** @defgroup gbsv gbsv
     * @{
     * @ingroup LU LU
     */
    /*! @brief  GBSV computes the solution to system of linear equations A * X = B for GB matrices
*
* @details
* \b Purpose:
* \verbatim
  GBSV computes the solution to a real system of linear equations
  A * X = B, where A is a band matrix of order N with KL subdiagonals
  and KU superdiagonals, and X and B are N-by-NRHS matrices.
  \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the matrix A in band storage, in rows KL+1 to
          2*KL+KU+1; rows 1 to KL of the array need not be set.
          The j-th column of A is stored in the j-th column of the
          array AB as follows: \n
          AB(KL+KU+1+i-j,j) = A(i,j) for fla_max(1,j-KU)<=i<=min(N,j+KL) \n
          On exit, details of the factorization: U is stored as an
          upper triangular band matrix with KL+KU superdiagonals in
          rows 1 to KL+KU+1, and the multipliers used during the
          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
          See below for further details. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices that define the permutation matrix P;
          row i of the matrix was interchanged with row IPIV(i). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B. \n
          On exit, if INFO = 0, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
                has been completed, but the factor U is exactly
                singular, and the solution has not been computed. \n

 *  * */
    template <typename T>
    void gbsv(integer *n, integer *kl, integer *ku, integer *nrhs, T *ab, integer *ldab,
              integer *ipiv, T *b, integer *ldb, integer *info)
    {
        gbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    }
    /** @}*/ // end of gbsv

    /** @defgroup gbsvx gbsvx
     * @{
     * @ingroup LU LU
     */
    /*! @brief It computes the solution to system of linear equations A * X = B for GB matrices
along with error bounds
*
*  @details
* \b Purpose:
  \verbatim
   GBSVX uses the LU factorization to compute the solution to a real
   system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
   where A is a band matrix of order N with KL subdiagonals and KU
   superdiagonals, and X and B are N-by-NRHS matrices.

   Error bounds on the solution and a condition estimate are also
   provided.
   \endverbatim

 * @param[in] FACT
          FACT is CHARACTER*1 \n
          Specifies whether or not the factored form of the matrix A is
          supplied on entry, and if not, whether the matrix A should be
          equilibrated before it is factored. \n
          = 'F':  On entry, AFB and IPIV contain the factored form of
                  A.  If EQUED is not 'N', the matrix A has been
                  equilibrated with scaling factors given by R and C.
                  AB, AFB, and IPIV are not modified. \n
          = 'N':  The matrix A will be copied to AFB and factored. \n
          = 'E':  The matrix A will be equilibrated if necessary, then
                  copied to AFB and factored. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations. \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
          The j-th column of A is stored in the j-th column of the
          array AB as follows: \n
          AB(KU+1+i-j,j) = A(i,j) for fla_max(1,j-KU)<=i<=min(N,j+kl) \n
 \n
          If FACT = 'F' and EQUED is not 'N', then A must have been
          equilibrated by the scaling factors in R and/or C.  AB is not
          modified if FACT = 'F' or 'N', or if FACT = 'E' and
          EQUED = 'N' on exit. \n
 \n
          On exit, if EQUED .ne. 'N', A is scaled as follows: \n
          EQUED = 'R':  A := diag(R) * A \n
          EQUED = 'C':  A := A * diag(C) \n
          EQUED = 'B':  A := diag(R) * A * diag(C). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[in,out] AFB
          AFB is REAL array, dimension (LDAFB,N) \n
          If FACT = 'F', then AFB is an input argument and on entry
          contains details of the LU factorization of the band matrix
          A, as computed by SGBTRF.  U is stored as an upper triangular
          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
          and the multipliers used during the factorization are stored
          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is
          the factored form of the equilibrated matrix A. \n
 \n
          If FACT = 'N', then AFB is an output argument and on exit
          returns details of the LU factorization of A. \n
 \n
          If FACT = 'E', then AFB is an output argument and on exit
          returns details of the LU factorization of the equilibrated
          matrix A (see the description of AB for the form of the
          equilibrated matrix). \n
 * @param[in] LDAFB
          LDAFB is INTEGER \n
          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          If FACT = 'F', then IPIV is an input argument and on entry
          contains the pivot indices from the factorization A = L*U
          as computed by SGBTRF; row i of the matrix was interchanged
          with row IPIV(i). \n
 \n
          If FACT = 'N', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = L*U
          of the original matrix A. \n
 \n
          If FACT = 'E', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = L*U
          of the equilibrated matrix A. \n
 * @param[in,out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done. \n
          = 'N':  No equilibration (always true if FACT = 'N'). \n
          = 'R':  Row equilibration, i.e., A has been premultiplied by
                  diag(R). \n
          = 'C':  Column equilibration, i.e., A has been postmultiplied
                  by diag(C). \n
          = 'B':  Both row and column equilibration, i.e., A has been
                  replaced by diag(R) * A * diag(C). \n
          EQUED is an input argument if FACT = 'F'; otherwise, it is an
          output argument. \n
 * @param[in,out] R
          R is REAL array, dimension (N) \n
          The row scale factors for A.  If EQUED = 'R' or 'B', A is
          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
          is not accessed.  R is an input argument if FACT = 'F';
          otherwise, R is an output argument.  If FACT = 'F' and
          EQUED = 'R' or 'B', each element of R must be positive. \n
 * @param[in,out] C
          C is REAL array, dimension (N) \n
          The column scale factors for A.  If EQUED = 'C' or 'B', A is
          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
          is not accessed.  C is an input argument if FACT = 'F';
          otherwise, C is an output argument.  If FACT = 'F' and
          EQUED = 'C' or 'B', each element of C must be positive. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B. \n
          On exit, \n
          if EQUED = 'N', B is not modified; \n
          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
          diag(R)*B; \n
          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
          overwritten by diag(C)*B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
          X is REAL array, dimension (LDX,NRHS) \n
          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
          to the original system of equations.  Note that A and B are
          modified on exit if EQUED .ne. 'N', and the solution to the
          equilibrated system is inv(diag(C))*X if TRANS = 'N' and
          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
          and EQUED = 'R' or 'B'. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          The estimate of the reciprocal condition number of the matrix
          A after equilibration (if done). If RCOND is less than the
          machine precision (in particular, if RCOND = 0), the matrix
          is singular to working precision.  This condition is
          indicated by a return code of INFO > 0. \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (N) \n
          On exit, RWORK(1) contains the reciprocal pivot growth
          factor norm(A)/norm(U). The "max absolute element" norm is
          used. If RWORK(1) is much less than 1, then the stability
          of the LU factorization of the (equilibrated) matrix A
          could be poor. This also means that the solution X, condition
          estimator RCOND, and forward error bound FERR could be
          unreliable. If factorization fails with 0<INFO<=N, then
          RWORK(1) contains the reciprocal pivot growth factor for the
          leading INFO columns of A. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, and i is \n
                <= N:  U(i,i) is exactly zero.  The factorization
                       has been completed, but the factor U is exactly
                       singular, so the solution and error bounds
                       could not be computed. RCOND = 0 is returned. \n
                = N+1: U is nonsingular, but RCOND is less than machine
                       precision, meaning that the matrix is singular
                       to working precision.  Nevertheless, the
                       solution and error bounds are computed because
                       there are a number of situations where the
                       computed solution can be more accurate than the
                       value of RCOND would suggest. \n

 *  * */
    template <typename T>
    void gbsvx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
               integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, T *r, T *c, T *b,
               integer *ldb, T *x, integer *ldx, T *rcond, T *ferr, T *berr, T *work,
               integer *iwork, integer *info)
    {
        gbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx,
              rcond, ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void gbsvx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
               integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, Ta *r, Ta *c,
               T *b, integer *ldb, T *x, integer *ldx, Ta *rcond, Ta *ferr, Ta *berr, T *work,
               Ta *rwork, integer *info)
    {
        gbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx,
              rcond, ferr, berr, work, rwork, info);
    }
    /** @}*/ // end of gbsvx

    /** @defgroup gbsvxx gbsvxx
     * @{
     * @ingroup LU LU
     */
    /*! @brief  GBSVXX computes the solution to system of linear equations A * X = B for GB matrices
*
* @details
* \b Purpose
* \verbatim
    GBSVXX uses the LU factorization to compute the solution to a
    real system of linear equations  A * X = B,  where A is an
    N-by-N matrix and X and B are N-by-NRHS matrices.

    If requested, both normwise and maximum componentwise error bounds
    are returned. SGBSVXX will return a solution with a tiny
    guaranteed error (O(eps) where eps is the working machine
    precision) unless the matrix is very ill-conditioned, in which
    case a warning is returned. Relevant condition numbers also are
    calculated and returned.

    SGBSVXX accepts user-provided factorizations and equilibration
    factors; see the definitions of the FACT and EQUED options.
    Solving with refinement and using a factorization from a previous
    SGBSVXX call will also produce a solution with either O(eps)
    errors or warnings, but we cannot make that claim for general
    user-provided factorizations and equilibration factors if they
    differ from what SGBSVXX would itself produce.
   \endverbatim

 * @param[in] FACT
        FACT is CHARACTER*1 \n
        Specifies whether or not the factored form of the matrix A is
        supplied on entry, and if not, whether the matrix A should be
        equilibrated before it is factored. \n
         = 'F':  On entry, AF and IPIV contain the factored form of A.
                 If EQUED is not 'N', the matrix A has been
                 equilibrated with scaling factors given by R and C.
                 A, AF, and IPIV are not modified. \n
         = 'N':  The matrix A will be copied to AF and factored. \n
         = 'E':  The matrix A will be equilibrated if necessary, then
                 copied to AF and factored. \n
 * @param[in] TRANS
        TRANS is CHARACTER*1 \n
        Specifies the form of the system of equations: \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) \n
 * @param[in] N
        N is INTEGER \n
        The number of linear equations, i.e., the order of the
        matrix A.  N >= 0. \n
 * @param[in] KL
        KL is INTEGER \n
        The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
        KU is INTEGER \n
        The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] NRHS
        NRHS is INTEGER \n
        The number of right hand sides, i.e., the number of columns
        of the matrices B and X.  NRHS >= 0. \n
 * @param[in,out] AB
        AB is REAL array, dimension (LDAB,N) \n
        On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
        The j-th column of A is stored in the j-th column of the
        array AB as follows: \n
        AB(KU+1+i-j,j) = A(i,j) for fla_max(1,j-KU)<=i<=min(N,j+kl) \n
 \n
        If FACT = 'F' and EQUED is not 'N', then AB must have been
        equilibrated by the scaling factors in R and/or C.  AB is not
        modified if FACT = 'F' or 'N', or if FACT = 'E' and
        EQUED = 'N' on exit. \n
 \n
        On exit, if EQUED .ne. 'N', A is scaled as follows: \n
        EQUED = 'R':  A := diag(R) * A \n
        EQUED = 'C':  A := A * diag(C) \n
        EQUED = 'B':  A := diag(R) * A * diag(C). \n
 * @param[in] LDAB
        LDAB is INTEGER \n
        The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[in,out] AFB
        AFB is REAL array, dimension (LDAFB,N) \n
        If FACT = 'F', then AFB is an input argument and on entry
        contains details of the LU factorization of the band matrix
        A, as computed by SGBTRF.  U is stored as an upper triangular
        band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
        and the multipliers used during the factorization are stored
        in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is
        the factored form of the equilibrated matrix A. \n
 \n
        If FACT = 'N', then AF is an output argument and on exit
        returns the factors L and U from the factorization A = P*L*U
        of the original matrix A. \n
 \n
        If FACT = 'E', then AF is an output argument and on exit
        returns the factors L and U from the factorization A = P*L*U
        of the equilibrated matrix A (see the description of A for
        the form of the equilibrated matrix). \n
 * @param[in] LDAFB
        LDAFB is INTEGER \n
        The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. \n
 * @param[in,out] IPIV
        IPIV is INTEGER array, dimension (N) \n
        If FACT = 'F', then IPIV is an input argument and on entry
        contains the pivot indices from the factorization A = P*L*U
        as computed by SGETRF; row i of the matrix was interchanged
        with row IPIV(i). \n
 \n
        If FACT = 'N', then IPIV is an output argument and on exit
        contains the pivot indices from the factorization A = P*L*U
        of the original matrix A. \n
 \n
        If FACT = 'E', then IPIV is an output argument and on exit
        contains the pivot indices from the factorization A = P*L*U
        of the equilibrated matrix A. \n
 * @param[in,out] EQUED
        EQUED is CHARACTER*1 \n
        Specifies the form of equilibration that was done. \n
         = 'N':  No equilibration (always true if FACT = 'N'). \n
         = 'R':  Row equilibration, i.e., A has been premultiplied by
                 diag(R). \n
         = 'C':  Column equilibration, i.e., A has been postmultiplied
                 by diag(C). \n
         = 'B':  Both row and column equilibration, i.e., A has been
                 replaced by diag(R) * A * diag(C). \n
        EQUED is an input argument if FACT = 'F'; otherwise, it is an
        output argument. \n
 * @param[in,out] R
        R is REAL array, dimension (N) \n
        The row scale factors for A.  If EQUED = 'R' or 'B', A is
        multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
        is not accessed.  R is an input argument if FACT = 'F';
        otherwise, R is an output argument.  If FACT = 'F' and
        EQUED = 'R' or 'B', each element of R must be positive. \n
        If R is output, each element of R is a power of the radix. \n
        If R is input, each element of R should be a power of the radix
        to ensure a reliable solution and error estimates. Scaling by
        powers of the radix does not cause rounding errors unless the
        result underflows or overflows. Rounding errors during scaling
        lead to refining with a matrix that is not equivalent to the
        input matrix, producing error estimates that may not be
        reliable. \n
 * @param[in,out] C
        C is REAL array, dimension (N) \n
        The column scale factors for A.  If EQUED = 'C' or 'B', A is
        multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
        is not accessed.  C is an input argument if FACT = 'F';
        otherwise, C is an output argument.  If FACT = 'F' and
        EQUED = 'C' or 'B', each element of C must be positive. \n
        If C is output, each element of C is a power of the radix. \n
        If C is input, each element of C should be a power of the radix
        to ensure a reliable solution and error estimates. Scaling by
        powers of the radix does not cause rounding errors unless the
        result underflows or overflows. Rounding errors during scaling
        lead to refining with a matrix that is not equivalent to the
        input matrix, producing error estimates that may not be
        reliable. \n
 * @param[in,out] B
        B is REAL array, dimension (LDB,NRHS) \n
        On entry, the N-by-NRHS right hand side matrix B. \n
        On exit, \n
        if EQUED = 'N', B is not modified; \n
        if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
          diag(R)*B; \n
        if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
          overwritten by diag(C)*B. \n
 * @param[in] LDB
        LDB is INTEGER \n
        The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
        X is REAL array, dimension (LDX,NRHS) \n
        If INFO = 0, the N-by-NRHS solution matrix X to the original
        system of equations.  Note that A and B are modified on exit
        if EQUED .ne. 'N', and the solution to the equilibrated system is
        inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or
        inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'. \n
 * @param[in] LDX
        LDX is INTEGER \n
        The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
        RCOND is REAL \n
        Reciprocal scaled condition number. This is an estimate of the
        reciprocal Skeel condition number of the matrix A after
        equilibration (if done).  If this is less than the machine
        precision (in particular, if it is zero), the matrix is singular
        to working precision.  Note that the error may still be small even
        if this number is very small and the matrix appears ill-
        conditioned. \n
 * @param[out] RPVGRW
        RPVGRW is REAL \n
        Reciprocal pivot growth.  On exit, this contains the reciprocal
        pivot growth factor norm(A)/norm(U). The "max absolute element"
        norm is used.  If this is much less than 1, then the stability of
        the LU factorization of the (equilibrated) matrix A could be poor.
        This also means that the solution X, estimated condition numbers,
        and error bounds could be unreliable. If factorization fails with
        0<INFO<=N, then this contains the reciprocal pivot growth factor
        for the leading INFO columns of A.  In SGESVX, this quantity is
        returned in WORK(1). \n
 * @param[out] BERR
        BERR is REAL array, dimension (NRHS) \n
        Componentwise relative backward error.  This is the
        componentwise relative backward error of each solution vector X(j)
        (i.e., the smallest relative change in any element of A or B that
        makes X(j) an exact solution). \n
 * @param[in] N_ERR_BNDS
        N_ERR_BNDS is INTEGER \n
        Number of error bounds to return for each right hand side
        and each type (normwise or componentwise).  See ERR_BNDS_NORM and
        ERR_BNDS_COMP below. \n
 * @param[out] ERR_BNDS_NORM
        ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS)
        For each right-hand side, this array contains information about
        various error bounds and condition numbers corresponding to the
        normwise relative error, which is defined as follows: \n
 \n
        Normwise relative error in the ith solution vector: \n

        \f[\frac{{max\_j}\  (abs(XTRUE(j,i) - X(j,i)))} {{max\_j}\ abs(X(j,i))} \f] \n

 \n
        The array is indexed by the type of error information as described
        below. There currently are up to three pieces of information
        returned. \n
 \n
        The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
        right-hand side. \n
 \n
        The second index in ERR_BNDS_NORM(:,err) contains the following
        three fields: \n
        err = 1 "Trust/don't trust" boolean. Trust the answer if the
                reciprocal condition number is less than the threshold
                sqrt(n) * slamch('Epsilon'). \n
 \n
        err = 2 "Guaranteed" error bound: The estimated forward error,
                almost certainly within a factor of 10 of the true error
                so long as the next entry is greater than the threshold
                sqrt(n) * slamch('Epsilon'). This error bound should only
                be trusted if the previous boolean is true. \n
 \n
        err = 3  Reciprocal condition number: Estimated normwise
                reciprocal condition number.  Compared with the threshold
                sqrt(n) * slamch('Epsilon') to determine if the error
                estimate is "guaranteed". These reciprocal condition
                numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                appropriately scaled matrix Z.
                Let Z = S*A, where S scales each row by a power of the
                radix so all absolute row sums of Z are approximately 1. \n
 \n
        See Lapack Working Note 165 for further details and extra
        cautions. \n
 * @param[out] ERR_BNDS_COMP
        ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
        For each right-hand side, this array contains information about
        various error bounds and condition numbers corresponding to the
        componentwise relative error, which is defined as follows: \n
 \n
        Componentwise relative error in the ith solution vector: \n

        \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
        The array is indexed by the right-hand side i (on which the
        componentwise relative error depends), and the type of error
        information as described below. There currently are up to three
        pieces of information returned for each right-hand side. If
        componentwise accuracy is not requested (PARAMS(3) = 0.0), then
        ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
        the first (:,N_ERR_BNDS) entries are returned. \n
 \n
        The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
        right-hand side. \n
 \n
        The second index in ERR_BNDS_COMP(:,err) contains the following
        three fields: \n
        err = 1 "Trust/don't trust" boolean. Trust the answer if the
                reciprocal condition number is less than the threshold
                sqrt(n) * slamch('Epsilon'). \n
 \n
        err = 2 "Guaranteed" error bound: The estimated forward error,
                almost certainly within a factor of 10 of the true error
                so long as the next entry is greater than the threshold
                sqrt(n) * slamch('Epsilon'). This error bound should only
                be trusted if the previous boolean is true. \n
 \n
        err = 3  Reciprocal condition number: Estimated componentwise
                reciprocal condition number.  Compared with the threshold
                sqrt(n) * slamch('Epsilon') to determine if the error
                estimate is "guaranteed". These reciprocal condition
                numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                appropriately scaled matrix Z.
                Let Z = S*(A*diag(x)), where x is the solution for the
                current right-hand side and S scales each row of
                A*diag(x) by a power of the radix so all absolute row
                sums of Z are approximately 1. \n
 \n
        See Lapack Working Note 165 for further details and extra
        cautions. \n
 * @param[in] NPARAMS
        NPARAMS is INTEGER \n
        Specifies the number of parameters set in PARAMS.  If <= 0, the
        PARAMS array is never referenced and default values are used. \n
 * @param[in,out] PARAMS
        PARAMS is REAL array, dimension NPARAMS \n
        Specifies algorithm parameters.  If an entry is < 0.0, then
        that entry will be filled with default value used for that
        parameter.  Only positions up to NPARAMS are accessed; defaults
        are used for higher-numbered parameters. \n
 \n
        PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
            refinement or not. \n
         Default: 1.0 \n
            = 0.0:  No refinement is performed, and no error bounds are
                    computed. \n
            = 1.0:  Use the double-precision refinement algorithm,
                    possibly with doubled-single computations if the
                    compilation environment does not support DOUBLE
                    PRECISION.
              (other values are reserved for future use) \n
 \n
        PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
            computations allowed for refinement. \n
         Default: 10 \n
         Aggressive: Set to 100 to permit convergence using approximate
                     factorizations or factorizations other than LU. If
                     the factorization uses a technique other than
                     Gaussian elimination, the guarantees in
                     err_bnds_norm and err_bnds_comp may no longer be
                     trustworthy. \n
 \n
        PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
            will attempt to find a solution with small componentwise
            relative error in the double-precision algorithm.  Positive
            is true, 0.0 is false. \n
         Default: 1.0 (attempt componentwise convergence) \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  Successful exit. The solution to every right-hand side is
           guaranteed. \n
          < 0:  If INFO = -i, the i-th argument had an illegal value \n
          > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
           has been completed, but the factor U is exactly singular, so
           the solution and error bounds could not be computed. RCOND = 0
           is returned. \n
          = N+J: The solution corresponding to the Jth right-hand side is
           not guaranteed. The solutions corresponding to other right-
           hand sides K with K > J may not be guaranteed as well, but
           only the first such right-hand side is reported. If a small
           componentwise error is not requested (PARAMS(3) = 0.0) then
           the Jth right-hand side is the first with a normwise error
           bound that is not guaranteed (the smallest J such
           that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
           the Jth right-hand side is the first with either a normwise or
           componentwise error bound that is not guaranteed (the smallest
           J such that either ERR_BNDS_NORM(J,1) = 0.0 or
           ERR_BNDS_COMP(J,1) = 0.0). See the definition of
           ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
           about all of the right-hand sides check ERR_BNDS_NORM or
           ERR_BNDS_COMP. \n

 *  * */
    template <typename T>
    void gbsvxx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
                integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, T *r, T *c, T *b,
                integer *ldb, T *x, integer *ldx, T *rcond, T *rpvgrw, T *berr, integer *n_err_bnds,
                T *err_bnds_norm, T *err_bnds_comp, integer *nparams, T *params, T *work,
                integer *iwork, integer *info)
    {
        gbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x,
               ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
               work, iwork, info);
    }
    template <typename T, typename Ta>
    void gbsvxx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
                integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, Ta *r, Ta *c,
                T *b, integer *ldb, T *x, integer *ldx, Ta *rcond, Ta *rpvgrw, Ta *berr,
                integer *n_err_bnds, Ta *err_bnds_norm, Ta *err_bnds_comp, integer *nparams,
                Ta *params, T *work, Ta *rwork, integer *info)
    {
        gbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x,
               ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
               work, rwork, info);
    }
    /** @}*/ // end of gbsvxx

    /** @defgroup gtsv gtsv
     * @{
     * @ingroup LU LU
     */
    /*! @brief GTSV computes the solution to system of linear equations A * X = B for GT matrices

 * @details
 * \b Purpose:
    \verbatim
     GTSV  solves the equation
        A*X = B,
     where A is an n by n tridiagonal matrix, by Gaussian elimination with
     partial pivoting.

     Note that the equation  A**T*X = B  may be solved by interchanging the
     order of the arguments DU and DL.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] DL
          DL is REAL array, dimension (N-1) \n
          On entry, DL must contain the (n-1) sub-diagonal elements of
          A. \n
 \n
          On exit, DL is overwritten by the (n-2) elements of the
          second super-diagonal of the upper triangular matrix U from
          the LU factorization of A, in DL(1), ..., DL(n-2). \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, D must contain the diagonal elements of A. \n
 \n
          On exit, D is overwritten by the n diagonal elements of U. \n
 * @param[in,out] DU
          DU is REAL array, dimension (N-1) \n
          On entry, DU must contain the (n-1) super-diagonal elements
          of A. \n
 \n
          On exit, DU is overwritten by the (n-1) elements of the first
          super-diagonal of U. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N by NRHS matrix of right hand side matrix B. \n
          On exit, if INFO = 0, the N by NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
               has not been computed.  The factorization has not been
               completed unless i = N. \n

 *  * */
    template <typename T>
    void gtsv(integer *n, integer *nrhs, T *dl, T *d, T *du, T *b, integer *ldb, integer *info)
    {
        gtsv(n, nrhs, dl, d, du, b, ldb, info);
    }
    /** @}*/ // end of gtsv

    /** @defgroup gtsvx gtsvx
     * @{
     * @ingroup LU LU
     */
    /*! @brief GTSVX uses the LU factorization to compute the solution to a real system of linear
 equations

 * @details
 * \b Purpose:
    \verbatim
     GTSVX uses the LU factorization to compute the solution to a real
     system of linear equations A * X = B or A**T * X = B,
     where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     matrices.

     Error bounds on the solution and a condition estimate are also
     provided.

 * \b Description:
     The following steps are performed:

     1. If FACT = 'N', the LU decomposition is used to factor the matrix A
        as A = L * U, where L is a product of permutation and unit lower
        bidiagonal matrices and U is upper triangular with nonzeros in
        only the main diagonal and first two superdiagonals.

     2. If some U(i,i)=0, so that U is exactly singular, then the routine
        returns with INFO = i. Otherwise, the factored form of A is used
        to estimate the condition number of the matrix A.  If the
        reciprocal of the condition number is less than machine precision,
        INFO = N+1 is returned as a warning, but the routine still goes on
        to solve for X and compute error bounds as described below.

     3. The system of equations is solved for X using the factored form
        of A.

     4. Iterative refinement is applied to improve the computed solution
        matrix and calculate error bounds and backward error estimates
        for it.
    \endverbatim

 * @param[in] FACT
          FACT is CHARACTER*1 \n
          Specifies whether or not the factored form of A has been
          supplied on entry. \n
          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored
                  form of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV
                  will not be modified. \n
          = 'N':  The matrix will be copied to DLF, DF, and DUF
                  and factored. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] DL
          DL is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of A. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of A. \n
 * @param[in] DU
          DU is REAL array, dimension (N-1) \n
          The (n-1) superdiagonal elements of A. \n
 * @param[in,out] DLF
          DLF is REAL array, dimension (N-1) \n
          If FACT = 'F', then DLF is an input argument and on entry
          contains the (n-1) multipliers that define the matrix L from
          the LU factorization of A as computed by SGTTRF. \n
 \n
          If FACT = 'N', then DLF is an output argument and on exit
          contains the (n-1) multipliers that define the matrix L from
          the LU factorization of A. \n
 * @param[in,out] DF
          DF is REAL array, dimension (N) \n
          If FACT = 'F', then DF is an input argument and on entry
          contains the n diagonal elements of the upper triangular
          matrix U from the LU factorization of A. \n
 \n
          If FACT = 'N', then DF is an output argument and on exit
          contains the n diagonal elements of the upper triangular
          matrix U from the LU factorization of A. \n
 * @param[in,out] DUF
          DUF is REAL array, dimension (N-1) \n
          If FACT = 'F', then DUF is an input argument and on entry
          contains the (n-1) elements of the first superdiagonal of U. \n
 \n
          If FACT = 'N', then DUF is an output argument and on exit
          contains the (n-1) elements of the first superdiagonal of U. \n
 * @param[in,out] DU2
          DU2 is REAL array, dimension (N-2) \n
          If FACT = 'F', then DU2 is an input argument and on entry
          contains the (n-2) elements of the second superdiagonal of
          U. \n
 \n
          If FACT = 'N', then DU2 is an output argument and on exit
          contains the (n-2) elements of the second superdiagonal of
          U. \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          If FACT = 'F', then IPIV is an input argument and on entry
          contains the pivot indices from the LU factorization of A as
          computed by SGTTRF. \n
 \n
          If FACT = 'N', then IPIV is an output argument and on exit
          contains the pivot indices from the LU factorization of A;
          row i of the matrix was interchanged with row IPIV(i).
          IPIV(i) will always be either i or i+1; IPIV(i) = i indicates
          a row interchange was not required. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The N-by-NRHS right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
          X is REAL array, dimension (LDX,NRHS) \n
          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          The estimate of the reciprocal condition number of the matrix
          A.  If RCOND is less than the machine precision (in
          particular, if RCOND = 0), the matrix is singular to working
          precision.  This condition is indicated by a return code of
          INFO > 0. \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X). \n
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, and i is \n
                <= N:  U(i,i) is exactly zero.  The factorization
                       has not been completed unless i = N, but the
                       factor U is exactly singular, so the solution
                       and error bounds could not be computed.
                       RCOND = 0 is returned. \n
                = N+1: U is nonsingular, but RCOND is less than machine
                       precision, meaning that the matrix is singular
                       to working precision.  Nevertheless, the
                       solution and error bounds are computed because
                       there are a number of situations where the
                       computed solution can be more accurate than the
                       value of RCOND would suggest. \n

 *  * */
    template <typename T>
    void gtsvx(char *fact, char *trans, integer *n, integer *nrhs, T *dl, T *d, T *du, T *dlf,
               T *df, T *duf, T *du2, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx,
               T *rcond, T *ferr, T *berr, T *work, integer *iwork, integer *info)
    {
        gtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr,
              berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void gtsvx(char *fact, char *trans, integer *n, integer *nrhs, T *dl, T *d, T *du, T *dlf,
               T *df, T *duf, T *du2, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx,
               Ta *rcond, Ta *ferr, Ta *berr, T *work, Ta *rwork, integer *info)
    {
        gtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr,
              berr, work, rwork, info);
    }
    /** @}*/ // end of gtsvx
    /** @}*/ // end of LU

    /** @defgroup LU_Computational LU Computational: computational routines (factor, cond, etc.)
     * @ingroup LinearSolve
     * @{
     */
    /** @defgroup gecon gecon
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GECON estimates the reciprocal of the condition number of a general real matrix A

 * @details
 * \b Purpose:
    \verbatim
    GECON estimates the reciprocal of the condition number of a general
    real matrix A, in either the 1-norm or the infinity-norm, using
    the LU factorization computed by SGETRF.

    An estimate is obtained for norm(inv(A)), and the reciprocal of the
    condition number is computed as
       RCOND = 1 / ( norm(A) * norm(inv(A))).
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required: \n
          = '1' or 'O':  1-norm; \n
          = 'I':         Infinity-norm. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The factors L and U from the factorization A = P*L*U
          as computed by SGETRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] ANORM
          ANORM is REAL \n
          If NORM = '1' or 'O', the 1-norm of the original matrix A. \n
          If NORM = 'I', the infinity-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(norm(A) * norm(inv(A))). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (4*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gecon(char *norm, integer *n, T *a, integer *lda, T *anorm, T *rcond, T *work,
               integer *iwork, integer *info)
    {
        gecon(norm, n, a, lda, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void gecon(char *norm, integer *n, T *a, integer *lda, Ta *anorm, Ta *rcond, T *work, Ta *rwork,
               integer *info)
    {
        gecon(norm, n, a, lda, anorm, rcond, work, rwork, info);
    }
    /** @}*/ // end of gecon

    /** @defgroup getrf getrf
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LU factorization of a general m-by-n matrix a
    *   using partial pivoting with row interchanges.
    *
    *  @details
    * \b Purpose:
    * \verbatim
        LU factorization of a general m-by-n matrix a using partial pivoting with row interchanges.
        The factorization has the form
            A = P * L * U
        where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower
        trapezoidal if M >  , and U is upper triangular (upper trapezoidal if M < N).

        This is the right-looking Level 3 BLAS version of the algorithm.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is REAL/DOUBLE PRECISION/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the m-by-n matrix to be factored. \n
              On exit, the factors L and U from the factorization
              A = P*L*U; the unit diagonal elements of L are not stored. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,m) \n
    * @param[out] ipiv
              ipiv is integer array, dimension (min(m,n)) \n
              The pivot indices; for 1 <= i <= min(m,n), row i of the
              matrix was interchanged with row ipiv(i). \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                    has been completed, but the factor U is exactly
                    singular, and division by zero will occur if it is used
                    to solve a system of equations. \n

    *     *  */
    template <typename T>
    void getrf(integer *m, integer *n, T *a, integer *lda, integer *ipiv, integer *info)
    {
        getrf(m, n, a, lda, ipiv, info);
    }

    /** @}*/ // end of getrf

    /** @defgroup getrf2 getrf2
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GETRF2 computes an LU factorization of a general M-by-N matrix A
 * @details
 * \b Purpose:
    \verbatim
    GETRF2 computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.
    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).
    This is the recursive version of the algorithm. It divides
    the matrix into four submatrices:
           [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
       A = [ -----|----- ]  with n1 = min(m,n)/2
           [  A21 | A22  ]       n2 = n-n1
                                          [ A11 ]
    The subroutine calls itself to factor [ --- ],
                                          [ A12 ]
                    [ A12 ]
    do the swaps on [ --- ], solve A12, update A22,
                    [ A22 ]
    then calls itself to factor A22 and do the swaps on A21.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix to be factored. \n
          On exit, the factors L and U from the factorization
          A = P*L*U; the unit diagonal elements of L are not stored. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (min(M,N)) \n
          The pivot indices; for 1 <= i <= min(M,N), row i of the
          matrix was interchanged with row IPIV(i). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                has been completed, but the factor U is exactly
                singular, and division by zero will occur if it is used
                to solve a system of equations. \n

 *  * */
    template <typename T>
    void getrf2(integer *m, integer *n, T *a, integer *lda, integer *ipiv, integer *info)
    {
        getrf2(m, n, a, lda, ipiv, info);
    }
    /** @}*/ // end of getrf2

    /** @defgroup getf2 getf2
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LU factorization of a general m-by-n matrix a
    *   using partial pivoting with row interchanges.
    *
    * @details
    * \b Purpose:
    * \verbatim
        LU factorization of a general m-by-n matrix a using partial pivoting with row interchanges.
        The factorization has the form
            A = P * L * U
        where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower
        trapezoidal if M > N), and U is upper triangular (upper trapezoidal if M < N).

        This is the right-looking Level 2 BLAS version of the algorithm.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is REAL/DOUBLE PRECISION/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the m-by-n matrix to be factored. \n
              On exit, the factors L and U from the factorization
              A = P*L*U; the unit diagonal elements of L are not stored. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,m) \n
    * @param[out] ipiv
              ipiv is integer array, dimension (min(m,n)) \n
              ipiv is integer array, dimension (min(m,n)) \n
              The pivot indices; for 1 <= i <= min(m,n), row i of the
              matrix was interchanged with row ipiv(i). \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n
              > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                    has been completed, but the factor U is exactly
                    singular, and division by zero will occur if it is used
                    to solve a system of equations. \n

    *     *  */
    template <typename T>
    void getf2(integer *m, integer *n, T *a, integer *lda, integer *ipiv, integer *info)
    {
        getf2(m, n, a, lda, ipiv, info);
    }
    /** @}*/ // end of getf2

    /** @defgroup getrs getrs
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GETRS solves a system of linear equations
 * @details
 * \b Purpose:
    \verbatim
     GETRS solves a system of linear equations
        A * X = B  or  A**T * X = B
     with a general N-by-N matrix A using the LU factorization computed
     by GETRF.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T* X = B  (Transpose) \n
          = 'C':  A**T* X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The factors L and U from the factorization A = P*L*U
          as computed by GETRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from GETRF; for 1<=i<=N, row i of the
          matrix was interchanged with row IPIV(i). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void getrs(char *trans, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
               integer *ldb, integer *info)
    {
        getrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
    }

    /** @}*/ // end of getrs

    /** @defgroup getri getri
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GETRI computes the inverse of a matrix using the LU factorization computed by GETRF
 * @details
 * \b Purpose:
    \verbatim
     GETRI computes the inverse of a matrix using the LU factorization
     computed by GETRF.

     This method inverts U and then computes inv(A) by solving the system
     inv(A)*L = inv(U) for inv(A).
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the factors L and U from the factorization
          A = P*L*U as computed by GETRF. \n
          On exit, if INFO = 0, the inverse of the original matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from GETRF; for 1<=i<=N, row i of the
          matrix was interchanged with row IPIV(i). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO=0, then WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.  LWORK >= fla_max(1,N).
          For optimal performance LWORK >= N*NB, where NB is
          the optimal blocksize returned by ILAENV. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                singular and its inverse could not be computed. \n

 *  * */
    template <typename T>
    void getri(integer *n, T *a, integer *lda, integer *ipiv, T *work, integer *lwork,
               integer *info)
    {
        getri(n, a, lda, ipiv, work, lwork, info);
    }
    /** @}*/ // end of getri

    /** @defgroup gerfs gerfs
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GERFS improves the computed solution to a system of linear equations \n
  and provides error bounds and backward error estimates for the solution
 * @details
 * \b Purpose:
    \verbatim
    GERFS improves the computed solution to a system of linear
    equations and provides error bounds and backward error estimates for
    the solution.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The original N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The factors L and U from the factorization A = P*L*U
          as computed by GETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from GETRF; for 1<=i<=N, row i of the
          matrix was interchanged with row IPIV(i). \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by GETRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gerfs(char *trans, integer *n, integer *nrhs, T *a, integer *lda, T *af, integer *ldaf,
               integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr, T *work,
               integer *iwork, integer *info)
    {
        gerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork,
              info);
    }
    template <typename T, typename Ta>
    void gerfs(char *trans, integer *n, integer *nrhs, T *a, integer *lda, T *af, integer *ldaf,
               integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work,
               Ta *rwork, integer *info)
    {
        gerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork,
              info);
    }
    /** @}*/ // end of gerfs

    /** @defgroup gerfsx gerfsx
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GERFSX improves the computed solution to a system of linear equations \n
    and provides error bounds and backward error estimates for the solution
 *  @details
 *  \b Purpose:
    \verbatim
     GERFSX improves the computed solution to a system of linear
     equations and provides error bounds and backward error estimates
     for the solution.  In addition to normwise error bound, the code
     provides maximum componentwise error bound if possible.  See
     comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the
     error bounds.

     The original system of linear equations may have been equilibrated
     before calling this routine, as described by arguments EQUED, R
     and C below. In this case, the solution and error bounds returned
     are for the original unequilibrated system.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
            = 'N':  A * X = B     (No transpose) \n
            = 'T':  A**T * X = B  (Transpose) \n
            = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done to A
          before calling this routine. This is needed to compute
          the solution and error bounds correctly. \n
           = 'N':  No equilibration \n
           = 'R':  Row equilibration, i.e., A has been premultiplied by
                   diag(R). \n
           = 'C':  Column equilibration, i.e., A has been postmultiplied
                   by diag(C). \n
           = 'B':  Both row and column equilibration, i.e., A has been
                   replaced by diag(R) * A * diag(C).
                   The right hand side B has been changed accordingly. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The original N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The factors L and U from the factorization A = P*L*U
          as computed by GETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from GETRF; for 1<=i<=N, row i of the
          matrix was interchanged with row IPIV(i). \n
 * @param[in] R
          R is REAL array, dimension (N) \n
          The row scale factors for A.  If EQUED = 'R' or 'B', A is
          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
          is not accessed. \n
          If R is accessed, each element of R should be a power of the radix
          to ensure a reliable solution and error estimates. Scaling by
          powers of the radix does not cause rounding errors unless the
          result underflows or overflows. Rounding errors during scaling
          lead to refining with a matrix that is not equivalent to the
          input matrix, producing error estimates that may not be
          reliable. \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The column scale factors for A.  If EQUED = 'C' or 'B', A is
          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
          is not accessed. \n
          If C is accessed, each element of C should be a power of the radix
          to ensure a reliable solution and error estimates. Scaling by
          powers of the radix does not cause rounding errors unless the
          result underflows or overflows. Rounding errors during scaling
          lead to refining with a matrix that is not equivalent to the
          input matrix, producing error estimates that may not be
          reliable. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by GETRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          Reciprocal scaled condition number. This is an estimate of the
          reciprocal Skeel condition number of the matrix A after
          equilibration (if done).  If this is less than the machine
          precision (in particular, if it is zero), the matrix is singular
          to working precision.  Note that the error may still be small even
          if this number is very small and the matrix appears ill-
          conditioned. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          Componentwise relative backward error. This is the
          componentwise relative backward error of each solution vector X(j)
          (i.e., the smallest relative change in any element of A or B that
          makes X(j) an exact solution). \n
 * @param[in] N_ERR_BNDS
          N_ERR_BNDS is INTEGER \n
          Number of error bounds to return for each right hand side
          and each type (normwise or componentwise).  See ERR_BNDS_NORM and
          ERR_BNDS_COMP below. \n
 * @param[out] ERR_BNDS_NORM
          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS)
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          normwise relative error, which is defined as follows: \n
 \n
          Normwise relative error in the ith solution vector: \n

          \f[ \frac{{max\_j}\ abs(XTRUE(j,i) - X(j,i))}{{max\_j}\ abs(X(j,i))}\f] \n

 \n
          The array is indexed by the type of error information as described
          below. There currently are up to three pieces of information
          returned. \n
 \n
          The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_NORM(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated normwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*A, where S scales each row by a power of the
                  radix so all absolute row sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[out] ERR_BNDS_COMP
          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          componentwise relative error, which is defined as follows: \n
 \n
          Componentwise relative error in the ith solution vector: \n

          \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
          The array is indexed by the right-hand side i (on which the
          componentwise relative error depends), and the type of error
          information as described below. There currently are up to three
          pieces of information returned for each right-hand side. If
          componentwise accuracy is not requested (PARAMS(3) = 0.0), then
          ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
          the first (:,N_ERR_BNDS) entries are returned. \n
 \n
          The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_COMP(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated componentwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*(A*diag(x)), where x is the solution for the
                  current right-hand side and S scales each row of
                  A*diag(x) by a power of the radix so all absolute row
                  sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in] NPARAMS
          NPARAMS is INTEGER \n
          Specifies the number of parameters set in PARAMS.  If <= 0, the
          PARAMS array is never referenced and default values are used. \n
 * @param[in,out] PARAMS
          PARAMS is REAL array, dimension NPARAMS \n
          Specifies algorithm parameters.  If an entry is < 0.0, then
          that entry will be filled with default value used for that
          parameter.  Only positions up to NPARAMS are accessed; defaults
          are used for higher-numbered parameters. \n
 \n
          PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
              refinement or not. \n
           Default: 1.0 \n
              = 0.0:  No refinement is performed, and no error bounds are
                      computed. \n
              = 1.0:  Use the double-precision refinement algorithm,
                      possibly with doubled-single computations if the
                      compilation environment does not support DOUBLE
                      PRECISION.
                (other values are reserved for future use) \n
 \n
          PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
              computations allowed for refinement. \n
           Default: 10 \n
           Aggressive: Set to 100 to permit convergence using approximate
                       factorizations or factorizations other than LU. If
                       the factorization uses a technique other than
                       Gaussian elimination, the guarantees in
                       err_bnds_norm and err_bnds_comp may no longer be
                       trustworthy. \n
 \n
          PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
              will attempt to find a solution with small componentwise
              relative error in the double-precision algorithm.  Positive
              is true, 0.0 is false. \n
           Default: 1.0 (attempt componentwise convergence) \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  Successful exit. The solution to every right-hand side is
           guaranteed. \n
          < 0:  If INFO = -i, the i-th argument had an illegal value \n
          > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
           has been completed, but the factor U is exactly singular, so
           the solution and error bounds could not be computed. RCOND = 0
           is returned. \n
          = N+J: The solution corresponding to the Jth right-hand side is
           not guaranteed. The solutions corresponding to other right-
           hand sides K with K > J may not be guaranteed as well, but
           only the first such right-hand side is reported. If a small
           componentwise error is not requested (PARAMS(3) = 0.0) then
           the Jth right-hand side is the first with a normwise error
           bound that is not guaranteed (the smallest J such
           that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
           the Jth right-hand side is the first with either a normwise or
           componentwise error bound that is not guaranteed (the smallest
           J such that either ERR_BNDS_NORM(J,1) = 0.0 or
           ERR_BNDS_COMP(J,1) = 0.0). See the definition of
           ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
           about all of the right-hand sides check ERR_BNDS_NORM or
           ERR_BNDS_COMP. \n

*  * */
    template <typename T>
    void gerfsx(char *trans, char *equed, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, T *r, T *c, T *b, integer *ldb, T *x, integer *ldx,
                T *rcond, T *berr, integer *n_err_bnds, T *err_bnds_norm, T *err_bnds_comp,
                integer *nparams, T *params, T *work, integer *iwork, integer *info)
    {
        gerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr,
               n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
    }
    template <typename T, typename Ta>
    void gerfsx(char *trans, char *equed, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, Ta *r, Ta *c, T *b, integer *ldb, T *x, integer *ldx,
                Ta *rcond, Ta *berr, integer *n_err_bnds, Ta *err_bnds_norm, Ta *err_bnds_comp,
                integer *nparams, Ta *params, T *work, Ta *rwork, integer *info)
    {
        gerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr,
               n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    }
    /** @}*/ // end of gerfsx

    /** @defgroup gerqf gerqf
     * @ingroup RQ_factorization RQ_factorization
     * @{
     */
    /*! @brief GERQF computes a RQ factorization of a M-by-N matrix
     * @details
     * \b Purpose:
        \verbatim
         GERQF computes an RQ factorization of a real M-by-N matrix A:
         A = R * Q.
        \endverbatim 

     * @param[in] M
              M is INTEGER \n
              The number of rows of the matrix A.  M >= 0. \n
     * @param[in] N
              N is INTEGER \n
              The number of columns of the matrix A.  N >= 0. \n
     * @param[in,out] A
              A is REAL array, dimension (LDA,N) \n
              On entry, the M-by-N matrix A. \n
              On exit,
              if m <= n, the upper triangle of the subarray
              A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R; \n
              if m >= n, the elements on and above the (m-n)-th subdiagonal
              contain the M-by-N upper trapezoidal matrix R;
              the remaining elements, with the array TAU, represent the
              orthogonal matrix Q as a product of min(m,n) elementary
              reflectors (see Further Details). \n
     * @param[in] LDA
              LDA is INTEGER \n
              The leading dimension of the array A.  LDA >= fla_max(1,M). \n
     * @param[out] TAU
              TAU is REAL array, dimension (min(M,N)) \n
              The scalar factors of the elementary reflectors (see Further
              Details). \n
     * @param[out]	WORK	
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
     * @param[in]	LWORK	
              LWORK is INTEGER \n
              The dimension of the array WORK.  LWORK >= fla_max(1,M).
              For optimum performance LWORK >= M*NB, where NB is
              the optimal blocksize. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
     * @param[out]	INFO	
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

     *  * */
    template< typename T >
    void gerqf(integer* m, integer* n, T* a, integer* lda, T* tau, T* work, integer* lwork, integer* info)
    {
        gerqf(m, n, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of gerqf

    /** @defgroup geequ geequ
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GEEQU computes row and column scalings intended to equilibrate an \n
     M-by-N matrix A and reduce its condition number

 * @details
 * \b Purpose:
    \verbatim
     GEEQU computes row and column scalings intended to equilibrate an
     M-by-N matrix A and reduce its condition number.  R returns the row
     scale factors and C the column scale factors, chosen to try to make
     the largest element in each row and column of the matrix B with
     elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.

     R(i) and C(j) are restricted to be between SMLNUM = smallest safe
     number and BIGNUM = largest safe number.  Use of these scaling
     factors is not guaranteed to reduce the condition number of A but
     works well in practice.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The M-by-N matrix whose equilibration factors are
          to be computed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] R
          R is REAL array, dimension (M) \n
          If INFO = 0 or INFO > M, R contains the row scale factors
          for A. \n
 * @param[out] C
          C is REAL array, dimension (N) \n
          If INFO = 0,  C contains the column scale factors for A. \n
 * @param[out] ROWCND
          ROWCND is REAL \n
          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
          AMAX is neither too large nor too small, it is not worth
          scaling by R. \n
 * @param[out] COLCND
          COLCND is REAL \n
          If INFO = 0, COLCND contains the ratio of the smallest
          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
          worth scaling by C. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i,  and i is \n
                <= M:  the i-th row of A is exactly zero \n
                >  M:  the (i-M)-th column of A is exactly zero \n

 *  * */
    template <typename T>
    void geequ(integer *m, integer *n, T *a, integer *lda, T *r, T *c, T *rowcnd, T *colcnd,
               T *amax, integer *info)
    {
        geequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    }
    template <typename T, typename Ta>
    void geequ(integer *m, integer *n, T *a, integer *lda, Ta *r, Ta *c, Ta *rowcnd, Ta *colcnd,
               Ta *amax, integer *info)
    {
        geequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    }
    /** @}*/ // end of geequ

    /** @defgroup geequb geequb
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief  GEEQUB - computes row and column scaling to reduce condition number of matrix
 *
 *  @details
 *  \b Purpose:
 *  \verbatim
    GEEQUB computes row and column scalings intended to equilibrate an
    M-by-N matrix A and reduce its condition number.  R returns the row
    scale factors and C the column scale factors, chosen to try to make
    the largest element in each row and column of the matrix B with
    elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most
    the radix.

    R(i) and C(j) are restricted to be a power of the radix between
    SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use
    of these scaling factors is not guaranteed to reduce the condition
    number of A but works well in practice.

    This routine differs from SGEEQU by restricting the scaling factors
    to a power of the radix.  Barring over- and underflow, scaling by
    these factors introduces no additional rounding errors.  However, the
    scaled entries' magnitudes are no longer approximately 1 but lie
    between sqrt(radix) and 1/sqrt(radix).
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The M-by-N matrix whose equilibration factors are
          to be computed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] R
          R is REAL array, dimension (M) \n
          If INFO = 0 or INFO > M, R contains the row scale factors
          for A. \n
 * @param[out] C
          C is REAL array, dimension (N) \n
          If INFO = 0,  C contains the column scale factors for A. \n
 * @param[out] ROWCND
          ROWCND is REAL \n
          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
          AMAX is neither too large nor too small, it is not worth
          scaling by R. \n
 * @param[out] COLCND
          COLCND is REAL \n
          If INFO = 0, COLCND contains the ratio of the smallest
          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
          worth scaling by C. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i,  and i is \n
                <= M:  the i-th row of A is exactly zero \n
                >  M:  the (i-M)-th column of A is exactly zero \n

 *  * */
    template <typename T>
    void geequb(integer *m, integer *n, T *a, integer *lda, T *r, T *c, T *rowcnd, T *colcnd,
                T *amax, integer *info)
    {
        geequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    }
    template <typename T, typename Ta>
    void geequb(integer *m, integer *n, T *a, integer *lda, Ta *r, Ta *c, Ta *rowcnd, Ta *colcnd,
                Ta *amax, integer *info)
    {
        geequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info);
    }
    /** @}*/ // end of geequb

    /** @defgroup laqge laqge
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LAQGE scales a general rectangular matrix, using row and column scaling factors
 computed by geequ

 * @details
 * \b Purpose:
    \verbatim

    LAQGE equilibrates a general M by N matrix A using the row and
    column scaling factors in the vectors R and C.

    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M by N matrix A.
          On exit, the equilibrated matrix.  See EQUED for the form of
          the equilibrated matrix. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(M,1). \n
 * @param[in] R
          R is REAL array, dimension (M) \n
          The row scale factors for A. \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The column scale factors for A. \n
 * @param[in] ROWCND
          ROWCND is REAL \n
          Ratio of the smallest R(i) to the largest R(i). \n
 * @param[in] COLCND
          COLCND is REAL \n
          Ratio of the smallest C(i) to the largest C(i). \n
 * @param[in] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix entry. \n
 * @param[out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done. \n
          = 'N':  No equilibration \n
          = 'R':  Row equilibration, i.e., A has been premultiplied by
                  diag(R). \n
          = 'C':  Column equilibration, i.e., A has been postmultiplied
                  by diag(C). \n
          = 'B':  Both row and column equilibration, i.e., A has been
                  replaced by diag(R) * A * diag(C).  \n

 * */
    template <typename T>
    void laqge(integer *m, integer *n, T *a, integer *lda, T *r, T *c, T *rowcnd, T *colcnd,
               T *amax, char *equed)
    {
        laqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
    }
    template <typename T, typename Ta>
    void laqge(integer *m, integer *n, T *a, integer *lda, Ta *r, Ta *c, Ta *rowcnd, Ta *colcnd,
               Ta *amax, char *equed)
    {
        laqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed);
    }
    /** @}*/ // end of laqge

    /** @defgroup laswp laswp
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief Swap rows

  * @details
  * \b Purpose:
    \verbatim
     Perform a series of row interchanges on the matrix a.
     One row interchange is initiated for each of rows k1 through k2 of A.
    \endverbatim

    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the matrix of column dimension n to which the row
              interchanges will be applied. \n
              On exit, the permuted matrix. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
    * @param[in] k1
              k1 is integer* \n
              The first element of ipiv for which a row interchange will
              be done. \n
    * @param[in] k2
              k2 is integer* \n
              (k2-k1+1) is the number of elements of ipiv for which a row
              interchange will be done. \n
    * @param[in] ipiv
              ipiv is integer array, dimension (k1+(k2-k1)*abs(incx)) \n
              The vector of pivot indices. Only the elements in positions
              k1 through k1+(k2-k1)*abs(incx) of ipiv are accessed. \n
              ipiv(k1+(K-k1)*abs(incx)) = L implies rows K and L are to be
              interchanged. \n
    * @param[in] incx
              incx is integer* \n
              The increment between successive values of ipiv. If incx
              is negative, the pivots are applied in reverse order. \n

    *     *  */
    template <typename T>
    void laswp(integer *n, T *a, integer *lda, integer *k1, integer *k2, integer *ipiv,
               integer *incx)
    {
        laswp(n, a, lda, k1, k2, ipiv, incx);
    }
    /** @}*/ // end of laswp

    /** @defgroup gesc2 gesc2
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GESC2 solves a system of linear equations using the LU factorization with complete
 pivoting computed by getc2

 * @details
 * \b Purpose:
    \verbatim
    GESC2 solves a system of linear equations

              A * X = scale* RHS

    with a general N-by-N matrix A using the LU factorization with
    complete pivoting computed by GETC2.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the  LU part of the factorization of the n-by-n
          matrix A computed by SGETC2:  A = P * L * U * Q \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1, N). \n
 * @param[in,out] RHS
          RHS is REAL array, dimension (N). \n
          On entry, the right hand side vector b.
          On exit, the solution vector X. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N). \n
          The pivot indices; for 1 <= i <= N, row i of the
          matrix has been interchanged with row IPIV(i). \n
 * @param[in] JPIV
          JPIV is INTEGER array, dimension (N). \n
          The pivot indices; for 1 <= j <= N, column j of the
          matrix has been interchanged with column JPIV(j). \n
 * @param[out] SCALE
          SCALE is REAL \n
          On exit, SCALE contains the scale factor. SCALE is chosen
          0 <= SCALE <= 1 to prevent overflow in the solution. \n

 *  * */
    template <typename T>
    void gesc2(integer *n, T *a, integer *lda, T *rhs, integer *ipiv, integer *jpiv, T *scale)
    {
        gesc2(n, a, lda, rhs, ipiv, jpiv, scale);
    }
    template <typename T, typename Ta>
    void gesc2(integer *n, T *a, integer *lda, T *rhs, integer *ipiv, integer *jpiv, Ta *scale)
    {
        gesc2(n, a, lda, rhs, ipiv, jpiv, scale);
    }
    /** @}*/ // end of gesc2

    /** @defgroup latdf latdf
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 \n
     and computes a contribution to the reciprocal Dif-estimate
 * @details
 * \b Purpose:
    \verbatim
    LATDF uses the LU factorization of the n-by-n matrix Z computed by
    SGETC2 and computes a contribution to the reciprocal Dif-estimate
    by solving Z * x = b for x, and choosing the r.h.s. b such that
    the norm of x is as large as possible. On entry RHS = b holds the
    contribution from earlier solved sub-systems, and on   return RHS = x.

    The factorization of Z   returned by SGETC2 has the form Z = P*L*U*Q,
    where P and Q are permutation matrices. L is lower triangular with
    unit diagonal elements and U is upper triangular.
    \endverbatim

 * @param[in] IJOB
          IJOB is INTEGER \n
          IJOB = 2: First compute an approximative null-vector e
              of Z using SGECON, e is normalized and solve for
              Zx = +-e - f with the sign giving the greater value
              of 2-norm(x). About 5 times as expensive as Default. \n
          IJOB .ne. 2: Local look ahead strategy where all entries of
              the r.h.s. b is chosen as either +1 or -1 (Default). \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix Z. \n
 * @param[in] Z
          Z is REAL array, dimension (LDZ, N) \n
          On entry, the LU part of the factorization of the n-by-n
          matrix Z computed by SGETC2:  Z = P * L * U * Q \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z.  LDA >= fla_max(1, N). \n
 * @param[in,out] RHS
          RHS is REAL array, dimension N. \n
          On entry, RHS contains contributions from other subsystems.
          On exit, RHS contains the solution of the subsystem with
          entries according to the value of IJOB (see above). \n
 * @param[in,out] RDSUM
          RDSUM is REAL \n
          On entry, the sum of squares of computed contributions to
          the Dif-estimate under computation by STGSYL, where the
          scaling factor RDSCAL (see below) has been factored out. \n
          On exit, the corresponding sum of squares updated with the
          contributions from the current sub-system.
          If TRANS = 'T' RDSUM is not touched.
          NOTE: RDSUM only makes sense when STGSY2 is called by STGSYL. \n
 * @param[in,out] RDSCAL
          RDSCAL is REAL \n
          On entry, scaling factor used to prevent overflow in RDSUM. \n
          On exit, RDSCAL is updated w.r.t. the current contributions
          in RDSUM. \n
          If TRANS = 'T', RDSCAL is not touched. \n
          NOTE: RDSCAL only makes sense when STGSY2 is called by
                STGSYL. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N). \n
          The pivot indices; for 1 <= i <= N, row i of the
          matrix has been interchanged with row IPIV(i). \n
 * @param[in] JPIV
          JPIV is INTEGER array, dimension (N). \n
          The pivot indices; for 1 <= j <= N, column j of the
          matrix has been interchanged with column JPIV(j).  \n

 *  * */
    template <typename T>
    void latdf(integer *ijob, integer *n, T *z, integer *ldz, T *rhs, T *rdsum, T *rdscal,
               integer *ipiv, integer *jpiv)
    {
        latdf(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
    }
    template <typename T, typename Ta>
    void latdf(integer *ijob, integer *n, T *z, integer *ldz, T *rhs, Ta *rdsum, Ta *rdscal,
               integer *ipiv, integer *jpiv)
    {
        latdf(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv);
    }
    /** @}*/ // end of latdf

    /** @defgroup la_gercond la_gercond
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LA_GERCOND estimates the Skeel condition number for a general matrix

 * @details
 * \b Purpose:
    \verbatim
    LA_GERCOND estimates the Skeel condition number of op(A) * op2(C)
    where op2 is determined by CMODE as follows
    CMODE =  1    op2(C) = C
    CMODE =  0    op2(C) = I
    CMODE = -1    op2(C) = inv(C)
    The Skeel condition number cond(A) = norminf(|inv(A)||A|)
    is computed by computing scaling factors R such that
    diag(R)*A*op2(C) is row equilibrated and computing the standard
    infinity-norm condition number.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
            = 'N':  A * X = B     (No transpose) \n
            = 'T':  A**T * X = B  (Transpose) \n
            = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The factors L and U from the factorization
          A = P*L*U as computed by SGETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from the factorization A = P*L*U
          as computed by SGETRF; row i of the matrix was interchanged
          with row IPIV(i). \n
 * @param[in] CMODE
          CMODE is INTEGER \n
          Determines op2(C) in the formula op(A) * op2(C) as follows: \n
          CMODE =  1    op2(C) = C \n
          CMODE =  0    op2(C) = I \n
          CMODE = -1    op2(C) = inv(C) \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The vector C in the formula op(A) * op2(C). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          i > 0:  The ith argument is invalid. \n
 * @param[out] WORK
          WORK is REAL array, dimension (3*N). \n
          Workspace. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (N). \n
          Workspace.2 \n

 *  * */
    template <typename T>
    T la_gercond(char *trans, integer *n, T *a, integer *lda, T *af, integer *ldaf, integer *ipiv,
                 integer *cmode, T *c, integer *info, T *work, integer *iwork)
    {
        return la_gercond(*trans, *n, a, *lda, af, *ldaf, *ipiv, *cmode, c, *info, work, *iwork);
    }
    template <typename T, typename Ta>
    Ta la_gercond_x(char *trans, integer *n, T *a, integer *lda, T *af, integer *ldaf,
                    integer *ipiv, T *x, integer *info, T *work, Ta *rwork)
    {
        return cla_gercond_x(*trans, *n, a, *lda, af, *ldaf, *ipiv, x, *info, work, *rwork);
    }
    template <typename T, typename Ta>
    Ta la_gercond_c(char *trans, integer *n, T *a, integer *lda, T *af, integer *ldaf,
                    integer *ipiv, Ta *c, logical *capply, integer *info, T *work, Ta *rwork)
    {
        return cla_gercond_c(*trans, *n, a, *lda, af, *ldaf, *ipiv, c, *capply, *info, work,
                             *rwork);
    }
    /** @}*/ // end of la_gercond

    /** @defgroup la_gerpvgrw la_gerpvgrw
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LA_GERPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U)

 * @details
 * \b Purpose:
    \verbatim
    LA_GERPVGRW computes the reciprocal pivot growth factor
    norm(A)/norm(U). The "max absolute element" norm is used. If this is
    much less than 1, the stability of the LU factorization of the
    (equilibrated) matrix A could be poor. This also means that the
    solution X, estimated condition numbers, and error bounds could be
    unreliable.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NCOLS
          NCOLS is INTEGER \n
          The number of columns of the matrix A. NCOLS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The factors L and U from the factorization
          A = P*L*U as computed by SGETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n

 *  * */
    template <typename T>
    T la_gerpvgrw(integer *n, integer *ncols, T *a, integer *lda, T *af, integer *ldaf)
    {
        return la_gerpvgrw(n, ncols, a, lda, af, ldaf);
    }

    /** @}*/ // end of la_gerpvgrw

    /** @defgroup la_gerfsx_extended la_gerfsx_extended
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LA_GERFSX_EXTENDED improves the computed solution to a system of      \n
     linear equations for general matrices by performing extra-precise iterative \n
     refinement and provides error bounds and backward error estimates           \n
     for the solution
 * @details
 * \b Purpose:
    \verbatim
    LA_GERFSX_EXTENDED improves the computed solution to a system of
    linear equations by performing extra-precise iterative refinement
    and provides error bounds and backward error estimates for the solution.
    This subroutine is called by SGERFSX to perform iterative refinement.
    In addition to normwise error bound, the code provides maximum
    componentwise error bound if possible. See comments for ERRS_N
    and ERRS_C for details of the error bounds. Note that this
    subroutine is only resonsible for setting the second fields of
    ERRS_N and ERRS_C.
    \endverbatim

 * @param[in] PREC_TYPE
          PREC_TYPE is INTEGER \n
          Specifies the intermediate precision to be used in refinement.
          The value is defined by ILAPREC(P) where P is a CHARACTER and P \n
              = 'S':  Single \n
              = 'D':  Double \n
              = 'I':  Indigenous \n
              = 'X' or 'E':  Extra \n
 * @param[in] TRANS_TYPE
          TRANS_TYPE is INTEGER \n
          Specifies the transposition operation on A. \n
          The value is defined by ILATRANS(T) where T is a CHARACTER and T \n
              = 'N':  No transpose \n
              = 'T':  Transpose \n
              = 'C':  Conjugate transpose \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right-hand-sides, i.e., the number of columns of the
          matrix B. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The factors L and U from the factorization
          A = P*L*U as computed by SGETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from the factorization A = P*L*U
          as computed by SGETRF; row i of the matrix was interchanged
          with row IPIV(i). \n
 * @param[in] COLEQU
          COLEQU is LOGICAL \n
          If .TRUE. then column equilibration was done to A before calling
          this routine. This is needed to compute the solution and error
          bounds correctly. \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The column scale factors for A. If COLEQU = .FALSE., C
          is not accessed. If C is input, each element of C should be a power
          of the radix to ensure a reliable solution and error estimates.
          Scaling by powers of the radix does not cause rounding errors unless
          the result underflows or overflows. Rounding errors during scaling
          lead to refining with a matrix that is not equivalent to the
          input matrix, producing error estimates that may not be
          reliable. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right-hand-side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] Y
          Y is REAL array, dimension (LDY,NRHS) \n
          On entry, the solution matrix X, as computed by SGETRS.
          On exit, the improved solution matrix Y. \n
 * @param[in] LDY
          LDY is INTEGER \n
          The leading dimension of the array Y.  LDY >= fla_max(1,N). \n
 * @param[out] BERR_OUT
          BERR_OUT is REAL array, dimension (NRHS) \n
          On exit, BERR_OUT(j) contains the componentwise relative backward
          error for right-hand-side j from the formula \n
             fla_max(i) (abs(RES(i)) / (abs(op(A_s))*abs(Y) + abs(B_s))(i)) \n
          where abs(Z) is the componentwise absolute value of the matrix
          or vector Z. This is computed by SLA_LIN_BERR. \n
 * @param[in] N_NORMS
          N_NORMS is INTEGER \n
          Determines which error bounds to   return (see ERRS_N
          and ERRS_C). \n
          If N_NORMS >= 1   return normwise error bounds. \n
          If N_NORMS >= 2   return componentwise error bounds. \n
 * @param[in,out] ERRS_N
          ERRS_N is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          normwise relative error, which is defined as follows: \n
 \n
          Normwise relative error in the ith solution vector: \n

          \f[ \frac{{max\_j}\  abs(XTRUE(j,i) - X(j,i))}{{max\_j}\ abs(X(j,i))} \f] \n

 \n
          The array is indexed by the type of error information as described
          below. There currently are up to three pieces of information
           returned. \n
 \n
          The first index in ERRS_N(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERRS_N(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon').
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true.
 \n
          err = 3  Reciprocal condition number: Estimated normwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*A, where S scales each row by a power of the
                  radix so all absolute row sums of Z are approximately 1.
 \n
          This subroutine is only responsible for setting the second field
          above. \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in,out] ERRS_C
          ERRS_C is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          componentwise relative error, which is defined as follows: \n
 \n
          Componentwise relative error in the ith solution vector: \n

          \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
          The array is indexed by the right-hand side i (on which the
          componentwise relative error depends), and the type of error
          information as described below. There currently are up to three
          pieces of information   returned for each right-hand side. If
          componentwise accuracy is not requested (PARAMS(3) = 0.0), then
          ERRS_C is not accessed.  If N_ERR_BNDS < 3, then at most
          the first (:,N_ERR_BNDS) entries are returned. \n
 \n
          The first index in ERRS_C(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERRS_C(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon').
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true.
 \n
          err = 3  Reciprocal condition number: Estimated componentwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*(A*diag(x)), where x is the solution for the
                  current right-hand side and S scales each row of
                  A*diag(x) by a power of the radix so all absolute row
                  sums of Z are approximately 1.
 \n
          This subroutine is only responsible for setting the second field
          above. \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in] RES
          RES is REAL array, dimension (N) \n
          Workspace to hold the intermediate residual. \n
 * @param[in] AYB
          AYB is REAL array, dimension (N) \n
          Workspace. This can be the same workspace passed for Y_TAIL. \n
 * @param[in] DY
          DY is REAL array, dimension (N) \n
          Workspace to hold the intermediate solution. \n
 * @param[in] Y_TAIL
          Y_TAIL is REAL array, dimension (N) \n
          Workspace to hold the trailing bits of the intermediate solution. \n
 * @param[in] RCOND
          RCOND is REAL \n
          Reciprocal scaled condition number.  This is an estimate of the
          reciprocal Skeel condition number of the matrix A after
          equilibration (if done).  If this is less than the machine
          precision (in particular, if it is zero), the matrix is singular
          to working precision.  Note that the error may still be small even
          if this number is very small and the matrix appears ill-
          conditioned. \n
 * @param[in] ITHRESH
          ITHRESH is INTEGER \n
          The maximum number of residual computations allowed for
          refinement. The default is 10. For 'aggressive' set to 100 to
          permit convergence using approximate factorizations or
          factorizations other than LU. If the factorization uses a
          technique other than Gaussian elimination, the guarantees in
          ERRS_N and ERRS_C may no longer be trustworthy. \n
 * @param[in] RTHRESH
          RTHRESH is REAL \n
          Determines when to stop refinement if the error estimate stops
          decreasing. Refinement will stop when the next solution no longer
          satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is
          the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The
          default value is 0.5. For 'aggressive' set to 0.9 to permit
          convergence on extremely ill-conditioned matrices. See LAWN 165
          for more details. \n
 * @param[in] DZ_UB
          DZ_UB is REAL \n
          Determines when to start considering componentwise convergence.
          Componentwise convergence is only considered after each component
          of the solution Y is stable, which we definte as the relative
          change in each component being less than DZ_UB. The default value
          is 0.25, requiring the first bit to be stable. See LAWN 165 for
          more details. \n
 * @param[in] IGNORE_CWISE
          IGNORE_CWISE is LOGICAL \n
          If .TRUE. then ignore componentwise convergence. Default value
          is .FALSE.. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          < 0:  if INFO = -i, the ith argument to SGETRS had an illegal
                value \n

 * */
    template <typename T>
    void la_gerfsx_extended(integer *prec_type, integer *trans_type, integer *n, integer *nrhs,
                            T *a, integer *lda, T *af, integer *ldaf, integer *ipiv,
                            logical *colequ, T *c, T *b, integer *ldb, T *y, integer *ldy,
                            T *berr_out, integer *n_norms, T *errs_n, T *errs_c, T *res, T *ayb,
                            T *dy, T *y_tail, T *rcond, integer *ithresh, T *rthresh, T *dz_ub,
                            logical *ignore_cwise, integer *info)
    {
        la_gerfsx_extended(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                           ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail,
                           rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
    }
    template <typename T, typename Ta>
    void la_gerfsx_extended(integer *prec_type, integer *trans_type, integer *n, integer *nrhs,
                            T *a, integer *lda, T *af, integer *ldaf, integer *ipiv,
                            logical *colequ, Ta *c, T *b, integer *ldb, T *y, integer *ldy,
                            Ta *berr_out, integer *n_norms, Ta *errs_n, Ta *errs_c, T *res, Ta *ayb,
                            T *dy, T *y_tail, Ta *rcond, integer *ithresh, Ta *rthresh, Ta *dz_ub,
                            logical *ignore_cwise, integer *info)
    {
        la_gerfsx_extended(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                           ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy, y_tail,
                           rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
    }
    /** @}*/ // end of la_gerfsx_extended

    /** @defgroup gbcon gbcon
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief estimates the reciprocal of the condition number of a real general band matrix A, in
     either the 1-norm or the infinity-norm,
     *
     * @details
     * \b Purpose :
       \verbatim
        GBCON estimates the reciprocal of the condition number of a real
        general band matrix A, in either the 1-norm or the infinity-norm,
        using the LU factorization computed by SGBTRF.

        An estimate is obtained for norm(inv(A)), and the reciprocal of the
        condition number is computed as
          RCOND = 1 / ( norm(A) * norm(inv(A))).
     \endverbatim

     * @param[in] NORM
              NORM is CHARACTER*1 \n
              Specifies whether the 1-norm condition number or the
              infinity-norm condition number is required: \n
              = '1' or 'O':  1-norm; \n
              = 'I':         Infinity-norm. \n
     * @param[in] N
              N is INTEGER \n
              The order of the matrix A.  N >= 0. \n
     * @param[in] KL
              KL is INTEGER \n
              The number of subdiagonals within the band of A.  KL >= 0. \n
     * @param[in] KU
              KU is INTEGER \n
              The number of superdiagonals within the band of A.  KU >= 0. \n
     * @param[in] AB
              AB is REAL array, dimension (LDAB,N) \n
              Details of the LU factorization of the band matrix A, as
              computed by SGBTRF.  U is stored as an upper triangular band
              matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
              the multipliers used during the factorization are stored in
              rows KL+KU+2 to 2*KL+KU+1. \n
     * @param[in] LDAB
              LDAB is INTEGER \n
              The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. \n
     * @param[in] IPIV
              IPIV is INTEGER array, dimension (N) \n
              The pivot indices; for 1 <= i <= N, row i of the matrix was
              interchanged with row IPIV(i). \n
     * @param[in] ANORM
              ANORM is REAL \n
              If NORM = '1' or 'O', the 1-norm of the original matrix A. \n
              If NORM = 'I', the infinity-norm of the original matrix A. \n
     * @param[out] RCOND
              RCOND is REAL \n
              The reciprocal of the condition number of the matrix A,
              computed as RCOND = 1/(norm(A) * norm(inv(A))). \n
     * @param[out]	WORK
              WORK is COMPLEX array, dimension (2*N) \n
     * @param[out]	RWORK
              RWORK is REAL array, dimension (N) \n
     * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0: if INFO = -i, the i-th argument had an illegal value \n

     * * */
    template <typename T>
    void gbcon(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab,
               integer *ipiv, T *anorm, T *rcond, T *work, integer *iwork, integer *info)
    {
        gbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void gbcon(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab,
               integer *ipiv, Ta *anorm, Ta *rcond, T *work, Ta *rwork, integer *info)
    {
        gbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, rwork, info);
    }
    /** @}*/ // end of gbcon

    /** @defgroup gbtrf gbtrf
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief \b GBTRF computes the LU factorization of a general band matrix using the pivot
 *  @details
 *  \b Purpose:
 *  \verbatim
    GBTRF computes an LU factorization of a real m-by-n band matrix A
    using partial pivoting with row interchanges.

    This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0.
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the matrix A in band storage, in rows KL+1 to
          2*KL+KU+1; rows 1 to KL of the array need not be set.
          The j-th column of A is stored in the j-th column of the
          array AB as follows: \n
          AB(kl+ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(m,j+kl) \n
 \n
          On exit, details of the factorization: U is stored as an
          upper triangular band matrix with KL+KU superdiagonals in
          rows 1 to KL+KU+1, and the multipliers used during the
          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
          See below for further details. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. \n
 * @param[in] IPIV
           IPIV is INTEGER array, dimension (min(M,N)) \n
           The pivot indices; for 1 <= i <= min(M,N), row i of the
           matrix was interchanged with row IPIV(i). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
               has been completed, but the factor U is exactly
               singular, and division by zero will occur if it is used
               to solve a system of equations. \n

 *   * */
    template <typename T>
    void gbtrf(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab,
               integer *ipiv, integer *info)
    {
        gbtrf(m, n, kl, ku, ab, ldab, ipiv, info);
    }

    /** @}*/ // end of gbtrf

    /** @defgroup gbtf2 gbtf2
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GBTF2 computes the LU factorization of a general band matrix using the unblocked
 version of the algorithm

 * @details
 * \b Purpose:
    \verbatim
    GBTF2 computes an LU factorization of a real m-by-n band matrix A
    using partial pivoting with row interchanges.

    This is the unblocked version of the algorithm, calling Level 2 BLAS.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the matrix A in band storage, in rows KL+1 to
          2*KL+KU+1; rows 1 to KL of the array need not be set.
          The j-th column of A is stored in the j-th column of the
          array AB as follows: \n
          AB(kl+ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(m,j+kl) \n
 \n
          On exit, details of the factorization: U is stored as an
          upper triangular band matrix with KL+KU superdiagonals in
          rows 1 to KL+KU+1, and the multipliers used during the
          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
          See below for further details. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (min(M,N)) \n
          The pivot indices; for 1 <= i <= min(M,N), row i of the
          matrix was interchanged with row IPIV(i). \n
 * @param[out] INFO
            INFO is INTEGER \n
            = 0: successful exit \n
            < 0: if INFO = -i, the i-th argument had an illegal value \n
            > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
                 has been completed, but the factor U is exactly
                 singular, and division by zero will occur if it is used
                 to solve a system of equations. \n

 *  * */
    template <typename T>
    void gbtf2(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab,
               integer *ipiv, integer *info)
    {
        gbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
    }
    /** @}*/ // end of gbtf2

    /** @defgroup gbtrs gbtrs
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief \b GBTRS solves a system of linear equationsA * X = B  or  A**T * X = B with a
general band matrix A using the LU factorization computed by GBTRF
 *  @details
    \b Purpose:
    \verbatim
    GBTRS solves a system of linear equations
    A * X = B  or  A**T * X = B
    with a general band matrix A using the LU factorization computed
    by GBTRF.
    \endverbatim
 * @param[in] TRANS
        TRANS is CHARACTER*1 \n
        Specifies the form of the system of equations. \n
        = 'N':  A * X = B  (No transpose) \n
        = 'T':  A**T* X = B  (Transpose) \n
        = 'C':  A**T* X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
        N is INTEGER \n
        The order of the matrix A.  N >= 0. \n
 * @param[in] KL
        KL is INTEGER \n
        The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
        KU is INTEGER \n
        The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] NRHS
        NRHS is INTEGER \n
        The number of right hand sides, i.e., the number of columns
        of the matrix B.  NRHS >= 0. \n
 * @param[in] AB
        AB is REAL array, dimension (LDAB,N) \n
        Details of the LU factorization of the band matrix A, as
        computed by SGBTRF.  U is stored as an upper triangular band
        matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
        the multipliers used during the factorization are stored in
        rows KL+KU+2 to 2*KL+KU+1. \n
 * @param[in] LDAB
        LDAB is INTEGER \n
        The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. \n
 * @param[in] IPIV
        IPIV is INTEGER \n array, dimension (N)
        The pivot indices; for 1 <= i <= N, row i of the matrix was
        interchanged with row IPIV(i). \n
 * @param[in,out] B
        B is REAL array, dimension (LDB,NRHS) \n
        On entry, the right hand side matrix B.
        On exit, the solution matrix X. \n
 * @param[in] LDB
        LDB is INTEGER \n
        The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
        INFO is INTEGER \n
        = 0:  successful exit \n
        < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void gbtrs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
               integer *ldab, integer *ipiv, T *b, integer *ldb, integer *info)
    {
        gbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    }
    /** @}*/ // end of gbtrs

    /** @defgroup gbrfs gbrfs
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief Improves computed solution to banded matrix and provides error bounds and backward
estimates
*
* @details

* \b Purpose:
  \verbatim
    GBRFS improves the computed solution to a system of linear
    equations when the coefficient matrix is banded, and provides
    error bounds and backward error estimates for the solution.
  \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns \n
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The original band matrix A, stored in rows 1 to KL+KU+1. \n
          The j-th column of A is stored in the j-th column of the
          array AB as follows: \n
          AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(n,j+kl). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[in] AFB
          AFB is REAL array, dimension (LDAFB,N) \n
          Details of the LU factorization of the band matrix A, as
          computed by SGBTRF.  U is stored as an upper triangular band
          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
          the multipliers used during the factorization are stored in
          rows KL+KU+2 to 2*KL+KU+1. \n
 * @param[in] LDAFB
          LDAFB is INTEGER \n
          The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from SGBTRF; for 1<=i<=N, row i of the
          matrix was interchanged with row IPIV(i). \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by SGBTRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X). \n
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gbrfs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
               integer *ldab, T *afb, integer *ldafb, integer *ipiv, T *b, integer *ldb, T *x,
               integer *ldx, T *ferr, T *berr, T *work, integer *iwork, integer *info)
    {
        gbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work,
              iwork, info);
    }
    template <typename T, typename Ta>
    void gbrfs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
               integer *ldab, T *afb, integer *ldafb, integer *ipiv, T *b, integer *ldb, T *x,
               integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork, integer *info)
    {
        gbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work,
              rwork, info);
    }
    /** @}*/ // end of gbrfs

    /** @defgroup gbrfsx gbrfsx
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief Improves computed solution to banded matrix and provides backward estimates,component
wise and normwise error bound
 *
 *  @details
 *  \b Purpose:
 *  \verbatim
    Improves computed solution to banded matrix and provides backward estimates,component wise and
normwise error bounds \endverbatim

 * @param[in] TRANS
        TRANS is CHARACTER*1 \n
        Specifies the form of the system of equations: \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] EQUED
        EQUED is CHARACTER*1 \n
        Specifies the form of equilibration that was done to A
        before calling this routine. This is needed to compute
        the solution and error bounds correctly. \n
           = 'N':  No equilibration \n
           = 'R':  Row equilibration, i.e., A has been premultiplied by
                   diag(R). \n
           = 'C':  Column equilibration, i.e., A has been postmultiplied
                   by diag(C). \n
           = 'B':  Both row and column equilibration, i.e., A has been
                   replaced by diag(R) * A * diag(C).
                   The right hand side B has been changed accordingly. \n
 * @param[in] N
        N is INTEGER \n
        The order of the matrix A.  N >= 0. \n
 * @param[in] KL
        KL is INTEGER  \n
        The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
        KU is INTEGER \n
        The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] NRHS
        NRHS is INTEGER \n
        The number of right hand sides, i.e., the number of columns
        of the matrices B and X.  NRHS >= 0. \n
 * @param[in] AB
        AB is REAL array, dimension (LDAB,N) \n
        The original band matrix A, stored in rows 1 to KL+KU+1. \n
        The j-th column of A is stored in the j-th column of the
        array AB as follows: \n
        AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(n,j+kl). \n
 * @param[in] LDAB
        LDAB is INTEGER \n
        The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[in] AFB
        AFB is REAL array, dimension (LDAFB,N) \n
        Details of the LU factorization of the band matrix A, as
        computed by DGBTRF.  U is stored as an upper triangular band
        matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
        the multipliers used during the factorization are stored in
        rows KL+KU+2 to 2*KL+KU+1. \n
 * @param[in] LDAFB
        LDAFB is INTEGER \n
        The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1. \n
 * @param[in] IPIV
        IPIV is INTEGER \n array, dimension (N) \n
        The pivot indices from SGETRF; for 1<=i<=N, row i of the
        matrix was interchanged with row IPIV(i). \n
 * @param[in,out] R
        R is REAL array, dimension (N) \n
        The row scale factors for A.  If EQUED = 'R' or 'B', A is
        multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
        is not accessed.  R is an input argument if FACT = 'F';
        otherwise, R is an output argument.  If FACT = 'F' and
        EQUED = 'R' or 'B', each element of R must be positive. \n
        If R is output, each element of R is a power of the radix.
        If R is input, each element of R should be a power of the radix
        to ensure a reliable solution and error estimates. Scaling by
        powers of the radix does not cause rounding errors unless the
        result underflows or overflows. Rounding errors during scaling
        lead to refining with a matrix that is not equivalent to the
        input matrix, producing error estimates that may not be
        reliable. \n
 * @param[in,out] C
        C is REAL array, dimension (N) \n
        The column scale factors for A.  If EQUED = 'C' or 'B', A is
        multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
        is not accessed.  C is an input argument if FACT = 'F';
        otherwise, C is an output argument.  If FACT = 'F' and
        EQUED = 'C' or 'B', each element of C must be positive. \n
        If C is output, each element of C is a power of the radix. \n
        If C is input, each element of C should be a power of the radix
        to ensure a reliable solution and error estimates. Scaling by
        powers of the radix does not cause rounding errors unless the
        result underflows or overflows. Rounding errors during scaling
        lead to refining with a matrix that is not equivalent to the
        input matrix, producing error estimates that may not be
        reliable. \n
 * @param[in] B
        B is REAL array, dimension (LDB,NRHS) \n
        The right hand side matrix B. \n
 * @param[in] LDB
        LDB is INTEGER \n
        The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
        X is REAL array, dimension (LDX,NRHS) \n
        On entry, the solution matrix X, as computed by SGETRS. \n
        On exit, the improved solution matrix X. \n
 * @param[in] LDX
        LDX is INTEGER \n
        The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
        RCOND is REAL \n
        Reciprocal scaled condition number. This is an estimate of the
        reciprocal Skeel condition number of the matrix A after
        equilibration (if done).  If this is less than the machine
        precision (in particular, if it is zero), the matrix is singular
        to working precision.  Note that the error may still be small even
        if this number is very small and the matrix appears ill-
        conditioned. \n
 * @param[out] BERR
        BERR is REAL array, dimension (NRHS) \n
        Componentwise relative backward error. This is the
        componentwise relative backward error of each solution vector X(j)
        (i.e., the smallest relative change in any element of A or B that
        makes X(j) an exact solution). \n
 * @param[in] N_ERR_BNDS
        N_ERR_BNDS is INTEGER \n
        Number of error bounds to return for each right hand side
        and each type (normwise or componentwise).  See ERR_BNDS_NORM and
        ERR_BNDS_COMP below. \n
 * @param[out] ERR_BNDS_NORM
        ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) \n
        For each right-hand side, this array contains information about
        various error bounds and condition numbers corresponding to the
        normwise relative error, which is defined as follows: \n
 \n
        Normwise relative error in the ith solution vector: \n

        \f[ \frac{{max\_j}\  abs(XTRUE(j,i) - X(j,i))}{{max\_j}\ abs(X(j,i))} \f] \n

 \n
        The array is indexed by the type of error information as described
        below. There currently are up to three pieces of information
        returned. \n
 \n
        The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
        right-hand side. \n
 \n
        The second index in ERR_BNDS_NORM(:,err) contains the following
        three fields: \n
        err = 1 "Trust/don't trust" boolean. Trust the answer if the
                reciprocal condition number is less than the threshold
                sqrt(n) * slamch('Epsilon'). \n
 \n
        err = 2 "Guaranteed" error bound: The estimated forward error,
                almost certainly within a factor of 10 of the true error
                so long as the next entry is greater than the threshold
                sqrt(n) * slamch('Epsilon'). This error bound should only
                be trusted if the previous boolean is true. \n
 \n
        err = 3  Reciprocal condition number: Estimated normwise
                reciprocal condition number.  Compared with the threshold
                sqrt(n) * slamch('Epsilon') to determine if the error
                estimate is "guaranteed". These reciprocal condition
                numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                appropriately scaled matrix Z.
                Let Z = S*A, where S scales each row by a power of the
                radix so all absolute row sums of Z are approximately 1. \n
 * @param[out] ERR_BNDS_COMP
        ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
        For each right-hand side, this array contains information about
        various error bounds and condition numbers corresponding to the
        componentwise relative error, which is defined as follows: \n
 \n
        Componentwise relative error in the ith solution vector: \n

        \f[{max\_j} \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
        The array is indexed by the right-hand side i (on which the
        componentwise relative error depends), and the type of error
        information as described below. There currently are up to three
        pieces of information returned for each right-hand side. If
        componentwise accuracy is not requested (PARAMS(3) = 0.0), then
        ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
        the first (:,N_ERR_BNDS) entries are returned. \n
 \n
        The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
        right-hand side. \n
 \n
        The second index in ERR_BNDS_COMP(:,err) contains the following
        three fields: \n
        err = 1 "Trust/don't trust" boolean. Trust the answer if the
                reciprocal condition number is less than the threshold
                sqrt(n) * slamch('Epsilon'). \n
 \n
        err = 2 "Guaranteed" error bound: The estimated forward error,
                almost certainly within a factor of 10 of the true error
                so long as the next entry is greater than the threshold
                sqrt(n) * slamch('Epsilon'). This error bound should only
                be trusted if the previous boolean is true. \n
 \n
        err = 3  Reciprocal condition number: Estimated componentwise
                reciprocal condition number.  Compared with the threshold
                sqrt(n) * slamch('Epsilon') to determine if the error
                estimate is "guaranteed". These reciprocal condition
                numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                appropriately scaled matrix Z.
                Let Z = S*(A*diag(x)), where x is the solution for the
                current right-hand side and S scales each row of
                A*diag(x) by a power of the radix so all absolute row
                sums of Z are approximately 1. \n
 \n
        See Lapack Working Note 165 for further details and extra
        cautions. \n
 * @param[in] NPARAMS
        NPARAMS is INTEGER \n
        Specifies the number of parameters set in PARAMS.  If <= 0, the
        PARAMS array is never referenced and default values are used. \n
 * @param[in,out] PARAMS
        PARAMS is REAL array, dimension NPARAMS \n
        Specifies algorithm parameters.  If an entry is < 0.0, then
        that entry will be filled with default value used for that
        parameter.  Only positions up to NPARAMS are accessed; defaults
        are used for higher-numbered parameters. \n
 \n
        PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
            refinement or not. \n
         Default: 1.0 \n
            = 0.0:  No refinement is performed, and no error bounds are
                    computed. \n
            = 1.0:  Use the double-precision refinement algorithm,
                    possibly with doubled-single computations if the
                    compilation environment does not support DOUBLE
                    PRECISION. \n
              (other values are reserved for future use) \n
 \n
        PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
            computations allowed for refinement. \n
         Default: 10 \n
         Aggressive: Set to 100 to permit convergence using approximate
                     factorizations or factorizations other than LU. If
                     the factorization uses a technique other than
                     Gaussian elimination, the guarantees in
                     err_bnds_norm and err_bnds_comp may no longer be
                     trustworthy. \n
 \n
        PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
            will attempt to find a solution with small componentwise
            relative error in the double-precision algorithm.  Positive
            is true, 0.0 is false. \n
         Default: 1.0 (attempt componentwise convergence) \n
 * @param[out]	WORK
        WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	RWORK
        RWORK is REAL array, dimension (2*N) \n
 * @param[out]	INFO
        INFO is INTEGER \n
        = 0:  Successful exit. The solution to every right-hand side is
         guaranteed. \n
        < 0:  If INFO = -i, the i-th argument had an illegal value \n
        > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
         has been completed, but the factor U is exactly singular, so
         the solution and error bounds could not be computed. RCOND = 0
         is returned. \n
        = N+J: The solution corresponding to the Jth right-hand side is
         not guaranteed. The solutions corresponding to other right-
         hand sides K with K > J may not be guaranteed as well, but
         only the first such right-hand side is reported. If a small
         componentwise error is not requested (PARAMS(3) = 0.0) then
         the Jth right-hand side is the first with a normwise error
         bound that is not guaranteed (the smallest J such
         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
         the Jth right-hand side is the first with either a normwise or
         componentwise error bound that is not guaranteed (the smallest
         J such that either ERR_BNDS_NORM(J,1) = 0.0 or
         ERR_BNDS_COMP(J,1) = 0.0). See the definition of
         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
         about all of the right-hand sides check ERR_BNDS_NORM or
         ERR_BNDS_COMP. \n
* */
    template <typename T>
    void gbrfsx(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs,
                T *ab, integer *ldab, T *afb, integer *ldafb, integer *ipiv, T *r, T *c, T *b,
                integer *ldb, T *x, integer *ldx, T *rcond, T *berr, integer *n_err_bnds,
                T *err_bnds_norm, T *err_bnds_comp, integer *nparams, T *params, T *work,
                integer *iwork, integer *info)
    {
        gbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx,
               rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork,
               info);
    }
    template <typename T, typename Ta>
    void gbrfsx(char *trans, char *equed, integer *n, integer *kl, integer *ku, integer *nrhs,
                T *ab, integer *ldab, T *afb, integer *ldafb, integer *ipiv, Ta *r, Ta *c, T *b,
                integer *ldb, T *x, integer *ldx, Ta *rcond, Ta *berr, integer *n_err_bnds,
                Ta *err_bnds_norm, Ta *err_bnds_comp, integer *nparams, Ta *params, T *work,
                Ta *rwork, integer *info)
    {
        gbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x, ldx,
               rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork,
               info);
    }
    /** @}*/ // end of gbrfsx

    /** @defgroup gbequ gbequ
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief Computes row and column scalings intended to equilibrate  band matrix and reduce its
 condition number
 *
 * @details
 * \b Purpose:
 * \verbatim
    GBEQU computes row and column scalings intended to equilibrate an
    M-by-N band matrix A and reduce its condition number.
   \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th
          column of A is stored in the j-th column of the array AB as
          follows: \n
          AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(m,j+kl). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[out] R
          R is REAL array, dimension (M) \n
          If INFO = 0, or INFO > M, R contains the row scale factors
          for A. \n
 * @param[out] C
          C is REAL array, dimension (N) \n
          If INFO = 0, C contains the column scale factors for A. \n
 * @param[out] ROWCND
          ROWCND is REAL \n
          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
          AMAX is neither too large nor too small, it is not worth
          scaling by R. \n
 * @param[out] COLCND
          COLCND is REAL \n
          If INFO = 0, COLCND contains the ratio of the smallest
          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
          worth scaling by C. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, and i is \n
                <= M:  the i-th row of A is exactly zero \n
                >  M:  the (i-M)-th column of A is exactly zero \n

 *  * */
    template <typename T>
    void gbequ(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *r, T *c,
               T *rowcnd, T *colcnd, T *amax, integer *info)
    {
        gbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    }
    template <typename T, typename Ta>
    void gbequ(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, Ta *r, Ta *c,
               Ta *rowcnd, Ta *colcnd, Ta *amax, integer *info)
    {
        gbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    }
    /** @}*/ // end of gbequ

    /** @defgroup gbequb gbequb
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief Computes row and column scalings intended to equilibrate band matrix and reduce its
 condition number(scaling factor:power of radix)
 *
 * @details
 * \b Purpose:
 * \verbatim
    GBEQUB computes row and column scalings intended to equilibrate an
    M-by-N matrix A and reduce its condition number.  R returns the row
    scale factors and C the column scale factors, chosen to try to make
    the largest element in each row and column of the matrix B with
    elements B(i,j)=R(i)*A(i,j)*C(j) have an absolute value of at most
    the radix.

    R(i) and C(j) are restricted to be a power of the radix between
    SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use
    of these scaling factors is not guaranteed to reduce the condition
    number of A but works well in practice.

    This routine differs from SGEEQU by restricting the scaling factors
    to a power of the radix.  Barring over- and underflow, scaling by
    these factors introduces no additional rounding errors.  However, the
    scaled entries' magnitudes are no longer approximately 1 but lie
    between sqrt(radix) and 1/sqrt(radix).
   \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th
          column of A is stored in the j-th column of the array AB as
          follows: \n
          AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(m,j+kl). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[out] R
          R is REAL array, dimension (M) \n
          If INFO = 0, or INFO > M, R contains the row scale factors
          for A. \n
 * @param[out] C
          C is REAL array, dimension (N) \n
          If INFO = 0, C contains the column scale factors for A. \n
 * @param[out] ROWCND
          ROWCND is REAL \n
          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
          AMAX is neither too large nor too small, it is not worth
          scaling by R. \n
 * @param[out] COLCND
          COLCND is REAL \n
          If INFO = 0, COLCND contains the ratio of the smallest
          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
          worth scaling by C. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix element.  If AMAX is very
          close to overflow or very close to underflow, the matrix
          should be scaled. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i,  and i is \n
                <= M:  the i-th row of A is exactly zero \n
                >  M:  the (i-M)-th column of A is exactly zero \n

 *  * */
    template <typename T>
    void gbequb(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *r, T *c,
                T *rowcnd, T *colcnd, T *amax, integer *info)
    {
        gbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    }
    template <typename T, typename Ta>
    void gbequb(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, Ta *r,
                Ta *c, Ta *rowcnd, Ta *colcnd, Ta *amax, integer *info)
    {
        gbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info);
    }
    /** @}*/ // end of gbequb

    /** @defgroup laqgb laqgb
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LAQGB scales a general band matrix, using row and column scaling factors computed by
 gbequ

 * @details
 * \b Purpose:
    \verbatim
    LAQGB equilibrates a general M by N band matrix A with KL
    subdiagonals and KU superdiagonals using the row and scaling factors
    in the vectors R and C.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
          The j-th column of A is stored in the j-th column of the
          array AB as follows:
          AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(m,j+kl)
 \n
          On exit, the equilibrated matrix, in the same storage format
          as A.  See EQUED for the form of the equilibrated matrix. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDA >= KL+KU+1. \n
 * @param[in] R
          R is REAL array, dimension (M) \n
          The row scale factors for A. \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The column scale factors for A. \n
 * @param[in] ROWCND
          ROWCND is REAL \n
          Ratio of the smallest R(i) to the largest R(i). \n
 * @param[in] COLCND
          COLCND is REAL \n
          Ratio of the smallest C(i) to the largest C(i). \n
 * @param[in] AMAX
          AMAX is REAL \n
          Absolute value of largest matrix entry. \n
 * @param[out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done. \n
          = 'N':  No equilibration \n
          = 'R':  Row equilibration, i.e., A has been premultiplied by
                  diag(R). \n
          = 'C':  Column equilibration, i.e., A has been postmultiplied
                  by diag(C). \n
          = 'B':  Both row and column equilibration, i.e., A has been
                  replaced by diag(R) * A * diag(C).  \n

 * */
    template <typename T>
    void laqgb(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *r, T *c,
               T *rowcnd, T *colcnd, T *amax, char *equed)
    {
        laqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
    }
    template <typename T, typename Ta>
    void laqgb(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, Ta *r, Ta *c,
               Ta *rowcnd, Ta *colcnd, Ta *amax, char *equed)
    {
        laqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed);
    }
    /** @}*/ // end of laqgb

    /** @defgroup gbrcond gbrcond
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LA_GBRCOND estimates the Skeel condition number for a general banded matrix

 * @details
 * \b Purpose:
    \verbatim
    LA_GBRCOND Estimates the Skeel condition number of  op(A) * op2(C)
    where op2 is determined by CMODE as follows
    CMODE =  1    op2(C) = C
    CMODE =  0    op2(C) = I
    CMODE = -1    op2(C) = inv(C)
    The Skeel condition number  cond(A) = norminf(|inv(A)||A|)
    is computed by computing scaling factors R such that
    diag(R)*A*op2(C) is row equilibrated and computing the standard
    infinity-norm condition number.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
           = 'N':  A * X = B     (No transpose) \n
           = 'T':  A**T * X = B  (Transpose) \n
           = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
          The j-th column of A is stored in the j-th column of the
          array AB as follows: \n
          AB(KU+1+i-j,j) = A(i,j) for fla_max(1,j-KU)<=i<=min(N,j+kl) \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[in] AFB
          AFB is REAL array, dimension (LDAFB,N) \n
          Details of the LU factorization of the band matrix A, as
          computed by SGBTRF.  U is stored as an upper triangular
          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
          and the multipliers used during the factorization are stored
          in rows KL+KU+2 to 2*KL+KU+1. \n
 * @param[in] LDAFB
          LDAFB is INTEGER \n
          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from the factorization A = P*L*U \n
          as computed by SGBTRF; row i of the matrix was interchanged
          with row IPIV(i). \n
 * @param[in] CMODE
          CMODE is INTEGER \n
          Determines op2(C) in the formula op(A) * op2(C) as follows: \n
          CMODE =  1    op2(C) = C \n
          CMODE =  0    op2(C) = I \n
          CMODE = -1    op2(C) = inv(C) \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The vector C in the formula op(A) * op2(C). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          i > 0:  The ith argument is invalid. \n
 * @param[out] WORK
          WORK is REAL array, dimension (5*N). \n
          Workspace. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (N). \n
          Workspace. \n

 *  * */
    template <typename T>
    T la_gbrcond(char *trans, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *afb,
                 integer *ldafb, integer *ipiv, integer *cmode, T *c, integer *info, T *work,
                 integer *iwork)
    {
        return la_gbrcond(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work,
                          iwork);
    }
    /** @}*/ // end of gbrcond

    /** @defgroup la_gbrpvgrw la_gbrpvgrw
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LA_GBRPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U) for a general
 banded matrix.

 * @details
 * \b Purpose:
    \verbatim
    LA_GBRPVGRW computes the reciprocal pivot growth factor
    norm(A)/norm(U). The "max absolute element" norm is used. If this is
    much less than 1, the stability of the LU factorization of the
    (equilibrated) matrix A could be poor. This also means that the
    solution X, estimated condition numbers, and error bounds could be
    unreliable.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] NCOLS
          NCOLS is INTEGER \n
          The number of columns of the matrix A.  NCOLS >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
          The j-th column of A is stored in the j-th column of the
          array AB as follows: \n
          AB(KU+1+i-j,j) = A(i,j) for fla_max(1,j-KU)<=i<=min(N,j+kl) \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[in] AFB
          AFB is REAL array, dimension (LDAFB,N) \n
          Details of the LU factorization of the band matrix A, as
          computed by SGBTRF.  U is stored as an upper triangular
          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
          and the multipliers used during the factorization are stored
          in rows KL+KU+2 to 2*KL+KU+1. \n
 * @param[in] LDAFB
          LDAFB is INTEGER \n
          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. \n

 *  * */
    template <typename T>
    T la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, T *ab, integer *ldab,
                  T *afb, integer *ldafb)
    {
        return la_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb);
    }
    template <typename T, typename Ta>
    Ta la_gbrpvgrw(integer *n, integer *kl, integer *ku, integer *ncols, T *ab, integer *ldab,
                   T *afb, integer *ldafb)
    {
        return la_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb);
    }
    /** @}*/ // end of la_gbrpvgrw

    /** @defgroup la_gbrfsx_extended
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief LA_GBRFSX_EXTENDED improves the computed solution to a system of     \n
     linear equations for general banded matrices by performing extra-precise   \n
     iterative refinement and provides error bounds and backward error          \n
     estimates for the solution
 * @details
 * \b Purpose:
    \verbatim
    LA_GBRFSX_EXTENDED improves the computed solution to a system of
    linear equations by performing extra-precise iterative refinement
    and provides error bounds and backward error estimates for the solution.
    This subroutine is called by SGBRFSX to perform iterative refinement.
    In addition to normwise error bound, the code provides maximum
    componentwise error bound if possible. See comments for ERR_BNDS_NORM
    and ERR_BNDS_COMP for details of the error bounds. Note that this
    subroutine is only resonsible for setting the second fields of
    ERR_BNDS_NORM and ERR_BNDS_COMP.
    \endverbatim

 * @param[in] PREC_TYPE
          PREC_TYPE is INTEGER \n
          Specifies the intermediate precision to be used in refinement.
          The value is defined by ILAPREC(P) where P is a CHARACTER and P \n
               = 'S':  Single \n
               = 'D':  Double \n
               = 'I':  Indigenous \n
               = 'X' or 'E':  Extra \n
 * @param[in] TRANS_TYPE
          TRANS_TYPE is INTEGER \n
          Specifies the transposition operation on A. \n
          The value is defined by ILATRANS(T) where T is a CHARACTER and T \n
               = 'N':  No transpose \n
               = 'T':  Transpose \n
               = 'C':  Conjugate transpose \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0 \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right-hand-sides, i.e., the number of columns of the
          matrix B. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the N-by-N matrix AB. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= fla_max(1,N). \n
 * @param[in] AFB
          AFB is REAL array, dimension (LDAFB,N) \n
          The factors L and U from the factorization
          A = P*L*U as computed by SGBTRF. \n
 * @param[in] LDAFB
          LDAFB is INTEGER \n
          The leading dimension of the array AF.  LDAFB >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices from the factorization A = P*L*U
          as computed by SGBTRF; row i of the matrix was interchanged
          with row IPIV(i). \n
 * @param[in] COLEQU
          COLEQU is LOGICAL \n
          If .TRUE. then column equilibration was done to A before calling
          this routine. This is needed to compute the solution and error
          bounds correctly. \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The column scale factors for A. If COLEQU = .FALSE., C
          is not accessed. If C is input, each element of C should be a power
          of the radix to ensure a reliable solution and error estimates.
          Scaling by powers of the radix does not cause rounding errors unless
          the result underflows or overflows. Rounding errors during scaling
          lead to refining with a matrix that is not equivalent to the
          input matrix, producing error estimates that may not be
          reliable. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right-hand-side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] Y
          Y is REAL array, dimension (LDY,NRHS) \n
          On entry, the solution matrix X, as computed by SGBTRS.
          On exit, the improved solution matrix Y. \n
 * @param[in] LDY
          LDY is INTEGER \n
          The leading dimension of the array Y.  LDY >= fla_max(1,N). \n
 * @param[out] BERR_OUT
          BERR_OUT is REAL array, dimension (NRHS) \n
          On exit, BERR_OUT(j) contains the componentwise relative backward
          error for right-hand-side j from the formula \n
             fla_max(i) (abs(RES(i)) / (abs(op(A_s))*abs(Y) + abs(B_s))(i)) \n
          where abs(Z) is the componentwise absolute value of the matrix
          or vector Z. This is computed by SLA_LIN_BERR. \n
 * @param[in] N_NORMS
          N_NORMS is INTEGER \n
          Determines which error bounds to   return (see ERR_BNDS_NORM
          and ERR_BNDS_COMP). \n
          If N_NORMS >= 1   return normwise error bounds. \n
          If N_NORMS >= 2   return componentwise error bounds. \n
 * @param[in,out] ERR_BNDS_NORM
          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          normwise relative error, which is defined as follows: \n
 \n
          Normwise relative error in the ith solution vector: \n

          \f[ \frac{{max\_j}\  abs(XTRUE(j,i) - X(j,i))}{{max\_j}\ abs(X(j,i))} \f] \n

 \n
          The array is indexed by the type of error information as described
          below. There currently are up to three pieces of information
           returned. \n
 \n
          The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_NORM(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon').
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true.
 \n
          err = 3  Reciprocal condition number: Estimated normwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*A, where S scales each row by a power of the
                  radix so all absolute row sums of Z are approximately 1.
 \n
          This subroutine is only responsible for setting the second field
          above. \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in,out] ERR_BNDS_COMP
          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          componentwise relative error, which is defined as follows: \n
 \n
          Componentwise relative error in the ith solution vector: \n

          \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
          The array is indexed by the right-hand side i (on which the
          componentwise relative error depends), and the type of error
          information as described below. There currently are up to three
          pieces of information   returned for each right-hand side. If
          componentwise accuracy is not requested (PARAMS(3) = 0.0), then
          ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
          the first (:,N_ERR_BNDS) entries are returned. \n
 \n
          The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
          right-hand side.
 \n
          The second index in ERR_BNDS_COMP(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon').
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true.
 \n
          err = 3  Reciprocal condition number: Estimated componentwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*(A*diag(x)), where x is the solution for the
                  current right-hand side and S scales each row of
                  A*diag(x) by a power of the radix so all absolute row
                  sums of Z are approximately 1.
 \n
          This subroutine is only responsible for setting the second field
          above. \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in] RES
          RES is REAL array, dimension (N) \n
          Workspace to hold the intermediate residual. \n
 * @param[in] AYB
          AYB is REAL array, dimension (N) \n
          Workspace. This can be the same workspace passed for Y_TAIL. \n
 * @param[in] DY
          DY is REAL array, dimension (N) \n
          Workspace to hold the intermediate solution. \n
 * @param[in] Y_TAIL
          Y_TAIL is REAL array, dimension (N) \n
          Workspace to hold the trailing bits of the intermediate solution. \n
 * @param[in] RCOND
          RCOND is REAL \n
          Reciprocal scaled condition number. This is an estimate of the
          reciprocal Skeel condition number of the matrix A after
          equilibration (if done). If this is less than the machine
          precision (in particular, if it is zero), the matrix is singular
          to working precision. Note that the error may still be small even
          if this number is very small and the matrix appears ill-
          conditioned. \n
 * @param[in] ITHRESH
          ITHRESH is INTEGER \n
          The maximum number of residual computations allowed for
          refinement. The default is 10. For 'aggressive' set to 100 to
          permit convergence using approximate factorizations or
          factorizations other than LU. If the factorization uses a
          technique other than Gaussian elimination, the guarantees in
          ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy. \n
 * @param[in] RTHRESH
          RTHRESH is REAL \n
          Determines when to stop refinement if the error estimate stops
          decreasing. Refinement will stop when the next solution no longer
          satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is
          the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The
          default value is 0.5. For 'aggressive' set to 0.9 to permit
          convergence on extremely ill-conditioned matrices. See LAWN 165
          for more details. \n
 * @param[in] DZ_UB
          DZ_UB is REAL \n
          Determines when to start considering componentwise convergence.
          Componentwise convergence is only considered after each component
          of the solution Y is stable, which we definte as the relative
          change in each component being less than DZ_UB. The default value
          is 0.25, requiring the first bit to be stable. See LAWN 165 for
          more details. \n
 * @param[in] IGNORE_CWISE
          IGNORE_CWISE is LOGICAL \n
          If .TRUE. then ignore componentwise convergence. Default value
          is .FALSE.. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          < 0:  if INFO = -i, the ith argument to SGBTRS had an illegal
                value  \n

 * */
    template <typename T>
    void la_gbrfsx_extended(integer *prec_type, integer *trans_type, integer *n, integer *kl,
                            integer *ku, integer *nrhs, T *ab, integer *ldab, T *afb,
                            integer *ldafb, integer *ipiv, logical *colequ, T *c, T *b,
                            integer *ldb, T *y, integer *ldy, T *berr_out, integer *n_norms,
                            T *err_bnds_norm, T *err_bnds_comp, T *res, T *ayb, T *dy, T *y_tail,
                            T *rcond, integer *ithresh, T *rthresh, T *dz_ub, logical *ignore_cwise,
                            integer *info)
    {
        la_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv,
                           colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm,
                           err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                           ignore_cwise, info);
    }
    template <typename T, typename Ta>
    void la_gbrfsx_extended(integer *prec_type, integer *trans_type, integer *n, integer *kl,
                            integer *ku, integer *nrhs, T *ab, integer *ldab, T *afb,
                            integer *ldafb, integer *ipiv, logical *colequ, Ta *c, T *b,
                            integer *ldb, T *y, integer *ldy, Ta *berr_out, integer *n_norms,
                            Ta *err_bnds_norm, Ta *err_bnds_comp, T *res, Ta *ayb, T *dy, T *y_tail,
                            Ta *rcond, integer *ithresh, Ta *rthresh, Ta *dz_ub,
                            logical *ignore_cwise, integer *info)
    {
        la_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv,
                           colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm,
                           err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                           ignore_cwise, info);
    }
    /** @}*/ // end of la_gbrfsx_extended

    /** @defgroup gtcon gtcon
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GTCON estimates the reciprocal of the condition number of a real tridiagonal matrix A
 * @details
 * \b Purpose:
    \verbatim
     GTCON estimates the reciprocal of the condition number of a real
     tridiagonal matrix A using the LU factorization as computed by
     SGTTRF.

     An estimate is obtained for norm(inv(A)), and the reciprocal of the
     condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required: \n
          = '1' or 'O':  1-norm; \n
          = 'I':         Infinity-norm. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] DL
          DL is REAL array, dimension (N-1) \n
          The (n-1) multipliers that define the matrix L from the
          LU factorization of A as computed by SGTTRF. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the upper triangular matrix U from
          the LU factorization of A. \n
 * @param[in] DU
          DU is REAL array, dimension (N-1) \n
          The (n-1) elements of the first superdiagonal of U. \n
 * @param[in] DU2
          DU2 is REAL array, dimension (N-2) \n
          The (n-2) elements of the second superdiagonal of U. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices; for 1 <= i <= n, row i of the matrix was
          interchanged with row IPIV(i).  IPIV(i) will always be either
          i or i+1; IPIV(i) = i indicates a row interchange was not
          required. \n
 * @param[in] ANORM
          ANORM is REAL \n
          If NORM = '1' or 'O', the 1-norm of the original matrix A. \n
          If NORM = 'I', the infinity-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (2*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

*  * */
    template <typename T>
    void gtcon(char *norm, integer *n, T *dl, T *d, T *du, T *du2, integer *ipiv, T *anorm,
               T *rcond, T *work, integer *iwork, integer *info)
    {
        gtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void gtcon(char *norm, integer *n, T *dl, T *d, T *du, T *du2, integer *ipiv, Ta *anorm,
               Ta *rcond, T *work, integer *info)
    {
        gtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info);
    }
    /** @}*/ // end of gtcon

    /** @defgroup gttrf gttrf
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GTTRF computes an LU factorization of a real tridiagonal matrix A
 *
 * @details
 * \b Purpose:
    \verbatim
     GTTRF computes an LU factorization of a real tridiagonal matrix A
     using elimination with partial pivoting and row interchanges.
     The factorization has the form
        A = L * U
     where L is a product of permutation and unit lower bidiagonal
     matrices and U is upper triangular with nonzeros in only the main
     diagonal and first two superdiagonals.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. \n
 * @param[in,out] DL
          DL is REAL array, dimension (N-1) \n
          On entry, DL must contain the (n-1) sub-diagonal elements of
          A. \n
 \n
          On exit, DL is overwritten by the (n-1) multipliers that
          define the matrix L from the LU factorization of A. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, D must contain the diagonal elements of A. \n
 \n
          On exit, D is overwritten by the n diagonal elements of the
          upper triangular matrix U from the LU factorization of A. \n
 * @param[in,out] DU
          DU is REAL array, dimension (N-1) \n
          On entry, DU must contain the (n-1) super-diagonal elements
          of A. \n
 \n
          On exit, DU is overwritten by the (n-1) elements of the first
          super-diagonal of U. \n
 * @param[out] DU2
          DU2 is REAL array, dimension (N-2) \n
          On exit, DU2 is overwritten by the (n-2) elements of the
          second super-diagonal of U. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices; for 1 <= i <= n, row i of the matrix was
          interchanged with row IPIV(i).  IPIV(i) will always be either
          i or i+1; IPIV(i) = i indicates a row interchange was not
          required. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -k, the k-th argument had an illegal value \n
          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization
                has been completed, but the factor U is exactly
                singular, and division by zero will occur if it is used
                to solve a system of equations. \n

 *  * */
    template <typename T>
    void gttrf(integer *n, T *dl, T *d, T *du, T *du2, integer *ipiv, integer *info)
    {
        gttrf(n, dl, d, du, du2, ipiv, info);
    }
    /** @}*/ // end of gttrf

    /** @defgroup gttrs gttrs
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GTTRS solves one of the systems of equations
 * @details
 * \b Purpose:
    \verbatim
     GTTRS solves one of the systems of equations
        A*X = B  or  A**T*X = B,
     with a tridiagonal matrix A using the LU factorization computed
     by SGTTRF.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations. \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T* X = B  (Transpose) \n
          = 'C':  A**T* X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] DL
          DL is REAL array, dimension (N-1) \n
          The (n-1) multipliers that define the matrix L from the
          LU factorization of A. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the upper triangular matrix U from
          the LU factorization of A. \n
 * @param[in] DU
          DU is REAL array, dimension (N-1) \n
          The (n-1) elements of the first super-diagonal of U. \n
 * @param[in] DU2
          DU2 is REAL array, dimension (N-2) \n
          The (n-2) elements of the second super-diagonal of U. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices; for 1 <= i <= n, row i of the matrix was
          interchanged with row IPIV(i).  IPIV(i) will always be either
          i or i+1; IPIV(i) = i indicates a row interchange was not
          required. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the matrix of right hand side vectors B. \n
          On exit, B is overwritten by the solution vectors X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gttrs(char *trans, integer *n, integer *nrhs, T *dl, T *d, T *du, T *du2, integer *ipiv,
               T *b, integer *ldb, integer *info)
    {
        gttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info);
    }
    /** @}*/ // end of gttrs

    /** @defgroup gtts2 gtts2
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GTTS2 solves a system of linear equations with a tridiagonal matrix using the LU
 factorization computed by sgttrf

 * @details
 * \b Purpose:
    \verbatim
    GTTS2 solves one of the systems of equations
       A*X = B  or  A**T*X = B,
    with a tridiagonal matrix A using the LU factorization computed
    by GTTRF.
    \endverbatim

 * @param[in] ITRANS
          ITRANS is INTEGER \n
          Specifies the form of the system of equations. \n
          = 0:  A * X = B  (No transpose) \n
          = 1:  A**T* X = B  (Transpose) \n
          = 2:  A**T* X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] DL
          DL is REAL array, dimension (N-1) \n
          The (n-1) multipliers that define the matrix L from the
          LU factorization of A. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the upper triangular matrix U from
          the LU factorization of A. \n
 * @param[in] DU
          DU is REAL array, dimension (N-1) \n
          The (n-1) elements of the first super-diagonal of U. \n
 * @param[in] DU2
          DU2 is REAL array, dimension (N-2) \n
          The (n-2) elements of the second super-diagonal of U. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices; for 1 <= i <= n, row i of the matrix was
          interchanged with row IPIV(i).  IPIV(i) will always be either
          i or i+1; IPIV(i) = i indicates a row interchange was not
          required. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the matrix of right hand side vectors B.
          On exit, B is overwritten by the solution vectors X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n

 *  * */
    template <typename T>
    void gtts2(integer *itrans, integer *n, integer *nrhs, T *dl, T *d, T *du, T *du2,
               integer *ipiv, T *b, integer *ldb)
    {
        gtts2(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
    }
    /** @}*/ // end of gtts2

    /** @defgroup gtrfs gtrfs
     * @ingroup LU_Computational LU_Computational
     * @{
     */
    /*! @brief GTRFS improves the computed solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
     GTRFS improves the computed solution to a system of linear
     equations when the coefficient matrix is tridiagonal, and provides
     error bounds and backward error estimates for the solution.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B     (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] DL
          DL is REAL array, dimension (N-1) \n
          The (n-1) subdiagonal elements of A. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The diagonal elements of A. \n
 * @param[in] DU
          DU is REAL array, dimension (N-1) \n
          The (n-1) superdiagonal elements of A. \n
 * @param[in] DLF
          DLF is REAL array, dimension (N-1) \n
          The (n-1) multipliers that define the matrix L from the
          LU factorization of A as computed by SGTTRF. \n
 * @param[in] DF
          DF is REAL array, dimension (N) \n
          The n diagonal elements of the upper triangular matrix U from
          the LU factorization of A. \n
 * @param[in] DUF
          DUF is REAL array, dimension (N-1) \n
          The (n-1) elements of the first superdiagonal of U. \n
 * @param[in] DU2
          DU2 is REAL array, dimension (N-2) \n
          The (n-2) elements of the second superdiagonal of U. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          The pivot indices; for 1 <= i <= n, row i of the matrix was
          interchanged with row IPIV(i).  IPIV(i) will always be either
          i or i+1; IPIV(i) = i indicates a row interchange was not
          required. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by SGTTRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gtrfs(char *trans, integer *n, integer *nrhs, T *dl, T *d, T *du, T *dlf, T *df, T *duf,
               T *du2, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr,
               T *work, T *iwork, integer *info)
    {
        gtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work,
              iwork, info);
    }
    template <typename T, typename Ta>
    void gtrfs(char *trans, integer *n, integer *nrhs, T *dl, T *d, T *du, T *dlf, T *df, T *duf,
               T *du2, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr,
               T *work, Ta *rwork, integer *info)
    {
        gtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work,
              rwork, info);
    }
    /** @}*/ // end of gtrfs

    /** @}*/ // end of LU_Computational

    /** @defgroup LDL LDL
     * @ingroup LinearSolve
     * @{
     */
    /** @defgroup  sysv sysv
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SYSV computes the solution to system of linear equations A * X = B for SY matrices

 * @details
 * \b Purpose:
    \verbatim
     SYSV computes the solution to a real system of linear equations
        A * X = B,
     where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     matrices.

     The diagonal pivoting method is used to factor A as
        A = U * D * U**T,  if UPLO = 'U', or
        A = L * D * L**T,  if UPLO = 'L',
     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and D is symmetric and block diagonal with
     1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     used to solve the system of equations A * X = B.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, if INFO = 0, the block diagonal matrix D and the
          multipliers used to obtain the factor U or L from the
          factorization A = U*D*U**T or A = L*D*L**T as computed by
          SSYTRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D, as
          determined by SSYTRF.  If IPIV(k) > 0, then rows and columns
          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
          then rows and columns k-1 and -IPIV(k) were interchanged and
          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
          diagonal block. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >= 1, and for best performance
          LWORK >= fla_max(1,N*NB), where NB is the optimal blocksize for
          SSYTRF. \n
          for LWORK < N, TRS will be done with Level BLAS 2 \n
          for LWORK >= N, TRS will be done with Level BLAS 3 \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular, so the solution could not be computed. \n

 *  * */
    template <typename T>
    void sysv(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
              integer *ldb, T *work, integer *lwork, integer *info)
    {
        sysv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of sysv

    /** @defgroup sysv_rook sysv_rook
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SSYSV_ROOK computes the solution to system of linear equations A * X = B for SY
 matrices

 * @details
 * \b Purpose:
    \verbatim
     SYSV_ROOK computes the solution to a real system of linear
     equations
        A * X = B,
     where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     matrices.

     The diagonal pivoting method is used to factor A as
        A = U * D * U**T,  if UPLO = 'U', or
        A = L * D * L**T,  if UPLO = 'L',
     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and D is symmetric and block diagonal with
     1-by-1 and 2-by-2 diagonal blocks.

     SSYTRF_ROOK is called to compute the factorization of a real
     symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     pivoting method.

     The factored form of A is then used to solve the system
     of equations A * X = B by calling SYTRS_ROOK.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, if INFO = 0, the block diagonal matrix D and the
          multipliers used to obtain the factor U or L from the
          factorization A = U*D*U**T or A = L*D*L**T as computed by
          SSYTRF_ROOK. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D,
          as determined by SSYTRF_ROOK. \n
 \n
          If UPLO = 'U': \n
               If IPIV(k) > 0, then rows and columns k and IPIV(k)
               were interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
               If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
               columns k and -IPIV(k) were interchanged and rows and
               columns k-1 and -IPIV(k-1) were inerchaged,
               D(k-1:k,k-1:k) is a 2-by-2 diagonal block. \n
 \n
          If UPLO = 'L': \n
               If IPIV(k) > 0, then rows and columns k and IPIV(k)
               were interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
               If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
               columns k and -IPIV(k) were interchanged and rows and
               columns k+1 and -IPIV(k+1) were inerchaged,
               D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B. \n
          On exit, if INFO = 0, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >= 1, and for best performance
          LWORK >= fla_max(1,N*NB), where NB is the optimal blocksize for
          SSYTRF_ROOK. \n
 \n
          TRS will be done with Level 2 BLAS \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular, so the solution could not be computed. \n

 *  * */
    template <typename T>
    void sysv_rook(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
                   integer *ldb, T *work, integer *lwork, integer *info)
    {
        sysv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of sysvrook

    /** @defgroup sysv_rk sysv_rk
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SSYSV_RK computes the solution to system of linear equations A * X = B for SY
 matrices

 * @details
 * \b Purpose:
    \verbatim
     SYSV_RK computes the solution to a real system of linear
     equations A * X = B, where A is an N-by-N symmetric matrix
     and X and B are N-by-NRHS matrices.

     The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     to factor A as
        A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
        A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     where U (or L) is unit upper (or lower) triangular matrix,
     U**T (or L**T) is the transpose of U (or L), P is a permutation
     matrix, P**T is the transpose of P, and D is symmetric and block
     diagonal with 1-by-1 and 2-by-2 diagonal blocks.

     SSYTRF_RK is called to compute the factorization of a real
     symmetric matrix.  The factored form of A is then used to solve
     the system of equations A * X = B by calling BLAS3 routine SYTRS_3.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A. \n
            If UPLO = 'U': the leading N-by-N upper triangular part
            of A contains the upper triangular part of the matrix A,
            and the strictly lower triangular part of A is not
            referenced. \n
 \n
            If UPLO = 'L': the leading N-by-N lower triangular part
            of A contains the lower triangular part of the matrix A,
            and the strictly upper triangular part of A is not
            referenced. \n
 \n
          On exit, if INFO = 0, diagonal of the block diagonal
          matrix D and factors U or L  as computed by SSYTRF_RK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                are stored on exit in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 \n
          For more info see the description of DSYTRF_RK routine. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] E
          E is REAL array, dimension (N) \n
          On exit, contains the output computed by the factorization
          routine DSYTRF_RK, i.e. the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where \n
          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; \n
          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. \n
 \n
          NOTE: For 1-by-1 diagonal block D(k), where
          1 <= k <= N, the element E(k) is set to 0 in both
          UPLO = 'U' or UPLO = 'L' cases. \n
 \n
          For more info see the description of DSYTRF_RK routine. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D,
          as determined by SSYTRF_RK. \n
 \n
          For more info see the description of DSYTRF_RK routine. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension ( MAX(1,LWORK) ). \n
          Work array used in the factorization stage.
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >= 1. For best performance
          of factorization stage LWORK >= fla_max(1,N*NB), where NB is
          the optimal blocksize for DSYTRF_RK. \n
 \n
          If LWORK = -1, then a workspace query is assumed;
          the routine only calculates the optimal size of the WORK
          array for factorization stage, returns this value as
          the first entry of the WORK array, and no error message
          related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: If INFO = -k, the k-th argument had an illegal value \n\
          > 0: If INFO = k, the matrix A is singular, because:
                 If UPLO = 'U': column k in the upper
                 triangular part of A contains all zeros.
                 If UPLO = 'L': column k in the lower
                 triangular part of A contains all zeros. \n
 \n
               Therefore D(k,k) is exactly zero, and superdiagonal
               elements of column k of U (or subdiagonal elements of
               column k of L ) are all zeros. The factorization has
               been completed, but the block diagonal matrix D is
               exactly singular, and division by zero will occur if
               it is used to solve a system of equations. \n
 \n
               NOTE: INFO only stores the first occurrence of
               a singularity, any subsequent occurrence of singularity
               is not stored in INFO even though the factorization
               always completes. \n

 *  * */
    template <typename T>
    void sysv_rk(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *e, integer *ipiv,
                 T *b, integer *ldb, T *work, integer *lwork, integer *info)
    {
        sysv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of sysv_rk

    /** @defgroup sysvx sysvx
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SYSVX computes the solution to system of linear equations A * X = B for SY matrices

 * @details
 * \b Purpose:
    \verbatim
     SYSVX uses the diagonal pivoting factorization to compute the
     solution to a real system of linear equations A * X = B,
     where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     matrices.

     Error bounds on the solution and a condition estimate are also
     provided.

     * \b Description:
         ==============
     The following steps are performed:

     1. If FACT = 'N', the diagonal pivoting method is used to factor A.
        The form of the factorization is
           A = U * D * U**T,  if UPLO = 'U', or
           A = L * D * L**T,  if UPLO = 'L',
        where U (or L) is a product of permutation and unit upper (lower)
        triangular matrices, and D is symmetric and block diagonal with
        1-by-1 and 2-by-2 diagonal blocks.

     2. If some D(i,i)=0, so that D is exactly singular, then the routine
        returns with INFO = i. Otherwise, the factored form of A is used
        to estimate the condition number of the matrix A.  If the
        reciprocal of the condition number is less than machine precision,
        INFO = N+1 is returned as a warning, but the routine still goes on
        to solve for X and compute error bounds as described below.

     3. The system of equations is solved for X using the factored form
        of A.

     4. Iterative refinement is applied to improve the computed solution
        matrix and calculate error bounds and backward error estimates
        for it.
     \endverbatim

 * @param[in] FACT
          FACT is CHARACTER*1 \n
          Specifies whether or not the factored form of A has been
          supplied on entry. \n
          = 'F':  On entry, AF and IPIV contain the factored form of
                  A.  AF and IPIV will not be modified. \n
          = 'N':  The matrix A will be copied to AF and factored. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of A contains the upper triangular part
          of the matrix A, and the strictly lower triangular part of A
          is not referenced.  If UPLO = 'L', the leading N-by-N lower
          triangular part of A contains the lower triangular part of
          the matrix A, and the strictly upper triangular part of A is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] AF
          AF is REAL array, dimension (LDAF,N) \n
          If FACT = 'F', then AF is an input argument and on entry
          contains the block diagonal matrix D and the multipliers used
          to obtain the factor U or L from the factorization
          A = U*D*U**T or A = L*D*L**T as computed by SSYTRF. \n
 \n
          If FACT = 'N', then AF is an output argument and on exit
          returns the block diagonal matrix D and the multipliers used
          to obtain the factor U or L from the factorization
          A = U*D*U**T or A = L*D*L**T. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          If FACT = 'F', then IPIV is an input argument and on entry
          contains details of the interchanges and the block structure
          of D, as determined by SSYTRF. \n
          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
          interchanged and D(k,k) is a 1-by-1 diagonal block.
          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 \n
          If FACT = 'N', then IPIV is an output argument and on exit
          contains details of the interchanges and the block structure
          of D, as determined by SSYTRF. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The N-by-NRHS right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
          X is REAL array, dimension (LDX,NRHS) \n
          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          The estimate of the reciprocal condition number of the matrix
          A.  If RCOND is less than the machine precision (in
          particular, if RCOND = 0), the matrix is singular to working
          precision.  This condition is indicated by a return code of
          INFO > 0. \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >= fla_max(1,3*N), and for best
          performance, when FACT = 'N', LWORK >= fla_max(1,3*N,N*NB), where
          NB is the optimal blocksize for SSYTRF. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, and i is \n
                <= N:  D(i,i) is exactly zero.  The factorization
                       has been completed but the factor D is exactly
                       singular, so the solution and error bounds could
                       not be computed. RCOND = 0 is returned. \n
                = N+1: D is nonsingular, but RCOND is less than machine
                       precision, meaning that the matrix is singular
                       to working precision.  Nevertheless, the
                       solution and error bounds are computed because
                       there are a number of situations where the
                       computed solution can be more accurate than the
                       value of RCOND would suggest. \n

 *  * */
    template <typename T>
    void sysvx(char *fact, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af,
               integer *ldaf, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, T *rcond,
               T *ferr, T *berr, T *work, integer *lwork, integer *iwork, integer *info)
    {
        sysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
              lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void sysvx(char *fact, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af,
               integer *ldaf, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, Ta *rcond,
               Ta *ferr, Ta *berr, T *work, integer *lwork, Ta *rwork, integer *info)
    {
        sysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
              lwork, rwork, info);
    }
    /** @}*/ // end of sysvx

    /** @defgroup sysvxx sysvxx
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief  SYSVXX uses the diagonal pivoting factorization to compute the \n
    solution to a real system of linear equations

 * @details
 * \b Purpose:
    \verbatim
    SYSVXX uses the diagonal pivoting factorization to compute the
    solution to a real system of linear equations A * X = B, where A
    is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.

    If requested, both normwise and maximum componentwise error bounds
    are returned. SSYSVXX will return a solution with a tiny
    guaranteed error (O(eps) where eps is the working machine
    precision) unless the matrix is very ill-conditioned, in which
    case a warning is returned. Relevant condition numbers also are
    calculated and returned.

    SSYSVXX accepts user-provided factorizations and equilibration
    factors; see the definitions of the FACT and EQUED options.
    Solving with refinement and using a factorization from a previous
    SSYSVXX call will also produce a solution with either O(eps)
    errors or warnings, but we cannot make that claim for general
    user-provided factorizations and equilibration factors if they
    differ from what SYSVXX would itself produce.

 * \b Description:
      ============
    The following steps are performed:

    1. If FACT = 'E', real scaling factors are computed to equilibrate
    the system:

      diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B

    Whether or not the system will be equilibrated depends on the
    scaling of the matrix A, but if equilibration is used, A is
    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.

    2. If FACT = 'N' or 'E', the LU decomposition is used to factor
    the matrix A (after equilibration if FACT = 'E') as

       A = U * D * U**T,  if UPLO = 'U', or
       A = L * D * L**T,  if UPLO = 'L',

    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, and D is symmetric and block diagonal with
    1-by-1 and 2-by-2 diagonal blocks.

    3. If some D(i,i)=0, so that D is exactly singular, then the
    routine returns with INFO = i. Otherwise, the factored form of A
    is used to estimate the condition number of the matrix A (see
    argument RCOND).  If the reciprocal of the condition number is
    less than machine precision, the routine still goes on to solve
    for X and compute error bounds as described below.

    4. The system of equations is solved for X using the factored form
    of A.

    5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),
    the routine will use iterative refinement to try to get a small
    error and error bounds.  Refinement calculates the residual to at
    least twice the working precision.

    6. If equilibration was used, the matrix X is premultiplied by
    diag(R) so that it solves the original system before
    equilibration.
    \endverbatim

 * @param[in] FACT
          FACT is CHARACTER*1 \n
          Specifies whether or not the factored form of the matrix A is
          supplied on entry, and if not, whether the matrix A should be
          equilibrated before it is factored. \n
           = 'F':  On entry, AF and IPIV contain the factored form of A.
                   If EQUED is not 'N', the matrix A has been
                   equilibrated with scaling factors given by S.
                   A, AF, and IPIV are not modified. \n
           = 'N':  The matrix A will be copied to AF and factored. \n
           = 'E':  The matrix A will be equilibrated if necessary, then
                     copied to AF and factored. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of A contains the upper triangular
          part of the matrix A, and the strictly lower triangular
          part of A is not referenced.  If UPLO = 'L', the leading
          N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
          diag(S)*A*diag(S). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] AF
          AF is REAL array, dimension (LDAF,N) \n
          If FACT = 'F', then AF is an input argument and on entry
          contains the block diagonal matrix D and the multipliers
          used to obtain the factor U or L from the factorization A =
          U*D*U**T or A = L*D*L**T as computed by SSYTRF. \n
 \n
          If FACT = 'N', then AF is an output argument and on exit
          returns the block diagonal matrix D and the multipliers
          used to obtain the factor U or L from the factorization A =
          U*D*U**T or A = L*D*L**T. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          If FACT = 'F', then IPIV is an input argument and on entry
          contains details of the interchanges and the block
          structure of D, as determined by SSYTRF.  If IPIV(k) > 0,
          then rows and columns k and IPIV(k) were interchanged and
          D(k,k) is a 1-by-1 diagonal block.  If UPLO = 'U' and
          IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and
          -IPIV(k) were interchanged and D(k-1:k,k-1:k) is a 2-by-2
          diagonal block.  If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0,
          then rows and columns k+1 and -IPIV(k) were interchanged
          and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 \n
          If FACT = 'N', then IPIV is an output argument and on exit
          contains details of the interchanges and the block
          structure of D, as determined by SSYTRF. \n
 * @param[in,out] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done. \n
           = 'N':  No equilibration (always true if FACT = 'N'). \n
           = 'Y':  Both row and column equilibration, i.e., A has been
                   replaced by diag(S) * A * diag(S). \n
          EQUED is an input argument if FACT = 'F'; otherwise, it is an
          output argument. \n
 * @param[in,out] S
          S is REAL array, dimension (N) \n
          The scale factors for A.  If EQUED = 'Y', A is multiplied on
          the left and right by diag(S).  S is an input argument if FACT =
          'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED
          = 'Y', each element of S must be positive.  If S is output, each
          element of S is a power of the radix. If S is input, each element
          of S should be a power of the radix to ensure a reliable solution
          and error estimates. Scaling by powers of the radix does not cause
          rounding errors unless the result underflows or overflows.
          Rounding errors during scaling lead to refining with a matrix that
          is not equivalent to the input matrix, producing error estimates
          that may not be reliable. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B. \n
          On exit,
          if EQUED = 'N', B is not modified; \n
          if EQUED = 'Y', B is overwritten by diag(S)*B; \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
          X is REAL array, dimension (LDX,NRHS) \n
          If INFO = 0, the N-by-NRHS solution matrix X to the original
          system of equations.  Note that A and B are modified on exit if
          EQUED .ne. 'N', and the solution to the equilibrated system is
          inv(diag(S))*X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          Reciprocal scaled condition number.  This is an estimate of the
          reciprocal Skeel condition number of the matrix A after
          equilibration (if done).  If this is less than the machine
          precision (in particular, if it is zero), the matrix is singular
          to working precision.  Note that the error may still be small even
          if this number is very small and the matrix appears ill-
          conditioned. \n
 * @param[out] RPVGRW
          RPVGRW is REAL \n
          Reciprocal pivot growth.  On exit, this contains the reciprocal
          pivot growth factor norm(A)/norm(U). The "max absolute element"
          norm is used.  If this is much less than 1, then the stability of
          the LU factorization of the (equilibrated) matrix A could be poor.
          This also means that the solution X, estimated condition numbers,
          and error bounds could be unreliable. If factorization fails with
          0<INFO<=N, then this contains the reciprocal pivot growth factor
          for the leading INFO columns of A. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          Componentwise relative backward error.  This is the
          componentwise relative backward error of each solution vector X(j)
          (i.e., the smallest relative change in any element of A or B that
          makes X(j) an exact solution). \n
 * @param[in] N_ERR_BNDS
          N_ERR_BNDS is INTEGER \n
          Number of error bounds to return for each right hand side
          and each type (normwise or componentwise).  See ERR_BNDS_NORM and
          ERR_BNDS_COMP below. \n
 * @param[out] ERR_BNDS_NORM
          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS)
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          normwise relative error, which is defined as follows: \n
 \n
          Normwise relative error in the ith solution vector: \n

          \f[ \frac{{max\_j}\  abs(XTRUE(j,i) - X(j,i))}{{max\_j}\ abs(X(j,i))} \f] \n

 \n
          The array is indexed by the type of error information as described
          below. There currently are up to three pieces of information
          returned. \n
 \n
          The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_NORM(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated normwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*A, where S scales each row by a power of the
                  radix so all absolute row sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[out] ERR_BNDS_COMP
          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS)
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          componentwise relative error, which is defined as follows: \n
 \n
          Componentwise relative error in the ith solution vector: \n

          \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
          The array is indexed by the right-hand side i (on which the
          componentwise relative error depends), and the type of error
          information as described below. There currently are up to three
          pieces of information returned for each right-hand side. If
          componentwise accuracy is not requested (PARAMS(3) = 0.0), then
          ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
          the first (:,N_ERR_BNDS) entries are returned. \n
 \n
          The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_COMP(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated componentwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*(A*diag(x)), where x is the solution for the
                  current right-hand side and S scales each row of
                  A*diag(x) by a power of the radix so all absolute row
                  sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in] NPARAMS
          NPARAMS is INTEGER \n
          Specifies the number of parameters set in PARAMS.  If <= 0, the
          PARAMS array is never referenced and default values are used. \n
 * @param[in,out] PARAMS
          PARAMS is REAL array, dimension NPARAMS \n
          Specifies algorithm parameters.  If an entry is < 0.0, then
          that entry will be filled with default value used for that
          parameter.  Only positions up to NPARAMS are accessed; defaults
          are used for higher-numbered parameters. \n
 \n
           PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
                refinement or not. \n
             Default: 1.0 \n
                = 0.0:  No refinement is performed, and no error bounds are
                        computed. \n
                = 1.0:  Use the double-precision refinement algorithm,
                        possibly with doubled-single computations if the
                        compilation environment does not support DOUBLE
                        PRECISION. \n
                  (other values are reserved for future use) \n
 \n
           PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
                computations allowed for refinement. \n
             Default: 10 \n
             Aggressive: Set to 100 to permit convergence using approximate
                         factorizations or factorizations other than LU. If
                         the factorization uses a technique other than
                         Gaussian elimination, the guarantees in
                         err_bnds_norm and err_bnds_comp may no longer be
                         trustworthy. \n
 \n
           PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
                will attempt to find a solution with small componentwise
                relative error in the double-precision algorithm.  Positive
                is true, 0.0 is false. \n
             Default: 1.0 (attempt componentwise convergence) \n
 * @param[out]	WORK
          WORK is REAL array, dimension (4*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  Successful exit. The solution to every right-hand side is
           guaranteed. \n
          < 0:  If INFO = -i, the i-th argument had an illegal value \n
          > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
           has been completed, but the factor U is exactly singular, so
           the solution and error bounds could not be computed. RCOND = 0
           is returned. \n
          = N+J: The solution corresponding to the Jth right-hand side is
           not guaranteed. The solutions corresponding to other right-
           hand sides K with K > J may not be guaranteed as well, but
           only the first such right-hand side is reported. If a small
           componentwise error is not requested (PARAMS(3) = 0.0) then
           the Jth right-hand side is the first with a normwise error
           bound that is not guaranteed (the smallest J such
           that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
           the Jth right-hand side is the first with either a normwise or
           componentwise error bound that is not guaranteed (the smallest
           J such that either ERR_BNDS_NORM(J,1) = 0.0 or
           ERR_BNDS_COMP(J,1) = 0.0). See the definition of
           ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
           about all of the right-hand sides check ERR_BNDS_NORM or
           ERR_BNDS_COMP. \n

 *  * */
    template <typename T>
    void sysvxx(char *fact, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, char *equed, T *s, T *b, integer *ldb, T *x,
                integer *ldx, T *rcond, T *rpvgrw, T *berr, integer *n_err_bnds, T *err_bnds_norm,
                T *err_bnds_comp, integer *nparams, T *params, T *work, integer *iwork,
                integer *info)
    {
        sysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw,
               berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
    }
    template <typename T, typename Ta>
    void sysvxx(char *fact, char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, char *equed, Ta *s, T *b, integer *ldb, T *x,
                integer *ldx, Ta *rcond, Ta *rpvgrw, Ta *berr, integer *n_err_bnds,
                Ta *err_bnds_norm, Ta *err_bnds_comp, integer *nparams, Ta *params, T *work,
                Ta *rwork, integer *info)
    {
        sysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx, rcond, rpvgrw,
               berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    }
    /** @}*/ // end of sysvxx

    /** @defgroup spsv spsv
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SPSV computes the solution to system of linear equations A * X = B for OTHER matrices

 * @details
 * \b Purpose:
    \verbatim
     SPSV computes the solution to a real system of linear equations
        A * X = B,
     where A is an N-by-N symmetric matrix stored in packed format and X
     and B are N-by-NRHS matrices.

     The diagonal pivoting method is used to factor A as
        A = U * D * U**T,  if UPLO = 'U', or
        A = L * D * L**T,  if UPLO = 'L',
     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, D is symmetric and block diagonal with 1-by-1
     and 2-by-2 diagonal blocks.  The factored form of A is then used to
     solve the system of equations A * X = B.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
          See below for further details. \n
 \n
          On exit, the block diagonal matrix D and the multipliers used
          to obtain the factor U or L from the factorization
          A = U*D*U**T or A = L*D*L**T as computed by SSPTRF, stored as
          a packed triangular matrix in the same storage format as A. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D, as
          determined by SSPTRF.  If IPIV(k) > 0, then rows and columns
          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
          then rows and columns k-1 and -IPIV(k) were interchanged and
          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
          diagonal block. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
                has been completed, but the block diagonal matrix D is
                exactly singular, so the solution could not be
                computed. \n

 *  * */
    template <typename T>
    void spsv(char *uplo, integer *n, integer *nrhs, T *ap, integer *ipiv, T *b, integer *ldb,
              integer *info)
    {
        spsv(uplo, n, nrhs, ap, ipiv, b, ldb, info);
    }
    template <typename T>
    void hpsv(char *uplo, integer *n, integer *nrhs, T *ap, integer *ipiv, T *b, integer *ldb,
              integer *info)
    {
        hpsv(uplo, n, nrhs, ap, ipiv, b, ldb, info);
    }
    /** @}*/ // end of spsv

    /** @defgroup spsvx spsvx
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SPSVX computes the solution to system of linear equations A * X = B for OTHER
 matrices

 * @details
 * \b Purpose:
    \verbatim
     SPSVX uses the diagonal pivoting factorization A = U*D*U**T or
     A = L*D*L**T to compute the solution to a real system of linear
     equations A * X = B, where A is an N-by-N symmetric matrix stored
     in packed format and X and B are N-by-NRHS matrices.

     Error bounds on the solution and a condition estimate are also
     provided.

     * \b Description:
       =================

     The following steps are performed:

     1. If FACT = 'N', the diagonal pivoting method is used to factor A as
           A = U * D * U**T,  if UPLO = 'U', or
           A = L * D * L**T,  if UPLO = 'L',
        where U (or L) is a product of permutation and unit upper (lower)
        triangular matrices and D is symmetric and block diagonal with
        1-by-1 and 2-by-2 diagonal blocks.

     2. If some D(i,i)=0, so that D is exactly singular, then the routine
        returns with INFO = i. Otherwise, the factored form of A is used
        to estimate the condition number of the matrix A.  If the
        reciprocal of the condition number is less than machine precision,
        INFO = N+1 is returned as a warning, but the routine still goes on
        to solve for X and compute error bounds as described below.

     3. The system of equations is solved for X using the factored form
        of A.

     4. Iterative refinement is applied to improve the computed solution
        matrix and calculate error bounds and backward error estimates
        for it.
    \endverbatim

 * @param[in] FACT
          FACT is CHARACTER*1 \n
          Specifies whether or not the factored form of A has been
          supplied on entry. \n
          = 'F':  On entry, AFP and IPIV contain the factored form of
                  A.  AP, AFP and IPIV will not be modified. \n
          = 'N':  The matrix A will be copied to AFP and factored. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangle of the symmetric matrix A, packed
          columnwise in a linear array.  The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
          See below for further details. \n
 * @param[in,out] AFP
          AFP is REAL array, dimension (N*(N+1)/2) \n
          If FACT = 'F', then AFP is an input argument and on entry
          contains the block diagonal matrix D and the multipliers used
          to obtain the factor U or L from the factorization
          A = U*D*U**T or A = L*D*L**T as computed by SSPTRF, stored as
          a packed triangular matrix in the same storage format as A. \n
 \n
          If FACT = 'N', then AFP is an output argument and on exit
          contains the block diagonal matrix D and the multipliers used
          to obtain the factor U or L from the factorization
          A = U*D*U**T or A = L*D*L**T as computed by SSPTRF, stored as
          a packed triangular matrix in the same storage format as A. \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          If FACT = 'F', then IPIV is an input argument and on entry
          contains details of the interchanges and the block structure
          of D, as determined by SSPTRF. \n
          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
          interchanged and D(k,k) is a 1-by-1 diagonal block. \n
          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 \n
          If FACT = 'N', then IPIV is an output argument and on exit
          contains details of the interchanges and the block structure
          of D, as determined by SSPTRF. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The N-by-NRHS right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out] X
          X is REAL array, dimension (LDX,NRHS) \n
          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          The estimate of the reciprocal condition number of the matrix
          A.  If RCOND is less than the machine precision (in
          particular, if RCOND = 0), the matrix is singular to working
          precision.  This condition is indicated by a return code of
          INFO > 0. \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, and i is \n
                <= N:  D(i,i) is exactly zero.  The factorization
                       has been completed but the factor D is exactly
                       singular, so the solution and error bounds could
                       not be computed. RCOND = 0 is returned. \n
                = N+1: D is nonsingular, but RCOND is less than machine
                       precision, meaning that the matrix is singular
                       to working precision.  Nevertheless, the
                       solution and error bounds are computed because
                       there are a number of situations where the
                       computed solution can be more accurate than the
                       value of RCOND would suggest. \n

 *  * */
    template <typename T>
    void spsvx(char *fact, char *uplo, integer *n, integer *nrhs, T *ap, T *afp, integer *ipiv,
               T *b, integer *ldb, T *x, integer *ldx, T *rcond, T *ferr, T *berr, T *work,
               integer *iwork, integer *info)
    {
        spsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, iwork,
              info);
    }
    template <typename T, typename Ta>
    void spsvx(char *fact, char *uplo, integer *n, integer *nrhs, T *ap, T *afp, integer *ipiv,
               T *b, integer *ldb, T *x, integer *ldx, Ta *rcond, Ta *ferr, Ta *berr, T *work,
               Ta *rwork, integer *info)
    {
        spsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork,
              info);
    }
    template <typename T, typename Ta>
    void hpsvx(char *fact, char *uplo, integer *n, integer *nrhs, T *ap, T *afp, integer *ipiv,
               T *b, integer *ldb, T *x, integer *ldx, Ta *rcond, Ta *ferr, Ta *berr, T *work,
               Ta *rwork, integer *info)
    {
        hpsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork,
              info);
    }
    /** @}*/ // end of spsvx

    /** @defgroup sysv_aa sysv_aa
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SYSV_AA computes the solution to system of linear equations A * X = B for SY matrices

 * @details
 * \b Purpose:
    \verbatim
    SYSV computes the solution to a real system of linear equations
       A * X = B,
    where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
    matrices.

    Aasen's algorithm is used to factor A as
       A = U**T * T * U,  if UPLO = 'U', or
       A = L * T * L**T,  if UPLO = 'L',
    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, and T is symmetric tridiagonal. The factored
    form of A is then used to solve the system of equations A * X = B.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, if INFO = 0, the tridiagonal matrix T and the
          multipliers used to obtain the factor U or L from the
          factorization A = U**T*T*U or A = L*T*L**T as computed by
          SSYTRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of A were interchanged with the
          row and column IPIV(k). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N-by-NRHS right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >= MAX(1,2*N,3*N-2), and for
          the best performance, LWORK >= MAX(1,N*NB), where NB is
          the optimal blocksize for SSYTRF_AA. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular, so the solution could not be computed. \n

 *  * */
    template <typename T>
    void sysv_aa(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
                 integer *ldb, T *work, integer *lwork, integer *info)
    {
        sysv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of sysv_Aa

    /** @defgroup sysv_aa_2stage sysv_aa_2stage
     * @ingroup LDL LDL
     * @{
     */
    /*! @brief SYSV_AA_2STAGE computes the solution to system of linear equations A * X = B for SY
 matrices

 * @details
 * \b Purpose:
    \verbatim
     SYSV_AA_2STAGE computes the solution to a real system of
     linear equations
        A * X = B,
     where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     matrices.

     Aasen's 2-stage algorithm is used to factor A as
        A = U**T * T * U,  if UPLO = 'U', or
        A = L * T * L**T,  if UPLO = 'L',
     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and T is symmetric and band. The matrix T is
     then LU-factored with partial pivoting. The factored form of A
     is then used to solve the system of equations A * X = B.

     This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, L is stored below (or above) the subdiaonal blocks,
          when UPLO  is 'L' (or 'U'). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] TB
          TB is REAL array, dimension (LTB) \n
          On exit, details of the LU factorization of the band matrix. \n
 * @param[in] LTB
          LTB is INTEGER \n
          The size of the array TB. LTB >= 4*N, internally
          used to select NB such that LTB >= (3*NB+1)*N. \n
 \n
          If LTB = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of LTB,
          returns this value as the first entry of TB, and
          no error message related to LTB is issued by XERBLA. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of A were interchanged with the
          row and column IPIV(k). \n
 * @param[out] IPIV2
          IPIV2 is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of T were interchanged with the
          row and column IPIV(k). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS)
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL workspace of size LWORK \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The size of WORK. LWORK >= N, internally used to select NB
          such that LWORK >= N*NB. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the WORK array,
          returns this value as the first entry of the WORK array, and
          no error message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, band LU factorization failed on i-th column \n

 *  * */
    template <typename T>
    void sysv_aa_2stage(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *tb,
                        integer *ltb, integer *ipiv, integer *ipiv2, T *b, integer *ldb, T *work,
                        integer *lwork, integer *info)
    {
        sysv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of sysv_aa_2stage
    /** @}*/ // end of LDL

    /** @defgroup LDL_computation LDL computation
     * @ingroup LinearSolve
     * @{
     */
    /** @defgroup hecon hecon
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief HECON estimates the reciprocal of the condition number of a complex Hermitian matrix
 A

 * @details
 * \b Purpose:
    \verbatim
    HECON estimates the reciprocal of the condition number of a complex
    Hermitian matrix A using the factorization A = U*D*U**H or
    A = L*D*L**H computed by CHETRF.

    An estimate is obtained for norm(inv(A)), and the reciprocal of the
    condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**H; \n
          = 'L':  Lower triangular, form is A = L*D*L**H. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by CHETRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by CHETRF. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T, typename Ta>
    void hecon(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, Ta *anorm, Ta *rcond,
               T *work, integer *info)
    {
        hecon(uplo, n, a, lda, ipiv, anorm, rcond, work, info);
    }
    /** @}*/ // end of hecon

    /** @defgroup sycon sycon
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYCON estimates the reciprocal of the condition number of a real symmetric matrix A

 * @details
 * \b Purpose:
    \verbatim
     SYCON estimates the reciprocal of the condition number (in the
     1-norm) of a real symmetric matrix A using the factorization
     A = U*D*U**T or A = L*D*L**T computed by SYTRF.

     An estimate is obtained for norm(inv(A)), and the reciprocal of the
     condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by SSYTRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (2*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sycon(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *anorm, T *rcond,
               T *work, integer *iwork, integer *info)
    {
        sycon(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void sycon(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, Ta *anorm, Ta *rcond,
               T *work, integer *info)
    {
        sycon(uplo, n, a, lda, ipiv, anorm, rcond, work, info);
    }
    /** @}*/ // end of sycon

    /** @defgroup sytrf sytrf
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRF computes the factorization of a real symmetric matrix A \n
     using the Bunch-Kaufman diagonal pivoting method
 * @details
 * \b Purpose:
    \verbatim
     SYTRF computes the factorization of a real symmetric matrix A using
     the Bunch-Kaufman diagonal pivoting method.  The form of the
     factorization is

        A = U**T*D*U  or  A = L*D*L**T

     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and D is symmetric and block diagonal with
     1-by-1 and 2-by-2 diagonal blocks.

     This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, the block diagonal matrix D and the multipliers used
          to obtain the factor U or L (see below for further details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D.
          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
          interchanged and D(k,k) is a 1-by-1 diagonal block. \n
          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >=1.  For best performance
          LWORK >= N*NB, where NB is the block size returned by ILAENV. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
                has been completed, but the block diagonal matrix D is
                exactly singular, and division by zero will occur if it
                is used to solve a system of equations. \n

 *  * */
    template <typename T>
    void sytrf(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *work, integer *lwork,
               integer *info)
    {
        sytrf(uplo, n, a, lda, ipiv, work, lwork, info);
    }
    /** @}*/ // end of sytrf

    /** @defgroup lasyf lasyf
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LASYF computes a partial factorization of a real symmetric matrix using the
 Bunch-Kaufman diagonal pivoting method

 * @details
 * \b Purpose:
    \verbatim
    LASYF computes a partial factorization of a real symmetric matrix A
    using the Bunch-Kaufman diagonal pivoting method. The partial
    factorization has the form:

    A  =  (I  U12) (A11  0 ) ( I       0   )  if UPLO = 'U', or:
          (0  U22) ( 0   D ) (U12**T U22**T)

    A  =  (L11  0) ( D   0 ) (L11**T L21**T)  if UPLO = 'L'
          (L21  I) ( 0  A22) ( 0       I   )

    where the order of D is at most NB. The actual order is   returned in
    the argument KB, and is either NB or NB-1, or N if N <= NB.

    SLASYF is an auxiliary routine called by SSYTRF. It uses blocked code
    (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
    A22 (if UPLO = 'L').
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NB
          NB is INTEGER \n
          The maximum number of columns of the matrix A that should be
          factored.  NB should be at least 2 to allow for 2-by-2 pivot
          blocks. \n
 * @param[out] KB
          KB is INTEGER \n
          The number of columns of A that were actually factored.
          KB is either NB-1 or NB, or N if N <= NB. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          n-by-n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n-by-n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
          On exit, A contains details of the partial factorization. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D.
 \n
          If UPLO = 'U': \n
             Only the last KB elements of IPIV are set.
 \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
             interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
             is a 2-by-2 diagonal block.
 \n
          If UPLO = 'L': \n
             Only the first KB elements of IPIV are set.
 \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
             interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
             is a 2-by-2 diagonal block. \n
 * @param[out] W
          W is REAL array, dimension (LDW,NB) \n
 * @param[in] LDW
          LDW is INTEGER \n
          The leading dimension of the array W.  LDW >= fla_max(1,N). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          > 0: if INFO = k, D(k,k) is exactly zero. The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular.  \n

 *  * */
    template <typename T>
    void lasyf(char *uplo, integer *n, integer *nb, integer *kb, T *a, integer *lda, integer *ipiv,
               T *w, integer *ldw, integer *info)
    {
        lasyf(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
    }
    /** @}*/ // end of lasyf

    /** @defgroup sytf2 sytf2
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTF2 computes the factorization of a real symmetric indefinite matrix, using the
 diagonal pivoting method

 * @details
 * \b Purpose:
    \verbatim
    SYTF2 computes the factorization of a real symmetric matrix A using
    the Bunch-Kaufman diagonal pivoting method:

       A = U*D*U**T  or  A = L*D*L**T

    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, U**T is the transpose of U, and D is symmetric and
    block diagonal with 1-by-1 and 2-by-2 diagonal blocks.

    This is the unblocked version of the algorithm, calling Level 2 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          n-by-n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n-by-n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, the block diagonal matrix D and the multipliers used
          to obtain the factor U or L (see below for further details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D.
 \n
          If UPLO = 'U': \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
             interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) = IPIV(k-1) < 0, then rows and columns
             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
             is a 2-by-2 diagonal block.
 \n
          If UPLO = 'L': \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
             interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) = IPIV(k+1) < 0, then rows and columns
             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)
             is a 2-by-2 diagonal block. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -k, the k-th argument had an illegal value \n
          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular, and division by zero will occur if it
               is used to solve a system of equations. \n

 *  * */
    template <typename T>
    void sytf2(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, integer *info)
    {
        sytf2(uplo, n, a, lda, ipiv, info);
    }
    /** @}*/ // end of sytf2

    /** @defgroup sytrs sytrs
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRS solves a system of linear equations A*X = B with a real \n
     symmetric matrix A
 * @details
 * \b Purpose:
    \verbatim
     SYTRS solves a system of linear equations A*X = B with a real
     symmetric matrix A using the factorization A = U*D*U**T or
     A = L*D*L**T computed by SYTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by SSYTRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sytrs(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
               integer *ldb, integer *info)
    {
        sytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info);
    }
    /** @}*/ // end of sytrs

    /** @defgroup sytri sytri
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRI computes the inverse of a real symmetric indefinite matrix A \n
    using the factorization A = U*D*U**T or A = L*D*L**T computed by SYTRF
 * @details
 * \b Purpose:
    \verbatim
    SYTRI computes the inverse of a real symmetric indefinite matrix
    A using the factorization A = U*D*U**T or A = L*D*L**T computed by
    SYTRF.
    \endverbatim
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the block diagonal matrix D and the multipliers
          used to obtain the factor U or L as computed by SSYTRF. \n
 \n
          On exit, if INFO = 0, the (symmetric) inverse of the original
          matrix.  If UPLO = 'U', the upper triangular part of the
          inverse is formed and the part of A below the diagonal is not
          referenced; if UPLO = 'L' the lower triangular part of the
          inverse is formed and the part of A above the diagonal is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SYTRF. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
               inverse could not be computed. \n

 *  * */
    template <typename T>
    void sytri(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *work, integer *info)
    {
        sytri(uplo, n, a, lda, ipiv, work, info);
    }
    /** @}*/ // end of sytri

    /** @defgroup syrfs syrfs
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYRFS improves the computed solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
     SYRFS improves the computed solution to a system of linear
     equations when the coefficient matrix is symmetric indefinite, and
     provides error bounds and backward error estimates for the solution.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of A contains the upper triangular part
          of the matrix A, and the strictly lower triangular part of A
          is not referenced.  If UPLO = 'L', the leading N-by-N lower
          triangular part of A contains the lower triangular part of
          the matrix A, and the strictly upper triangular part of A is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The factored form of the matrix A.  AF contains the block
          diagonal matrix D and the multipliers used to obtain the
          factor U or L from the factorization A = U*D*U**T or
          A = L*D*L**T as computed by SSYTRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by SSYTRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void syrfs(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af, integer *ldaf,
               integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr, T *work,
               integer *iwork, integer *info)
    {
        syrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void syrfs(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *af, integer *ldaf,
               integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work,
               Ta *rwork, integer *info)
    {
        syrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }
    /** @}*/ // end of syrfs

    /** @defgroup syrfsx syrfsx
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYRFSX improves the computed solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
    SYRFSX improves the computed solution to a system of linear
    equations when the coefficient matrix is symmetric indefinite, and
    provides error bounds and backward error estimates for the
    solution.  In addition to normwise error bound, the code provides
    maximum componentwise error bound if possible.  See comments for
    ERR_BNDS_NORM and ERR_BNDS_COMP for details of the error bounds.

    The original system of linear equations may have been equilibrated
    before calling this routine, as described by arguments EQUED and S
    below. In this case, the solution and error bounds returned are
    for the original unequilibrated system.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] EQUED
          EQUED is CHARACTER*1 \n
          Specifies the form of equilibration that was done to A
          before calling this routine. This is needed to compute
          the solution and error bounds correctly. \n
            = 'N':  No equilibration \n
            = 'Y':  Both row and column equilibration, i.e., A has been
                    replaced by diag(S) * A * diag(S).
                    The right hand side B has been changed accordingly. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of A contains the upper triangular
          part of the matrix A, and the strictly lower triangular
          part of A is not referenced.  If UPLO = 'L', the leading
          N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is REAL array, dimension (LDAF,N) \n
          The factored form of the matrix A.  AF contains the block
          diagonal matrix D and the multipliers used to obtain the
          factor U or L from the factorization A = U*D*U**T or A =
          L*D*L**T as computed by SSYTRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF. \n
 * @param[in,out] S
          S is REAL array, dimension (N) \n
          The scale factors for A.  If EQUED = 'Y', A is multiplied on
          the left and right by diag(S).  S is an input argument if FACT =
          'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED
          = 'Y', each element of S must be positive.  If S is output, each
          element of S is a power of the radix. If S is input, each element
          of S should be a power of the radix to ensure a reliable solution
          and error estimates. Scaling by powers of the radix does not cause
          rounding errors unless the result underflows or overflows.
          Rounding errors during scaling lead to refining with a matrix that
          is not equivalent to the input matrix, producing error estimates
          that may not be reliable. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by SGETRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          Reciprocal scaled condition number.  This is an estimate of the
          reciprocal Skeel condition number of the matrix A after
          equilibration (if done).  If this is less than the machine
          precision (in particular, if it is zero), the matrix is singular
          to working precision.  Note that the error may still be small even
          if this number is very small and the matrix appears ill-
          conditioned. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          Componentwise relative backward error.  This is the
          componentwise relative backward error of each solution vector X(j)
          (i.e., the smallest relative change in any element of A or B that
          makes X(j) an exact solution). \n
 * @param[in] N_ERR_BNDS
          N_ERR_BNDS is INTEGER \n
          Number of error bounds to return for each right hand side
          and each type (normwise or componentwise).  See ERR_BNDS_NORM and
          ERR_BNDS_COMP below. \n
 * @param[out] ERR_BNDS_NORM
          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          normwise relative error, which is defined as follows: \n
 \n
          Normwise relative error in the ith solution vector: \n

          \f[ \frac{{max\_j}\  abs(XTRUE(j,i) - X(j,i))}{{max\_j}\ abs(X(j,i))} \f] \n

 \n
          The array is indexed by the type of error information as described
          below. There currently are up to three pieces of information
          returned. \n
 \n
          The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_NORM(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated normwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*A, where S scales each row by a power of the
                  radix so all absolute row sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[out] ERR_BNDS_COMP
          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          componentwise relative error, which is defined as follows: \n
 \n
          Componentwise relative error in the ith solution vector: \n

          \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
          The array is indexed by the right-hand side i (on which the
          componentwise relative error depends), and the type of error
          information as described below. There currently are up to three
          pieces of information returned for each right-hand side. If
          componentwise accuracy is not requested (PARAMS(3) = 0.0), then
          ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
          the first (:,N_ERR_BNDS) entries are returned. \n
 \n
          The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_COMP(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon'). \n
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true. \n
 \n
          err = 3  Reciprocal condition number: Estimated componentwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*(A*diag(x)), where x is the solution for the
                  current right-hand side and S scales each row of
                  A*diag(x) by a power of the radix so all absolute row
                  sums of Z are approximately 1. \n
 \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in] NPARAMS
          NPARAMS is INTEGER \n
          Specifies the number of parameters set in PARAMS.  If <= 0, the
          PARAMS array is never referenced and default values are used. \n
 * @param[in,out] PARAMS
          PARAMS is REAL array, dimension NPARAMS \n
          Specifies algorithm parameters.  If an entry is < 0.0, then
          that entry will be filled with default value used for that
          parameter.  Only positions up to NPARAMS are accessed; defaults
          are used for higher-numbered parameters. \n
 \n
          PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
              refinement or not. \n
           Default: 1.0 \n
              = 0.0:  No refinement is performed, and no error bounds are
                      computed. \n
              = 1.0:  Use the double-precision refinement algorithm,
                      possibly with doubled-single computations if the
                      compilation environment does not support DOUBLE
                      PRECISION. \n
                (other values are reserved for future use) \n
 \n
          PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
              computations allowed for refinement. \n
           Default: 10 \n
           Aggressive: Set to 100 to permit convergence using approximate
                       factorizations or factorizations other than LU. If
                       the factorization uses a technique other than
                       Gaussian elimination, the guarantees in
                       err_bnds_norm and err_bnds_comp may no longer be
                       trustworthy. \n
 \n
          PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
              will attempt to find a solution with small componentwise
              relative error in the double-precision algorithm.  Positive
              is true, 0.0 is false. \n
           Default: 1.0 (attempt componentwise convergence) \n
 * @param[out]	WORK
          WORK is REAL array, dimension (4*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
           = 0:  Successful exit. The solution to every right-hand side is
             guaranteed. \n
           < 0:  If INFO = -i, the i-th argument had an illegal value \n
           > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
             has been completed, but the factor U is exactly singular, so
             the solution and error bounds could not be computed. RCOND = 0
             is returned. \n
           = N+J: The solution corresponding to the Jth right-hand side is
             not guaranteed. The solutions corresponding to other right-
             hand sides K with K > J may not be guaranteed as well, but
             only the first such right-hand side is reported. If a small
             componentwise error is not requested (PARAMS(3) = 0.0) then
             the Jth right-hand side is the first with a normwise error
             bound that is not guaranteed (the smallest J such
             that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
             the Jth right-hand side is the first with either a normwise or
             componentwise error bound that is not guaranteed (the smallest
             J such that either ERR_BNDS_NORM(J,1) = 0.0 or
             ERR_BNDS_COMP(J,1) = 0.0). See the definition of
             ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
             about all of the right-hand sides check ERR_BNDS_NORM or
             ERR_BNDS_COMP. \n

 *  * */
    template <typename T>
    void syrfsx(char *uplo, char *equed, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, T *s, T *b, integer *ldb, T *x, integer *ldx,
                T *rcond, T *berr, integer *n_err_bnds, T *err_bnds_norm, T *err_bnds_comp,
                integer *nparams, T *params, T *work, integer *iwork, integer *info)
    {
        syrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr,
               n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info);
    }
    template <typename T, typename Ta>
    void syrfsx(char *uplo, char *equed, integer *n, integer *nrhs, T *a, integer *lda, T *af,
                integer *ldaf, integer *ipiv, Ta *s, T *b, integer *ldb, T *x, integer *ldx,
                Ta *rcond, Ta *berr, integer *n_err_bnds, Ta *err_bnds_norm, Ta *err_bnds_comp,
                integer *nparams, Ta *params, T *work, Ta *rwork, integer *info)
    {
        syrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr,
               n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info);
    }
    /** @}*/ // end of syrfsx

    /** @defgroup syequb syequb
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYEQUB computes row and column scalings intended to equilibrate a symmetric matrix A

 * @details
 * \b Purpose:
    \verbatim
     SYEQUB computes row and column scalings intended to equilibrate a
     symmetric matrix A (with respect to the Euclidean norm) and reduce
     its condition number. The scale factors S are computed by the BIN
     algorithm (see references) so that the scaled matrix B with elements
     B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     the smallest possible condition number over all possible diagonal
     scalings.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The N-by-N symmetric matrix whose scaling factors are to be
          computed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[out] S
          S is REAL array, dimension (N) \n
          If INFO = 0, S contains the scale factors for A. \n
 * @param[out] SCOND
          SCOND is REAL \n
          If INFO = 0, S contains the ratio of the smallest S(i) to
          the largest S(i). If SCOND >= 0.1 and AMAX is neither too
          large nor too small, it is not worth scaling by S. \n
 * @param[out] AMAX
          AMAX is REAL \n
          Largest absolute value of any matrix element. If AMAX is
          very close to overflow or very close to underflow, the
          matrix should be scaled. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (2*N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the i-th diagonal element is nonpositive. \n

 *  * */
    template <typename T>
    void syequb(char *uplo, integer *n, T *a, integer *lda, T *s, T *scond, T *amax, T *work,
                integer *info)
    {
        syequb(uplo, n, a, lda, s, scond, amax, work, info);
    }
    template <typename T, typename Ta>
    void syequb(char *uplo, integer *n, T *a, integer *lda, Ta *s, Ta *scond, Ta *amax, T *work,
                integer *info)
    {
        syequb(uplo, n, a, lda, s, scond, amax, work, info);
    }
    /** @}*/ // end of syequb

    /** @defgroup syconv syconv
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYCONV convert A given by TRF into L and D and vice-versa

 * @details
 * \b Purpose:
    \verbatim
     SYCONV convert A given by TRF into L and D and vice-versa.
     Get Non-diag elements of D (returned in workspace) and
     apply or reverse permutation done in TRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] WAY
          WAY is CHARACTER*1 \n
          = 'C': Convert \n
          = 'R': Revert \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by SSYTRF. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF. \n
 * @param[out] E
          E is REAL array, dimension (N) \n
          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1
          or 2-by-2 block diagonal matrix D in LDLT. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void syconv(char *uplo, char *way, integer *n, T *a, integer *lda, integer *ipiv, T *work,
                integer *info)
    {
        syconv(uplo, way, n, a, lda, ipiv, work, info);
    }
    /** @}*/ // end of syconv

    /** @defgroup sycon_3 sycon_3
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYCON_3 estimates the reciprocal of the condition number of a real symmetric matrix A

 * @details
 * \b Purpose:
    \verbatim
     SYCON_3 estimates the reciprocal of the condition number (in the
     1-norm) of a real symmetric matrix A using the factorization
     computed by DSYTRF_RK or DSYTRF_BK:

        A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),

     where U (or L) is unit upper (or lower) triangular matrix,
     U**T (or L**T) is the transpose of U (or L), P is a permutation
     matrix, P**T is the transpose of P, and D is symmetric and block
     diagonal with 1-by-1 and 2-by-2 diagonal blocks.

     An estimate is obtained for norm(inv(A)), and the reciprocal of the
     condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
     This routine uses BLAS3 solver SSYTRS_3.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are
          stored as an upper or lower triangular matrix: \n
          = 'U':  Upper triangular, form is A = P*U*D*(U**T)*(P**T); \n
          = 'L':  Lower triangular, form is A = P*L*D*(L**T)*(P**T). \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          Diagonal of the block diagonal matrix D and factors U or L
          as computed by SSYTRF_RK and SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                should be provided on entry in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] E
          E is REAL array, dimension (N) \n
          On entry, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where
          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. \n

          NOTE: For 1-by-1 diagonal block D(k), where
          1 <= k <= N, the element E(k) is not referenced in both
          UPLO = 'U' or UPLO = 'L' cases. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF_RK or SSYTRF_BK. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (2*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sycon_3(char *uplo, integer *n, T *a, integer *lda, T *e, integer *ipiv, T *anorm,
                 T *rcond, T *work, integer *iwork, integer *info)
    {
        sycon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void sycon_3(char *uplo, integer *n, T *a, integer *lda, T *e, integer *ipiv, Ta *anorm,
                 Ta *rcond, T *work, integer *info)
    {
        sycon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info);
    }
    /** @}*/ // end of sycon_3

    /** @defgroup sytri2 sytri2
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRI2 computes the inverse of a REAL symmetric indefinite matrix \n
     A using the factorization
 * @details
 * \b Purpose:
    \verbatim
     SYTRI2 computes the inverse of a REAL symmetric indefinite matrix
     A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     SYTRF. SYTRI2 sets the LEADING DIMENSION of the workspace
     before calling SYTRI2X that actually computes the inverse.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the block diagonal matrix D and the multipliers
          used to obtain the factor U or L as computed by SSYTRF. \n
 \n
          On exit, if INFO = 0, the (symmetric) inverse of the original
          matrix.  If UPLO = 'U', the upper triangular part of the
          inverse is formed and the part of A below the diagonal is not
          referenced; if UPLO = 'L' the lower triangular part of the
          inverse is formed and the part of A above the diagonal is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SYTRF. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N+NB+1)*(NB+3) \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.
          WORK is size >= (N+NB+1)*(NB+3) \n
          If LWORK = -1, then a workspace query is assumed; the routine
           calculates: \n
              - the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, \n
              - and no error message related to LWORK is issued by XERBLA.
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
               inverse could not be computed. \n

 *  * */
    template <typename T>
    void sytri2(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *work, integer *lwork,
                integer *info)
    {
        sytri2(uplo, n, a, lda, ipiv, work, lwork, info);
    }
    /** @}*/ // end of sytri2

    /** @defgroup sytri2x sytri2x
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRI2X computes the inverse of a real symmetric indefinite matrix
     A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     SYTRF.
 * @details
 * \b Purpose:
    \verbatim
     SYTRI2X computes the inverse of a real symmetric indefinite matrix
     A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     SYTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the NNB diagonal matrix D and the multipliers
          used to obtain the factor U or L as computed by SSYTRF. \n
 \n
          On exit, if INFO = 0, the (symmetric) inverse of the original
          matrix.  If UPLO = 'U', the upper triangular part of the
          inverse is formed and the part of A below the diagonal is not
          referenced; if UPLO = 'L' the lower triangular part of the
          inverse is formed and the part of A above the diagonal is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the NNB structure of D
          as determined by SYTRF. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
               inverse could not be computed. \n

 *  * */
    template <typename T>
    void sytri2x(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, float *work,
                 integer *nb, integer *info)
    {
        sytri2x(uplo, n, a, lda, ipiv, work, nb, info);
    }
    /** @}*/ // end of sytri2x

    /** @defgroup sytri_3 sytri_3
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRI_3 computes the inverse of a real symmetric indefinite  \n
     matrix A using the factorization computed by SYTRF_RK or SYTRF_BK
 * @details
 * \b Purpose:
    \verbatim
     SYTRI_3 computes the inverse of a real symmetric indefinite
     matrix A using the factorization computed by SYTRF_RK or SYTRF_BK:

         A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),

     where U (or L) is unit upper (or lower) triangular matrix,
     U**T (or L**T) is the transpose of U (or L), P is a permutation
     matrix, P**T is the transpose of P, and D is symmetric and block
     diagonal with 1-by-1 and 2-by-2 diagonal blocks.

     SYTRI_3 sets the leading dimension of the workspace  before calling
     SYTRI_3X that actually computes the inverse.  This is the blocked
     version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are
          stored as an upper or lower triangular matrix. \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, diagonal of the block diagonal matrix D and
          factors U or L as computed by SSYTRF_RK and SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                should be provided on entry in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 \n
          On exit, if INFO = 0, the symmetric inverse of the original
          matrix. \n
             If UPLO = 'U': the upper triangular part of the inverse
             is formed and the part of A below the diagonal is not
             referenced; \n
             If UPLO = 'L': the lower triangular part of the inverse
             is formed and the part of A above the diagonal is not
             referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] E
          E is REAL array, dimension (N) \n
          On entry, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where \n
          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; \n
          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. \n
 \n
          NOTE: For 1-by-1 diagonal block D(k), where
          1 <= k <= N, the element E(k) is not referenced in both
          UPLO = 'U' or UPLO = 'L' cases. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF_RK or SSYTRF_BK. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N+NB+1)*(NB+3). \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK. LWORK >= (N+NB+1)*(NB+3). \n
 \n
          If LDWORK = -1, then a workspace query is assumed;
          the routine only calculates the optimal size of the optimal
          size of the WORK array, returns this value as the first
          entry of the WORK array, and no error message related to
          LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
               inverse could not be computed. \n

 *  * */
    template <typename T>
    void sytri_3(char *uplo, integer *n, T *a, integer *lda, T *e, integer *ipiv, T *work,
                 integer *lwork, integer *info)
    {
        sytri_3(uplo, n, a, lda, e, ipiv, work, lwork, info);
    }
    /** @}*/ // end of sytri_3

    /** @defgroup sytri_3x sytri_3x
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRI_3X computes the inverse of a real symmetric indefinite matrix A

 * @details
 * \b Purpose:
    \verbatim
    SYTRI_3X computes the inverse of a real symmetric indefinite
    matrix A using the factorization computed by SSYTRF_RK or SSYTRF_BK:

        A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),

    where U (or L) is unit upper (or lower) triangular matrix,
    U**T (or L**T) is the transpose of U (or L), P is a permutation
    matrix, P**T is the transpose of P, and D is symmetric and block
    diagonal with 1-by-1 and 2-by-2 diagonal blocks.

    This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are
          stored as an upper or lower triangular matrix. \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, diagonal of the block diagonal matrix D and
          factors U or L as computed by SYTRF_RK and SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                should be provided on entry in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A.
 \n
          On exit, if INFO = 0, the symmetric inverse of the original
          matrix. \n
             If UPLO = 'U': the upper triangular part of the inverse
             is formed and the part of A below the diagonal is not
             referenced; \n
             If UPLO = 'L': the lower triangular part of the inverse
             is formed and the part of A above the diagonal is not
             referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] E
          E is REAL array, dimension (N) \n
          On entry, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where \n
          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) not referenced; \n
          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) not referenced. \n
 \n
          NOTE: For 1-by-1 diagonal block D(k), where
          1 <= k <= N, the element E(k) is not referenced in both
          UPLO = 'U' or UPLO = 'L' cases. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF_RK or SSYTRF_BK. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N+NB+1,NB+3). \n
 * @param[in] NB
          NB is INTEGER \n
          Block size. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
               inverse could not be computed. \n

 *  * */
    template <typename T>
    void sytri_3x(char *uplo, integer *n, T *a, integer *lda, T *e, integer *ipiv, T *work,
                  integer *nb, integer *info)
    {
        sytri_3x(uplo, n, a, lda, e, ipiv, work, nb, info);
    }
    /** @}*/ // end of sytri_3x

    /** @defgroup sytrs2 sytrs2
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRS2 solves a system of linear equations A*X = B with a real \n
     symmetric matrix A
 * @details
 * \b Purpose:
    \verbatim
     SYTRS2 solves a system of linear equations A*X = B with a real
     symmetric matrix A using the factorization A = U*D*U**T or
     A = L*D*L**T computed by SSYTRF and converted by SSYCONV.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by SSYTRF.
          Note that A is input / output. This might be counter-intuitive,
          and one may think that A is input only. A is input / output. This
          is because, at the start of the subroutine, we permute A in a
          "better" form and then we permute A back to its original form at
          the end. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sytrs2(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
                integer *ldb, T *work, integer *info)
    {
        sytrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info);
    }
    /** @}*/ // end of sytrs2

    /** @defgroup sytrs_3 sytrs_3
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRS_3 solves a system of linear equations A * X = B with a real \n
     symmetric matrix A
 * @details
 * \b Purpose:
    \verbatim
     SYTRS_3 solves a system of linear equations A * X = B with a real
     symmetric matrix A using the factorization computed
     by SSYTRF_RK or SSYTRF_BK:

        A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),

     where U (or L) is unit upper (or lower) triangular matrix,
     U**T (or L**T) is the transpose of U (or L), P is a permutation
     matrix, P**T is the transpose of P, and D is symmetric and block
     diagonal with 1-by-1 and 2-by-2 diagonal blocks.

     This algorithm is using Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are
          stored as an upper or lower triangular matrix: \n
          = 'U':  Upper triangular, form is A = P*U*D*(U**T)*(P**T); \n
          = 'L':  Lower triangular, form is A = P*L*D*(L**T)*(P**T). \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          Diagonal of the block diagonal matrix D and factors U or L
          as computed by SSYTRF_RK and SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                should be provided on entry in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] E
          E is REAL array, dimension (N) \n
          On entry, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where \n
          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; \n
          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. \n
 \n
          NOTE: For 1-by-1 diagonal block D(k), where
          1 <= k <= N, the element E(k) is not referenced in both
          UPLO = 'U' or UPLO = 'L' cases. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF_RK or SSYTRF_BK. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sytrs_3(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *e, integer *ipiv,
                 T *b, integer *ldb, integer *info)
    {
        sytrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info);
    }
    /** @}*/ // end of sytrs_3

    /** @defgroup syswapr syswapr
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief  SYSWAPR applies an elementary permutation on the rows and columns of a symmetric
 matrix.

 * @details
 * \b Purpose:
    \verbatim
     SYSWAPR applies an elementary permutation on the rows and the columns of
     a symmetric matrix.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the NB diagonal matrix D and the multipliers
          used to obtain the factor U or L as computed by SSYTRF. \n
 \n
          On exit, if INFO = 0, the (symmetric) inverse of the original
          matrix.  If UPLO = 'U', the upper triangular part of the
          inverse is formed and the part of A below the diagonal is not
          referenced; if UPLO = 'L' the lower triangular part of the
          inverse is formed and the part of A above the diagonal is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] I1
          I1 is INTEGER \n
          Index of the first row to swap \n
 * @param[in] I2
          I2 is INTEGER \n
          Index of the second row to swap \n

 *  * */
    template <typename T>
    void syswapr(char *uplo, integer *n, T *a, integer *lda, integer *i1, integer *i2)
    {
        syswapr(uplo, n, a, lda, i1, i2);
    }
    /** @}*/ // end of syswapr

    /** @defgroup la_hercond_c la_hercond_c
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LA_HERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for
 Hermitian indefinite matrices

 * @details
 * \b Purpose:
    \verbatim
    LA_HERCOND_C computes the infinity norm condition number of
    op(A) * inv(diag(C)) where C is a REAL vector.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N)
          On entry, the N-by-N matrix A \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is COMPLEX array, dimension (LDAF,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by CHETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by CHETRF. \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The vector C in the formula op(A) * inv(diag(C)). \n
 * @param[in] CAPPLY
          CAPPLY is LOGICAL \n
          If .TRUE. then access the vector C in the formula above. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          i > 0:  The ith argument is invalid. \n
 * @param[out] WORK
          WORK is COMPLEX array, dimension (2*N). \n
          Workspace. \n
 * @param[out] RWORK
          RWORK is REAL array, dimension (N). \n
          Workspace. \n

 *  * */
    template <typename T, typename Ta>
    Ta la_hercond_c(char *uplo, integer *n, T *a, integer *lda, T *af, integer *ldaf, integer *ipiv,
                    Ta *c, logical *capply, integer *info, T *work, Ta *rwork)
    {
        return la_hercond_c(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork);
    }
    /** @}*/ // end of la_hercond_c

    /** @defgroup la_hercond_x la_hercond_x
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LA_HERCOND_X computes the infinity norm condition number of op(A)*diag(x) for
 Hermitian indefinite matrices

 * @details
 * \b Purpose:
    \verbatim
    LA_HERCOND_X computes the infinity norm condition number of
    op(A) * diag(X) where X is a COMPLEX vector.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is COMPLEX array, dimension (LDAF,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by CHETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by CHETRF. \n
 * @param[in] X
          X is COMPLEX array, dimension (N) \n
          The vector X in the formula op(A) * diag(X). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          i > 0:  The ith argument is invalid. \n
 * @param[out] WORK
          WORK is COMPLEX array, dimension (2*N). \n
          Workspace. \n
 * @param[out] RWORK
          RWORK is REAL array, dimension (N). \n
          Workspace.  \n

 *  * */
    template <typename T, typename Ta>
    Ta la_hercond_x(char *uplo, integer *n, T *a, integer *lda, T *af, integer *ldaf, integer *ipiv,
                    T *x, integer *info, T *work, Ta *rwork)
    {
        return la_hercond_x(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork);
    }
    /** @}*/ // end of la_hercond_x

    /** @defgroup la_herfsx_extended la_herfsx_extended
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LA_HERFSX_EXTENDED improves the computed solution to a system \n
     of linear equations for Hermitian indefinite matrices by performing \n
     extra-precise iterative refinement and provides error bounds and    \n
     backward error estimates for the solution
 * @details
 * \b Purpose:
    \verbatim
    LA_HERFSX_EXTENDED improves the computed solution to a system of
    linear equations by performing extra-precise iterative refinement
    and provides error bounds and backward error estimates for the solution.
    This subroutine is called by CHERFSX to perform iterative refinement.
    In addition to normwise error bound, the code provides maximum
    componentwise error bound if possible. See comments for ERR_BNDS_NORM
    and ERR_BNDS_COMP for details of the error bounds. Note that this
    subroutine is only resonsible for setting the second fields of
    ERR_BNDS_NORM and ERR_BNDS_COMP.
    \endverbatim

 * @param[in] PREC_TYPE
          PREC_TYPE is INTEGER \n
          Specifies the intermediate precision to be used in refinement. \n
          The value is defined by ILAPREC(P) where P is a CHARACTER and P \n
              = 'S':  Single \n
              = 'D':  Double \n
              = 'I':  Indigenous \n
              = 'X' or 'E':  Extra \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right-hand-sides, i.e., the number of columns of the
          matrix B. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is COMPLEX array, dimension (LDAF,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by CHETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by CHETRF. \n
 * @param[in] COLEQU
          COLEQU is LOGICAL \n
          If .TRUE. then column equilibration was done to A before calling
          this routine. This is needed to compute the solution and error
          bounds correctly. \n
 * @param[in] C
          C is REAL array, dimension (N) \n
          The column scale factors for A. If COLEQU = .FALSE., C
          is not accessed. If C is input, each element of C should be a power
          of the radix to ensure a reliable solution and error estimates.
          Scaling by powers of the radix does not cause rounding errors unless
          the result underflows or overflows. Rounding errors during scaling
          lead to refining with a matrix that is not equivalent to the
          input matrix, producing error estimates that may not be
          reliable. \n
 * @param[in] B
          B is COMPLEX array, dimension (LDB,NRHS) \n
          The right-hand-side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] Y
          Y is COMPLEX array, dimension (LDY,NRHS) \n
          On entry, the solution matrix X, as computed by CHETRS.
          On exit, the improved solution matrix Y. \n
 * @param[in] LDY
          LDY is INTEGER \n
          The leading dimension of the array Y.  LDY >= fla_max(1,N). \n
 * @param[out] BERR_OUT
          BERR_OUT is REAL array, dimension (NRHS) \n
          On exit, BERR_OUT(j) contains the componentwise relative backward
          error for right-hand-side j from the formula \n
             fla_max(i) (abs(RES(i)) / (abs(op(A_s))*abs(Y) + abs(B_s))(i)) \n
          where abs(Z) is the componentwise absolute value of the matrix
          or vector Z. This is computed by CLA_LIN_BERR. \n
 * @param[in] N_NORMS
          N_NORMS is INTEGER
          Determines which error bounds to   return (see ERR_BNDS_NORM
          and ERR_BNDS_COMP). \n
          If N_NORMS >= 1   return normwise error bounds. \n
          If N_NORMS >= 2   return componentwise error bounds. \n
 * @param[in,out] ERR_BNDS_NORM
          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          normwise relative error, which is defined as follows: \n
 \n
          Normwise relative error in the ith solution vector: \n

          \f[ \frac{{max\_j}\  abs(XTRUE(j,i) - X(j,i))}{{max\_j}\ abs(X(j,i))} \f] \n

 \n
          The array is indexed by the type of error information as described
          below. There currently are up to three pieces of information
          returned. \n
 \n
          The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_NORM(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon').
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true.
 \n
          err = 3  Reciprocal condition number: Estimated normwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*A, where S scales each row by a power of the
                  radix so all absolute row sums of Z are approximately 1.
 \n
          This subroutine is only responsible for setting the second field
          above. \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in,out] ERR_BNDS_COMP
          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) \n
          For each right-hand side, this array contains information about
          various error bounds and condition numbers corresponding to the
          componentwise relative error, which is defined as follows: \n
 \n
          Componentwise relative error in the ith solution vector: \n

          \f[{max\_j}\  \frac{abs(XTRUE(j,i) - X(j,i))}{abs(X(j,i))} \f] \n

 \n
          The array is indexed by the right-hand side i (on which the
          componentwise relative error depends), and the type of error
          information as described below. There currently are up to three
          pieces of information   returned for each right-hand side. If
          componentwise accuracy is not requested (PARAMS(3) = 0.0), then
          ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
          the first (:,N_ERR_BNDS) entries are   returned. \n
 \n
          The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
          right-hand side. \n
 \n
          The second index in ERR_BNDS_COMP(:,err) contains the following
          three fields: \n
          err = 1 "Trust/don't trust" boolean. Trust the answer if the
                  reciprocal condition number is less than the threshold
                  sqrt(n) * slamch('Epsilon').
 \n
          err = 2 "Guaranteed" error bound: The estimated forward error,
                  almost certainly within a factor of 10 of the true error
                  so long as the next entry is greater than the threshold
                  sqrt(n) * slamch('Epsilon'). This error bound should only
                  be trusted if the previous boolean is true.
 \n
          err = 3  Reciprocal condition number: Estimated componentwise
                  reciprocal condition number.  Compared with the threshold
                  sqrt(n) * slamch('Epsilon') to determine if the error
                  estimate is "guaranteed". These reciprocal condition
                  numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
                  appropriately scaled matrix Z.
                  Let Z = S*(A*diag(x)), where x is the solution for the
                  current right-hand side and S scales each row of
                  A*diag(x) by a power of the radix so all absolute row
                  sums of Z are approximately 1.
 \n
          This subroutine is only responsible for setting the second field
          above. \n
          See Lapack Working Note 165 for further details and extra
          cautions. \n
 * @param[in] RES
          RES is COMPLEX array, dimension (N) \n
          Workspace to hold the intermediate residual. \n
 * @param[in] AYB
          AYB is REAL array, dimension (N) \n
          Workspace. \n
 * @param[in] DY
          DY is COMPLEX array, dimension (N) \n
          Workspace to hold the intermediate solution. \n
 * @param[in] Y_TAIL
          Y_TAIL is COMPLEX array, dimension (N) \n
          Workspace to hold the trailing bits of the intermediate solution. \n
 * @param[in] RCOND
          RCOND is REAL \n
          Reciprocal scaled condition number.  This is an estimate of the
          reciprocal Skeel condition number of the matrix A after
          equilibration (if done).  If this is less than the machine
          precision (in particular, if it is zero), the matrix is singular
          to working precision.  Note that the error may still be small even
          if this number is very small and the matrix appears ill-
          conditioned. \n
 * @param[in] ITHRESH
          ITHRESH is INTEGER \n
          The maximum number of residual computations allowed for
          refinement. The default is 10. For 'aggressive' set to 100 to
          permit convergence using approximate factorizations or
          factorizations other than LU. If the factorization uses a
          technique other than Gaussian elimination, the guarantees in
          ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy. \n
 * @param[in] RTHRESH
          RTHRESH is REAL \n
          Determines when to stop refinement if the error estimate stops
          decreasing. Refinement will stop when the next solution no longer
          satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is
          the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The
          default value is 0.5. For 'aggressive' set to 0.9 to permit
          convergence on extremely ill-conditioned matrices. See LAWN 165
          for more details. \n
 * @param[in] DZ_UB
          DZ_UB is REAL \n
          Determines when to start considering componentwise convergence.
          Componentwise convergence is only considered after each component
          of the solution Y is stable, which we definte as the relative
          change in each component being less than DZ_UB. The default value
          is 0.25, requiring the first bit to be stable. See LAWN 165 for
          more details. \n
 * @param[in] IGNORE_CWISE
          IGNORE_CWISE is LOGICAL \n
          If .TRUE. then ignore componentwise convergence. Default value
          is .FALSE.. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  Successful exit. \n
          < 0:  if INFO = -i, the ith argument to CLA_HERFSX_EXTENDED had an illegal
                value \n

 * */
    template <typename T, typename Ta>
    void la_herfsx_extended(integer *prec_type, char *uplo, integer *n, integer *nrhs, T *a,
                            integer *lda, T *af, integer *ldaf, integer *ipiv, logical *colequ,
                            Ta *c, T *b, integer *ldb, T *y, integer *ldy, Ta *berr_out,
                            integer *n_norms, Ta *err_bnds_norm, Ta *err_bnds_comp, T *res, Ta *ayb,
                            T *dy, T *y_tail, Ta *rcond, integer *ithresh, Ta *rthresh, Ta *dz_ub,
                            logical *ignore_cwise, integer *info)
    {
        la_herfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b, ldb, y,
                           ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy,
                           y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info);
    }
    /** @}*/ // end of la_herfsx_extended

    /** @defgroup la_herpvgrw la_herpvgrw
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LA_HERPVGRW computes the reciprocal pivot growth factor norm(A)/norm(U)

 * @details
 * \b Purpose:
    \verbatim
    LA_HERPVGRW computes the reciprocal pivot growth factor
    norm(A)/norm(U). The "max absolute element" norm is used. If this is
    much less than 1, the stability of the LU factorization of the
    (equilibrated) matrix A could be poor. This also means that the
    solution X, estimated condition numbers, and error bounds could be
    unreliable.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] INFO
          INFO is INTEGER \n
          The value of INFO returned from SSYTRF, .i.e., the pivot in
          column INFO is exactly 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the N-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] AF
          AF is COMPLEX array, dimension (LDAF,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by CHETRF. \n
 * @param[in] LDAF
          LDAF is INTEGER \n
          The leading dimension of the array AF.  LDAF >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by CHETRF. \n
 * @param[out] WORK
          WORK is REAL array, dimension (2*N) \n

 *  * */
    template <typename T, typename Ta>
    Ta la_herpvgrw(char *uplo, integer *n, integer *info, T *a, integer *lda, T *af, integer *ldaf,
                   integer *ipiv, Ta *work)
    {
        return la_herpvgrw(uplo, n, info, a, lda, af, ldaf, ipiv, work);
    }
    /** @}*/ // end of la_herpvgrw

    /** @defgroup hpcon hpcon
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SPCON estimates the reciprocal of the condition number

 * @details
 * \b Purpose:
    \verbatim
     SPCON estimates the reciprocal of the condition number (in the
     1-norm) of a real symmetric packed matrix A using the factorization
     A = U*D*U**T or A = L*D*L**T computed by SSPTRF.

     An estimate is obtained for norm(inv(A)), and the reciprocal of the
     condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by SSPTRF, stored as a
          packed triangular matrix. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSPTRF. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (2*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void spcon(char *uplo, integer *n, T *ap, integer *ipiv, T *anorm, T *rcond, T *work,
               integer *iwork, integer *info)
    {
        spcon(uplo, n, ap, ipiv, anorm, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void spcon(char *uplo, integer *n, T *ap, integer *ipiv, Ta *anorm, Ta *rcond, T *work,
               integer *info)
    {
        spcon(uplo, n, ap, ipiv, anorm, rcond, work, info);
    }
    template <typename T, typename Ta>
    void hpcon(char *uplo, integer *n, T *ap, integer *ipiv, Ta *anorm, Ta *rcond, T *work,
               integer *info)
    {
        hpcon(uplo, n, ap, ipiv, anorm, rcond, work, info);
    }
    /** @}*/ // end of hpcon
    /** @defgroup sptrf sptrf
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SPTRF computes the factorization of a real symmetric matrix A

 * @details
 * \b Purpose:
    \verbatim
     SPTRF computes the factorization of a real symmetric matrix A stored
     in packed format using the Bunch-Kaufman diagonal pivoting method:

        A = U*D*U**T  or  A = L*D*L**T

     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and D is symmetric and block diagonal with
     1-by-1 and 2-by-2 diagonal blocks.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangle of the symmetric matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n

          On exit, the block diagonal matrix D and the multipliers used
          to obtain the factor U or L, stored as a packed triangular
          matrix overwriting A (see below for further details). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D.
          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
          interchanged and D(k,k) is a 1-by-1 diagonal block. \n
          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular, and division by zero will occur if it
               is used to solve a system of equations. \n

 *  * */
    template <typename T>
    void sptrf(char *uplo, integer *n, T *ap, integer *ipiv, integer *info)
    {
        sptrf(uplo, n, ap, ipiv, info);
    }
    template <typename T>
    void hptrf(char *uplo, integer *n, T *ap, integer *ipiv, integer *info)
    {
        hptrf(uplo, n, ap, ipiv, info);
    }
    /** @}*/ // end of sptrf

    /** @defgroup sptrs sptrs
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SPTRS solves a system of linear equations A*X = B

 * @details
 * \b Purpose:
    \verbatim
     SPTRS solves a system of linear equations A*X = B with a real
     symmetric matrix A stored in packed format using the factorization
     A = U*D*U**T or A = L*D*L**T computed by SSPTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by SSPTRF, stored as a
          packed triangular matrix. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSPTRF. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sptrs(char *uplo, integer *n, integer *nrhs, T *ap, integer *ipiv, T *b, integer *ldb,
               integer *info)
    {
        sptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info);
    }
    template <typename T>
    void hptrs(char *uplo, integer *n, integer *nrhs, T *ap, integer *ipiv, T *b, integer *ldb,
               integer *info)
    {
        hptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info);
    }
    /** @}*/ // end of sptrs

    /** @defgroup sptri sptri
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SPTRI computes the inverse of a real symmetric indefinite matrix A

 * @details
 * \b Purpose:
    \verbatim
     SPTRI computes the inverse of a real symmetric indefinite matrix
     A in packed storage using the factorization A = U*D*U**T or
     A = L*D*L**T computed by SPTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the block diagonal matrix D and the multipliers
          used to obtain the factor U or L as computed by SSPTRF,
          stored as a packed triangular matrix. \n
 \n
          On exit, if INFO = 0, the (symmetric) inverse of the original
          matrix, stored as a packed triangular matrix. The j-th column
          of inv(A) is stored in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j; \n
          if UPLO = 'L',
             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SPTRF. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
               inverse could not be computed. \n

 *  * */
    template <typename T>
    void sptri(char *uplo, integer *n, T *ap, integer *ipiv, T *work, integer *info)
    {
        sptri(uplo, n, ap, ipiv, work, info);
    }
    template <typename T>
    void hptri(char *uplo, integer *n, T *ap, integer *ipiv, T *work, integer *info)
    {
        hptri(uplo, n, ap, ipiv, work, info);
    }
    /** @}*/ // end of sptri

    /** @defgroup sprfs sprfs
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SPRFS improves the computed solution to a system of linear equations \n
     when the coefficient matrix is symmetric indefinite and packed
 * @details
 * \b Purpose:
    \verbatim
     SPRFS improves the computed solution to a system of linear
     equations when the coefficient matrix is symmetric indefinite
     and packed, and provides error bounds and backward error estimates
     for the solution.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangle of the symmetric matrix A, packed
          columnwise in a linear array.  The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[in] AFP
          AFP is REAL array, dimension (N*(N+1)/2) \n
          The factored form of the matrix A. AFP contains the block
          diagonal matrix D and the multipliers used to obtain the
          factor U or L from the factorization A = U*D*U**T or
          A = L*D*L**T as computed by SSPTRF, stored as a packed
          triangular matrix. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSPTRF. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,NRHS) \n
          On entry, the solution matrix X, as computed by SSPTRS.
          On exit, the improved solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sprfs(char *uplo, integer *n, integer *nrhs, T *ap, T *afp, integer *ipiv, T *b,
               integer *ldb, T *x, integer *ldx, T *ferr, T *berr, T *work, integer *iwork,
               integer *info)
    {
        sprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void sprfs(char *uplo, integer *n, integer *nrhs, T *ap, T *afp, integer *ipiv, T *b,
               integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork,
               integer *info)
    {
        sprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }
    template <typename T, typename Ta>
    void hprfs(char *uplo, integer *n, integer *nrhs, T *ap, T *afp, integer *ipiv, T *b,
               integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork,
               integer *info)
    {
        hprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }
    /** @}*/ // end of sprfs

    /** @defgroup hecon_rook hecon_rook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief HECON_ROOK estimates the reciprocal of the condition number fort HE matrices \n
     using factorization obtained with one of the bounded diagonal pivoting methods (max 2
 interchanges)

 * @details
 * \b Purpose:
    \verbatim
    HECON_ROOK estimates the reciprocal of the condition number of a complex
    Hermitian matrix A using the factorization A = U*D*U**H or
    A = L*D*L**H computed by CHETRF_ROOK.

    An estimate is obtained for norm(inv(A)), and the reciprocal of the
    condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**H; \n
          = 'L':  Lower triangular, form is A = L*D*L**H. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by CHETRF_ROOK. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by CHETRF_ROOK. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out] WORK
          WORK is COMPLEX array, dimension (2*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T, typename Ta>
    void hecon_rook(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, Ta *anorm, Ta *rcond,
                    T *work, integer *info)
    {
        hecon_rook(uplo, n, a, lda, ipiv, anorm, rcond, work, info);
    }
    /** @}*/ // end of hecon_rook hecon_rook

    /** @defgroup  sycon_rook sycon_rook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYCON_ROOK estimates the reciprocal of the condition number of a real symmetric
 matrix A

 * @details
 * \b Purpose:
    \verbatim
    SYCON_ROOK estimates the reciprocal of the condition number (in the
    1-norm) of a real symmetric matrix A using the factorization
    A = U*D*U**T or A = L*D*L**T computed by SSYTRF_ROOK.

    An estimate is obtained for norm(inv(A)), and the reciprocal of the
    condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The block diagonal matrix D and the multipliers used to
          obtain the factor U or L as computed by SSYTRF_ROOK. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF_ROOK. \n
 * @param[in] ANORM
          ANORM is REAL \n
          The 1-norm of the original matrix A. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
          estimate of the 1-norm of inv(A) computed in this routine. \n
 * @param[out] WORK
          WORK is REAL array, dimension (2*N) \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sycon_rook(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *anorm, T *rcond,
                    T *work, integer *iwork, integer *info)
    {
        sycon_rook(*uplo, *n, a, *lda, *ipiv, anorm, rcond, work, *iwork, *info);
    }
    template <typename T, typename Ta>
    void sycon_rook(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, Ta *anorm, Ta *rcond,
                    T *work, integer *iwork, integer *info)
    {
        sycon_rook(*uplo, *n, a, *lda, *ipiv, anorm, rcond, work, *iwork, *info);
    }

    /** @}*/ // end of sycon_rook

    /** @defgroup sytrf_hook sytrf_hook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRF_ROOK computes the factorization of a real symmetric matrix A \n
     using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
 * @details
 * \b Purpose:
    \verbatim
     SYTRF_ROOK computes the factorization of a real symmetric matrix A
     using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
     The form of the factorization is

        A = U*D*U**T  or  A = L*D*L**T

     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and D is symmetric and block diagonal with
     1-by-1 and 2-by-2 diagonal blocks.

     This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, the block diagonal matrix D and the multipliers used
          to obtain the factor U or L (see below for further details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D. \n
 \n
          If UPLO = 'U': \n
               If IPIV(k) > 0, then rows and columns k and IPIV(k)
               were interchanged and D(k,k) is a 1-by-1 diagonal block. \n
 \n
               If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
               columns k and -IPIV(k) were interchanged and rows and
               columns k-1 and -IPIV(k-1) were inerchaged,
               D(k-1:k,k-1:k) is a 2-by-2 diagonal block. \n
 \n
          If UPLO = 'L': \n
               If IPIV(k) > 0, then rows and columns k and IPIV(k)
               were interchanged and D(k,k) is a 1-by-1 diagonal block. \n
 \n
               If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
               columns k and -IPIV(k) were interchanged and rows and
               columns k+1 and -IPIV(k+1) were inerchaged,
               D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)). \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >=1.  For best performance
          LWORK >= N*NB, where NB is the block size returned by ILAENV. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
                has been completed, but the block diagonal matrix D is
                exactly singular, and division by zero will occur if it
                is used to solve a system of equations. \n

 *  * */
    template <typename T>
    void sytrf_rook(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *work,
                    integer *lwork, integer *info)
    {
        sytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
    }
    /** @}*/ // end of sytrf_hook

    /** @defgroup lasyf_rook lasyf_rook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LASYF_ROOK computes a partial factorization of a real symmetric    \n
     matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method
 * @details
 * \b Purpose:
    \verbatim
    LASYF_ROOK computes a partial factorization of a real symmetric
    matrix A using the bounded Bunch-Kaufman ("rook") diagonal
    pivoting method. The partial factorization has the form:

    A  =  (I  U12) (A11  0 ) ( I       0   )  if UPLO = 'U', or:
          (0  U22) ( 0   D ) (U12**T U22**T)

    A  =  (L11  0) ( D   0 ) (L11**T L21**T)  if UPLO = 'L'
          (L21  I) ( 0  A22) ( 0       I   )

    where the order of D is at most NB. The actual order is   returned in
    the argument KB, and is either NB or NB-1, or N if N <= NB.

    SLASYF_ROOK is an auxiliary routine called by SSYTRF_ROOK. It uses
    blocked code (calling Level 3 BLAS) to update the submatrix
    A11 (if UPLO = 'U') or A22 (if UPLO = 'L').
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NB
          NB is INTEGER \n
          The maximum number of columns of the matrix A that should be
          factored.  NB should be at least 2 to allow for 2-by-2 pivot
          blocks. \n
 * @param[out] KB
          KB is INTEGER \n
          The number of columns of A that were actually factored.
          KB is either NB-1 or NB, or N if N <= NB. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          n-by-n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n-by-n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
          On exit, A contains details of the partial factorization. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D.
 \n
          If UPLO = 'U': \n
             Only the last KB elements of IPIV are set.
 \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
             interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
             columns k and -IPIV(k) were interchanged and rows and
             columns k-1 and -IPIV(k-1) were inerchaged,
             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
 \n
          If UPLO = 'L': \n
             Only the first KB elements of IPIV are set.
 \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k)
             were interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
             columns k and -IPIV(k) were interchanged and rows and
             columns k+1 and -IPIV(k+1) were inerchaged,
             D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 * @param[out] W
          W is REAL array, dimension (LDW,NB) \n
 * @param[in] LDW
          LDW is INTEGER \n
          The leading dimension of the array W.  LDW >= fla_max(1,N). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular.  \n

 *  * */
    template <typename T>
    void lasyf_rook(char *uplo, integer *n, integer *nb, integer *kb, T *a, integer *lda,
                    integer *ipiv, T *w, integer *ldw, integer *info)
    {
        lasyf_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info);
    }
    /** @}*/ // end of lasyf_rook

    /** @defgroup sytf2_rook sytf2_rook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTF2_ROOK computes the factorization of a real symmetric indefinite \n
     matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method
 * @details
 * \b Purpose:
    \verbatim
    SYTF2_ROOK computes the factorization of a real symmetric matrix A
    using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:

       A = U*D*U**T  or  A = L*D*L**T

    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, U**T is the transpose of U, and D is symmetric and
    block diagonal with 1-by-1 and 2-by-2 diagonal blocks.

    This is the unblocked version of the algorithm, calling Level 2 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          n-by-n upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading n-by-n lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, the block diagonal matrix D and the multipliers used
          to obtain the factor U or L (see below for further details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D. \n
 \n
          If UPLO = 'U': \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k)
             were interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
             columns k and -IPIV(k) were interchanged and rows and
             columns k-1 and -IPIV(k-1) were inerchaged,
             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
 \n
          If UPLO = 'L': \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k)
             were interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
             columns k and -IPIV(k) were interchanged and rows and
             columns k+1 and -IPIV(k+1) were inerchaged,
             D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -k, the k-th argument had an illegal value \n
          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
               has been completed, but the block diagonal matrix D is
               exactly singular, and division by zero will occur if it
               is used to solve a system of equations. \n

 *  * */
    template <typename T>
    void sytf2_rook(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, integer *info)
    {
        sytf2_rook(uplo, n, a, lda, ipiv, info);
    }
    /** @}*/ // end of sytf2_rook

    /** @defgroup sytri_rook sytri_rook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRI_ROOK computes the inverse of a real symmetric matrix A

 * @details
 * \b Purpose:
    \verbatim
    SYTRI_ROOK computes the inverse of a real symmetric
    matrix A using the factorization A = U*D*U**T or A = L*D*L**T
    computed by SSYTRF_ROOK.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U*D*U**T; \n
          = 'L':  Lower triangular, form is A = L*D*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the block diagonal matrix D and the multipliers
          used to obtain the factor U or L as computed by SSYTRF_ROOK.
 \n
          On exit, if INFO = 0, the (symmetric) inverse of the original
          matrix.  If UPLO = 'U', the upper triangular part of the
          inverse is formed and the part of A below the diagonal is not
          referenced; if UPLO = 'L' the lower triangular part of the
          inverse is formed and the part of A above the diagonal is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D
          as determined by SSYTRF_ROOK. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
               inverse could not be computed. \n

 *  * */
    template <typename T>
    void sytri_rook(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *work,
                    integer *info)
    {
        sytri_rook(uplo, n, a, lda, ipiv, work, info);
    }
    /** @}*/ // end of sytri_rook

    /** @defgroup sytrf_rk sytrf_rk
     * @ingroup LDL_computation LDL_computation
     * @{
     */

    /*! @brief SYTRF_RK computes the factorization of a real symmetric indefinite matrix \n
         using the bounded Bunch-Kaufman (rook) diagonal pivoting method (BLAS3 blocked algorithm).
     * @details
     * \b Purpose:
        \verbatim
         SYTRF_RK computes the factorization of a real symmetric matrix A
         using the bounded Bunch-Kaufman (rook) diagonal pivoting method:

            A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),

         where U (or L) is unit upper (or lower) triangular matrix,
         U**T (or L**T) is the transpose of U (or L), P is a permutation
         matrix, P**T is the transpose of P, and D is symmetric and block
         diagonal with 1-by-1 and 2-by-2 diagonal blocks.

         This is the blocked version of the algorithm, calling Level 3 BLAS.
         For more information see Further Details section.
        \endverbatim

     * @param[in] UPLO
              UPLO is CHARACTER*1 \n
              Specifies whether the upper or lower triangular part of the
              symmetric matrix A is stored: \n
              = 'U':  Upper triangular \n
              = 'L':  Lower triangular \n
     * @param[in] N
              N is INTEGER \n
              The order of the matrix A.  N >= 0. \n
     * @param[in,out] A
              A is REAL array, dimension (LDA,N) \n
              On entry, the symmetric matrix A. \n
                If UPLO = 'U': the leading N-by-N upper triangular part
                of A contains the upper triangular part of the matrix A,
                and the strictly lower triangular part of A is not
                referenced. \n
     \n
                If UPLO = 'L': the leading N-by-N lower triangular part
                of A contains the lower triangular part of the matrix A,
                and the strictly upper triangular part of A is not
                referenced. \n
     \n
              On exit, contains: \n
                a) ONLY diagonal elements of the symmetric block diagonal
                   matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
                   (superdiagonal (or subdiagonal) elements of D
                    are stored on exit in array E), and \n
                b) If UPLO = 'U': factor U in the superdiagonal part of A.
                   If UPLO = 'L': factor L in the subdiagonal part of A. \n
     * @param[in] LDA
              LDA is INTEGER \n
              The leading dimension of the array A.  LDA >= fla_max(1,N). \n
     * @param[out] E
              E is REAL array, dimension (N) \n
              On exit, contains the superdiagonal (or subdiagonal)
              elements of the symmetric block diagonal matrix D
              with 1-by-1 or 2-by-2 diagonal blocks, where \n
              If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; \n
              If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. \n
     \n
              NOTE: For 1-by-1 diagonal block D(k), where
              1 <= k <= N, the element E(k) is set to 0 in both
              UPLO = 'U' or UPLO = 'L' cases. \n
     * @param[out] IPIV
              IPIV is INTEGER array, dimension (N) \n
              IPIV describes the permutation matrix P in the factorization
              of matrix A as follows. The absolute value of IPIV(k)
              represents the index of row and column that were
              interchanged with the k-th row and column. The value of UPLO
              describes the order in which the interchanges were applied.
              Also, the sign of IPIV represents the block structure of
              the symmetric block diagonal matrix D with 1-by-1 or 2-by-2
              diagonal blocks which correspond to 1 or 2 interchanges
              at each factorization step. For more info see Further
              Details section. \n
     \n
              If UPLO = 'U',
              ( in factorization order, k decreases from N to 1): \n
                a) A single positive entry IPIV(k) > 0 means:
                   D(k,k) is a 1-by-1 diagonal block.
                   If IPIV(k) != k, rows and columns k and IPIV(k) were
                   interchanged in the matrix A(1:N,1:N);
                   If IPIV(k) = k, no interchange occurred. \n
     \n
                b) A pair of consecutive negative entries
                   IPIV(k) < 0 and IPIV(k-1) < 0 means:
                   D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
                   (NOTE: negative entries in IPIV appear ONLY in pairs). \n
                   1) If -IPIV(k) != k, rows and columns
                      k and -IPIV(k) were interchanged
                      in the matrix A(1:N,1:N).
                      If -IPIV(k) = k, no interchange occurred. \n
                   2) If -IPIV(k-1) != k-1, rows and columns
                      k-1 and -IPIV(k-1) were interchanged
                      in the matrix A(1:N,1:N).
                      If -IPIV(k-1) = k-1, no interchange occurred. \n
     \n
                c) In both cases a) and b), always ABS( IPIV(k)) <= k. \n
     \n
                d) NOTE: Any entry IPIV(k) is always NONZERO on output. \n
     \n
              If UPLO = 'L',
              ( in factorization order, k increases from 1 to N): \n
                a) A single positive entry IPIV(k) > 0 means:
                   D(k,k) is a 1-by-1 diagonal block.
                   If IPIV(k) != k, rows and columns k and IPIV(k) were
                   interchanged in the matrix A(1:N,1:N).
                   If IPIV(k) = k, no interchange occurred. \n
     \n
                b) A pair of consecutive negative entries
                   IPIV(k) < 0 and IPIV(k+1) < 0 means:
                   D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
                   (NOTE: negative entries in IPIV appear ONLY in pairs). \n
                   1) If -IPIV(k) != k, rows and columns
                      k and -IPIV(k) were interchanged
                      in the matrix A(1:N,1:N).
                      If -IPIV(k) = k, no interchange occurred. \n
                   2) If -IPIV(k+1) != k+1, rows and columns
                      k-1 and -IPIV(k-1) were interchanged
                      in the matrix A(1:N,1:N).
                      If -IPIV(k+1) = k+1, no interchange occurred. \n
     \n
                c) In both cases a) and b), always ABS( IPIV(k)) >= k. \n
     \n
                d) NOTE: Any entry IPIV(k) is always NONZERO on output. \n
     * @param[out]	WORK
              WORK is REAL array, dimension ( MAX(1,LWORK) ). \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
     * @param[in]	LWORK
              LWORK is INTEGER \n
              The length of WORK.  LWORK >=1.  For best performance
              LWORK >= N*NB, where NB is the block size returned
              by ILAENV. \n
     \n
              If LWORK = -1, then a workspace query is assumed;
              the routine only calculates the optimal size of the WORK
              array, returns this value as the first entry of the WORK
              array, and no error message related to LWORK is issued
              by XERBLA. \n
     * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: If INFO = -k, the k-th argument had an illegal value \n
              > 0: If INFO = k, the matrix A is singular, because:
                     If UPLO = 'U': column k in the upper
                     triangular part of A contains all zeros.
                     If UPLO = 'L': column k in the lower
                     triangular part of A contains all zeros. \n
     \n
                   Therefore D(k,k) is exactly zero, and superdiagonal
                   elements of column k of U (or subdiagonal elements of
                   column k of L ) are all zeros. The factorization has
                   been completed, but the block diagonal matrix D is
                   exactly singular, and division by zero will occur if
                   it is used to solve a system of equations. \n
     \n
                   NOTE: INFO only stores the first occurrence of
                   a singularity, any subsequent occurrence of singularity
                   is not stored in INFO even though the factorization
                   always completes. \n

     *  * */
    template <typename T>
    void sytrf_rk(char *uplo, integer *n, T *a, integer *lda, T *e, integer *ipiv, T *work,
                  integer *lwork, integer *info)
    {
        sytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info);
    }
    /** @}*/ // end of sytrf_rk

    /** @defgroup lasyf_rk lasyf_rk
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LASYF_RK computes a partial factorization of a real symmetric indefinite \n
     matrix using bounded Bunch-Kaufman (rook) diagonal pivoting method
 * @details
 * \b Purpose:
    \verbatim
    LASYF_RK computes a partial factorization of a real symmetric
    matrix A using the bounded Bunch-Kaufman (rook) diagonal
    pivoting method. The partial factorization has the form:

    A  =  (I  U12) (A11  0 ) ( I       0   )  if UPLO = 'U', or:
          (0  U22) ( 0   D ) (U12**T U22**T)

    A  =  (L11  0) ( D   0 ) (L11**T L21**T)  if UPLO = 'L',
          (L21  I) ( 0  A22) ( 0       I   )

    where the order of D is at most NB. The actual order is   returned in
    the argument KB, and is either NB or NB-1, or N if N <= NB.

    LASYF_RK is an auxiliary routine called by SSYTRF_RK. It uses
    blocked code (calling Level 3 BLAS) to update the submatrix
    A11 (if UPLO = 'U') or A22 (if UPLO = 'L').
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NB
          NB is INTEGER \n
          The maximum number of columns of the matrix A that should be
          factored.  NB should be at least 2 to allow for 2-by-2 pivot
          blocks. \n
 * @param[out] KB
          KB is INTEGER \n
          The number of columns of A that were actually factored.
          KB is either NB-1 or NB, or N if N <= NB. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A. \n
            If UPLO = 'U': the leading N-by-N upper triangular part
            of A contains the upper triangular part of the matrix A,
            and the strictly lower triangular part of A is not
            referenced.
 \n
            If UPLO = 'L': the leading N-by-N lower triangular part
            of A contains the lower triangular part of the matrix A,
            and the strictly upper triangular part of A is not
            referenced.
 \n
          On exit, contains: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                are stored on exit in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] E
          E is REAL array, dimension (N) \n
          On exit, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where \n
          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0; \n
          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0. \n
 \n
          NOTE: For 1-by-1 diagonal block D(k), where
          1 <= k <= N, the element E(k) is set to 0 in both
          UPLO = 'U' or UPLO = 'L' cases. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          IPIV describes the permutation matrix P in the factorization
          of matrix A as follows. The absolute value of IPIV(k)
          represents the index of row and column that were
          interchanged with the k-th row and column. The value of UPLO
          describes the order in which the interchanges were applied.
          Also, the sign of IPIV represents the block structure of
          the symmetric block diagonal matrix D with 1-by-1 or 2-by-2
          diagonal blocks which correspond to 1 or 2 interchanges
          at each factorization step. \n
 \n
          If UPLO = 'U',
          (in factorization order, k decreases from N to 1): \n
            a) A single positive entry IPIV(k) > 0 means:
               D(k,k) is a 1-by-1 diagonal block.
               If IPIV(k) != k, rows and columns k and IPIV(k) were
               interchanged in the submatrix A(1:N,N-KB+1:N);
               If IPIV(k) = k, no interchange occurred. \n
 \n
            b) A pair of consecutive negative entries
               IPIV(k) < 0 and IPIV(k-1) < 0 means:
               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
               (NOTE: negative entries in IPIV appear ONLY in pairs). \n
               1) If -IPIV(k) != k, rows and columns
                  k and -IPIV(k) were interchanged
                  in the matrix A(1:N,N-KB+1:N).
                  If -IPIV(k) = k, no interchange occurred. \n
               2) If -IPIV(k-1) != k-1, rows and columns
                  k-1 and -IPIV(k-1) were interchanged
                  in the submatrix A(1:N,N-KB+1:N).
                  If -IPIV(k-1) = k-1, no interchange occurred. \n
 \n
            c) In both cases a) and b) is always ABS(IPIV(k)) <= k. \n
 \n
            d) NOTE: Any entry IPIV(k) is always NONZERO on output. \n
 \n
          If UPLO = 'L',
          (in factorization order, k increases from 1 to N): \n
            a) A single positive entry IPIV(k) > 0 means:
               D(k,k) is a 1-by-1 diagonal block.
               If IPIV(k) != k, rows and columns k and IPIV(k) were
               interchanged in the submatrix A(1:N,1:KB).
               If IPIV(k) = k, no interchange occurred. \n
 \n
            b) A pair of consecutive negative entries
               IPIV(k) < 0 and IPIV(k+1) < 0 means:
               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
               (NOTE: negative entries in IPIV appear ONLY in pairs). \n
               1) If -IPIV(k) != k, rows and columns
                  k and -IPIV(k) were interchanged
                  in the submatrix A(1:N,1:KB).
                  If -IPIV(k) = k, no interchange occurred. \n
               2) If -IPIV(k+1) != k+1, rows and columns
                  k-1 and -IPIV(k-1) were interchanged
                  in the submatrix A(1:N,1:KB).
                  If -IPIV(k+1) = k+1, no interchange occurred. \n
 \n
            c) In both cases a) and b) is always ABS(IPIV(k)) >= k. \n
 \n
            d) NOTE: Any entry IPIV(k) is always NONZERO on output. \n
 * @param[out] W
          W is REAL array, dimension (LDW,NB) \n
 * @param[in] LDW
          LDW is INTEGER \n
          The leading dimension of the array W.  LDW >= fla_max(1,N). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: If INFO = -k, the k-th argument had an illegal value \n
          > 0: If INFO = k, the matrix A is singular, because: \n
                 If UPLO = 'U': column k in the upper
                 triangular part of A contains all zeros. \n
                 If UPLO = 'L': column k in the lower
                 triangular part of A contains all zeros. \n
 \n
               Therefore D(k,k) is exactly zero, and superdiagonal
               elements of column k of U (or subdiagonal elements of
               column k of L) are all zeros. The factorization has
               been completed, but the block diagonal matrix D is
               exactly singular, and division by zero will occur if
               it is used to solve a system of equations.
 \n
               NOTE: INFO only stores the first occurrence of
               a singularity, any subsequent occurrence of singularity
               is not stored in INFO even though the factorization
               always completes.  \n

 *  * */
    template <typename T>
    void lasyf_rk(char *uplo, integer *n, integer *nb, integer *kb, T *a, integer *lda, T *e,
                  integer *ipiv, T *w, integer *ldw, integer *info)
    {
        lasyf_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info);
    }
    /** @}*/ // end of lasyf_rk

    /** @defgroup sytf2_rk sytf2_rk
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTF2_RK computes the factorization of a real symmetric indefinite \n
     matrix using the bounded Bunch-Kaufman (rook) diagonal pivoting method

 * @details
 * \b Purpose:
    \verbatim
    SYTF2_RK computes the factorization of a real symmetric matrix A
    using the bounded Bunch-Kaufman (rook) diagonal pivoting method:

       A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),

    where U (or L) is unit upper (or lower) triangular matrix,
    U**T (or L**T) is the transpose of U (or L), P is a permutation
    matrix, P**T is the transpose of P, and D is symmetric and block
    diagonal with 1-by-1 and 2-by-2 diagonal blocks.

    This is the unblocked version of the algorithm, calling Level 2 BLAS.
    For more information see Further Details section.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is stored: \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A. \n
            If UPLO = 'U': the leading N-by-N upper triangular part
            of A contains the upper triangular part of the matrix A,
            and the strictly lower triangular part of A is not
            referenced.
 \n
            If UPLO = 'L': the leading N-by-N lower triangular part
            of A contains the lower triangular part of the matrix A,
            and the strictly upper triangular part of A is not
            referenced.
 \n
          On exit, contains: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                are stored on exit in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A. \n
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] E
          E is REAL array, dimension (N) \n
          On exit, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where
          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
 \n
          NOTE: For 1-by-1 diagonal block D(k), where
          1 <= k <= N, the element E(k) is set to 0 in both
          UPLO = 'U' or UPLO = 'L' cases. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          IPIV describes the permutation matrix P in the factorization
          of matrix A as follows. The absolute value of IPIV(k)
          represents the index of row and column that were
          interchanged with the k-th row and column. The value of UPLO
          describes the order in which the interchanges were applied.
          Also, the sign of IPIV represents the block structure of
          the symmetric block diagonal matrix D with 1-by-1 or 2-by-2
          diagonal blocks which correspond to 1 or 2 interchanges
          at each factorization step. For more info see Further
          Details section. \n
 \n
          If UPLO = 'U',
          (in factorization order, k decreases from N to 1): \n
            a) A single positive entry IPIV(k) > 0 means:
               D(k,k) is a 1-by-1 diagonal block.
               If IPIV(k) != k, rows and columns k and IPIV(k) were
               interchanged in the matrix A(1:N,1:N);
               If IPIV(k) = k, no interchange occurred. \n
 \n
            b) A pair of consecutive negative entries
               IPIV(k) < 0 and IPIV(k-1) < 0 means:
               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
               (NOTE: negative entries in IPIV appear ONLY in pairs). \n
               1) If -IPIV(k) != k, rows and columns
                  k and -IPIV(k) were interchanged
                  in the matrix A(1:N,1:N).
                  If -IPIV(k) = k, no interchange occurred. \n
               2) If -IPIV(k-1) != k-1, rows and columns
                  k-1 and -IPIV(k-1) were interchanged
                  in the matrix A(1:N,1:N).
                  If -IPIV(k-1) = k-1, no interchange occurred. \n
 \n
            c) In both cases a) and b), always ABS(IPIV(k)) <= k. \n
 \n
            d) NOTE: Any entry IPIV(k) is always NONZERO on output. \n
 \n
          If UPLO = 'L',
          (in factorization order, k increases from 1 to N): \n
            a) A single positive entry IPIV(k) > 0 means:
               D(k,k) is a 1-by-1 diagonal block.
               If IPIV(k) != k, rows and columns k and IPIV(k) were
               interchanged in the matrix A(1:N,1:N).
               If IPIV(k) = k, no interchange occurred.
 \n
            b) A pair of consecutive negative entries
               IPIV(k) < 0 and IPIV(k+1) < 0 means:
               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
               (NOTE: negative entries in IPIV appear ONLY in pairs). \n
               1) If -IPIV(k) != k, rows and columns
                  k and -IPIV(k) were interchanged
                  in the matrix A(1:N,1:N).
                  If -IPIV(k) = k, no interchange occurred. \n
               2) If -IPIV(k+1) != k+1, rows and columns
                  k-1 and -IPIV(k-1) were interchanged
                  in the matrix A(1:N,1:N).
                  If -IPIV(k+1) = k+1, no interchange occurred. \n
 \n
            c) In both cases a) and b), always ABS(IPIV(k)) >= k. \n
 \n
            d) NOTE: Any entry IPIV(k) is always NONZERO on output. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: If INFO = -k, the k-th argument had an illegal value \n
          > 0: If INFO = k, the matrix A is singular, because: \n
                 If UPLO = 'U': column k in the upper
                 triangular part of A contains all zeros. \n
                 If UPLO = 'L': column k in the lower
                 triangular part of A contains all zeros. \n
 \n
               Therefore D(k,k) is exactly zero, and superdiagonal
               elements of column k of U (or subdiagonal elements of
               column k of L) are all zeros. The factorization has
               been completed, but the block diagonal matrix D is
               exactly singular, and division by zero will occur if
               it is used to solve a system of equations.
 \n
               NOTE: INFO only stores the first occurrence of
               a singularity, any subsequent occurrence of singularity
               is not stored in INFO even though the factorization
               always completes. \n

 *  * */
    template <typename T>
    void sytf2_rk(char *uplo, integer *n, T *a, integer *lda, T *e, integer *ipiv, integer *info)
    {
        sytf2_rk(uplo, n, a, lda, e, ipiv, info);
    }
    /** @}*/ // end of systf2_rk

    /** @defgroup syconvf syconvf
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYCONVF converts the factorization output format

 * @details
 * \b Purpose:
    \verbatim
    If parameter WAY = 'C':
    SYCONVF converts the factorization output format used in
    SYTRF provided on entry in parameter A into the factorization
    output format used in SYTRF_RK (or SYTRF_BK) that is stored
    on exit in parameters A and E. It also coverts in place details of
    the intechanges stored in IPIV from the format used in SYTRF into
    the format used in SYTRF_RK (or SYTRF_BK).

    If parameter WAY = 'R':
    SYCONVF performs the conversion in reverse direction, i.e.
    converts the factorization output format used in SYTRF_RK
    (or SYTRF_BK) provided on entry in parameters A and E into
    the factorization output format used in SSYTRF that is stored
    on exit in parameter A. It also coverts in place details of
    the intechanges stored in IPIV from the format used in SYTRF_RK
    (or SYTRF_BK) into the format used in SYTRF.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are
          stored as an upper or lower triangular matrix A. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] WAY
          WAY is CHARACTER*1 \n
          = 'C': Convert \n
          = 'R': Revert \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
 \n
          1) If WAY ='C':
 \n
          On entry, contains factorization details in format used in
          SSYTRF: \n
            a) all elements of the symmetric block diagonal
               matrix D on the diagonal of A and on superdiagonal
               (or subdiagonal) of A, and \n
            b) If UPLO = 'U': multipliers used to obtain factor U
               in the superdiagonal part of A.
               If UPLO = 'L': multipliers used to obtain factor L
               in the superdiagonal part of A. \n
 \n
          On exit, contains factorization details in format used in
          SSYTRF_RK or SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                are stored on exit in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A. \n
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 \n
          2) If WAY = 'R':
 \n
          On entry, contains factorization details in format used in
          SSYTRF_RK or SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                are stored on exit in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A. \n
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 \n
          On exit, contains factorization details in format used in
          SSYTRF: \n
            a) all elements of the symmetric block diagonal
               matrix D on the diagonal of A and on superdiagonal
               (or subdiagonal) of A, and \n
            b) If UPLO = 'U': multipliers used to obtain factor U
               in the superdiagonal part of A.
               If UPLO = 'L': multipliers used to obtain factor L
               in the superdiagonal part of A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] E
          E is REAL array, dimension (N) \n
 \n
          1) If WAY ='C':
 \n
          On entry, just a workspace.
 \n
          On exit, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where
          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
 \n
          2) If WAY = 'R':
 \n
          On entry, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where
          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
 \n
          On exit, is not changed \n
 * @param[in,out] IPIV
          IPIV is INTEGER array, dimension (N) \n
 \n
          1) If WAY ='C': \n
          On entry, details of the interchanges and the block
          structure of D in the format used in SSYTRF. \n
          On exit, details of the interchanges and the block
          structure of D in the format used in SSYTRF_RK
          (or SSYTRF_BK). \n
 \n
          1) If WAY ='R': \n
          On entry, details of the interchanges and the block
          structure of D in the format used in SSYTRF_RK
          (or SSYTRF_BK). \n
          On exit, details of the interchanges and the block
          structure of D in the format used in SSYTRF. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void syconvf(char *uplo, char *way, integer *n, T *a, integer *lda, T *e, integer *ipiv,
                 integer *info)
    {
        syconvf(uplo, way, n, a, lda, e, ipiv, info);
    }
    /** @}*/ // end of syconvf

    /** @defgroup syconvf_rook syconvf_rook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYCONVF_ROOK converts the factorization output format

 * @details
 * \b Purpose:
    \verbatim
    If parameter WAY = 'C':
    SSYCONVF_ROOK converts the factorization output format used in
    SSYTRF_ROOK provided on entry in parameter A into the factorization
    output format used in SSYTRF_RK (or SSYTRF_BK) that is stored
    on exit in parameters A and E. IPIV format for SSYTRF_ROOK and
    SSYTRF_RK (or SSYTRF_BK) is the same and is not converted.

    If parameter WAY = 'R':
    SSYCONVF_ROOK performs the conversion in reverse direction, i.e.
    converts the factorization output format used in SSYTRF_RK
    (or SSYTRF_BK) provided on entry in parameters A and E into
    the factorization output format used in SSYTRF_ROOK that is stored
    on exit in parameter A. IPIV format for SSYTRF_ROOK and
    SSYTRF_RK (or SSYTRF_BK) is the same and is not converted.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are
          stored as an upper or lower triangular matrix A. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] WAY
          WAY is CHARACTER*1 \n
          = 'C': Convert \n
          = 'R': Revert \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
 \n
          1) If WAY ='C':
 \n
          On entry, contains factorization details in format used in
          SSYTRF_ROOK: \n
            a) all elements of the symmetric block diagonal
               matrix D on the diagonal of A and on superdiagonal
               (or subdiagonal) of A, and \n
            b) If UPLO = 'U': multipliers used to obtain factor U
               in the superdiagonal part of A.
               If UPLO = 'L': multipliers used to obtain factor L
               in the superdiagonal part of A. \n
 \n
          On exit, contains factorization details in format used in
          SSYTRF_RK or SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                are stored on exit in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A. \n
 \n
          2) If WAY = 'R':
 \n
          On entry, contains factorization details in format used in
          SSYTRF_RK or SSYTRF_BK: \n
            a) ONLY diagonal elements of the symmetric block diagonal
               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
               (superdiagonal (or subdiagonal) elements of D
                are stored on exit in array E), and \n
            b) If UPLO = 'U': factor U in the superdiagonal part of A.
               If UPLO = 'L': factor L in the subdiagonal part of A.
 \n
          On exit, contains factorization details in format used in
          SSYTRF_ROOK: \n
            a) all elements of the symmetric block diagonal
               matrix D on the diagonal of A and on superdiagonal
               (or subdiagonal) of A, and \n
            b) If UPLO = 'U': multipliers used to obtain factor U
               in the superdiagonal part of A.
               If UPLO = 'L': multipliers used to obtain factor L
               in the superdiagonal part of A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] E
          E is REAL array, dimension (N) \n
 \n
          1) If WAY ='C':
 \n
          On entry, just a workspace.
 \n
          On exit, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where
          If UPLO = 'U': E(i) = D(i-1,i), i=2:N, E(1) is set to 0;
          If UPLO = 'L': E(i) = D(i+1,i), i=1:N-1, E(N) is set to 0.
 \n
          2) If WAY = 'R':
 \n
          On entry, contains the superdiagonal (or subdiagonal)
          elements of the symmetric block diagonal matrix D
          with 1-by-1 or 2-by-2 diagonal blocks, where \n
          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced; \n
          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced. \n
 \n
          On exit, is not changed \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          On entry, details of the interchanges and the block
          structure of D as determined: \n
          1) by SSYTRF_ROOK, if WAY ='C'; \n
          2) by SSYTRF_RK (or SSYTRF_BK), if WAY ='R'.
          The IPIV format is the same for all these routines. \n
 \n
          On exit, is not changed. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void syconvf_rook(char *uplo, char *way, integer *n, T *a, integer *lda, T *e, integer *ipiv,
                      integer *info)
    {
        syconvf_rook(uplo, way, n, a, lda, e, ipiv, info);
    }
    /** @}*/ // end of sycovf_rook

    /** @defgroup sytrf_aa sytrf_aa
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRF_AA computes the factorization of a real symmetric matrix A using the Aasen's
 algorithm

 * @details
 * \b Purpose:
    \verbatim
     SYTRF_AA computes the factorization of a real symmetric matrix A
     using the Aasen's algorithm.  The form of the factorization is

        A = U**T*T*U  or  A = L*T*L**T

     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and T is a symmetric tridiagonal matrix.

     This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, the tridiagonal matrix is stored in the diagonals
          and the subdiagonals of A just below (or above) the diagonals,
          and L is stored below (or above) the subdiaonals, when UPLO
          is 'L' (or 'U'). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of A were interchanged with the
          row and column IPIV(k). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >= MAX(1,2*N). For optimum performance
          LWORK >= N*(1+NB), where NB is the optimal blocksize. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void sytrf_aa(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *work,
                  integer *lwork, integer *info)
    {
        sytrf_aa(uplo, n, a, lda, ipiv, work, lwork, info);
    }
    /** @}*/ // end of sytrf_aa

    /** @defgroup lasyf_aa lasyf_aa
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LASYF_AA factorizes a panel of a real symmetric matrix A using the Aasen's algorithm

 * @details
 * \b Purpose:
    \verbatim
    LASYF_AA factorizes a panel of a real symmetric matrix A using
    the Aasen's algorithm. The panel consists of a set of NB rows of A
    when UPLO is U, or a set of NB columns when UPLO is L.

    In order to factorize the panel, the Aasen's algorithm requires the
    last row, or column, of the previous panel. The first row, or column,
    of A is set to be the first row, or column, of an identity matrix,
    which is used to factorize the first panel.

    The resulting J-th row of U, or J-th column of L, is stored in the
    (J-1)-th row, or column, of A (without the unit diagonals), while
    the diagonal and subdiagonal of A are overwritten by those of T.
    \endverbatim

  * @param[in] UPLO
           UPLO is CHARACTER*1 \n
           = 'U':  Upper triangle of A is stored; \n
           = 'L':  Lower triangle of A is stored. \n
  * @param[in] J1
           J1 is INTEGER \n
           The location of the first row, or column, of the panel
           within the submatrix of A, passed to this routine, e.g.,
           when called by SSYTRF_AA, for the first panel, J1 is 1,
           while for the remaining panels, J1 is 2. \n
  * @param[in] M
           M is INTEGER \n
           The dimension of the submatrix. M >= 0. \n
  * @param[in] NB
           NB is INTEGER \n
           The dimension of the panel to be facotorized. \n
  * @param[in,out] A
           A is REAL array, dimension (LDA,M) for \n
           the first panel, while dimension (LDA,M+1) for the
           remaining panels. \n
  \n
           On entry, A contains the last row, or column, of
           the previous panel, and the trailing submatrix of A
           to be factorized, except for the first panel, only
           the panel is passed.
  \n
           On exit, the leading panel is factorized. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,M). \n
  * @param[out] IPIV
           IPIV is INTEGER array, dimension (M) \n
           Details of the row and column interchanges,
           the row and column k were interchanged with the row and
           column IPIV(k). \n
  * @param[in,out] H
           H is REAL workspace, dimension (LDH,NB). \n
  * @param[in] LDH
           LDH is INTEGER \n
           The leading dimension of the workspace H. LDH >= fla_max(1,M). \n
  * @param[out] WORK
           WORK is REAL workspace, dimension (M).  \n

 *  * */
    template <typename T>
    void lasyf_aa(char *uplo, integer *j1, integer *m, integer *nb, T *a, integer *lda,
                  integer *ipiv, T *h, integer *ldh, T *work)
    {
        lasyf_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
    }
    /** @}*/ // end of lasyf_aa

    /** @defgroup sytrs_aa sytrs_aa
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRS_AA solves a system of linear equations A*X = B with a real  \n
     symmetric matrix A
 * @details
 * \b Purpose:
    \verbatim
     SYTRS_AA solves a system of linear equations A*X = B with a real
     symmetric matrix A using the factorization A = U**T*T*U or
     A = L*T*L**T computed by SYTRF_AA.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U**T*T*U; \n
          = 'L':  Lower triangular, form is A = L*T*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          Details of factors computed by SSYTRF_AA. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges as computed by SSYTRF_AA. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,3*N-2). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sytrs_aa(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
                  integer *ldb, T *work, integer *lwork, integer *info)
    {
        sytrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of sytrs_aa

    /** @defgroup hetrf_aa hetrf_aa
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief HETRF_AA computes the factorization of a complex hermitian matrix A

 * @details
 * \b Purpose:
    \verbatim
    HETRF_AA computes the factorization of a complex hermitian matrix A
    using the Aasen's algorithm.  The form of the factorization is

       A = U**H*T*U  or  A = L*T*L**H

    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, and T is a hermitian tridiagonal matrix.

    This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the hermitian matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, the tridiagonal matrix is stored in the diagonals
          and the subdiagonals of A just below (or above) the diagonals,
          and L is stored below (or above) the subdiaonals, when UPLO
          is 'L' (or 'U'). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of A were interchanged with the
          row and column IPIV(k). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >= 2*N. For optimum performance
          LWORK >= N*(1+NB), where NB is the optimal blocksize. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void hetrf_aa(char *uplo, integer *n, T *a, integer *lda, integer *ipiv, T *work,
                  integer *lwork, integer *info)
    {
        hetrf_aa(uplo, n, a, lda, ipiv, work, lwork, info);
    }
    /** @}*/ // end of hetrf_aa

    /** @defgroup lahef_aa lahef_aa
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief LAHEF_AA factorizes a panel of a complex hermitian matrix A using the Aasen's
 algorithm

 * @details
 * \b Purpose:
    \verbatim
    LAHEF_AA factorizes a panel of a complex hermitian matrix A using
    the Aasen's algorithm. The panel consists of a set of NB rows of A
    when UPLO is U, or a set of NB columns when UPLO is L.

    In order to factorize the panel, the Aasen's algorithm requires the
    last row, or column, of the previous panel. The first row, or column,
    of A is set to be the first row, or column, of an identity matrix,
    which is used to factorize the first panel.

    The resulting J-th row of U, or J-th column of L, is stored in the
    (J-1)-th row, or column, of A (without the unit diagonals), while
    the diagonal and subdiagonal of A are overwritten by those of T.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] J1
          J1 is INTEGER \n
          The location of the first row, or column, of the panel
          within the submatrix of A, passed to this routine, e.g.,
          when called by CHETRF_AA, for the first panel, J1 is 1,
          while for the remaining panels, J1 is 2. \n
 * @param[in] M
          M is INTEGER \n
          The dimension of the submatrix. M >= 0. \n
 * @param[in] NB
          NB is INTEGER \n
          The dimension of the panel to be facotorized. \n
 * @param[in,out] A
          A is COMPLEX array, dimension (LDA,M) for
          the first panel, while dimension (LDA,M+1) for the
          remaining panels. \n
 \n
          On entry, A contains the last row, or column, of
          the previous panel, and the trailing submatrix of A
          to be factorized, except for the first panel, only
          the panel is passed. \n
 \n
          On exit, the leading panel is factorized. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the row and column interchanges,
          the row and column k were interchanged with the row and
          column IPIV(k). \n
 * @param[in,out] H
          H is COMPLEX workspace, dimension (LDH,NB). \n
 * @param[in] LDH
          LDH is INTEGER \n
          The leading dimension of the workspace H. LDH >= fla_max(1,M). \n
 * @param[out] WORK
          WORK is COMPLEX workspace, dimension (M). \n

 * */
    template <typename T>
    void lahef_aa(char *uplo, integer *j1, integer *m, integer *nb, T *a, integer *lda,
                  integer *ipiv, T *h, integer *ldh, T *work)
    {
        lahef_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work);
    }
    /** @}*/ // end of lahef_aa

    /** @defgroup hetrs_aa hetrs_aa
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief HETRS_AA solves a system of linear equations A*X = B with a complex hermitian matrix
 A

 * @details
 * \b Purpose:
    \verbatim
     HETRS_AA solves a system of linear equations A*X = B with a complex
     hermitian matrix A using the factorization A = U**H*T*U or
     A = L*T*L**H computed by HETRF_AA.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U**H*T*U; \n
          = 'L':  Lower triangular, form is A = L*T*L**H. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          Details of factors computed by CHETRF_AA. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges as computed by CHETRF_AA. \n
 * @param[in,out] B
          B is COMPLEX array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,3*N-2). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void hetrs_aa(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, integer *ipiv, T *b,
                  integer *ldb, T *work, integer *lwork, integer *info)
    {
        hetrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of hetrs__aa

    /** @defgroup hetrs_aa_2stage hetrs_aa_2stage
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief HETRS_AA_2STAGE solves a system of linear equations A*X = B with a real hermitian
 matrix A

 * @details
 * \b Purpose:
    \verbatim
     HETRS_AA_2STAGE solves a system of linear equations A*X = B with a real
     hermitian matrix A using the factorization A = U**T*T*U or
     A = L*T*L**T computed by HETRF_AA_2STAGE.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U**T*T*U; \n
          = 'L':  Lower triangular, form is A = L*T*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          Details of factors computed by CHETRF_AA_2STAGE. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] TB
          TB is COMPLEX array, dimension (LTB) \n
          Details of factors computed by CHETRF_AA_2STAGE. \n
 * @param[in] LTB
          LTB is INTEGER \n
          The size of the array TB. LTB >= 4*N. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges as computed by
          CHETRF_AA_2STAGE. \n
 * @param[in] IPIV2
          IPIV2 is INTEGER array, dimension (N) \n
          Details of the interchanges as computed by
          CHETRF_AA_2STAGE. \n
 * @param[in,out] B
          B is COMPLEX array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void hetrs_aa_2stage(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *tb,
                         integer *ltb, integer *ipiv, integer *ipiv2, T *b, integer *ldb,
                         integer *info)
    {
        hetrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
    }
    /** @}*/ // end of hetrs_aa_2stage

    /** @defgroup sytrs_aa_2stage sytrs_aa_2stage
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief SYTRS_AA_2STAGE solves a system of linear equations A*X = B with a real \n
     symmetric matrix A
 * @details
 * \b Purpose:
    \verbatim
     SYTRS_AA_2STAGE solves a system of linear equations A*X = B with a real
     symmetric matrix A using the factorization A = U**T*T*U or
     A = L*T*L**T computed by SYTRF_AA_2STAGE.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the details of the factorization are stored
          as an upper or lower triangular matrix. \n
          = 'U':  Upper triangular, form is A = U**T*T*U; \n
          = 'L':  Lower triangular, form is A = L*T*L**T. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          Details of factors computed by SSYTRF_AA_2STAGE. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] TB
          TB is REAL array, dimension (LTB) \n
          Details of factors computed by SSYTRF_AA_2STAGE. \n
 * @param[in] LTB
          LTB is INTEGER \n
          The size of the array TB. LTB >= 4*N. \n
 * @param[in] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges as computed by
          SSYTRF_AA_2STAGE. \n
 * @param[in] IPIV2
          IPIV2 is INTEGER array, dimension (N) \n
          Details of the interchanges as computed by
          SSYTRF_AA_2STAGE. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void sytrs_aa_2stage(char *uplo, integer *n, integer *nrhs, T *a, integer *lda, T *tb,
                         integer *ltb, integer *ipiv, integer *ipiv2, T *b, integer *ldb,
                         integer *info)
    {
        sytrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info);
    }
    /** @}*/ // end of sytrs_aa_2stage

    /** @defgroup sytrf_aa_2stage sytrf_aa_2stage
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief  SYTRF_AA_2STAGE computes the factorization of a real symmetric matrix A using the
 Aasen's algorithm

 * @details
 * \b Purpose:
    \verbatim
     SYTRF_AA_2STAGE computes the factorization of a real symmetric matrix A
     using the Aasen's algorithm.  The form of the factorization is

        A = U**T*T*U  or  A = L*T*L**T

     where U (or L) is a product of permutation and unit upper (lower)
     triangular matrices, and T is a symmetric band matrix with the
     bandwidth of NB (NB is internally selected and stored in TB( 1), and T is
     LU factorized with partial pivoting).

     This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, L is stored below (or above) the subdiaonal blocks,
          when UPLO  is 'L' (or 'U'). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] TB
          TB is REAL array, dimension (LTB) \n
          On exit, details of the LU factorization of the band matrix. \n
 * @param[in] LTB
          LTB is INTEGER \n
          The size of the array TB. LTB >= 4*N, internally
          used to select NB such that LTB >= (3*NB+1)*N. \n
 \n
          If LTB = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of LTB,
          returns this value as the first entry of TB, and
          no error message related to LTB is issued by XERBLA. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of A were interchanged with the
          row and column IPIV(k). \n
 * @param[out] IPIV2
          IPIV2 is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of T were interchanged with the
          row and column IPIV(k). \n
 * @param[out]	WORK
          WORK is REAL workspace of size LWORK \n
 * @param[in]	LWORK
          LWORK is INTEGER
          The size of WORK. LWORK >= N, internally used to select NB
          such that LWORK >= N*NB. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the WORK array,
          returns this value as the first entry of the WORK array, and
          no error message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, band LU factorization failed on i-th column \n

 *  * */
    template <typename T>
    void sytrf_aa_2stage(char *uplo, integer *n, T *a, integer *lda, T *tb, integer *ltb,
                         integer *ipiv, integer *ipiv2, T *work, integer *lwork, integer *info)
    {
        sytrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
    }
    /** @}*/ // end of sytrf_aa_2stage

    /** @defgroup hetrf_aa_2stage
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief HETRF_AA_2STAGE computes the factorization of a real hermitian matrix A

 * @details
 * \b Purpose:
    \verbatim
    HETRF_AA_2STAGE computes the factorization of a real hermitian matrix A
    using the Aasen's algorithm.  The form of the factorization is

       A = U**T*T*U  or  A = L*T*L**T

    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, and T is a hermitian band matrix with the
    bandwidth of NB (NB is internally selected and stored in TB( 1), and T is
    LU factorized with partial pivoting).

    This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the hermitian matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, L is stored below (or above) the subdiaonal blocks,
          when UPLO  is 'L' (or 'U'). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] TB
          TB is COMPLEX array, dimension (LTB) \n
          On exit, details of the LU factorization of the band matrix. \n
 * @param[in] LTB
          LTB is INTEGER \n
          The size of the array TB. LTB >= 4*N, internally
          used to select NB such that LTB >= (3*NB+1)*N. \n
 \n
          If LTB = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of LTB,
          returns this value as the first entry of TB, and
          no error message related to LTB is issued by XERBLA. \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of A were interchanged with the
          row and column IPIV(k). \n
 * @param[out] IPIV2
          IPIV2 is INTEGER array, dimension (N) \n
          On exit, it contains the details of the interchanges, i.e.,
          the row and column k of T were interchanged with the
          row and column IPIV(k). \n
 * @param[out]	WORK
          WORK is COMPLEX workspace of size LWORK \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The size of WORK. LWORK >= N, internally used to select NB
          such that LWORK >= N*NB. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the
          routine only calculates the optimal size of the WORK array,
          returns this value as the first entry of the WORK array, and
          no error message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = i, band LU factorization failed on i-th column \n

 *  * */
    template <typename T>
    void hetrf_aa_2stage(char *uplo, integer *n, T *a, integer *lda, T *tb, integer *ltb,
                         integer *ipiv, integer *ipiv2, T *work, integer *lwork, integer *info)
    {
        hetrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info);
    }
    /** @}*/ // end of hetrf_aa_2stage

    /** @}*/ // end of LDL_computation

    /** @defgroup hetrf_rook
     * @ingroup LDL_computation LDL_computation
     * @{
     */
    /*! @brief HETRF_ROOK computes the factorization of a complex hermitian matrix A

 * @details
 * \b Purpose:
    \verbatim
    HETRF_ROOK computes the factorization of a complex Hermitian matrix A
    using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
   T he form of the factorization is

      A = U*D*U**T  or  A = L*D*L**T

    where U (or L) is a product of permutation and unit upper (lower)
    triangular matrices, and D is Hermitian and block diagonal with
    1-by-1 and 2-by-2 diagonal blocks.

    This is the blocked version of the algorithm, calling Level 3 BLAS.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is COMPLEX array, dimension (LDA,N) \n
          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 \n
          On exit, the block diagonal matrix D and the multipliers used
          to obtain the factor U or L (see below for further details). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] IPIV
          IPIV is INTEGER array, dimension (N) \n
          Details of the interchanges and the block structure of D.
 \n
          If UPLO = 'U': \n
             Only the last KB elements of IPIV are set.
 \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
             interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
             columns k and -IPIV(k) were interchanged and rows and
             columns k-1 and -IPIV(k-1) were inerchaged,
             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
 \n
          If UPLO = 'L': \n
             Only the first KB elements of IPIV are set.
 \n
             If IPIV(k) > 0, then rows and columns k and IPIV(k)
             were interchanged and D(k,k) is a 1-by-1 diagonal block.
 \n
             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
             columns k and -IPIV(k) were interchanged and rows and
             columns k+1 and -IPIV(k+1) were inerchaged,
             D(k:k+1,k:k+1) is a 2-by-2 diagonal block. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)). \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The length of WORK.  LWORK >=1.  For best performance
          LWORK >= N*NB, where NB is the block size returned by ILAENV. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
                has been completed, but the block diagonal matrix D is
                exactly singular, and division by zero will occur if it
                is used to solve a system of equations. \n

 *  * */
    template <typename T>
    void hetrf_rook(char* uplo, integer* n, T* a, integer* lda, integer* ipiv,
                    T* work, integer* lwork, integer* info)
    {
        hetrf_rook(uplo, n, a, lda, ipiv, work, lwork, info);
    }
    /** @}*/ // end of hetrf_rook

    /** @}*/ // end of LDL_computation

    /** @defgroup Triangular Triangular Computational
     * @ingroup LinearSolve
     * @{
     */
    /** @defgroup trcon trcon
     * @ingroup Triangular
     * @{
     */
    /*! @brief TRCON estimates the reciprocal of the condition number of a triangular matrix A

 * @details
 * \b Purpose:
    \verbatim
     TRCON estimates the reciprocal of the condition number of a
     triangular matrix A, in either the 1-norm or the infinity-norm.

     The norm of A is computed and an estimate is obtained for
     norm(inv(A)), then the reciprocal of the condition number is
     computed as
        RCOND = 1 / ( norm(A) * norm(inv(A))).
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required: \n
          = '1' or 'O':  1-norm; \n
          = 'I':         Infinity-norm. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of the array A contains the upper
          triangular matrix, and the strictly lower triangular part of
          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
          triangular part of the array A contains the lower triangular
          matrix, and the strictly upper triangular part of A is not
          referenced.  If DIAG = 'U', the diagonal elements of A are
          also not referenced and are assumed to be 1. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(norm(A) * norm(inv(A))). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void trcon(char *norm, char *uplo, char *diag, integer *n, T *a, integer *lda, T *rcond,
               T *work, integer *iwork, integer *info)
    {
        trcon(norm, uplo, diag, n, a, lda, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void trcon(char *norm, char *uplo, char *diag, integer *n, T *a, integer *lda, Ta *rcond,
               T *work, Ta *rwork, integer *info)
    {
        trcon(norm, uplo, diag, n, a, lda, rcond, work, rwork, info);
    }
    /** @}*/ // end of trcon

    /** @defgroup trtrs trtrs
     * @ingroup Triangular
     * @{
     */
    /*! @brief TRTRS solves a triangular system of the form A * X = B  or  A**T * X = B

 * @details
 * \b Purpose:
    \verbatim
    TRTRS solves a triangular system of the form

       A * X = B  or  A**T * X = B,

    where A is a triangular matrix of order N, and B is an N-by-NRHS
    matrix.  A check is made to verify that A is nonsingular.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of the array A contains the upper
          triangular matrix, and the strictly lower triangular part of
          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
          triangular part of the array A contains the lower triangular
          matrix, and the strictly upper triangular part of A is not
          referenced.  If DIAG = 'U', the diagonal elements of A are
          also not referenced and are assumed to be 1. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, if INFO = 0, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, the i-th diagonal element of A is zero,
               indicating that the matrix is singular and the solutions
               X have not been computed. \n

 *  * */
    template <typename T>
    void trtrs(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, T *a, integer *lda,
               T *b, integer *ldb, integer *info)
    {
        trtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
    }
    /** @}*/ // end of trtrs

    /** @defgroup latrs latrs
     * @ingroup Triangular
     * @{
     */
    /*! @brief LATRS solves a triangular system of equations with the scale factor set to prevent
 overflow

 * @details
 * \b Purpose:
    \verbatim
    LATRS solves one of the triangular systems

       A *x = s*b  or  A**T*x = s*b

    with scaling to prevent overflow.  Here A is an upper or lower
    triangular matrix, A**T denotes the transpose of A, x and b are
    n-element vectors, and s is a scaling factor, usually less than
    or equal to 1, chosen so that the components of x will be less than
    the overflow threshold.  If the unscaled problem will not cause
    overflow, the Level 2 BLAS routine STRSV is called.  If the matrix A
    is singular (A(j,j) = 0 for some j), then s is set to 0 and a
    non-trivial solution to A*x = 0 is   returned.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the matrix A is upper or lower triangular. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the operation applied to A. \n
          = 'N':  Solve A * x = s*b  (No transpose) \n
          = 'T':  Solve A**T* x = s*b  (Transpose) \n
          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          Specifies whether or not the matrix A is unit triangular. \n
          = 'N':  Non-unit triangular \n
          = 'U':  Unit triangular \n
 * @param[in] NORMIN
          NORMIN is CHARACTER*1 \n
          Specifies whether CNORM has been set or not. \n
          = 'Y':  CNORM contains the column norms on entry \n
          = 'N':  CNORM is not set on entry.  On exit, the norms will
                  be computed and stored in CNORM. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The triangular matrix A.  If UPLO = 'U', the leading n by n
          upper triangular part of the array A contains the upper
          triangular matrix, and the strictly lower triangular part of
          A is not referenced.  If UPLO = 'L', the leading n by n lower
          triangular part of the array A contains the lower triangular
          matrix, and the strictly upper triangular part of A is not
          referenced.  If DIAG = 'U', the diagonal elements of A are
          also not referenced and are assumed to be 1. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= max (1,N). \n
 * @param[in,out] X
          X is REAL array, dimension (N) \n
          On entry, the right hand side b of the triangular system.
          On exit, X is overwritten by the solution vector x. \n
 * @param[out] SCALE
          SCALE is REAL \n
          The scaling factor s for the triangular system \n
             A * x = s*b  or  A**T* x = s*b. \n
          If SCALE = 0, the matrix A is singular or badly scaled, and
          the vector x is an exact or approximate solution to A*x = 0. \n
 * @param[in,out] CNORM
          CNORM is REAL array, dimension (N) \n
 \n
          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
          contains the norm of the off-diagonal part of the j-th column
          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
          must be greater than or equal to the 1-norm. \n
 \n
          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
          returns the 1-norm of the offdiagonal part of the j-th column
          of A. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -k, the k-th argument had an illegal value  \n

 *  * */
    template <typename T>
    void latrs(char *uplo, char *trans, char *diag, char *normin, integer *n, T *a, integer *lda,
               T *x, T *scale, T *cnorm, integer *info)
    {
        latrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    }
    template <typename T, typename Ta>
    void latrs(char *uplo, char *trans, char *diag, char *normin, integer *n, T *a, integer *lda,
               T *x, Ta *scale, Ta *cnorm, integer *info)
    {
        latrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info);
    }
    /** @}*/ // end of latrs

    /** @defgroup trtri trtri
     * @ingroup Triangular
     * @{
     */
    /*! @brief Inverse of a real upper or lower triangular matrix.
        *
        * @details
        * \b Purpose:
        * \verbatim
            Computation of inverse of a real upper or lower triangular matrix

            This is the Level 3 BLAS version of the algorithm.
        \endverbatim

        * @param[in] uplo
                  uplo is char* \n
                  = 'U':  A is upper triangular; \n
                  = 'L':  A is lower triangular. \n
        * @param[in] diag
                  diag is char* \n
                  = 'N':  A is non-unit triangular; \n
                  = 'U':  A is unit triangular. \n
        * @param[in] n
                  n is integer*  \n
                  The order of the matrix a.  n >= 0. \n
        * @param[in,out] a
                  a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
                  On entry, the triangular matrix a.  If uplo = 'U', the
                  leading n-by-n upper triangular part of the array a contains
                  the upper triangular matrix, and the strictly lower
                  triangular part of a is not referenced.  If uplo = 'L', the
                  leading n-by-n lower triangular part of the array a contains
                  the lower triangular matrix, and the strictly upper
                  triangular part of a is not referenced.  If diag = 'U', the
                  diagonal elements of a are also not referenced and are
                  assumed to be 1. \n
                  On exit, the (triangular) inverse of the original matrix, in
                  the same storage format. \n
        * @param[in] lda
                  lda is integer* \n
                  The leading dimension of the matrix a, lda >= fla_max(1,n) \n
        * @param[out]	INFO
                  INFO is INTEGER \n
                  = 0: successful exit \n
                  < 0: if INFO = -i, the i-th argument had an illegal value \n
                  > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
                       matrix is singular and its inverse can not be computed. \n

        *     *  */
    template <typename T>
    void trtri(char *uplo, char *diag, integer *n, T *a, integer *lda, integer *info)
    {
        trtri(uplo, diag, n, a, lda, info);
    }
    /** @}*/ // end of trtri

    /** @defgroup trti2 trti2
     * @ingroup Triangular
     * @{
     */
    /*! @brief Inverse of a triangular matrix (unblocked algorithm).
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of inverse of a real upper or lower triangular matrix (unblocked algorithm)

        This is the Level 2 BLAS version of the algorithm.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  A is upper triangular; \n
              = 'L':  A is lower triangular. \n
    * @param[in] diag
              diag is char* \n
              = 'N':  A is non-unit triangular; \n
              = 'U':  A is unit triangular. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the triangular matrix a.  If uplo = 'U', the
              leading n by n upper triangular part of the array a contains
              the upper triangular matrix, and the strictly lower
              triangular part of a is not referenced.  If uplo = 'L', the
              leading n by n lower triangular part of the array a contains
              the lower triangular matrix, and the strictly upper
              triangular part of a is not referenced.  If diag = 'U', the
              diagonal elements of a are also not referenced and are
              assumed to be 1. \n
              On exit, the (triangular) inverse of the original matrix, in
              the same storage format. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,n) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -k, the k-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void trti2(char *uplo, char *diag, integer *n, T *a, integer *lda, integer *info)
    {
        trti2(uplo, diag, n, a, lda, info);
    }
    /** @}*/ // end of trti2

    /** @defgroup trrfs trrfs
     * @ingroup Triangular
     * @{
     */
    /*! @brief TRRFS provides error bounds and backward error estimates for the \n
     solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
     TRRFS provides error bounds and backward error estimates for the
     solution to a system of linear equations with a triangular
     coefficient matrix.

     The solution matrix X must be computed by STRTRS or some other
     means before entering this routine.  STRRFS does not do iterative
     refinement because doing so cannot improve the backward error.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of the array A contains the upper
          triangular matrix, and the strictly lower triangular part of
          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
          triangular part of the array A contains the lower triangular
          matrix, and the strictly upper triangular part of A is not
          referenced.  If DIAG = 'U', the diagonal elements of A are
          also not referenced and are assumed to be 1. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in] X
          X is REAL array, dimension (LDX,NRHS) \n
          The solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void trrfs(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, T *a, integer *lda,
               T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr, T *work, integer *iwork,
               integer *info)
    {
        trrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void trrfs(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, T *a, integer *lda,
               T *b, integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork,
               integer *info)
    {
        trrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }
    /** @}*/ // end of trrfs

    /** @defgroup lauum lauum
     * @ingroup Triangular
     * @{
     */
    /*! @brief Product UUH or LHL, where U and L are upper or lower triangular matrices (blocked
    algorithm).
    *
    * @details
    * \b Purpose:
    * \verbatim
        Product UUH or LHL, where U and L are upper or lower triangular matrices (blocked
    algorithm). Computation of the product U * U**T or L**T * L, where the triangular factor U or L
    is stored in the upper or lower triangular part of the array a.

        If uplo = 'U' or 'u' then the upper triangle of the result is stored, overwriting the factor
    U in A. If uplo = 'L' or 'l' then the lower triangle of the result is stored, overwriting the
    factor L in A.

        This is the blocked form of the algorithm, calling Level 3 BLAS.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              Specifies whether the triangular factor stored in the array a
              is upper or lower triangular: \n
              = 'U':  Upper triangular \n
              = 'L':  Lower triangular \n
    * @param[in] n
              n is integer* \n
              The order of the triangular factor U or L.  n >= 0. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the triangular factor U or L. \n
              On exit, if uplo = 'U', the upper triangle of a is
              overwritten with the upper triangle of the product U * U**T; \n
              if uplo = 'L', the lower triangle of a is overwritten with
              the lower triangle of the product L**T * L. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,n) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -k, the k-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void lauum(char *uplo, integer *n, T *a, integer *lda, integer *info)
    {
        lauum(uplo, n, a, lda, info);
    }
    /** @}*/ // end of lauum

    /** @defgroup lauu2 lauu2
     * @ingroup Triangular
     * @{
     */
    /*! @brief Product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked
    algorithm).
    *
    * @details
    * \b Purpose:
    * \verbatim
        Product UUH or LHL, where U and L are upper or lower triangular matrices (unblocked
    algorithm). Computation of the product U * U**T or L**T * L, where the triangular factor U or L
    is stored in the upper or lower triangular part of the array a.

        If uplo = 'U' or 'u' then the upper triangle of the result is stored, overwriting the factor
    U in A. If uplo = 'L' or 'l' then the lower triangle of the result is stored, overwriting the
    factor L in A.

        This is the unblocked form of the algorithm, calling Level 2 BLAS.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              Specifies whether the triangular factor stored in the array a
              is upper or lower triangular: \n
              = 'U':  Upper triangular \n
              = 'L':  Lower triangular \n
    * @param[in] n
              n is integer* \n
              The order of the triangular factor U or L.  n >= 0. \n
    * @param[in,out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the triangular factor U or L. \n
              On exit, if uplo = 'U', the Upper triangle of a is
              overwritten with the upper triangle of the product U * U**T; \n
              if uplo = 'L', the lower triangle of a is overwritten with
              the lower triangle of the product L**T * L. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,n) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -k, the k-th argument had an illegal value \n

    *  */
    template <typename T>
    void lauu2(char *uplo, integer *n, T *a, integer *lda, integer *info)
    {
        lauu2(uplo, n, a, lda, info);
    }
    /** @}*/ // end of lauu2

    /** @defgroup tpcon tpcon
     * @ingroup Triangular
     * @{
     */
    /*! @brief TPCON estimates the reciprocal of the condition number of a packed
    triangular matrix A
 * @details
 * \b Purpose:
    \verbatim
     TPCON estimates the reciprocal of the condition number of a packed
     triangular matrix A, in either the 1-norm or the infinity-norm.

     The norm of A is computed and an estimate is obtained for
     norm(inv(A)), then the reciprocal of the condition number is
     computed as
        RCOND = 1 / ( norm(A) * norm(inv(A))).
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required: \n
          = '1' or 'O':  1-norm; \n
          = 'I':         Infinity-norm. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangular matrix A, packed columnwise in
          a linear array.  The j-th column of A is stored in the array
          AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
          If DIAG = 'U', the diagonal elements of A are not referenced
          and are assumed to be 1. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(norm(A) * norm(inv(A))). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tpcon(char *norm, char *uplo, char *diag, integer *n, T *ap, T *rcond, T *work,
               integer *iwork, integer *info)
    {
        tpcon(norm, uplo, diag, n, ap, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void tpcon(char *norm, char *uplo, char *diag, integer *n, T *ap, Ta *rcond, T *work, Ta *rwork,
               integer *info)
    {
        tpcon(norm, uplo, diag, n, ap, rcond, work, rwork, info);
    }
    /** @}*/ // end of tpcon

    /** @defgroup tptrs tptrs
     * @ingroup Triangular
     * @{
     */
    /*! @brief TPTRS solves a triangular system of the form A * X = B  or  A**T * X = B

 * @details
 * \b Purpose:
    \verbatim
     TPTRS solves a triangular system of the form

        A * X = B  or  A**T * X = B,

     where A is a triangular matrix of order N stored in packed format,
     and B is an N-by-NRHS matrix.  A check is made to verify that A is
     nonsingular.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangular matrix A, packed columnwise in
          a linear array.  The j-th column of A is stored in the array
          AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, if INFO = 0, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the i-th diagonal element of A is zero,
                indicating that the matrix is singular and the
                solutions X have not been computed. \n

 *  * */
    template <typename T>
    void tptrs(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, T *ap, T *b,
               integer *ldb, integer *info)
    {
        tptrs(uplo, trans, diag, n, nrhs, ap, b, ldb, info);
    }
    /** @}*/ // end of trptrs

    /** @defgroup latps latps
     * @ingroup Triangular
     * @{
     */
    /*! @brief LATPS solves a triangular system of equations with the matrix held in packed storage

 * @details
 * \b Purpose:
    \verbatim
    LATPS solves one of the triangular systems

       A *x = s*b  or  A**T*x = s*b

    with scaling to prevent overflow, where A is an upper or lower
    triangular matrix stored in packed form.  Here A**T denotes the
    transpose of A, x and b are n-element vectors, and s is a scaling
    factor, usually less than or equal to 1, chosen so that the
    components of x will be less than the overflow threshold.  If the
    unscaled problem will not cause overflow, the Level 2 BLAS routine
    STPSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
    then s is set to 0 and a non-trivial solution to A*x = 0 is   returned.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the matrix A is upper or lower triangular. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the operation applied to A. \n
          = 'N':  Solve A * x = s*b  (No transpose) \n
          = 'T':  Solve A**T* x = s*b  (Transpose) \n
          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          Specifies whether or not the matrix A is unit triangular. \n
          = 'N':  Non-unit triangular \n
          = 'U':  Unit triangular \n
 * @param[in] NORMIN
          NORMIN is CHARACTER*1 \n
          Specifies whether CNORM has been set or not. \n
          = 'Y':  CNORM contains the column norms on entry \n
          = 'N':  CNORM is not set on entry.  On exit, the norms will
                  be computed and stored in CNORM. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangular matrix A, packed columnwise in
          a linear array.  The j-th column of A is stored in the array
          AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[in,out] X
          X is REAL array, dimension (N) \n
          On entry, the right hand side b of the triangular system.
          On exit, X is overwritten by the solution vector x. \n
 * @param[out] SCALE
          SCALE is REAL \n
          The scaling factor s for the triangular system \n
             A * x = s*b  or  A**T* x = s*b. \n
          If SCALE = 0, the matrix A is singular or badly scaled, and
          the vector x is an exact or approximate solution to A*x = 0. \n
 * @param[in,out] CNORM
          CNORM is REAL array, dimension (N) \n
 \n
          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
          contains the norm of the off-diagonal part of the j-th column
          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
          must be greater than or equal to the 1-norm.
 \n
          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
          returns the 1-norm of the offdiagonal part of the j-th column
          of A. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -k, the k-th argument had an illegal value  \n

 *  * */
    template <typename T>
    void latps(char *uplo, char *trans, char *diag, char *normin, integer *n, T *ap, T *x, T *scale,
               T *cnorm, integer *info)
    {
        latps(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
    }
    template <typename T, typename Ta>
    void latps(char *uplo, char *trans, char *diag, char *normin, integer *n, T *ap, T *x,
               Ta *scale, Ta *cnorm, integer *info)
    {
        latps(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info);
    }

    /** @}*/ // end of latps

    /** @defgroup tptri tptri
     * @ingroup Triangular
     * @{
     */
    /*! @brief TPTRI computes the inverse of a real upper or lower triangular \n
     matrix A
 * @details
 * \b Purpose:
    \verbatim
     TPTRI computes the inverse of a real upper or lower triangular
     matrix A stored in packed format.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On entry, the upper or lower triangular matrix A, stored
          columnwise in a linear array.  The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n. \n
          See below for further details.
          On exit, the (triangular) inverse of the original matrix, in
          the same packed storage format. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular
                matrix is singular and its inverse can not be computed. \n

 *  * */
    template <typename T>
    void tptri(char *uplo, char *diag, integer *n, T *ap, integer *info)
    {
        tptri(uplo, diag, n, ap, info);
    }
    /** @}*/ // end of tptri

    /** @defgroup tprfs tprfs
     * @ingroup Triangular
     * @{
     */
    /*! @brief TPRFS provides error bounds and backward error estimates \n
     for the solution to a system of linear equations
 * @details
 * \b Purpose:
    \verbatim
     TPRFS provides error bounds and backward error estimates for the
     solution to a system of linear equations with a triangular packed
     coefficient matrix.

     The solution matrix X must be computed by STPTRS or some other
     means before entering this routine.  STPRFS does not do iterative
     refinement because doing so cannot improve the backward error.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangular matrix A, packed columnwise in
          a linear array.  The j-th column of A is stored in the array
          AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. \n
          If DIAG = 'U', the diagonal elements of A are not referenced
          and are assumed to be 1. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in] X
          X is REAL array, dimension (LDX,NRHS) \n
          The solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X). \n
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tprfs(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, T *ap, T *b,
               integer *ldb, T *x, integer *ldx, T *ferr, T *berr, T *work, integer *iwork,
               integer *info)
    {
        tprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, iwork, info);
    }
    template <typename T, typename Ta>
    void tprfs(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, T *ap, T *b,
               integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work, Ta *rwork,
               integer *info)
    {
        tprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, rwork, info);
    }

    /** @}*/ // end of tprfs

    /** @defgroup tftri tftri
     * @ingroup Triangular
     * @{
     */
    /*! @brief TFTRI computes the inverse of a triangular matrix A stored in RFP format

 * @details
 * \b Purpose:
    \verbatim
     TFTRI computes the inverse of a triangular matrix A stored in RFP
     format.

     This is a Level 3 BLAS version of the algorithm.
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  The Normal TRANSR of RFP A is stored; \n
          = 'T':  The Transpose TRANSR of RFP A is stored. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (NT); \n
          NT=N*(N+1)/2. On entry, the triangular factor of a Hermitian
          Positive Definite matrix A in RFP format. RFP format is
          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'
          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is
          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'T' then RFP is
          the transpose of RFP A as defined when
          TRANSR = 'N'. The contents of RFP A are defined by UPLO as
          follows: If UPLO = 'U' the RFP A contains the nt elements of
          upper packed A; If UPLO = 'L' the RFP A contains the nt
          elements of lower packed A. The LDA of RFP A is (N+1)/2 when
          TRANSR = 'T'. When TRANSR is 'N' the LDA is N+1 when N is
          even and N is odd. See the Note below for more details. \n
 \n
          On exit, the (triangular) inverse of the original matrix, in
          the same storage format. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
               matrix is singular and its inverse can not be computed. \n

 *  * */
    template <typename T>
    void tftri(char *transr, char *uplo, char *diag, integer *n, T *a, integer *info)
    {
        tftri(transr, uplo, diag, n, a, info);
    }
    /** @}*/ // end of tftri

    /** @defgroup tbcon tbcon
     * @ingroup Triangular
     * @{
     */
    /*! @brief TBCON estimates the reciprocal of the condition number of a  \n
     triangular band matrix A

 * @details
 * \b Purpose:
    \verbatim
     TBCON estimates the reciprocal of the condition number of a
     triangular band matrix A, in either the 1-norm or the infinity-norm.

     The norm of A is computed and an estimate is obtained for
     norm(inv(A)), then the reciprocal of the condition number is
     computed as
        RCOND = 1 / ( norm(A) * norm(inv(A))).
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies whether the 1-norm condition number or the
          infinity-norm condition number is required: \n
          = '1' or 'O':  1-norm; \n
          = 'I':         Infinity-norm. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of superdiagonals or subdiagonals of the
          triangular band matrix A.  KD >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The upper or lower triangular band matrix A, stored in the
          first kd+1 rows of the array. The j-th column of A is stored
          in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
          If DIAG = 'U', the diagonal elements of A are not referenced
          and are assumed to be 1. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[out] RCOND
          RCOND is REAL \n
          The reciprocal of the condition number of the matrix A,
          computed as RCOND = 1/(norm(A) * norm(inv(A))). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tbcon(char *norm, char *uplo, char *diag, integer *n, integer *kd, T *ab, integer *ldab,
               T *rcond, T *work, integer *iwork, integer *info)
    {
        tbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, iwork, info);
    }
    template <typename T, typename Ta>
    void tbcon(char *norm, char *uplo, char *diag, integer *n, integer *kd, T *ab, integer *ldab,
               Ta *rcond, T *work, Ta *rwork, integer *info)
    {
        tbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, rwork, info);
    }
    /** @}*/ // end of tbcon

    /** @defgroup tbtrs tbtrs
     * @ingroup Triangular
     * @{
     */
    /*! @brief TBTRS solves a triangular system of the form A * X = B  or  A**T * X = B

 * @details
 * \b Purpose:
    \verbatim
     TBTRS solves a triangular system of the form

        A * X = B  or  A**T * X = B,

     where A is a triangular band matrix of order N, and B is an
     N-by NRHS matrix.  A check is made to verify that A is nonsingular.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form the system of equations: \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of superdiagonals or subdiagonals of the
          triangular band matrix A.  KD >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The upper or lower triangular band matrix A, stored in the
          first kd+1 rows of AB.  The j-th column of A is stored
          in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
          If DIAG = 'U', the diagonal elements of A are not referenced
          and are assumed to be 1. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the right hand side matrix B.
          On exit, if INFO = 0, the solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, the i-th diagonal element of A is zero,
                indicating that the matrix is singular and the
                solutions X have not been computed. \n

 *  * */
    template <typename T>
    void tbtrs(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, T *ab,
               integer *ldab, T *b, integer *ldb, integer *info)
    {
        tbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info);
    }
    /** @}*/ // end of tbtrs

    /** @defgroup latbs latbs
     * @ingroup Triangular
     * @{
     */
    /*! @brief LATBS solves a triangular banded system of equations

 * @details
 * \b Purpose:
    \verbatim
    LATBS solves one of the triangular systems

       A *x = s*b  or  A**T*x = s*b

    with scaling to prevent overflow, where A is an upper or lower
    triangular band matrix.  Here A**T denotes the transpose of A, x and b
    are n-element vectors, and s is a scaling factor, usually less than
    or equal to 1, chosen so that the components of x will be less than
    the overflow threshold.  If the unscaled problem will not cause
    overflow, the Level 2 BLAS routine STBSV is called.  If the matrix A
    is singular (A(j,j) = 0 for some j), then s is set to 0 and a
    non-trivial solution to A*x = 0 is   returned.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the matrix A is upper or lower triangular. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the operation applied to A. \n
          = 'N':  Solve A * x = s*b  (No transpose) \n
          = 'T':  Solve A**T* x = s*b  (Transpose) \n
          = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          Specifies whether or not the matrix A is unit triangular. \n
          = 'N':  Non-unit triangular \n
          = 'U':  Unit triangular \n
 * @param[in] NORMIN
          NORMIN is CHARACTER*1 \n
          Specifies whether CNORM has been set or not. \n
          = 'Y':  CNORM contains the column norms on entry \n
          = 'N':  CNORM is not set on entry.  On exit, the norms will
                  be computed and stored in CNORM. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of subdiagonals or superdiagonals in the
          triangular matrix A.  KD >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The upper or lower triangular band matrix A, stored in the
          first KD+1 rows of the array. The j-th column of A is stored
          in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[in,out] X
          X is REAL array, dimension (N) \n
          On entry, the right hand side b of the triangular system.
          On exit, X is overwritten by the solution vector x. \n
 * @param[out] SCALE
          SCALE is REAL \n
          The scaling factor s for the triangular system \n
             A * x = s*b  or  A**T* x = s*b. \n
          If SCALE = 0, the matrix A is singular or badly scaled, and
          the vector x is an exact or approximate solution to A*x = 0. \n
 * @param[in,out] CNORM
          CNORM is REAL array, dimension (N) \n
 \n
          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
          contains the norm of the off-diagonal part of the j-th column
          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
          must be greater than or equal to the 1-norm.
 \n
          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
          returns the 1-norm of the offdiagonal part of the j-th column
          of A. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -k, the k-th argument had an illegal value  \n

 *  * */
    template <typename T>
    void latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, T *ab,
               integer *ldab, T *x, T *scale, T *cnorm, integer *info)
    {
        latbs(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
    }
    template <typename T, typename Ta>
    void latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, T *ab,
               integer *ldab, T *x, Ta *scale, Ta *cnorm, integer *info)
    {
        latbs(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info);
    }
    /** @}*/ // end of latbs

    /** @defgroup tbrfs tbrfs
     * @ingroup Triangular
     * @{
     */
    /*! @brief TBRFS provides error bounds and backward error estimates for the  \n
     solution to a system of linear equations

 * @details
 * \b Purpose:
    \verbatim
     TBRFS provides error bounds and backward error estimates for the
     solution to a system of linear equations with a triangular band
     coefficient matrix.

     The solution matrix X must be computed by STBTRS or some other
     means before entering this routine.  STBRFS does not do iterative
     refinement because doing so cannot improve the backward error.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the form of the system of equations: \n
          = 'N':  A * X = B  (No transpose) \n
          = 'T':  A**T * X = B  (Transpose) \n
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          = 'N':  A is non-unit triangular; \n
          = 'U':  A is unit triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] KD
          KD is INTEGER \n
          The number of superdiagonals or subdiagonals of the
          triangular band matrix A.  KD >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The upper or lower triangular band matrix A, stored in the
          first kd+1 rows of the array. The j-th column of A is stored
          in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for fla_max(1,j-kd)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). \n
          If DIAG = 'U', the diagonal elements of A are not referenced
          and are assumed to be 1. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KD+1. \n
 * @param[in] B
          B is REAL array, dimension (LDB,NRHS) \n
          The right hand side matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,N). \n
 * @param[in] X
          X is REAL array, dimension (LDX,NRHS) \n
          The solution matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(1,N). \n
 * @param[out] FERR
          FERR is REAL array, dimension (NRHS) \n
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X). \n
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error. \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (3*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tbrfs(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, T *ab,
               integer *ldab, T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr, T *work,
               integer *iwork, integer *info)
    {
        tbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr, work, iwork,
              info);
    }
    template <typename T, typename Ta>
    void tbrfs(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, T *ab,
               integer *ldab, T *b, integer *ldb, T *x, integer *ldx, Ta *ferr, Ta *berr, T *work,
               Ta *rwork, integer *info)
    {
        tbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr, work, rwork,
              info);
    }
    /** @}*/ // end of tbrfs

    /** @}*/ // end of Triangular

    /** @defgroup Auxilliary Auxilliary routines
     * @ingroup LinearSolve
     * @{
     */

    /** @defgroup lacn2 lacn2
     * @ingroup Auxilliary
     * @{
     */
    /*! @brief LACN2 estimates the 1-norm of a square matrix

 * @details
 * \b Purpose:
    \verbatim
     LACN2 estimates the 1-norm of a square, real matrix A.
     Reverse communication is used for evaluating matrix-vector products.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 1. \n
 * @param[out] V
          V is REAL array, dimension (N) \n
          On the final return, V = A*W,  where  EST = norm(V)/norm(W)
          (W is not returned). \n
 * @param[in,out] X
          X is REAL array, dimension (N) \n
          On an intermediate return, X should be overwritten by \n
                A * X,   if KASE=1, \n
                A**T * X,  if KASE=2, \n
          and SLACN2 must be re-called with all the other parameters
          unchanged. \n
 * @param[out] ISGN
          ISGN is INTEGER array, dimension (N) \n
 * @param[in,out] EST
          EST is REAL \n
          On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
          unchanged from the previous call to SLACN2. \n
          On exit, EST is an estimate (a lower bound) for norm(A). \n
 * @param[in,out] KASE
          KASE is INTEGER \n
          On the initial call to SLACN2, KASE should be 0. \n
          On an intermediate return, KASE will be 1 or 2, indicating
          whether X should be overwritten by A * X  or A**T * X.
          On the final return from SLACN2, KASE will again be 0. \n
 * @param[in,out] ISAVE
          ISAVE is INTEGER array, dimension (3) \n
          ISAVE is used to save variables between calls to SLACN2 \n

 * */
    template <typename T>
    void lacn2(integer *n, T *v, T *x, integer *isgn, T *est, integer *kase, integer *isave)
    {
        lacn2(n, v, x, isgn, est, kase, isave);
    }
    template <typename T, typename Ta>
    void lacn2(integer *n, T *v, T *x, Ta *est, integer *kase, integer *isave)
    {
        lacn2(n, v, x, est, kase, isave);
    }
    /** @}*/ // end of lacn2

    /** @defgroup lacon lacon
     * @ingroup Auxilliary
     * @{
     */
    /*! @brief LACON estimates the 1-norm of a square matrix, using reverse communication for
 evaluating matrix-vector products

 * @details
 * \b Purpose:
    \verbatim
    LACON estimates the 1-norm of a square, real matrix A.
    Reverse communication is used for evaluating matrix-vector products.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The order of the matrix.  N >= 1. \n
 * @param[out] V
          V is REAL array, dimension (N) \n
          On the final   return, V = A*W, where  EST = norm(V)/norm(W)
          (W is not   returned). \n
 * @param[in,out] X
          X is REAL array, dimension (N) \n
          On an intermediate   return, X should be overwritten by \n
               A * X,  if KASE=1, \n
               A**T * X, if KASE=2, \n
          and SLACON must be re-called with all the other parameters
          unchanged. \n
 * @param[out] ISGN
          ISGN is INTEGER array, dimension (N) \n
 * @param[in,out] EST
          EST is REAL \n
          On entry with KASE = 1 or 2 and JUMP = 3, EST should be
          unchanged from the previous call to SLACON. \n
          On exit, EST is an estimate (a lower bound) for norm(A). \n
 * @param[in,out] KASE
          KASE is INTEGER \n
          On the initial call to SLACON, KASE should be 0.
          On an intermediate return, KASE will be 1 or 2, indicating
          whether X should be overwritten by A * X  or A**T * X.
          On the final return from SLACON, KASE will again be 0.  \n

 * */
    template <typename T>
    void lacon(integer *n, T *v, T *x, integer *isgn, T *est, integer *kase)
    {
        lacon(n, v, x, isgn, est, kase);
    }
    template <typename T, typename Ta>
    void lacon(integer *n, T *v, T *x, Ta *est, integer *kase)
    {
        lacon(n, v, x, est, kase);
    }
    /** @}*/ // end of lacon

    /** @defgroup la_lin_berr la_lin_berr
     * @ingroup Auxilliary
     * @{
     */
    /*! @brief LA_LIN_BERR computes a component-wise relative backward error

 * @details
 * \b Purpose:
    \verbatim
    LA_LIN_BERR computes componentwise relative backward error from
    the formula
        fla_max(i) (abs(R(i)) / (abs(op(A_s))*abs(Y) + abs(B_s))(i))
    where abs(Z) is the componentwise absolute value of the matrix
    or vector Z.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0. \n
 * @param[in] NZ
          NZ is INTEGER \n
          We add (NZ+1)*SLAMCH('Safe minimum') to R(i) in the numerator to
          guard against spuriously zero residuals. Default value is N. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices AYB, RES, and BERR.  NRHS >= 0. \n
 * @param[in] RES
          RES is REAL array, dimension (N,NRHS) \n
          The residual matrix, i.e., the matrix R in the relative backward
          error formula above. \n
 * @param[in] AYB
          AYB is REAL array, dimension (N, NRHS) \n
          The denominator in the relative backward error formula above, i.e.,
          the matrix abs(op(A_s))*abs(Y) + abs(B_s). The matrices A, Y, and B
          are from iterative refinement (see sla_gerfsx_extended.f). \n
 * @param[out] BERR
          BERR is REAL array, dimension (NRHS) \n
          The componentwise relative backward error from the formula above. \n

 * */
    template <typename T>
    void la_lin_berr(integer *n, integer *nz, integer *nrhs, T *res, T *ayb, T *berr)
    {
        la_lin_berr(n, nz, nrhs, res, ayb, berr);
    }
    template <typename T, typename Ta>
    void la_lin_berr(integer *n, integer *nz, integer *nrhs, T *res, Ta *ayb, Ta *berr)
    {
        la_lin_berr(n, nz, nrhs, res, ayb, berr);
    }
    /** @}*/ // end of la_lin_bierr
    /** @}*/ // end of Auxilliary
    /** @}*/ // end of LS
}
#endif