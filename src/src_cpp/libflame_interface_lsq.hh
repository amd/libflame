/******************************************************************************
 * Copyright (C) 2021-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface_lsq.hh
 *  libflame_interface.hh defines all the Least square routines for Libflame CPP templated public
 *  interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_LSQ_HH
#define LIBFLAME_INTERFACE_LSQ_HH

#include "libflame.hh"

namespace libflame
{
     /** @defgroup Least Least square
     * @ingroup LAPACK
     * @{
     */
    /** @defgroup std Standard Least square
     * @ingroup Least
     * @{
     */
    /** @defgroup gels gels
     * @ingroup std
     * @{
     */
    /*! @brief GELS solves overdetermined or underdetermined systems for GE matrices
* @details
* \b Purpose:
 \verbatim
    GELS solves overdetermined or underdetermined real linear systems
    involving an M-by-N matrix A, or its transpose, using a QR or LQ
    factorization of A.  It is assumed that A has full rank.

    The following options are provided:

    1. If TRANS = 'N' and m >= n:  find the least squares solution of
      an overdetermined system, i.e., solve the least squares problem
                   minimize || B - A*X ||.

    2. If TRANS = 'N' and m < n:  find the minimum norm solution of
      an underdetermined system A * X = B.

    3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
      an underdetermined system A**T * X = B.

    4. If TRANS = 'T' and m < n:  find the least squares solution of
      an overdetermined system, i.e., solve the least squares problem
                   minimize || B - A**T * X ||.

    Several right hand side vectors b and solution vectors x can be
    handled in a single call; they are stored as the columns of the
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution
    matrix X.
 \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          = 'N': the linear system involves A; \n
          = 'T': the linear system involves A**T. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of
          columns of the matrices B and X. NRHS >=0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, \n
            if M >= N, A is overwritten by details of its QR
                       factorization as returned by SGEQRF; \n
            if M <  N, A is overwritten by details of its LQ
                       factorization as returned by SGELQF. \n
 * @param[in] LDA
         LDA is INTEGER \n
         Computes row and column scaling to reduce condition number of matrix  \n
         The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the matrix B of right hand side vectors, stored
          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
          if TRANS = 'T'. \n
          On exit, if INFO = 0, B is overwritten by the solution
          vectors, stored columnwise: \n
          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
          squares solution vectors; the residual sum of squares for the
          solution in each column is given by the sum of squares of
          elements N+1 to M in that column; \n
          if TRANS = 'N' and m < n, rows 1 to N of B contain the
          minimum norm solution vectors; \n
          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
          minimum norm solution vectors; \n
          if TRANS = 'T' and m < n, rows 1 to M of B contain the
          least squares solution vectors; the residual sum of squares
          for the solution in each column is given by the sum of
          squares of elements M+1 to N in that column. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= MAX(1,M,N). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.
          LWORK >= fla_max( 1, MN + fla_max( MN, NRHS ) ). \n
          For optimal performance,
          LWORK >= fla_max( 1, MN + fla_max( MN, NRHS )*NB ).
          where MN = min(M,N) and NB is the optimum block size. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO =  i, the i-th diagonal element of the
                triangular factor of A is zero, so that A does not have
                full rank; the least squares solution could not be
                computed. \n

 *  * */
    template <typename T>
    void gels(char *trans, integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b,
              integer *ldb, T *work, integer *lwork, integer *info)
    {
        gels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
    }
    /** @}*/ // end of gels

    /** @defgroup gelss gelss
     * @ingroup std
     * @{
     */
    /*! @brief The minimum-norm solution to a real linear least squares problem
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of the minimum-norm solution to a real linear least squares problem:
        minimize 2-norm(| B - A*X |)
        using the singular value decomposition (SVD) of A. A is an m-by-n matrix which may be
        rank-deficient.

        Several right hand side vectors B and solution vectors X can be handled in a single call;
        they are stored as the columns of the m-by-nrhs  right hand side matrix b and the n-by-nrhs
        solution matrix X.

        The effective rank of A is determined by treating as zero those singular values which are
        less than rcond times the largest singular value.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in] nrhs
              nrhs is integer* \n
              The number of right hand sides, i.e., the number of columns
              of the matrices b and X. nrhs >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the m-by-n matrix. \n
              On exit, the first min(m,n) rows of A are overwritten with
              its right singular vectors, stored rowwise. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,m) \n
    * @param[in,out] b
              b is float/double array, dimension (ldb,nrhs) \n
              On entry, the m-by-nrhs  right hand side matrix b. \n
              On exit, b is overwritten by the n-by-nrhs solution
              matrix X.  If m >= n and rank = n, the residual
              sum-of-squares for the solution in the i-th column is given
              by the sum of squares of elements n+1:m in that column. \n
    * @param[in] ldb
              ldb is integer* \n
              The leading dimension of the matrix b, ldb >= fla_max(1,max(m,n)). \n
    * @param[out] s
              s is float/double array, dimension (min(m,n)) \n
              The singular values of A in decreasing order.
              The condition number of A in the 2-norm = S(1)/S(min(m,n)). \n
    * @param[in] rcond
              rcond is float/double* \n
              rcond is used to determine the effective rank of a.
              Singular values S(i) <= rcond*S(1) are treated as zero. \n
              If rcond < 0, machine precision is used instead. \n
    * @param[out] rank
              rank is integer* \n
              The effective rank of a, i.e., the number of singular values
              which are greater than rcond*S(1). \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= 1, and also: \n
              LWORK >=  2*min(M,N) + fla_max(M,N,NRHS) \n
              For good performance, LWORK should generally be larger. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	RWORK
              RWORK is REAL array, dimension (5*min(M,N)) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n
              > 0:  the algorithm for computing the SVD failed to converge;
                    if INFO = i, i off-diagonal elements of an intermediate
                    bidiagonal form did not converge to zero. \n

    *     *  */
    template <typename T>
    void gelss(integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b, integer *ldb, T *s,
               T *rcond, integer *rank, T *work, integer *lwork, integer *info)
    {
        gelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
    }
    template <typename T, typename Ta>
    void gelss(integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b, integer *ldb, Ta *s,
               Ta *rcond, integer *rank, T *work, integer *lwork, Ta *rwork, integer *info)
    {
        gelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info);
    }
    /** @} */ // end of gelss

    /** @defgroup gelsd gelsd
     * @ingroup std
     * @{
     */
    /*! @brief The minimum-norm solution to a real linear least squares problem
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of the minimum-norm solution to a real linear least squares problem:
        minimize 2-norm(| B - A*X |)
        using the singular value decomposition (SVD) of A. A is an m-by-n matrix which may be
        rank-deficient.

        Several right hand side vectors B and solution vectors X can be handled in a single call;
        they are stored as the columns of the m-by-nrhs  right hand side matrix b and the n-by-nrhs
        solution matrix X.

        The problem is solved in three steps:
        (1) Reduce the coefficient matrix a to bidiagonal form with Householder transformations,
        reducing the original problem into a "bidiagonal least squares problem" (BLS)
        (2) Solve the BLS using a divide and conquer approach.
        (3) Apply back all the Householder transformations to solve the original least squares
        problem.

        The effective rank of A is determined by treating as zero those singular values which are
        less than rcond times the largest singular value.

        The divide and conquer algorithm makes very mild assumptions about floating point
        arithmetic. It will work on machines with a guard digit in add/subtract, or on those
        binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP,
        Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
        guard digits, but we know of none.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in] nrhs
              nrhs is integer* \n
              The number of right hand sides, i.e., the number of columns
              of the matrices b and X. nrhs >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the m-by-n matrix. \n
              On exit, A has been destroyed. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the matrix a, lda >= fla_max(1,m) \n
    * @param[in,out] b
              b is float/double array, dimension (ldb,nrhs) \n
              On entry, the m-by-n matrix. \n
              On exit, b is overwritten by the n-by-nrhs solution
              matrix X.  If m >= n and rank = n, the residual
              sum-of-squares for the solution in the i-th column is given
              by the sum of squares of elements n+1:m in that column. \n
    * @param[in] ldb
              ldb is integer* \n
              The leading dimension of the matrix b, ldb >= fla_max(1,max(m,n)). \n
    * @param[out] s
              s is float/double array, dimension (min(m,n)) \n
              The singular values of A in decreasing order.
              The condition number of A in the 2-norm = S(1)/S(min(m,n)). \n
    * @param[in] rcond
              rcond is float/double* \n
              rcond is used to determine the effective rank of a.
              Singular values S(i) <= rcond*S(1) are treated as zero.
              If rcond < 0, machine precision is used instead. \n
    * @param[out] rank
              rank is integer* \n
              The effective rank of a, i.e., the number of singular values
              which are greater than rcond*S(1). \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK must be at least 1. \n
              The exact minimum amount of workspace needed depends on M,
              N and NRHS. As long as LWORK is at least \n
                  2 * N + N * NRHS \n
              if M is greater than or equal to N or \n
                  2 * M + M * NRHS \n
              if M is less than N, the code will execute correctly. \n
              For good performance, LWORK should generally be larger. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the array WORK and the
              minimum sizes of the arrays RWORK and IWORK, and returns
              these values as the first entries of the WORK, RWORK and
              IWORK arrays, and no error message related to LWORK is issued
              by XERBLA. \n
    * @param[out]	RWORK
              RWORK is REAL array, dimension (MAX(1,LRWORK)) \n
              LRWORK >= \n
                 10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + \n
                 MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ) \n
              if M is greater than or equal to N or \n
                 10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS + \n
                 MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS ) \n
              if M is less than N, the code will execute correctly. \n
              SMLSIZ is returned by ILAENV and is equal to the maximum
              size of the subproblems at the bottom of the computation
              tree (usually about 25), and
                 NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 ) \n
              On exit, if INFO = 0, RWORK(1) returns the minimum LRWORK. \n
    * @param[out]	IWORK
              IWORK is INTEGER array, dimension (MAX(1,LIWORK)) \n
              LIWORK >= fla_max(1, 3*MINMN*NLVL + 11*MINMN),
              where MINMN = MIN( M,N ). \n
              On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -i, the i-th argument had an illegal value. \n
              > 0:  the algorithm for computing the SVD failed to converge;
                    if INFO = i, i off-diagonal elements of an intermediate
                    bidiagonal form did not converge to zero. \n

    *     *  */
    template <typename T>
    void gelsd(integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b, integer *ldb, T *s,
               T *rcond, integer *rank, T *work, integer *lwork, integer *iwork, integer *info)
    {
        gelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void gelsd(integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b, integer *ldb, Ta *s,
               Ta *rcond, integer *rank, T *work, integer *lwork, Ta *rwork, integer *iwork,
               integer *info)
    {
        gelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info);
    }
    /** @} */ // end of gelsd

    /** @defgroup gelsy gelsy
     * @ingroup std
     * @{
     */
    /*! @brief GELSY solves overdetermined or underdetermined systems for GE matrices
* @details
* \b Purpose:
 \verbatim
    GELSY computes the minimum-norm solution to a real linear least
    squares problem:
       minimize || A * X - B ||
    using a complete orthogonal factorization of A.  A is an M-by-N
    matrix which may be rank-deficient.

    Several right hand side vectors b and solution vectors x can be
    handled in a single call; they are stored as the columns of the
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution
    matrix X.

    The routine first computes a QR factorization with column pivoting:
       A * P = Q * [ R11 R12 ]
                   [  0  R22 ]
    with R11 defined as the largest leading submatrix whose estimated
    condition number is less than 1/RCOND.  The order of R11, RANK,
    is the effective rank of A.

    Then, R22 is considered to be negligible, and R12 is annihilated
    by orthogonal transformations from the right, arriving at the
    complete orthogonal factorization:
      A * P = Q * [ T11 0 ] * Z
                  [  0  0 ]
    The minimum-norm solution is then
      X = P * Z**T [ inv(T11)*Q1**T*B ]
                   [        0         ]
    where Q1 consists of the first RANK columns of Q.

    This routine is basically identical to the original xGELSX except
    three differences:
     o The call to the subroutine xGEQPF has been substituted by the
       the call to the subroutine xGEQP3. This subroutine is a Blas-3
       version of the QR factorization with column pivoting.
     o Matrix B (the right hand side) is updated with Blas-3.
     o The permutation of matrix B (the right hand side) is faster and
       more simple.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of
          columns of matrices B and X. NRHS >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, A has been overwritten by details of its
          complete orthogonal factorization. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the M-by-NRHS right hand side matrix B. \n
          On exit, the N-by-NRHS solution matrix X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          Computes row and column scaling to reduce condition number of matrix B \n
          The leading dimension of the array B. LDB >= fla_max(1,M,N). \n
 * @param[in,out] JPVT
          JPVT is INTEGER array, dimension (N) \n
          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
          to the front of AP, otherwise column i is a free column. \n
          On exit, if JPVT(i) = k, then the i-th column of AP
          was the k-th column of A. \n
 * @param[in] RCOND
          RCOND is REAL \n
          RCOND is used to determine the effective rank of A, which
          is defined as the order of the largest leading triangular
          submatrix R11 in the QR factorization with pivoting of A,
          whose estimated condition number < 1/RCOND. \n
 * @param[out] RANK
          RANK is INTEGER \n
          The effective rank of A, i.e., the order of the submatrix
          R11.  This is the same as the order of the submatrix T11
          in the complete orthogonal factorization of A. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK.
          The unblocked strategy requires that:
             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
          where MN = min( M, N ). \n
          The block algorithm requires that:
             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
          where NB is an upper bound on the blocksize returned
          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR,
          and SORMRZ. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: If INFO = -i, the i-th argument had an illegal value. \n

 *  * */
    template <typename T>
    void gelsy(integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b, integer *ldb,
               integer *jpvt, T rcond, integer *rank, T *work, integer *lwork, integer *info)
    {
        gelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info);
    }
    template <typename T, typename Ta>
    void gelsy(integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b, integer *ldb,
               integer *jpvt, Ta rcond, integer *rank, T *work, integer *lwork, Ta *rwork,
               integer *info)
    {
        gelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info);
    }
    /** @} */ // end of gelsy

    /** @defgroup getsls getsls
     * @ingroup std
     * @{
     */
    /*! @brief GETSLS solves overdetermined or underdetermined real linear systems
 * @details
 * \b Purpose:
    \verbatim
     GETSLS solves overdetermined or underdetermined real linear systems
     involving an M-by-N matrix A, using a tall skinny QR or short wide LQ
     factorization of A.  It is assumed that A has full rank.

     The following options are provided:

     1. If TRANS = 'N' and m >= n:  find the least squares solution of
        an overdetermined system, i.e., solve the least squares problem
                     minimize || B - A*X ||.

     2. If TRANS = 'N' and m < n:  find the minimum norm solution of
        an underdetermined system A * X = B.

     3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
        an undetermined system A**T * X = B.

     4. If TRANS = 'T' and m < n:  find the least squares solution of
        an overdetermined system, i.e., solve the least squares problem
                     minimize || B - A**T * X ||.

     Several right hand side vectors b and solution vectors x can be
     handled in a single call; they are stored as the columns of the
     M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     matrix X.
    \endverbatim

  * @param[in] TRANS
           TRANS is CHARACTER*1 \n
           = 'N': the linear system involves A; \n
           = 'T': the linear system involves A**T. \n
  * @param[in] M
           M is INTEGER \n
           The number of rows of the matrix A.  M >= 0. \n
  * @param[in] N
           N is INTEGER \n
           The number of columns of the matrix A.  N >= 0. \n
  * @param[in] NRHS
           NRHS is INTEGER \n
           The number of right hand sides, i.e., the number of
           columns of the matrices B and X. NRHS >=0. \n
  * @param[in,out] A
           A is REAL array, dimension (LDA,N) \n
           On entry, the M-by-N matrix A. \n
           On exit,
           A is overwritten by details of its QR or LQ
           factorization as returned by SGEQR or SGELQ. \n
  * @param[in] LDA
           LDA is INTEGER \n
           The leading dimension of the array A.  LDA >= fla_max(1,M). \n
  * @param[in,out] B
           B is REAL array, dimension (LDB,NRHS) \n
           On entry, the matrix B of right hand side vectors, stored
           columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
           if TRANS = 'T'. \n
           On exit, if INFO = 0, B is overwritten by the solution
           vectors, stored columnwise: \n
           if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
           squares solution vectors. \n
           if TRANS = 'N' and m < n, rows 1 to N of B contain the
           minimum norm solution vectors; \n
           if TRANS = 'T' and m >= n, rows 1 to M of B contain the
           minimum norm solution vectors; \n
           if TRANS = 'T' and m < n, rows 1 to M of B contain the
           least squares solution vectors. \n
  * @param[in] LDB
           LDB is INTEGER \n
           The leading dimension of the array B. LDB >= MAX(1,M,N). \n
  * @param[out]	WORK
          (workspace) COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) contains optimal (or either minimal
          or optimal, if query was assumed) LWORK.
          See LWORK for details. \n
  * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          If LWORK = -1 or -2, then a workspace query is assumed. \n
          If LWORK = -1, the routine calculates optimal size of WORK for the
          optimal performance and returns this value in WORK(1). \n
          If LWORK = -2, the routine calculates minimal size of WORK and
          returns this value in WORK(1). \n
  * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO =  i, the i-th diagonal element of the
                triangular factor of A is zero, so that A does not have
                full rank; the least squares solution could not be
                computed. \n

 *  * */
    template <typename T>
    void getsls(char *trans, integer *m, integer *n, integer *nrhs, T *a, integer *lda, T *b,
                integer *ldb, T *work, integer *lwork, integer *info)
    {
        getsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
    }
    /** @} */ // end of getsls
    /** @} */ // end of std

    /** @defgroup cnstr Constrained least squares
     *@ingroup Least
     * @{
     */

    /** @defgroup gglse gglse
     * @ingroup cnstr
     * @{
     */
    /*! @brief GGLSE solves overdetermined or underdetermined systems for OTHER matrices
 * @details
 * \b Purpose:
    \verbatim
     GGLSE solves the linear equality-constrained least squares (LSE)
     problem:
             minimize || c - A*x ||_2   subject to   B*x = d
     where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
     M-vector, and d is a given P-vector. It is assumed that
     P <= N <= M+P, and
              rank(B) = P and  rank( (A)) = N.
                                   ( (B))
     These conditions ensure that the LSE problem has a unique solution,
     which is obtained using a generalized RQ factorization of the
     matrices (B, A) given by
        B = (0 R)*Q,   A = Z*T*Q.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrices A and B. N >= 0. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows of the matrix B. 0 <= P <= N <= M+P. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, the elements on and above the diagonal of the array
          contain the min(M,N)-by-N upper trapezoidal matrix T. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the P-by-N matrix B. \n
          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)
          contains the P-by-P upper triangular matrix R. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,P). \n
 * @param[in,out] C
          C is REAL array, dimension (M) \n
          On entry, C contains the right hand side vector for the
          least squares part of the LSE problem. \n
          On exit, the residual sum of squares for the solution
          is given by the sum of squares of elements N-P+1 to M of
          vector C. \n
 * @param[in,out] D
          D is REAL array, dimension (P) \n
          On entry, D contains the right hand side vector for the
          constrained equation. \n
          On exit, D is destroyed. \n
 * @param[out] X
          X is REAL array, dimension (N) \n
          On exit, X is the solution of the LSE problem. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,M+N+P). \n
          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,
          where NB is an upper bound for the optimal blocksizes for
          SGEQRF, SGERQF, SORMQR and SORMRQ. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1:  the upper triangular factor R associated with B in the
                generalized RQ factorization of the pair (B, A) is
                singular, so that rank(B) < P; the least squares
                solution could not be computed. \n
          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor
                T associated with A in the generalized RQ factorization
                of the pair (B, A) is singular, so that
                rank( (A) ) < N; the least squares solution could not
                    ( (B) )
                be computed. \n

 *  * */
    template <typename T>
    void gglse(integer *m, integer *n, integer *p, T *a, integer *lda, T *b, integer *ldb, T *c,
               T *d, T *x, T *work, integer *lwork, integer *info)
    {
        gglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info);
    }
    /** @}*/ // end of gglse

    /** @defgroup ggglm ggglm
     * @ingroup cnstr
     * @{
     */

    /*! @brief GGGLM solves a general Gauss-Markov linear model (GLM) problem
 * @details
 * \b Purpose:
    \verbatim
     GGGLM solves a general Gauss-Markov linear model (GLM) problem:
             minimize || y ||_2   subject to   d = A*x + B*y
                 x
     where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
     given N-vector. It is assumed that M <= N <= M+P, and
                rank(A) = M    and    rank( A B) = N.
     Under these assumptions, the constrained equation is always
     consistent, and there is a unique solution x and a minimal 2-norm
     solution y, which is obtained using a generalized QR factorization
     of the matrices (A, B) given by
        A = Q*(R),   B = Q*T*Z.
              (0)
     In particular, if matrix B is square nonsingular, then the problem
     GLM is equivalent to the following weighted linear least squares
     problem
                  minimize || inv(B)*(d-A*x) ||_2
                      x
     where inv(B) denotes the inverse of B.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of rows of the matrices A and B.  N >= 0. \n
 * @param[in] M
          M is INTEGER \n
          The number of columns of the matrix A.  0 <= M <= N. \n
 * @param[in] P
          P is INTEGER \n
          The number of columns of the matrix B.  P >= N-M. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,M) \n
          On entry, the N-by-M matrix A. \n
          On exit, the upper triangular part of the array A contains
          the M-by-M upper triangular matrix R. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,N). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,P) \n
          On entry, the N-by-P matrix B. \n
          On exit, if N <= P, the upper triangle of the subarray
          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T;
          if N > P, the elements on and above the (N-P)th subdiagonal
          contain the N-by-P upper trapezoidal matrix T. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,N). \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, D is the left hand side of the GLM equation.
          On exit, D is destroyed. \n
 * @param[out] X
          X is REAL array, dimension (M) \n
 * @param[out] Y
          Y is REAL array, dimension (P) \n
          On exit, X and Y are the solutions of the GLM problem. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. LWORK >= fla_max(1,N+M+P).
          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,
          where NB is an upper bound for the optimal blocksizes for
          SGEQRF, SGERQF, SORMQR and SORMRQ. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value.
          = 1:  the upper triangular factor R associated with A in the
                generalized QR factorization of the pair (A, B) is
                singular, so that rank(A) < M; the least squares
                solution could not be computed. \n
          = 2:  the bottom (N-M) by (N-M) part of the upper trapezoidal
                factor T associated with B in the generalized QR
                factorization of the pair (A, B) is singular, so that
                rank( A B ) < N; the least squares solution could not
                be computed. \n

 *  * */
    template <typename T>
    void ggglm(integer *n, integer *m, integer *p, T *a, integer *lda, T *b, integer *ldb, T *d,
               T *x, T *y, T *work, integer *lwork, integer *info)
    {
        ggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
    }
    /** @}*/ // end of ggglm
    /** @}*/ // end of cnstr

    /** @defgroup Aux Auxiliary routines
     * @ingroup Least
     * @{
     */

    /** @defgroup laic1 laic1
     *@ingroup Aux
     * @{
     */
    /*! @brief LAIC1 applies one step of incremental condition estimation

 * @details
 * \b Purpose:
    \verbatim
    LAIC1 applies one step of incremental condition estimation in
    its simplest version:

    Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
    lower triangular matrix L, such that
             twonorm(L*x) = sest
    Then SLAIC1 computes sestpr, s, c such that
    the vector
                    [ s*x ]
             xhat = [  c  ]
    is an approximate singular vector of
                    [ L      0  ]
             Lhat = [ w**T gamma ]
    in the sense that
             twonorm(Lhat*xhat) = sestpr.

    Depending on JOB, an estimate for the largest or smallest singular
    value is computed.

    Note that [s c]**T and sestpr**2 is an eigenpair of the system

        diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]
                                              [ gamma ]

    where  alpha =  x**T*w.
    \endverbatim

 * @param[in] JOB
          JOB is INTEGER \n
          = 1: an estimate for the largest singular value is computed. \n
          = 2: an estimate for the smallest singular value is computed. \n
 * @param[in] J
          J is INTEGER \n
          Length of X and W \n
 * @param[in] X
          X is REAL array, dimension (J) \n
          The j-vector x. \n
 * @param[in] SEST
          SEST is REAL \n
          Estimated singular value of j by j matrix L \n
 * @param[in] W
          W is REAL array, dimension (J) \n
          The j-vector w. \n
 * @param[in] GAMMA
          GAMMA is REAL \n
          The diagonal element gamma. \n
 * @param[out] SESTPR
          SESTPR is REAL \n
          Estimated singular value of (j+1) by (j+1) matrix Lhat. \n
 * @param[out] S
          S is REAL \n
          Sine needed in forming xhat. \n
 * @param[out] C
          C is REAL \n
          Cosine needed in forming xhat. \n

 * */
    template <typename T>
    void laic1(integer *job, integer *j, T *x, T *sest, T *w, T *gamma, T *sestpr, T *s, T *c__)
    {
        laic1(job, j, x, sest, w, gamma, sestpr, s, c__);
    }
    template <typename T, typename Ta>
    void laic1(integer *job, integer *j, T *x, Ta *sest, T *w, T *gamma, Ta *sestpr, T *s, T *c__)
    {
        laic1(job, j, x, sest, w, gamma, sestpr, s, c__);
    }
    /** @}*/ // end of laic1

    /** @defgroup lals0 lals0
     *@ingroup Aux
     * @{
     */

    /*! @brief LALS0 applies back multiplying factors in solving the least squares problem \n
     using divide and conquer SVD approach. Used by sgelsd
 * @details
 * \b Purpose:
    \verbatim
    LALS0 applies back the multiplying factors of either the left or the
    right singular vector matrix of a diagonal matrix appended by a row
    to the right hand side matrix B in solving the least squares problem
    using the divide-and-conquer SVD approach.

    For the left singular vector matrix, three types of orthogonal
    matrices are involved:

    (1L) Givens rotations: the number of such rotations is GIVPTR; the
         pairs of columns/rows they were applied to are stored in GIVCOL;
         and the C- and S-values of these rotations are stored in GIVNUM.

    (2L) Permutation. The (NL+1)-st row of B is to be moved to the first
         row, and for J=2:N, PERM(J)-th row of B is to be moved to the
         J-th row.

    (3L) The left singular vector matrix of the remaining matrix.

    For the right singular vector matrix, four types of orthogonal
    matrices are involved:

    (1R) The right singular vector matrix of the remaining matrix.

    (2R) If SQRE = 1, one extra Givens rotation to generate the right
         null space.

    (3R) The inverse transformation of (2L).

    (4R) The inverse transformation of (1L).
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          Specifies whether singular vectors are to be computed in
          factored form: \n
          = 0: Left singular vector matrix. \n
          = 1: Right singular vector matrix. \n
 * @param[in] NL
          NL is INTEGER \n
          The row dimension of the upper block. NL >= 1. \n
 * @param[in] NR
          NR is INTEGER \n
          The row dimension of the lower block. NR >= 1. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          = 0: the lower block is an NR-by-NR square matrix. \n
          = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
 \n
          The bidiagonal matrix has row dimension N = NL + NR + 1,
          and column dimension M = N + SQRE. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of columns of B and BX. NRHS must be at least 1. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, NRHS) \n
          On input, B contains the right hand sides of the least
          squares problem in rows 1 through M. On output, B contains
          the solution X in rows 1 through N. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B. LDB must be at least
          fla_max(1,MAX(M, N)). \n
 * @param[out] BX
          BX is REAL array, dimension (LDBX, NRHS) \n
 * @param[in] LDBX
          LDBX is INTEGER \n
          The leading dimension of BX. \n
 * @param[in] PERM
          PERM is INTEGER array, dimension (N) \n
          The permutations (from deflation and sorting) applied
          to the two blocks. \n
 * @param[in] GIVPTR
          GIVPTR is INTEGER \n
          The number of Givens rotations which took place in this
          subproblem. \n
 * @param[in] GIVCOL
          GIVCOL is INTEGER array, dimension (LDGCOL, 2) \n
          Each pair of numbers indicates a pair of rows/columns
          involved in a Givens rotation. \n
 * @param[in] LDGCOL
          LDGCOL is INTEGER \n
          The leading dimension of GIVCOL, must be at least N. \n
 * @param[in] GIVNUM
          GIVNUM is REAL array, dimension (LDGNUM, 2) \n
          Each number indicates the C or S value used in the
          corresponding Givens rotation. \n
 * @param[in] LDGNUM
          LDGNUM is INTEGER \n
          The leading dimension of arrays DIFR, POLES and
          GIVNUM, must be at least K. \n
 * @param[in] POLES
          POLES is REAL array, dimension (LDGNUM, 2) \n
          On entry, POLES(1:K, 1) contains the new singular
          values obtained from solving the secular equation, and
          POLES(1:K, 2) is an array containing the poles in the secular
          equation. \n
 * @param[in] DIFL
          DIFL is REAL array, dimension (K). \n
          On entry, DIFL(I) is the distance between I-th updated
          (undeflated) singular value and the I-th (undeflated) old
          singular value. \n
 * @param[in] DIFR
          DIFR is REAL array, dimension (LDGNUM, 2). \n
          On entry, DIFR(I, 1) contains the distances between I-th
          updated (undeflated) singular value and the I+1-th
          (undeflated) old singular value. And DIFR(I, 2) is the
          normalizing factor for the I-th right singular vector. \n
 * @param[in] Z
          Z is REAL array, dimension (K) \n
          Contain the components of the deflation-adjusted updating row
          vector. \n
 * @param[in] K
          K is INTEGER \n
          Contains the dimension of the non-deflated matrix,
          This is the order of the related secular equation. 1 <= K <=N. \n
 * @param[in] C
          C is REAL \n
          C contains garbage if SQRE =0 and the C-value of a Givens
          rotation related to the right null space if SQRE = 1. \n
 * @param[in] S
          S is REAL \n
          S contains garbage if SQRE =0 and the S-value of a Givens
          rotation related to the right null space if SQRE = 1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (K) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value.  \n

 * */
    template <typename T>
    void lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, T *b,
               integer *ldb, T *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
               integer *ldgcol, T *givnum, integer *ldgnum, T *poles, T *difl, T *difr, T *z,
               integer *k, T *c, T *s, T *work, integer *info)
    {
        lals0(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum,
              ldgnum, poles, difl, difr, z, k, c, s, work, info);
    }
    template <typename T, typename Ta>
    void lals0(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *nrhs, T *b,
               integer *ldb, T *bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
               integer *ldgcol, Ta *givnum, integer *ldgnum, Ta *poles, Ta *difl, Ta *difr, Ta *z,
               integer *k, Ta *c, Ta *s, Ta *work, integer *info)
    {
        lals0(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol, givnum,
              ldgnum, poles, difl, difr, z, k, c, s, work, info);
    }
    /** @}*/ // end of lals0

    /** @defgroup lalsa lalsa
     *@ingroup Aux
     * @{
     */
    /*! @brief LALSA computes the SVD of the coefficient matrix in compact form. Used by sgelsd

 * @details
 * \b Purpose:
    \verbatim
    LALSA is an itermediate step in solving the least squares problem
    by computing the SVD of the coefficient matrix in compact form (The
    singular vectors are computed as products of simple orthorgonal
    matrices.).

    If ICOMPQ = 0, SLALSA applies the inverse of the left singular vector
    matrix of an upper bidiagonal matrix to the right hand side; and if
    ICOMPQ = 1, SLALSA applies the right singular vector matrix to the
    right hand side. The singular vector matrices were generated in
    compact form by SLALSA.
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          Specifies whether the left or the right singular vector
          matrix is involved. \n
          = 0: Left singular vector matrix \n
          = 1: Right singular vector matrix \n
 * @param[in] SMLSIZ
          SMLSIZ is INTEGER \n
          The maximum size of the subproblems at the bottom of the
          computation tree. \n
 * @param[in] N
          N is INTEGER \n
          The row and column dimensions of the upper bidiagonal matrix. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of columns of B and BX. NRHS must be at least 1. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB, NRHS) \n
          On input, B contains the right hand sides of the least
          squares problem in rows 1 through M. \n
          On output, B contains the solution X in rows 1 through N. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B in the calling subprogram.
          LDB must be at least fla_max(1,MAX(M, N)). \n
 * @param[out] BX
          BX is REAL array, dimension (LDBX, NRHS) \n
          On exit, the result of applying the left or right singular
          vector matrix to B. \n
 * @param[in] LDBX
          LDBX is INTEGER \n
          The leading dimension of BX. \n
 * @param[in] U
          U is REAL array, dimension (LDU, SMLSIZ). \n
          On entry, U contains the left singular vector matrices of all
          subproblems at the bottom level. \n
 * @param[in] LDU
          LDU is INTEGER, LDU = > N. \n
         The leading dimension of arrays U, VT, DIFL, DIFR,
         POLES, GIVNUM, and Z. \n
 * @param[in] VT
          VT is REAL array, dimension (LDU, SMLSIZ+1). \n
         On entry, VT**T contains the right singular vector matrices of
         all subproblems at the bottom level. \n
 * @param[in] K
          K is INTEGER array, dimension (N). \n
 * @param[in] DIFL
          DIFL is REAL array, dimension (LDU, NLVL). \n
         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1. \n
 * @param[in] DIFR
          DIFR is REAL array, dimension (LDU, 2 * NLVL). \n
          On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record
          distances between singular values on the I-th level and
          singular values on the (I -1)-th level, and DIFR(*, 2 * I)
          record the normalizing factors of the right singular vectors
          matrices of subproblems on I-th level. \n
 * @param[in] Z
          Z is REAL array, dimension (LDU, NLVL). \n
          On entry, Z(1, I) contains the components of the deflation-
          adjusted updating row vector for subproblems on the I-th
          level. \n
 * @param[in] POLES
          POLES is REAL array, dimension (LDU, 2 * NLVL). \n
          On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old
          singular values involved in the secular equations on the I-th
          level. \n
 * @param[in] GIVPTR
          GIVPTR is INTEGER array, dimension (N). \n
          On entry, GIVPTR(I) records the number of Givens
          rotations performed on the I-th problem on the computation
          tree. \n
 * @param[in] GIVCOL
          GIVCOL is INTEGER array, dimension (LDGCOL, 2 * NLVL). \n
          On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the
          locations of Givens rotations performed on the I-th level on
          the computation tree. \n
 * @param[in] LDGCOL
          LDGCOL is INTEGER, LDGCOL = > N. \n
          The leading dimension of arrays GIVCOL and PERM. \n
 * @param[in] PERM
          PERM is INTEGER array, dimension (LDGCOL, NLVL). \n
          On entry, PERM(*, I) records permutations done on the I-th
          level of the computation tree. \n
 * @param[in] GIVNUM
          GIVNUM is REAL array, dimension (LDU, 2 * NLVL). \n
          On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-
          values of Givens rotations performed on the I-th level on the
          computation tree. \n
 * @param[in] C
          C is REAL array, dimension (N). \n
          On entry, if the I-th subproblem is not square,
          C(I) contains the C-value of a Givens rotation related to
          the right null space of the I-th subproblem. \n
 * @param[in] S
          S is REAL array, dimension (N). \n
          On entry, if the I-th subproblem is not square,
          S(I) contains the S-value of a Givens rotation related to
          the right null space of the I-th subproblem. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (3*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 * */
    template <typename T>
    void lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, T *b, integer *ldb,
               T *bx, integer *ldbx, T *u, integer *ldu, T *vt, integer *k, T *difl, T *difr, T *z,
               T *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm,
               T *givnum, T *c, T *s, T *work, integer *iwork, integer *info)
    {
        lalsa(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles,
              givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
    }
    template <typename T, typename Ta>
    void lalsa(integer *icompq, integer *smlsiz, integer *n, integer *nrhs, T *b, integer *ldb,
               T *bx, integer *ldbx, Ta *u, integer *ldu, Ta *vt, integer *k, Ta *difl, Ta *difr,
               Ta *z, Ta *poles, integer *givptr, integer *givcol, integer *ldgcol, integer *perm,
               Ta *givnum, Ta *c, Ta *s, Ta *work, integer *iwork, integer *info)
    {
        lalsa(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles,
              givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info);
    }
    /** @}*/ // end of lalsa

    /** @defgroup lalsd lalsd
     *@ingroup Aux
     * @{
     */
    /*! @brief LALSD uses the singular value decomposition of A to solve the least squares problem

 * @details
 * \b Purpose:
    \verbatim
    LALSD uses the singular value decomposition of A to solve the least
    squares problem of finding X to minimize the Euclidean norm of each
    column of A*X-B, where A is N-by-N upper bidiagonal, and X and B
    are N-by-NRHS. The solution X overwrites B.

    The singular values of A smaller than RCOND times the largest
    singular value are treated as zero in solving the least squares
    problem; in this case a minimum norm solution is   returned.
    The actual singular values are   returned in D in ascending order.

    This code makes very mild assumptions about floating point
    arithmetic. It will work on machines with a guard digit in
    add/subtract, or on those binary machines without guard digits
    which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
    It could conceivably fail on hexadecimal or decimal machines
    without guard digits, but we know of none.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U': D and E define an upper bidiagonal matrix. \n
          = 'L': D and E define a  lower bidiagonal matrix. \n
 * @param[in] SMLSIZ
          SMLSIZ is INTEGER \n
          The maximum size of the subproblems at the bottom of the
          computation tree. \n
 * @param[in] N
          N is INTEGER \n
          The dimension of the  bidiagonal matrix.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of columns of B. NRHS must be at least 1. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry D contains the main diagonal of the bidiagonal
          matrix. On exit, if INFO = 0, D contains its singular values. \n
 * @param[in,out] E
          E is REAL array, dimension (N-1) \n
          Contains the super-diagonal entries of the bidiagonal matrix.
          On exit, E has been destroyed. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On input, B contains the right hand sides of the least
          squares problem. On output, B contains the solution X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of B in the calling subprogram.
          LDB must be at least fla_max(1,N). \n
 * @param[in] RCOND
          RCOND is REAL \n
          The singular values of A less than or equal to RCOND times
          the largest singular value are treated as zero in solving
          the least squares problem. If RCOND is negative,
          machine precision is used instead. \n
          For example, if diag(S)*X=B were the least squares problem,
          where diag(S) is a diagonal matrix of singular values, the
          solution would be X(i) = B(i) / S(i) if S(i) is greater than
          RCOND*max(S), and X(i) = 0 if S(i) is less than or equal to
          RCOND*max(S). \n
 * @param[out] RANK
          RANK is INTEGER \n
          The number of singular values of A greater than RCOND times
          the largest singular value. \n
 * @param[out] WORK
          WORK is REAL array, dimension at least \n
          (9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2), \n
          where NLVL = fla_max(0, INT(log_2 (N/(SMLSIZ+1))) + 1). \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension at least
          (3*N*NLVL + 11*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  The algorithm failed to compute a singular value while
                working on the submatrix lying in rows and columns
                INFO/(N+1) through MOD(INFO,N+1). \n

 * */
    template <typename T>
    void lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, T *d, T *e, T *b,
               integer *ldb, T *rcond, integer *rank, T *work, integer *iwork, integer *info)
    {
        lalsd(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, iwork, info);
    }
    template <typename T, typename Ta>
    void lalsd(char *uplo, integer *smlsiz, integer *n, integer *nrhs, Ta *d, Ta *e, T *b,
               integer *ldb, Ta *rcond, integer *rank, T *work, Ta *rwork, integer *iwork,
               integer *info)
    {
        lalsd(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, rwork, iwork, info);
    }
    /** @}*/ // end of lalsd
    /** @}*/ // end of Aux
    /**@} */ // end of least square
}
#endif