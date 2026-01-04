/******************************************************************************
 * Copyright (C) 2021-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface_svd.hh
 *  libflame_interface.hh defines all the svd routines for Libflame CPP templated public
 *  interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_SVD_HH
#define LIBFLAME_INTERFACE_SVD_HH


namespace libflame
{
     /** @defgroup Singular Singular Value Decomposition (SVD)
     * @ingroup LAPACK
     * @{
     */
     /** @defgroup SVD SVD Computational Routines
     * @ingroup Singular
     * @{
     */

    /** @defgroup gesvd gesvd
     * @ingroup SVD
     * @{
     */
    /*! @brief General matrix singular value decomposition (QR algorithm)
   *
   * @details
   * \b Purpose:
     \verbatim
     Computation of singular value decomposition (SVD) of a real m-by-n matrix a, optionally
     computing the left and/or right singular vectors. The SVD is written

     A = U * SIGMA * transpose(V)

     where SIGMA is an m-by-n matrix which is zero except for its min(m,n) diagonal elements,
     U is an M-by-M orthogonal matrix, and V is an n-by-n orthogonal matrix.  The diagonal
     elements of SIGMA are the singular values of A; they are real and non-negative, and are
     returned in descending order. The first min(m,n) columns of U and V are the left and right
     singular vectors of A.

     Note that the routine returns V**T, not V.
    \endverbatim

    * @param[in] jobu
              jobu is char* \n
              Specifies options for computing all or part of the matrix U: \n
              = 'A':  all M columns of U are returned in array U: \n
              = 'S':  the first min(m,n) columns of U (the left singular
              vectors) are returned in the array U; \n
              = 'O':  the first min(m,n) columns of U (the left singular
              vectors) are overwritten on the array a; \n
              = 'N':  no columns of U (no left singular vectors) are
              computed. \n
    * @param[in] jobv
              jobv is char* \n
              Specifies options for computing all or part of the matrix
              V**T: \n
              = 'A':  all N rows of V**T are returned in the array vt; \n
              = 'S':  the first min(m,n) rows of V**T (the right singular
              vectors) are returned in the array vt; \n
              = 'O':  the first min(m,n) rows of V**T (the right singular
              vectors) are overwritten on the array a; \n
              = 'N':  no rows of V**T (no right singular vectors) are
              computed. \n \n
              jobv and jobu cannot both be 'O'. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the input matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the input matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the m-by-n matrix a. \n
              On exit, \n
              if jobu = 'O',  A is overwritten with the first min(m,n)
              columns of U (the left singular vectors,
              stored columnwise); \n
              if jobv = 'O', A is overwritten with the first min(m,n)
              rows of V**T (the right singular vectors,
              stored rowwise); \n
              if jobu .ne. 'O' and jobv .ne. 'O', the contents of A
              are destroyed. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,m). \n
    * @param[out] s
              s is float/double array, dimension (min(m,n)) \n
              The singular values of A, sorted so that S(i) >= S(i+1). \n
    * @param[out] u
              U is float/double array, dimension (ldu,UCOL) \n
              (ldu,M) if jobu = 'A' or (ldu,min(m,n)) if jobu = 'S'. \n
              If jobu = 'A', U contains the m-by-m orthogonal matrix U; \n
              if jobu = 'S', U contains the first min(m,n) columns of U
              (the left singular vectors, stored columnwise); \n
              if jobu = 'N' or 'O', U is not referenced. \n
    * @param[in] ldu
              ldu is integer* \n
              The leading dimension of the array U.  ldu >= 1; if
              jobu = 'S' or 'A', ldu >= M. \n
    * @param[out] vt
              vt is float/double array, dimension (ldvt,N) \n
              If jobv = 'A', vt contains the n-by-n orthogonal matrix
              V**T; \n
              if jobv = 'S', vt contains the first min(m,n) rows of
              V**T (the right singular vectors, stored rowwise); \n
              if jobv = 'N' or 'O', vt is not referenced. \n
    * @param[in] ldvt
              ldvt is integer* \n
              The leading dimension of the array vt.  ldvt >= 1; if
              jobv = 'A', ldvt >= N; if jobv = 'S', ldvt >= min(m,n). \n
    * @param[out] superb
              superb is float/double array, dimension (min(m,n)) \n
              Backup of data from working array. \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N)). \n
              For good performance, LWORK should generally be larger. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	RWORK
              RWORK is REAL array, dimension (5*min(M,N)) \n
              On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the
              unconverged superdiagonal elements of an upper bidiagonal
              matrix B whose diagonal is in S (not necessarily sorted).
              B satisfies A = U * B * VT, so it has the same singular
              values as A, and singular vectors related by U and VT. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit. \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n
              > 0:  if CBDSQR did not converge, INFO specifies how many
                    superdiagonals of an intermediate bidiagonal form B
                    did not converge to zero. See the description of RWORK
                    above for details. \n

    *     *  */
    template <typename T>
    void gesvd(char *jobu, char *jobv, integer *m, integer *n, T *a, integer *lda, T *s, T *u,
               integer *ldu, T *vt, integer *ldvt, T *work, integer *lwork, integer *info)
    {
        gesvd(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
    }
    template <typename T, typename Ta>
    void gesvd(char *jobu, char *jobv, integer *m, integer *n, T *a, integer *lda, Ta *s, T *u,
               integer *ldu, T *vt, integer *ldvt, T *work, integer *lwork, Ta *rwork,
               integer *info)
    {
        gesvd(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
    }
    /** @}*/ // end of gesvd

    /** @defgroup gesvdq gesvdq
     * @ingroup SVD
     * @{
     */
    /*! @brief GESVDQ computes the singular value decomposition (SVD) with a QR-Preconditioned QR
 SVD Method for GE matrices
 * @details
 * \b Purpose:
    \verbatim
     GESVDQ computes the singular value decomposition (SVD) of a real
     M-by-N matrix A, where M >= N. The SVD of A is written as
                                        [++]   [xx]   [x0]   [xx]
                  A = U * SIGMA * V^*,  [++] = [xx] * [ox] * [xx]
                                        [++]   [xx]
     where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal
     matrix, and V is an N-by-N orthogonal matrix. The diagonal elements
     of SIGMA are the singular values of A. The columns of U and V are the
     left and the right singular vectors of A, respectively.
    \endverbatim

 * @param[in] JOBP
          JOBP is CHARACTER*1 \n
          = 'P' The rows of A are ordered in decreasing order with respect to
                ||A(i,:)||_\infty. This enhances numerical accuracy at the cost
                of extra data movement. Recommended for numerical robustness. \n
          = 'N' No row pivoting. \n
 * @param[in] JOBR
          JOBR is CHARACTER*1 \n
          = 'T' After the initial pivoted QR factorization, GESVD is applied to
          the transposed R**T of the computed triangular factor R. This involves
          some extra data movement (matrix transpositions). Useful for
          experiments, research and development. \n
          = 'N' The triangular factor R is given as input to GESVD. This may be
          preferred as it involves less data movement. \n
 * @param[in] JOBU
          JOBU is CHARACTER*1 \n
          = 'A' All M left singular vectors are computed and returned in the
          matrix U. See the description of U. \n
          = 'S' or 'U' N = min(M,N) left singular vectors are computed and returned
          in the matrix U. See the description of U. \n
          = 'R' Numerical rank NUMRANK is determined and only NUMRANK left singular
          vectors are computed and returned in the matrix U. \n
          = 'F' The N left singular vectors are returned in factored form as the
          product of the Q factor from the initial QR factorization and the
          N left singular vectors of (R**T , 0)**T. If row pivoting is used,
          then the necessary information on the row pivoting is stored in
          IWORK(N+1:N+M-1). \n
          = 'N' The left singular vectors are not computed. \n
 * @param[in] JOBV
          JOBV is CHARACTER*1 \n
          = 'A', 'V' All N right singular vectors are computed and returned in
          the matrix V. \n
          = 'R' Numerical rank NUMRANK is determined and only NUMRANK right singular
          vectors are computed and returned in the matrix V. This option is
          allowed only if JOBU = 'R' or JOBU = 'N'; otherwise it is illegal. \n
          = 'N' The right singular vectors are not computed. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the input matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the input matrix A.  M >= N >= 0. \n
 * @param[in,out] A
          A is REAL array of dimensions LDA x N \n
          On entry, the input matrix A. \n
          On exit, if JOBU .NE. 'N' or JOBV .NE. 'N', the lower triangle of A contains
          the Householder vectors as stored by GEQP3. If JOBU = 'F', these Householder
          vectors together with WORK(1:N) can be used to restore the Q factors from
          the initial pivoted QR factorization of A. See the description of U. \n
 * @param[in] LDA
          LDA is INTEGER. \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] S
          S is REAL array of dimension N. \n
          The singular values of A, ordered so that S(i) >= S(i+1). \n
 * @param[out] U
          U is REAL array, dimension \n
          LDU x M if JOBU = 'A'; see the description of LDU. In this case,
          on exit, U contains the M left singular vectors. \n
          LDU x N if JOBU = 'S', 'U', 'R' ; see the description of LDU. In this
          case, U contains the leading N or the leading NUMRANK left singular vectors. \n
          LDU x N if JOBU = 'F' ; see the description of LDU. In this case U
          contains N x N orthogonal matrix that can be used to form the left
          singular vectors. \n
          If JOBU = 'N', U is not referenced. \n
 * @param[in] LDU
          LDU is INTEGER. \n
          The leading dimension of the array U. \n
          If JOBU = 'A', 'S', 'U', 'R',  LDU >= fla_max(1,M). \n
          If JOBU = 'F',                 LDU >= fla_max(1,N). \n
          Otherwise,                     LDU >= 1. \n
 * @param[out] V
          V is REAL array, dimension \n
          LDV x N if JOBV = 'A', 'V', 'R' or if JOBA = 'E' . \n
          If JOBV = 'A', or 'V',  V contains the N-by-N orthogonal matrix  V**T; \n
          If JOBV = 'R', V contains the first NUMRANK rows of V**T (the right
          singular vectors, stored rowwise, of the NUMRANK largest singular values). \n
          If JOBV = 'N' and JOBA = 'E', V is used as a workspace. \n
          If JOBV = 'N', and JOBA.NE.'E', V is not referenced. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. \n
          If JOBV = 'A', 'V', 'R',  or JOBA = 'E', LDV >= fla_max(1,N).
          Otherwise,                               LDV >= 1. \n
 * @param[out] NUMRANK
          NUMRANK is INTEGER \n
          NUMRANK is the numerical rank first determined after the rank
          revealing QR factorization, following the strategy specified by the
          value of JOBA. If JOBV = 'R' and JOBU = 'R', only NUMRANK
          leading singular values and vectors are then requested in the call
          of GESVD. The final value of NUMRANK might be further reduced if
          some singular values are computed as zeros. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (fla_max(1, LIWORK)). \n
          On exit, IWORK(1:N) contains column pivoting permutation of the
          rank revealing QR factorization. \n
          If JOBP = 'P', IWORK(N+1:N+M-1) contains the indices of the sequence
          of row swaps used in row pivoting. These can be used to restore the
          left singular vectors in the case JOBU = 'F'. \n
 \n
          If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0,
          IWORK(1) returns the minimal LIWORK. \n
 * @param[in]	LIWORK
          LIWORK is INTEGER \n
          The dimension of the array IWORK. \n
          LIWORK >= N + M - 1,  if JOBP = 'P'; \n
          LIWORK >= N           if JOBP = 'N'. \n
 \n
          If LIWORK = -1, then a workspace query is assumed; the routine
          only calculates and returns the optimal and minimal sizes
          for the CWORK, IWORK, and RWORK arrays, and no error
          message related to LCWORK is issued by XERBLA. \n
 * @param[out]	CWORK
          CWORK is COMPLEX array, dimension (fla_max(2, LCWORK)), used as a workspace.
          On exit, if, on entry, LCWORK.NE.-1, CWORK(1:N) contains parameters
          needed to recover the Q factor from the QR factorization computed by
          CGEQP3. \n
 \n
          If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0,
          CWORK(1) returns the optimal LCWORK, and
          CWORK(2) returns the minimal LCWORK. \n
 * @param[in,out]	LCWORK
          LCWORK is INTEGER \n
          The dimension of the array CWORK. It is determined as follows:
          Let  LWQP3 = N+1,  LWCON = 2*N, and let \n
          LWUNQ = { MAX( N, 1 ),  if JOBU = 'R', 'S', or 'U' \n
                  { MAX( M, 1 ),  if JOBU = 'A' \n
          LWSVD = MAX( 3*N, 1 ) \n
          LWLQF = MAX( N/2, 1 ), LWSVD2 = MAX( 3*(N/2), 1 ), LWUNLQ = MAX( N, 1 ), \n
          LWQRF = MAX( N/2, 1 ), LWUNQ2 = MAX( N, 1 ) \n
          Then the minimal value of LCWORK is: \n
          = MAX( N + LWQP3, LWSVD )        if only the singular values are needed; \n
          = MAX( N + LWQP3, LWCON, LWSVD ) if only the singular values are needed,
                                   and a scaled condition estimate requested; \n
 \n
          = N + MAX( LWQP3, LWSVD, LWUNQ ) if the singular values and the left
                                   singular vectors are requested; \n
          = N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ) if the singular values and the left
                                   singular vectors are requested, and also
                                   a scaled condition estimate requested; \n
 \n
          = N + MAX( LWQP3, LWSVD )        if the singular values and the right
                                   singular vectors are requested; \n
          = N + MAX( LWQP3, LWCON, LWSVD ) if the singular values and the right
                                   singular vectors are requested, and also
                                   a scaled condition etimate requested; \n
 \n
          = N + MAX( LWQP3, LWSVD, LWUNQ ) if the full SVD is requested with JOBV = 'R';
                                   independent of JOBR; \n
          = N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ) if the full SVD is requested,
                                   JOBV = 'R' and, also a scaled condition
                                   estimate requested; independent of JOBR; \n
          = MAX( N + MAX( LWQP3, LWSVD, LWUNQ ),
         N + MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ) ) if the
                         full SVD is requested with JOBV = 'A' or 'V', and
                         JOBR ='N' \n
          = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ),
         N + MAX( LWQP3, LWCON, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ ) )
                         if the full SVD is requested with JOBV = 'A' or 'V', and
                         JOBR ='N', and also a scaled condition number estimate
                         requested. \n
          = MAX( N + MAX( LWQP3, LWSVD, LWUNQ ),
         N + MAX( LWQP3, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ ) ) if the
                         full SVD is requested with JOBV = 'A', 'V', and JOBR ='T'
          = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ), \n
         N + MAX( LWQP3, LWCON, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ ) )
                         if the full SVD is requested with JOBV = 'A', 'V' and
                         JOBR ='T', and also a scaled condition number estimate
                         requested. \n
          Finally, LCWORK must be at least two: LCWORK = MAX( 2, LCWORK ). \n
 \n
          If LCWORK = -1, then a workspace query is assumed; the routine
          only calculates and returns the optimal and minimal sizes
          for the CWORK, IWORK, and RWORK arrays, and no error
          message related to LCWORK is issued by XERBLA. \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (fla_max(1, LRWORK)). \n
          On exit, \n
          1. If JOBA = 'E', RWORK(1) contains an estimate of the condition
          number of column scaled A. If A = C * D where D is diagonal and C
          has unit columns in the Euclidean norm, then, assuming full column rank,
          N^(-1/4) * RWORK(1) <= ||pinv(C)||_2 <= N^(1/4) * RWORK(1).
          Otherwise, RWORK(1) = -1. \n
          2. RWORK(2) contains the number of singular values computed as
          exact zeros in CGESVD applied to the upper triangular or trapezoidal
          R (from the initial QR factorization). In case of early exit (no call to
          CGESVD, such as in the case of zero matrix) RWORK(2) = -1. \n
 \n
          If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0,
          RWORK(1) returns the minimal LRWORK. \n
 * @param[in]	LRWORK
          LRWORK is INTEGER. \n
          The dimension of the array RWORK.
          If JOBP ='P', then LRWORK >= MAX(2, M, 5*N);
          Otherwise, LRWORK >= MAX(2, 5*N). \n
 \n
          If LRWORK = -1, then a workspace query is assumed; the routine
          only calculates and returns the optimal and minimal sizes
          for the CWORK, IWORK, and RWORK arrays, and no error
          message related to LCWORK is issued by XERBLA. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if CBDSQR did not converge, INFO specifies how many superdiagonals
          of an intermediate bidiagonal form B (computed in CGESVD) did not
          converge to zero. \n

 *  * */
    template <typename T>
    void gesvdq(char *joba, char *jobp, char *jobr, char *jobu, char *jobv, integer *m, integer *n,
                T *a, integer *lda, T *s, T *u, integer *ldu, T *v, integer *ldv, integer *numrank,
                integer *iwork, integer *liwork, T *work, integer *lwork, T *rwork, integer *lrwork,
                integer *info)
    {
        gesvdq(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork,
               liwork, work, lwork, rwork, lrwork, info);
    }
    template <typename T, typename Ta>
    void gesvdq(char *joba, char *jobp, char *jobr, char *jobu, char *jobv, integer *m, integer *n,
                T *a, integer *lda, Ta *s, T *u, integer *ldu, T *v, integer *ldv, integer *numrank,
                integer *iwork, integer *liwork, T *cwork, integer *lcwork, Ta *rwork,
                integer *lrwork, integer *info)
    {
        gesvdq(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank, iwork,
               liwork, cwork, lcwork, rwork, lrwork, info);
    }
    /** @}*/ // end of gesvdq

    /** @defgroup gesdd gesdd
     * @ingroup SVD
     * @{
     */
    /*! @brief General matrix singular value decomposition (divide-and-conquer)
  * @details
  * \b Purpose:
    \verbatim
    Computation of singular value decomposition (SVD) of a real m-by-n matrix a, optionally
    computing the left and right singular vectors.  If singular vectors are desired, it uses a
    divide-and-conquer algorithm.

    The SVD is written

    A = U * SIGMA * transpose(V)

    where SIGMA is an m-by-n matrix which is zero except for its min(m,n) diagonal elements,
    U is an M-by-M orthogonal matrix, and V is an n-by-n orthogonal matrix.  The diagonal
    elements of SIGMA are the singular values of A; they are real and non-negative, and are
    returned in descending order.  The first min(m,n) columns of U and V are the left and right
    singular vectors of A.

    Note that the routine returns vt = V**T, not V.

    The divide and conquer algorithm makes very mild assumptions about floating point
    arithmetic. It will work on machines with a guard digit in add/subtract, or on those
    binary machines without guard digits which subtract like the Cray X-MP, Cray Y-MP, Cray
    C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without
    guard digits, but we know of none.
   \endverbatim

    * @param[in] jobz
              jobz is char* \n
              Specifies options for computing all or part of the matrix U: \n
              = 'A':  all M columns of U and all N rows of V**T are
              returned in the arrays U and vt; \n
              = 'S':  the first min(m,n) columns of U and the first
              min(m,n) rows of V**T are returned in the arrays U
              and vt; \n
              = 'O':  If m >= n, the first N columns of U are overwritten
              on the array a and all rows of V**T are returned in
              the array vt; \n
              otherwise, all columns of U are returned in the
              array U and the first M rows of V**T are overwritten
              in the array a; \n
              = 'N':  no columns of U or rows of V**T are computed. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the input matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the input matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the m-by-n matrix a. \n
              On exit, \n
              if jobz = 'O',  A is overwritten with the first N columns
              of U (the left singular vectors, stored
              columnwise) if m >= n;
              A is overwritten with the first M rows
              of V**T (the right singular vectors, stored
              rowwise) otherwise. \n
              if jobz .ne. 'O', the contents of A are destroyed. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,m). \n
    * @param[out] s
              s is float/double array, dimension (min(m,n)) \n
              The singular values of A, sorted so that S(i) >= S(i+1). \n
    * @param[out] u
              U is float/double array, dimension (ldu,UCOL) \n
              UCOL = m if jobz = 'A' or jobz = 'O' and M < N; \n
              UCOL = min(m,n) if jobz = 'S'. \n
              If jobz = 'A' or jobz = 'O' and m < m, U contains the m-by-m
              orthogonal matrix U; \n
              if jobz = 'S', U contains the first min(m,n) columns of U
              (the left singular vectors, stored columnwise); \n
              if jobz = 'O' and m >= n, or jobz = 'N', U is not referenced. \n
    * @param[in] ldu
              ldu is integer* \n
              The leading dimension of the array U.  ldu >= 1; if
              jobz = 'S' or 'A' or jobz = 'O' and m < n, ldu >= m. \n
    * @param[out] vt
              vt is float/double array, dimension (ldvt,n) \n
              If jobz = 'A' or jobz = 'O' and m >= n, vt contains the
              n-by-n orthogonal matrix V**T; \n
              if jobz = 'S', vt contains the first min(m,n) rows of
              V**T (the right singular vectors, stored rowwise); \n
              if jobz = 'O' and M < N, or jobz = 'N', vt is not referenced. \n
    * @param[in] ldvt
              ldvt is integer* \n
              The leading dimension of the array vt.  ldvt >= 1; \n
              if jobz = 'A' or jobz = 'O' and m >= n, ldvt >= n; \n
              if jobz = 'S', ldvt >= min(m,n). \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= 1. \n
              If LWORK = -1, a workspace query is assumed.  The optimal
              size for the WORK array is calculated and stored in WORK(1),
              and no other work except argument checking is performed. \n
 \n
              Let mx = fla_max(M,N) and mn = min(M,N). \n
              If JOBZ = 'N', LWORK >= 2*mn + mx. \n
              If JOBZ = 'O', LWORK >= 2*mn*mn + 2*mn + mx. \n
              If JOBZ = 'S', LWORK >=   mn*mn + 3*mn. \n
              If JOBZ = 'A', LWORK >=   mn*mn + 2*mn + mx. \n
              These are not tight minimums in all cases; see comments inside code.
              For good performance, LWORK should generally be larger;
              a query is recommended. \n
    * @param[out]	RWORK
              RWORK is REAL array, dimension (MAX(1,LRWORK)) \n
              Let mx = fla_max(M,N) and mn = min(M,N). \n
              If JOBZ = 'N',    LRWORK >= 5*mn (LAPACK <= 3.6 needs 7*mn); \n
              else if mx >> mn, LRWORK >= 5*mn*mn + 5*mn; \n
              else              LRWORK >= fla_max( 5*mn*mn + 5*mn,
                                               2*mx*mn + 2*mn*mn + mn ). \n
    * @param[out]	IWORK
              IWORK is INTEGER array, dimension (8*min(M,N)) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit. \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n
              > 0:  The updating process of SBDSDC did not converge. \n

    *     *  */
    template <typename T>
    void gesdd(char *jobz, integer *m, integer *n, T *a, integer *lda, T *s, T *u, integer *ldu,
               T *vt, integer *ldvt, T *work, integer *lwork, integer *iwork, integer *info)
    {
        gesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void gesdd(char *jobz, integer *m, integer *n, T *a, integer *lda, Ta *s, T *u, integer *ldu,
               T *vt, integer *ldvt, T *work, integer *lwork, Ta *rwork, integer *iwork,
               integer *info)
    {
        gesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
    }
    /** @}*/ // end of gesdd

    /** @defgroup gesvdx gesvdx
     * @ingroup SVD
     * @{
     */
    /*! @brief GESVDX computes the singular value decomposition (SVD) for GE matrices
 * @details
 * \b Purpose:
    \verbatim
     GESVDX computes the singular value decomposition (SVD) of a real
     M-by-N matrix A, optionally computing the left and/or right singular
     vectors. The SVD is written

         A = U * SIGMA * transpose(V)

     where SIGMA is an M-by-N matrix which is zero except for its
     min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
     V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
     are the singular values of A; they are real and non-negative, and
     are returned in descending order.  The first min(m,n) columns of
     U and V are the left and right singular vectors of A.

     GESVDX uses an eigenvalue problem for obtaining the SVD, which
     allows for the computation of a subset of singular values and
     vectors. See BDSVDX for details.

     Note that the routine returns V**T, not V.
    \endverbatim

 * @param[in] JOBU
          JOBU is CHARACTER*1 \n
          Specifies options for computing all or part of the matrix U:
          = 'V':  the first min(m,n) columns of U (the left singular
                  vectors) or as specified by RANGE are returned in
                  the array U; \n
          = 'N':  no columns of U (no left singular vectors) are
                  computed. \n
 * @param[in] JOBVT
          JOBVT is CHARACTER*1 \n
          Specifies options for computing all or part of the matrix
          V**T: \n
          = 'V':  the first min(m,n) rows of V**T (the right singular
                  vectors) or as specified by RANGE are returned in
                  the array VT; \n
          = 'N':  no rows of V**T (no right singular vectors) are
                  computed. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all singular values will be found. \n
          = 'V': all singular values in the half-open interval (VL,VU]
                 will be found. \n
          = 'I': the IL-th through IU-th singular values will be found. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the input matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the input matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A.
          On exit, the contents of A are destroyed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[in] VL
          VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for singular values. VU > VL.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for singular values. VU > VL.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest singular value to be returned.
          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER \n
          If RANGE='I', the index of the
          largest singular value to be returned.
          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[out] NS
          NS is INTEGER \n
          The total number of singular values found,
          0 <= NS <= min(M,N). \n
          If RANGE = 'A', NS = min(M,N); if RANGE = 'I', NS = IU-IL+1. \n
 * @param[out] S
          S is REAL array, dimension (min(M,N)) \n
          The singular values of A, sorted so that S(i) >= S(i+1). \n
 * @param[out] U
          U is REAL array, dimension (LDU,UCOL) \n
          If JOBU = 'V', U contains columns of U (the left singular
          vectors, stored columnwise) as specified by RANGE; if
          JOBU = 'N', U is not referenced. \n
          Note: The user must ensure that UCOL >= NS; if RANGE = 'V',
          the exact value of NS is not known in advance and an upper
          bound must be used. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U.  LDU >= 1; if
          JOBU = 'V', LDU >= M. \n
 * @param[out] VT
          VT is REAL array, dimension (LDVT,N) \n
          If JOBVT = 'V', VT contains the rows of V**T (the right singular
          vectors, stored rowwise) as specified by RANGE; if JOBVT = 'N',
          VT is not referenced. \n
          Note: The user must ensure that LDVT >= NS; if RANGE = 'V',
          the exact value of NS is not known in advance and an upper
          bound must be used. \n
 * @param[in] LDVT
          LDVT is INTEGER \n
          The leading dimension of the array VT.  LDVT >= 1; if
          JOBVT = 'V', LDVT >= NS (see above). \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK; \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
          LWORK >= MAX(1,MIN(M,N)*(MIN(M,N)+4)) for the paths (see
          comments inside the code): \n
             - PATH 1  (M much larger than N) \n
             - PATH 1t (N much larger than M) \n
          LWORK >= MAX(1,MIN(M,N)*2+MAX(M,N)) for the other paths.
          For good performance, LWORK should generally be larger. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (MAX(1,LRWORK)) \n
          LRWORK >= MIN(M,N)*(MIN(M,N)*2+15*MIN(M,N)). \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (12*MIN(M,N)) \n
          If INFO = 0, the first NS elements of IWORK are zero. If INFO > 0,
          then IWORK contains the indices of the eigenvectors that failed
          to converge in SBDSVDX/SSTEVX. \n
 * @param[out]	INFO
          INFO is INTEGER \n
           = 0:  successful exit \n
           < 0:  if INFO = -i, the i-th argument had an illegal value \n
           > 0:  if INFO = i, then i eigenvectors failed to converge
                 in SBDSVDX/SSTEVX.
                 if INFO = N*2 + 1, an internal error occurred in
                 SBDSVDX \n

 *  * */
    template <typename T>
    void gesvdx(char *jobu, char *jobvt, char *range, integer *m, integer *n, T *a, integer *lda,
                T *vl, T *vu, integer *il, integer *iu, integer *ns, T *s, T *u, integer *ldu,
                T *vt, integer *ldvt, T *work, integer *lwork, integer *iwork, integer *info)
    {
        gesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work,
               lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void gesvdx(char *jobu, char *jobvt, char *range, integer *m, integer *n, T *a, integer *lda,
                Ta *vl, Ta *vu, integer *il, integer *iu, integer *ns, Ta *s, T *u, integer *ldu,
                T *vt, integer *ldvt, T *work, integer *lwork, Ta *rwork, integer *iwork,
                integer *info)
    {
        gesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work,
               lwork, rwork, iwork, info);
    }
    /** @}*/ // end of gesvdx

    /** @defgroup gejsv gejsv
     * @ingroup SVD
     * @{
     */
    /*! @brief Computes the singular value decomposition (SVD) of a real matrix
* @details
* \b Purpose:
* \verbatim
    GEJSV computes the singular value decomposition (SVD) of a real M-by-N
    matrix [A], where M >= N. The SVD of [A] is written as

                [A] = [U] * [SIGMA] * [V]^t,

    where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N
    diagonal elements, [U] is an M-by-N (or M-by-M) orthonormal matrix, and
    [V] is an N-by-N orthogonal matrix. The diagonal elements of [SIGMA] are
    the singular values of [A]. The columns of [U] and [V] are the left and
    the right singular vectors of [A], respectively. The matrices [U] and [V]
    are computed and stored in the arrays U and V, respectively. The diagonal
    of [SIGMA] is computed and stored in the array SVA.
    GEJSV can sometimes compute tiny singular values and their singular vectors much
    more accurately than other SVD routines, see below under Further Details.
 \endverbatim

 * @param[in] JOBA
         JOBA is CHARACTER*1 \n
         Specifies the level of accuracy: \n
       = 'C': This option works well (high relative accuracy) if A = B * D,
              with well-conditioned B and arbitrary diagonal matrix D.
              The accuracy cannot be spoiled by COLUMN scaling. The
              accuracy of the computed output depends on the condition of
              B, and the procedure aims at the best theoretical accuracy.
              The relative error max_{i=1:N}|d sigma_i| / sigma_i is
              bounded by f(M,N)*epsilon* cond(B), independent of D.
              The input matrix is preprocessed with the QRF with column
              pivoting. This initial preprocessing and preconditioning by
              a rank revealing QR factorization is common for all values of
              JOBA. Additional actions are specified as follows: \n
       = 'E': Computation as with 'C' with an additional estimate of the
              condition number of B. It provides a realistic error bound. \n
       = 'F': If A = D1 * C * D2 with ill-conditioned diagonal scalings
              D1, D2, and well-conditioned matrix C, this option gives
              higher accuracy than the 'C' option. If the structure of the
              input matrix is not known, and relative accuracy is
              desirable, then this option is advisable. The input matrix A
              is preprocessed with QR factorization with FULL (row and
              column) pivoting. \n
       = 'G': Computation as with 'F' with an additional estimate of the
              condition number of B, where A=D*B. If A has heavily weighted
              rows, then using this condition number gives too pessimistic
              error bound. \n
       = 'A': Small singular values are the noise and the matrix is treated
              as numerically rank deficient. The error in the computed
              singular values is bounded by f(m,n)*epsilon*||A||.
              The computed SVD A = U * S * V^t restores A up to
              f(m,n)*epsilon*||A||. \n
              This gives the procedure the licence to discard (set to zero)
              all singular values below N*epsilon*||A||. \n
       = 'R': Similar as in 'A'. Rank revealing property of the initial
              QR factorization is used do reveal (using triangular factor)
              a gap sigma_{r+1} < epsilon * sigma_r in which case the
              numerical RANK is declared to be r. The SVD is computed with
              absolute error bounds, but more accurately than with 'A'. \n
 * @param[in] JOBU
         JOBU is CHARACTER*1 \n
         Specifies whether to compute the columns of U: \n
       = 'U': N columns of U are returned in the array U. \n
       = 'F': full set of M left sing. vectors is returned in the array U. \n
       = 'W': U may be used as workspace of length M*N. See the description
              of U. \n
       = 'N': U is not computed. \n
 * @param[in] JOBV
         JOBV is CHARACTER*1 \n
         Specifies whether to compute the matrix V: \n
       = 'V': N columns of V are returned in the array V; Jacobi rotations
              are not explicitly accumulated. \n
       = 'J': N columns of V are returned in the array V, but they are
              computed as the product of Jacobi rotations. This option is
              allowed only if JOBU .NE. 'N', i.e. in computing the full SVD. \n
       = 'W': V may be used as workspace of length N*N. See the description
              of V. \n
       = 'N': V is not computed. \n
 * @param[in] JOBR
         JOBR is CHARACTER*1 \n
         Specifies the RANGE for the singular values. Issues the licence to
         set to zero small positive singular values if they are outside
         specified range. If A .NE. 0 is scaled so that the largest singular
         value of c*A is around SQRT(BIG), BIG=SLAMCH('O'), then JOBR issues
         the licence to kill columns of A whose norm in c*A is less than
         SQRT(SFMIN) (for JOBR = 'R'), or less than SMALL=SFMIN/EPSLN,
         where SFMIN=SLAMCH('S'), EPSLN=SLAMCH('E'). \n
       = 'N': Do not kill small columns of c*A. This option assumes that
              BLAS and QR factorizations and triangular solvers are
              implemented to work in that range. If the condition of A
              is greater than BIG, use SGESVJ. \n
       = 'R': RESTRICTED range for sigma(c*A) is [SQRT(SFMIN), SQRT(BIG)]
              (roughly, as described above). This option is recommended.
                                             ===========================
         For computing the singular values in the FULL range [SFMIN,BIG]
         use SGESVJ. \n
 * @param[in] JOBT
         JOBT is CHARACTER*1 \n
         If the matrix is square then the procedure may determine to use
         transposed A if A^t seems to be better with respect to convergence.
         If the matrix is not square, JOBT is ignored. This is subject to
         changes in the future. \n
         The decision is based on two values of entropy over the adjoint
         orbit of A^t * A. See the descriptions of WORK(6) and WORK(7). \n
       = 'T': transpose if entropy test indicates possibly faster
         convergence of Jacobi process if A^t is taken as input. If A is
         replaced with A^t, then the row pivoting is included automatically. \n
       = 'N': do not speculate.
         This option can be used to compute only the singular values, or the
         full SVD (U, SIGMA and V). For only one set of singular vectors
         (U or V), the caller should provide both U and V, as one of the
         matrices is used as workspace if the matrix A is transposed.
         The implementer can easily remove this constraint and make the
         code more complicated. See the descriptions of U and V. \n
 * @param[in] JOBP
         JOBP is CHARACTER*1 \n
         Issues the licence to introduce structured perturbations to drown
         denormalized numbers. This licence should be active if the
         denormals are poorly implemented, causing slow computation,
         especially in cases of fast convergence (!). For details see [1,2].
         For the sake of simplicity, this perturbations are included only
         when the full SVD or only the singular values are requested. The
         implementer/user can easily add the perturbation for the cases of
         computing one set of singular vectors. \n
       = 'P': introduce perturbation \n
       = 'N': do not perturb \n
 * @param[in] M
         M is INTEGER \n
         The number of rows of the input matrix A.  M >= 0. \n
 * @param[in] N
         N is INTEGER \n
         The number of columns of the input matrix A. M >= N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] SVA
          SVA is REAL array, dimension (N) \n
          On exit, \n
          - For WORK(1)/WORK(2) = ONE: The singular values of A. During the
            computation SVA contains Euclidean column norms of the
            iterated matrices in the array A. \n
          - For WORK(1) .NE. WORK(2): The singular values of A are
            (WORK(1)/WORK(2)) * SVA(1:N). This factored form is used if
            sigma_max(A) overflows or if small singular values have been
            saved from underflow by scaling the input matrix A. \n
          - If JOBR='R' then some of the singular values may be returned
            as exact zeros obtained by "set to zero" because they are
            below the numerical rank threshold or are denormalized numbers.
 * @param[out] U
          U is REAL array, dimension ( LDU, N) \n
          If JOBU = 'U', then U contains on exit the M-by-N matrix of
                         the left singular vectors. \n
          If JOBU = 'F', then U contains on exit the M-by-M matrix of
                         the left singular vectors, including an ONB
                         of the orthogonal complement of the Range(A). \n
          If JOBU = 'W'  .AND. (JOBV = 'V' .AND. JOBT = 'T' .AND. M = N),
                         then U is used as workspace if the procedure
                         replaces A with A^t. In that case, [V] is computed
                         in U as left singular vectors of A^t and then
                         copied back to the V array. This 'W' option is just
                         a reminder to the caller that in this case U is
                         reserved as workspace of length N*N. \n
          If JOBU = 'N'  U is not referenced, unless JOBT= - computes row and
                         column scaling to reduce condition number of matrix. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U,  LDU >= 1. \n
          IF  JOBU = 'U' or 'F' or 'W',  then LDU >= M. \n
 * @param[out] V
          V is REAL array, dimension ( LDV, N) \n
          If JOBV = 'V', 'J' then V contains on exit the N-by-N matrix of
                         the right singular vectors; \n
          If JOBV = 'W', AND (JOBU = 'U' AND JOBT = 'T' AND M = N),
                         then V is used as workspace if the pprocedure
                         replaces A with A^t. In that case, [U] is computed
                         in V as right singular vectors of A^t and then
                         copied back to the U array. This 'W' option is just
                         a reminder to the caller that in this case V is
                         reserved as workspace of length N*N. \n
          If JOBV = 'N'  V is not referenced, unless JOBT='T'. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V,  LDV >= 1. \n
          If JOBV = 'V' or 'J' or 'W', then LDV >= N. \n
 * @param[out]	CWORK
          CWORK is COMPLEX array, dimension (MAX(2,LWORK)) \n
          If the call to CGEJSV is a workspace query (indicated by LWORK=-1 or
          LRWORK=-1), then on exit CWORK(1) contains the required length of
          CWORK for the job parameters used in the call. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          Length of CWORK to confirm proper allocation of workspace. \n
          LWORK depends on the job: \n
 \n
          1. If only SIGMA is needed ( JOBU = 'N', JOBV = 'N' ) and \n
            1.1 .. no scaled condition estimate required (JOBA.NE.'E'.AND.JOBA.NE.'G'): \n
               LWORK >= 2*N+1. This is the minimal requirement.
               ->> For optimal performance (blocked code) the optimal value
               is LWORK >= N + (N+1)*NB. Here NB is the optimal
               block size for CGEQP3 and CGEQRF.
               In general, optimal LWORK is computed as
               LWORK >= fla_max(N+LWORK(CGEQP3),N+LWORK(CGEQRF), LWORK(CGESVJ)). \n
            1.2. .. an estimate of the scaled condition number of A is
               required (JOBA='E', or 'G'). In this case, LWORK the minimal
               requirement is LWORK >= N*N + 2*N.
               ->> For optimal performance (blocked code) the optimal value
               is LWORK >= fla_max(N+(N+1)*NB, N*N+2*N)=N**2+2*N.
               In general, the optimal length LWORK is computed as
               LWORK >= fla_max(N+LWORK(CGEQP3),N+LWORK(CGEQRF), LWORK(CGESVJ),
                            N*N+LWORK(CPOCON)). \n
          2. If SIGMA and the right singular vectors are needed (JOBV = 'V'),
             (JOBU = 'N') \n
            2.1   .. no scaled condition estimate requested (JOBE = 'N'): \n
            -> the minimal requirement is LWORK >= 3*N. \n
            -> For optimal performance,
               LWORK >= fla_max(N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
               where NB is the optimal block size for CGEQP3, CGEQRF, CGELQ,
               CUNMLQ. In general, the optimal length LWORK is computed as
               LWORK >= fla_max(N+LWORK(CGEQP3), N+LWORK(CGESVJ),
                       N+LWORK(CGELQF), 2*N+LWORK(CGEQRF), N+LWORK(CUNMLQ)). \n
            2.2 .. an estimate of the scaled condition number of A is
               required (JOBA='E', or 'G'). \n
            -> the minimal requirement is LWORK >= 3*N. \n
            -> For optimal performance,
               LWORK >= fla_max(N+(N+1)*NB, 2*N,2*N+N*NB)=2*N+N*NB,
               where NB is the optimal block size for CGEQP3, CGEQRF, CGELQ,
               CUNMLQ. In general, the optimal length LWORK is computed as
               LWORK >= fla_max(N+LWORK(CGEQP3), LWORK(CPOCON), N+LWORK(CGESVJ),
                       N+LWORK(CGELQF), 2*N+LWORK(CGEQRF), N+LWORK(CUNMLQ)). \n
          3. If SIGMA and the left singular vectors are needed \n
            3.1  .. no scaled condition estimate requested (JOBE = 'N'): \n
            -> the minimal requirement is LWORK >= 3*N. \n
            -> For optimal performance:
               if JOBU = 'U' :: LWORK >= fla_max(3*N, N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
               where NB is the optimal block size for CGEQP3, CGEQRF, CUNMQR.
               In general, the optimal length LWORK is computed as
               LWORK >= fla_max(N+LWORK(CGEQP3), 2*N+LWORK(CGEQRF), N+LWORK(CUNMQR)). \n
            3.2  .. an estimate of the scaled condition number of A is
               required (JOBA='E', or 'G'). \n
            -> the minimal requirement is LWORK >= 3*N. \n
            -> For optimal performance:
               if JOBU = 'U' :: LWORK >= fla_max(3*N, N+(N+1)*NB, 2*N+N*NB)=2*N+N*NB,
               where NB is the optimal block size for CGEQP3, CGEQRF, CUNMQR.
               In general, the optimal length LWORK is computed as \n
               LWORK >= fla_max(N+LWORK(CGEQP3),N+LWORK(CPOCON), \n
                        2*N+LWORK(CGEQRF), N+LWORK(CUNMQR)). \n
 \n
          4. If the full SVD is needed: (JOBU = 'U' or JOBU = 'F') and \n
            4.1. if JOBV = 'V' \n
               the minimal requirement is LWORK >= 5*N+2*N*N. \n
            4.2. if JOBV = 'J' the minimal requirement is
               LWORK >= 4*N+N*N. \n
            In both cases, the allocated CWORK can accommodate blocked runs
            of CGEQP3, CGEQRF, CGELQF, CUNMQR, CUNMLQ. \n
 \n
          If the call to CGEJSV is a workspace query (indicated by LWORK=-1 or
          LRWORK=-1), then on exit CWORK(1) contains the optimal and CWORK(2) contains the
          minimal length of CWORK for the job parameters used in the call. \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (MAX(7,LWORK)) \n
          On exit, \n
          RWORK(1) = Determines the scaling factor SCALE = RWORK(2) / RWORK(1)
                    such that SCALE*SVA(1:N) are the computed singular values
                    of A. (See the description of SVA().) \n
          RWORK(2) = See the description of RWORK(1). \n
          RWORK(3) = SCONDA is an estimate for the condition number of
                    column equilibrated A. (If JOBA = 'E' or 'G')
                    SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1). \n
                    It is computed using SPOCON. It holds
                    N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
                    where R is the triangular factor from the QRF of A.
                    However, if R is truncated and the numerical rank is
                    determined to be strictly smaller than N, SCONDA is
                    returned as -1, thus indicating that the smallest
                    singular values might be lost. \n
 \n
          If full SVD is needed, the following two condition numbers are
          useful for the analysis of the algorithm. They are provided for
          a developer/implementer who is familiar with the details of
          the method. \n
 \n
          RWORK(4) = an estimate of the scaled condition number of the
                    triangular factor in the first QR factorization. \n
          RWORK(5) = an estimate of the scaled condition number of the
                    triangular factor in the second QR factorization. \n
          The following two parameters are computed if JOBT = 'T'.
          They are provided for a developer/implementer who is familiar
          with the details of the method. \n
          RWORK(6) = the entropy of A^* * A :: this is the Shannon entropy
                    of diag(A^* * A) / Trace(A^* * A) taken as point in the
                    probability simplex. \n
          RWORK(7) = the entropy of A * A^*. (See the description of RWORK(6).)
          If the call to CGEJSV is a workspace query (indicated by LWORK=-1 or
          LRWORK=-1), then on exit RWORK(1) contains the required length of
          RWORK for the job parameters used in the call. \n
 * @param[in]	LRWORK
          LRWORK is INTEGER \n
          Length of RWORK to confirm proper allocation of workspace. \n
          LRWORK depends on the job: \n
 \n
       1. If only the singular values are requested i.e. if
          LSAME(JOBU,'N') .AND. LSAME(JOBV,'N')
          then:
          1.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
               then: LRWORK = fla_max( 7, 2 * M ). \n
          1.2. Otherwise, LRWORK  = fla_max( 7,  N ). \n
       2. If singular values with the right singular vectors are requested \n
          i.e. if
          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J')) .AND.
          .NOT.(LSAME(JOBU,'U').OR.LSAME(JOBU,'F'))
          then: \n
          2.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
          then LRWORK = fla_max( 7, 2 * M ). \n
          2.2. Otherwise, LRWORK  = fla_max( 7,  N ). \n
       3. If singular values with the left singular vectors are requested, i.e. if
          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND.
          .NOT.(LSAME(JOBV,'V').OR.LSAME(JOBV,'J'))
          then: \n
          3.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
          then LRWORK = fla_max( 7, 2 * M ). \n
          3.2. Otherwise, LRWORK  = fla_max( 7,  N ). \n
       4. If singular values with both the left and the right singular vectors
          are requested, i.e. if
          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND.
          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J'))
          then: \n
          4.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'),
          then LRWORK = fla_max( 7, 2 * M ). \n
          4.2. Otherwise, LRWORK  = fla_max( 7, N ). \n
  \n
          If, on entry, LRWORK = -1 or LWORK=-1, a workspace query is assumed and
          the length of RWORK is returned in RWORK(1).  \n
 * @param[out]	IWORK
          IWORK is INTEGER array, of dimension at least 4, that further depends
          on the job: \n
  \n
          1. If only the singular values are requested then: \n
             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
             then the length of IWORK is N+M; otherwise the length of IWORK is N. \n
          2. If the singular values and the right singular vectors are requested then: \n
             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
             then the length of IWORK is N+M; otherwise the length of IWORK is N.  \n
          3. If the singular values and the left singular vectors are requested then: \n
             If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
             then the length of IWORK is N+M; otherwise the length of IWORK is N.  \n
          4. If the singular values with both the left and the right singular vectors \n
             are requested, then: \n
             4.1. If LSAME(JOBV,'J') the length of IWORK is determined as follows: \n
                  If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
                  then the length of IWORK is N+M; otherwise the length of IWORK is N.  \n
             4.2. If LSAME(JOBV,'V') the length of IWORK is determined as follows: \n
                  If ( LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G') )
                  then the length of IWORK is 2*N+M; otherwise the length of IWORK is 2*N. \n
         \n
          On exit, \n
          IWORK(1) = the numerical rank determined after the initial
                     QR factorization with pivoting. See the descriptions
                     of JOBA and JOBR. \n
          IWORK(2) = the number of the computed nonzero singular values \n
          IWORK(3) = if nonzero, a warning message: \n
                     If IWORK(3) = 1 then some of the column norms of A
                     were denormalized floats. The requested high accuracy
                     is not warranted by the data. \n
          IWORK(4) = 1 or -1. If IWORK(4) = 1, then the procedure used A^* to
                     do the job as specified by the JOB parameters. \n
          If the call to CGEJSV is a workspace query (indicated by LWORK = -1 and
          LRWORK = -1), then on exit IWORK(1) contains the required length of
          IWORK for the job parameters used in the call. \n
 * @param[out]	INFO
           INFO is INTEGER \n
           < 0:  if INFO = -i, then the i-th argument had an illegal value. \n
           = 0:  successful exit; \n
           > 0:  CGEJSV  did not converge in the maximal allowed number
                 of sweeps. The computed values may be inaccurate. \n

 *  * */
    template <typename T>
    void gejsv(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, integer *m,
               integer *n, T *a, integer *lda, T *sva, T *u, integer *ldu, T *v, integer *ldv,
               T *stat, integer *istat, T *work, integer *lwork, integer *iwork, integer *info)
    {
        gejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork,
              iwork, info);
    }
    template <typename T, typename Ta>
    void gejsv(char *joba, char *jobu, char *jobv, char *jobr, char *jobt, char *jobp, integer *m,
               integer *n, T *a, integer *lda, Ta *sva, T *u, integer *ldu, T *v, integer *ldv,
               T *cwork, integer *lwork, Ta *rwork, integer *lrwork, integer *iwork, integer *info)
    {
        gejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, cwork, lwork,
              rwork, lrwork, iwork, info);
    }
    /** @}*/ // end of gejsv

    /** @defgroup gesvj gesvj
     * @ingroup SVD
     * @{
     */
    /*! @brief GESVJ computes the singular value decomposition (SVD) of a M-by-N matrix A
 * @details
 * \b Purpose:
    \verbatim
     GESVJ computes the singular value decomposition (SVD) of a real
     M-by-N matrix A, where M >= N. The SVD of A is written as
                                        [++]   [xx]   [x0]   [xx]
                  A = U * SIGMA * V^t,  [++] = [xx] * [ox] * [xx]
                                        [++]   [xx]
     where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal
     matrix, and V is an N-by-N orthogonal matrix. The diagonal elements
     of SIGMA are the singular values of A. The columns of U and V are the
     left and the right singular vectors of A, respectively.
     GESVJ can sometimes compute tiny singular values and their singular vectors much
     more accurately than other SVD routines, see below under Further Details.
    \endverbatim

 * @param[in] JOBA
          JOBA is CHARACTER*1 \n
          Specifies the structure of A. \n
          = 'L': The input matrix A is lower triangular; \n
          = 'U': The input matrix A is upper triangular; \n
          = 'G': The input matrix A is general M-by-N matrix, M >= N. \n
 * @param[in] JOBU
          JOBU is CHARACTER*1 \n
          Specifies whether to compute the left singular vectors
          (columns of U): \n
          = 'U': The left singular vectors corresponding to the nonzero
                 singular values are computed and returned in the leading
                 columns of A. See more details in the description of A.
                 The default numerical orthogonality threshold is set to
                 approximately TOL=CTOL*EPS, CTOL=SQRT(M), EPS=SLAMCH('E'). \n
          = 'C': Analogous to JOBU='U', except that user can control the
                 level of numerical orthogonality of the computed left
                 singular vectors. TOL can be set to TOL = CTOL*EPS, where
                 CTOL is given on input in the array WORK.
                 No CTOL smaller than ONE is allowed. CTOL greater
                 than 1 / EPS is meaningless. The option 'C'
                 can be used if M*EPS is satisfactory orthogonality
                 of the computed left singular vectors, so CTOL=M could
                 save few sweeps of Jacobi rotations.
                 See the descriptions of A and WORK(1). \n
          = 'N': The matrix U is not computed. However, see the
                 description of A. \n
 * @param[in] JOBV
          JOBV is CHARACTER*1 \n
          Specifies whether to compute the right singular vectors, that
          is, the matrix V: \n
          = 'V':  the matrix V is computed and returned in the array V \n
          = 'A':  the Jacobi rotations are applied to the MV-by-N
                  array V. In other words, the right singular vector
                  matrix V is not computed explicitly; instead it is
                  applied to an MV-by-N matrix initially stored in the
                  first MV rows of V. \n
          = 'N':  the matrix V is not computed and the array V is not
                  referenced \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the input matrix A. 1/SLAMCH('E') > M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the input matrix A.
          M >= N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit,
          If JOBU = 'U' .OR. JOBU = 'C': \n
                 If INFO = 0: \n
                 RANKA orthonormal columns of U are returned in the
                 leading RANKA columns of the array A. Here RANKA <= N
                 is the number of computed singular values of A that are
                 above the underflow threshold SLAMCH('S'). The singular
                 vectors corresponding to underflowed or zero singular
                 values are not computed. The value of RANKA is returned
                 in the array WORK as RANKA=NINT(WORK(2)). Also see the
                 descriptions of SVA and WORK. The computed columns of U
                 are mutually numerically orthogonal up to approximately
                 TOL=SQRT(M)*EPS (default); or TOL=CTOL*EPS (JOBU = 'C'),
                 see the description of JOBU. \n
                 If INFO > 0, \n
                 the procedure GESVJ did not converge in the given number
                 of iterations (sweeps). In that case, the computed
                 columns of U may not be orthogonal up to TOL. The output
                 U (stored in A), SIGMA (given by the computed singular
                 values in SVA(1:N)) and V is still a decomposition of the
                 input matrix A in the sense that the residual
                 ||A-SCALE*U*SIGMA*V^T||_2 / ||A||_2 is small. \n
          If JOBU = 'N': \n
                 If INFO = 0: \n
                 Note that the left singular vectors are 'for free' in the
                 one-sided Jacobi SVD algorithm. However, if only the
                 singular values are needed, the level of numerical
                 orthogonality of U is not an issue and iterations are
                 stopped when the columns of the iterated matrix are
                 numerically orthogonal up to approximately M*EPS. Thus,
                 on exit, A contains the columns of U scaled with the
                 corresponding singular values. \n
                 If INFO > 0: \n
                 the procedure GESVJ did not converge in the given number
                 of iterations (sweeps). \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] SVA
          SVA is REAL array, dimension (N) \n
          On exit, \n
          If INFO = 0 : \n
          depending on the value SCALE = WORK(1), we have: \n
                 If SCALE = ONE: \n
                 SVA(1:N) contains the computed singular values of A.
                 During the computation SVA contains the Euclidean column
                 norms of the iterated matrices in the array A.
                 If SCALE .NE. ONE:
                 The singular values of A are SCALE*SVA(1:N), and this
                 factored representation is due to the fact that some of the
                 singular values of A might underflow or overflow. \n
 \n
          If INFO > 0 : \n
          the procedure GESVJ did not converge in the given number of
          iterations (sweeps) and SCALE*SVA(1:N) may not be accurate. \n
 * @param[in] MV
          MV is INTEGER \n
          If JOBV = 'A', then the product of Jacobi rotations in GESVJ
          is applied to the first MV rows of V. See the description of JOBV. \n
 * @param[in,out] V
          V is REAL array, dimension (LDV,N) \n
          If JOBV = 'V', then V contains on exit the N-by-N matrix of
                         the right singular vectors; \n
          If JOBV = 'A', then V contains the product of the computed right
                         singular vector matrix and the initial matrix in
                         the array V. \n
          If JOBV = 'N', then V is not referenced. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V, LDV >= 1. \n
          If JOBV = 'V', then LDV >= fla_max(1,N). \n
          If JOBV = 'A', then LDV >= fla_max(1,MV) . \n
 * @param[in,out]	CWORK
          CWORK is COMPLEX array, dimension (fla_max(1,LWORK)) \n
          Used as workspace. \n
          If on entry LWORK = -1, then a workspace query is assumed and
          no computation is done; CWORK(1) is set to the minial (and optimal)
          length of CWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER. \n
          Length of CWORK, LWORK >= M+N. \n
 * @param[in,out]	RWORK
          RWORK is REAL array, dimension (fla_max(6,LRWORK)) \n
          On entry, \n
          If JOBU = 'C' : \n
          RWORK(1) = CTOL, where CTOL defines the threshold for convergence.
                    The process stops if all columns of A are mutually
                    orthogonal up to CTOL*EPS, EPS=SLAMCH('E').
                    It is required that CTOL >= ONE, i.e. it is not
                    allowed to force the routine to obtain orthogonality
                    below EPSILON. \n
          On exit, \n
          RWORK(1) = SCALE is the scaling factor such that SCALE*SVA(1:N)
                    are the computed singular values of A.
                    (See description of SVA().) \n
          RWORK(2) = NINT(RWORK(2)) is the number of the computed nonzero
                    singular values. \n
          RWORK(3) = NINT(RWORK(3)) is the number of the computed singular
                    values that are larger than the underflow threshold. \n
          RWORK(4) = NINT(RWORK(4)) is the number of sweeps of Jacobi
                    rotations needed for numerical convergence. \n
          RWORK(5) = max_{i.NE.j} |COS(A(:,i),A(:,j))| in the last sweep.
                    This is useful information in cases when CGESVJ did
                    not converge, as it can be used to estimate whether
                    the output is still useful and for post festum analysis. \n
          RWORK(6) = the largest absolute value over all sines of the
                    Jacobi rotation angles in the last sweep. It can be
                    useful for a post festum analysis. \n
         If on entry LRWORK = -1, then a workspace query is assumed and
         no computation is done; RWORK(1) is set to the minial (and optimal)
         length of RWORK. \n
 * @param[in]	LRWORK
         LRWORK is INTEGER \n
         Length of RWORK, LRWORK >= MAX(6,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, then the i-th argument had an illegal value \n
          > 0:  CGESVJ did not converge in the maximal allowed number
                (NSWEEP=30) of sweeps. The output may still be useful.
                See the description of RWORK. \n

 *  * */
    template <typename T>
    void gesvj(char *joba, char *jobu, char *jobv, integer *m, integer *n, T *a, integer *lda,
               T *sva, integer *mv, T *v, integer *ldv, T *work, integer *lwork, integer *info)
    {
        gesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info);
    }
    template <typename T, typename Ta>
    void gesvj(char *joba, char *jobu, char *jobv, integer *m, integer *n, T *a, integer *lda,
               Ta *sva, integer *mv, T *v, integer *ldv, T *cwork, integer *lwork, Ta *rwork,
               integer *lrwork, integer *info)
    {
        gesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork, lwork, rwork, lrwork, info);
    }
    /** @}*/ // end of gesvj

    /** @defgroup bdsqr bdsqr
     * @ingroup SVD
     * @{
     */
    /*! @brief Bidiagonal QR algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Bidiagonal QR algorithm.
        Computation of singular values and, optionally, the right and/or left singular vectors
        from the singular value decomposition (SVD) of a real n-by-n (upper or lower) bidiagonal
        matrix b using the implicit zero-shift QR algorithm.  The SVD of B has the form

        B = Q * S * P**T

        where S is the diagonal matrix of singular values, Q is an orthogonal matrix of left
        singular vectors, and P is an orthogonal matrix of right singular vectors.  If left
        singular vectors are requested, this subroutine actually returns U*Q instead of Q, and,
        if right singular vectors are requested, this subroutine returns P**T*vt instead of P**T,
        for given real input matrices U and vt.  When U and vt are the orthogonal matrices that
        reduce a general matrix a to bidiagonal form:

        A = U*B*vt, as computed by SGEBRD, then

        A = (U*Q) * S * (P**T*vt)

        is the SVD of A.  Optionally, the subroutine may also compute Q**T*C for a given real input
    matrix c.

        See,
        "Computing Small Singular Values of Bidiagonal Matrices With Guaranteed High Relative
        Accuracy," by J. Demmel and W. Kahan, LAPACK Working Note #3 (or SIAM J. Sci. Statist.
        Comput. vol. 11, no. 5, pp. 873-912, Sept 1990) and
        "Accurate singular values and differential qd algorithms," by B. Parlett and V. Fernando,
        Technical Report CPAM-554, Mathematics Department, University of California at Berkeley,
        July 1992 for a detailed description of the algorithm.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  B is upper bidiagonal; \n
              = 'L':  B is lower bidiagonal. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix B.  n >= 0. \n
    * @param[in] ncvt
              ncvt is integer* \n
              The number of columns of the matrix vt. ncvt >= 0. \n
    * @param[in] nru
              nru is integer* \n
              The number of rows of the matrix u. nru >= 0. \n
    * @param[in] ncc
              ncc is integer* \n
              The number of columns of the matrix c. ncc >= 0. \n
    * @param[in,out] d
              d is float/double array, dimension (n) \n
              On entry, the n diagonal elements of the bidiagonal matrix b. \n
              On exit, if info=0, the singular values of B in decreasing
              order. \n
    * @param[in,out] e
              e is float/double array, dimension (n-1) \n
              On entry, the N-1 offdiagonal elements of the bidiagonal
              matrix b. \n
              On exit, if info = 0, E is destroyed; if info > 0, D and E
              will contain the diagonal and superdiagonal elements of a
              bidiagonal matrix orthogonally equivalent to the one given
              as input. \n
    * @param[in,out] vt
              vt is float/double array, dimension (ldvt, ncvt) \n
              On entry, an n-by-ncvt matrix vt. \n
              On exit, vt is overwritten by P**T * vt. \n
              Not referenced if ncvt = 0. \n
    * @param[in] ldvt
              ldvt is integer* \n
              The leading dimension of the array vt. \n
              ldvt >= fla_max(1,n) if ncvt > 0; ldvt >= 1 if ncvt = 0. \n
    * @param[in,out] u
              u is float/double array, dimension (ldu, n) \n
              On entry, an nru-by-n matrix U. \n
              On exit, U is overwritten by U * Q. \n
              Not referenced if nru = 0. \n
    * @param[in] ldu
              ldu is integer* \n
              The leading dimension of the array u.  ldu >= fla_max(1,nru). \n
    * @param[in,out] c
              c is float/doublearray, dimension (ldc, ncc) \n
              On entry, an n-by-ncc matrix c. \n
              On exit, c is overwritten by Q**T * C. \n
              Not referenced if ncc = 0. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. \n
              ldc >= fla_max(1,n) if ncc > 0; ldc >=1 if ncc = 0. \n
              * \par
              * \verbatim
                  Internal Parameters:
                  ====================
                  TOLMUL  T*, default = fla_max(10,min(100,EPS**(-1/8)))
                  TOLMUL controls the convergence criterion of the QR loop.
                  If it is positive, TOLMUL*EPS is the desired relative
                  precision in the computed singular values.
                  If it is negative, abs(TOLMUL*EPS*sigma_max) is the
                  desired absolute accuracy in the computed singular
                  values (corresponds to relative accuracy
                  abs(TOLMUL*EPS) in the largest singular value.
                  abs(TOLMUL) should be between 1 and 1/EPS, and preferably
                  between 10 (for fast convergence) and .1/EPS
                  (for there to be some accuracy in the results).
                  Default is to lose at either one eighth or 2 of the
                  available decimal digits in each computed singular value
                  (whichever is smaller).

                  MAXITR  integer*, default = 6
                  MAXITR controls the maximum number of passes of the
                  algorithm through its inner loop. The algorithms stops
                  (and so fails to converge) if the number of passes
                  through the inner loop exceeds MAXITR*N**2.
              \endverbatim
    * @param[out]	WORK
              WORK is REAL array, dimension (4*N) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  If INFO = -i, the i-th argument had an illegal value \n
              > 0:
                 if NCVT = NRU = NCC = 0, \n
                    = 1, a split was marked by a positive value in E \n
                    = 2, current block of Z not diagonalized after 30*N \n
                         iterations (in inner while loop) \n
                    = 3, termination criterion of outer while loop not met
                         (program created more than N unreduced blocks) \n
                 else NCVT = NRU = NCC = 0, \n
                       the algorithm did not converge; D and E contain the
                       elements of a bidiagonal matrix which is orthogonally
                       similar to the input matrix B;  if INFO = i, i
                       elements of E have not converged to zero. \n

    *     *  */
    template <typename T>
    void bdsqr(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, T *d, T *e, T *vt,
               integer *ldvt, T *u, integer *ldu, T *c, integer *ldc, T *rwork, integer *info)
    {
        bdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
    }
    template <typename T, typename Ta>
    void bdsqr(char *uplo, integer *n, integer *ncvt, integer *nru, integer *ncc, Ta *d, Ta *e,
               T *vt, integer *ldvt, T *u, integer *ldu, T *c, integer *ldc, Ta *rwork,
               integer *info)
    {
        bdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info);
    }
    /** @}*/ // end of bdsqr

    /** @defgroup bdsdc bdsdc
     * @ingroup SVD
     * @{
     */
    /*! @brief BDSDC Bidiagonal divide-and-conquer algorithm
    *
    * @details
    * \b Purpose:
    * \verbatim
        Computation of singular value decomposition (SVD) of a real n-by-n (upper or lower)
        bidiagonal matrix
        B:  B = U * S * vt,
        using a divide and conquer method, where S is a diagonal matrix with non-negative diagonal
        elements (the singular values of B), and U and vt are orthogonal matrices of left and
        right singular vectors, respectively. SBDSDC can be used to compute all singular values,
        and optionally, singular vectors or singular vectors in compact form.

        This code makes very mild assumptions about floating point arithmetic. It will work on
        machines with a guard digit in add/subtract, or on those binary machines without guard
        digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could
        conceivably fail on hexadecimal or decimal machines without guard digits, but we know of
        none.  See SLASD3 for details.

        The code currently calls SLASDQ if singular values only are desired. However, it can be
        slightly modified to compute singular values using the divide and conquer method.
    \endverbatim

    * @param[in] uplo
              uplo is char* \n
              = 'U':  B is upper bidiagonal. \n
              = 'L':  B is lower bidiagonal. \n
    * @param[in] compq
              compq is char* \n
              Specifies whether singular vectors are to be computed
              as follows: \n
              = 'N':  Compute singular values only; \n
              = 'P':  Compute singular values and compute singular
              vectors in compact form; \n
              = 'I':  Compute singular values and singular vectors. \n
    * @param[in] n
              n is integer* \n
              The order of the matrix B.  n >= 0. \n
    * @param[in,out] d
              d is float/double array, dimension (n) \n
              On entry, the n diagonal elements of the bidiagonal matrix b. \n
              On exit, if info=0, the singular values of B. \n
    * @param[in,out] e
              e is float/double array, dimension (n-1) \n
              On entry, the elements of E contain the offdiagonal
              elements of the bidiagonal matrix whose SVD is desired. \n
              On exit, E has been destroyed. \n
    * @param[out] u
              U is float/double array, dimension (ldu,n) \n
              If  compq = 'I', then: \n
              On exit, if info = 0, U contains the left singular vectors
              of the bidiagonal matrix. \n
              For other values of compq, U is not referenced. \n
    * @param[in] ldu
              ldu is integer* \n
              The leading dimension of the array U.  ldu >= 1. \n
              If singular vectors are desired, then ldu >= fla_max( 1, n ). \n
    * @param[out] vt
              vt is float/double array, dimension (ldvt,n) \n
              If  compq = 'I', then: \n
              On exit, if info = 0, vt**T contains the right singular
              vectors of the bidiagonal matrix. \n
              For other values of compq, vt is not referenced. \n
    * @param[in] ldvt
              ldvt is integer* \n
              The leading dimension of the array vt.  ldvt >= 1. \n
              If singular vectors are desired, then ldvt >= fla_max( 1, n ). \n
    * @param[out] q
              q is float/double array, dimension (LDQ) \n
              If  compq = 'P', then: \n
              On exit, if info = 0, Q and iq contain the left
              and right singular vectors in a compact form,
              requiring O(N log N) space instead of 2*n**2. \n
              In particular, Q contains all the float data in
              LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
              words of memory, where SMLSIZ is returned by ILAENV and
              is equal to the maximum size of the subproblems at the
              bottom of the computation tree (usually about 25). \n
              For other values of compq, Q is not referenced. \n
    * @param[out] iq
              iq is integer array, dimension (LDIQ) \n
              If  compq = 'P', then: \n
              On exit, if info = 0, Q and iq contain the left
              and right singular vectors in a compact form,
              requiring O(N log N) space instead of 2*n**2. \n
              In particular, iq contains all integer data in
              LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
              words of memory, where SMLSIZ is returned by ILAENV and
              is equal to the maximum size of the subproblems at the
              bottom of the computation tree (usually about 25). \n
              For other values of compq, iq is not referenced. \n
    * @param[out]	WORK
              WORK is REAL array, dimension (MAX(1,LWORK)) \n
              If COMPQ = 'N' then LWORK >= (4 * N). \n
              If COMPQ = 'P' then LWORK >= (6 * N). \n
              If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N). \n
    * @param[out]	IWORK
              IWORK is INTEGER array, dimension (8*N) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit. \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n
              > 0:  The algorithm failed to compute a singular value.
                    The update process of divide and conquer failed. \n

    *     *  */
    template <typename T>
    void bdsdc(char *uplo, char *compq, integer *n, T *d, T *e, T *u, integer *ldu, T *vt,
               integer *ldvt, T *q, T *iq, T *work, integer *iwork, integer *info)
    {
        bdsdc(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info);
    }
    /** @}*/ // end of bdsdc

    /** @defgroup bdsvdx bdsvdx
     * @ingroup SVD
     * @{
     */
    /*! @brief BDSVDX computes the singular value decomposition (SVD) of a real  \n
     N-by-N (upper or lower) bidiagonal matrix B
 * @details
 * \b Purpose:
    \verbatim
     BDSVDX computes the singular value decomposition (SVD) of a real
     N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT,
     where S is a diagonal matrix with non-negative diagonal elements
     (the singular values of B), and U and VT are orthogonal matrices
     of left and right singular vectors, respectively.

     Given an upper bidiagonal B with diagonal D = [ d_1 d_2 ... d_N ]
     and superdiagonal E = [ e_1 e_2 ... e_N-1 ], SBDSVDX computes the
     singular value decompositon of B through the eigenvalues and
     eigenvectors of the N*2-by-N*2 tridiagonal matrix

           |  0  d_1                |
           | d_1  0  e_1            |
     TGK = |     e_1  0  d_2        |
           |         d_2  .   .     |
           |              .   .   . |

     If (s,u,v) is a singular triplet of B with ||u|| = ||v|| = 1, then
     (+/-s,q), ||q|| = 1, are eigenpairs of TGK, with q = P * ( u' +/-v') /
     sqrt(2) = ( v_1 u_1 v_2 u_2 ... v_n u_n) / sqrt(2), and
     P = [ e_{n+1} e_{1} e_{n+2} e_{2} ... ].

     Given a TGK matrix, one can either a) compute -s,-v and change signs
     so that the singular values (and corresponding vectors) are already in
     descending order (as in SGESVD/SGESDD) or b) compute s,v and reorder
     the values (and corresponding vectors). SBDSVDX implements a) by
     calling SSTEVX (bisection plus inverse iteration, to be replaced
     with a version of the Multiple Relative Robust Representation
     algorithm. (See P. Willems and B. Lang, A framework for the MR^3
     algorithm: theory and implementation, SIAM J. Sci. Comput.,
     35:740-766, 2013.)
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  B is upper bidiagonal; \n
          = 'L':  B is lower bidiagonal. \n
 * @param[in] JOBZ
          JOBZ is CHARACTER*1 \n
          = 'N':  Compute singular values only; \n
          = 'V':  Compute singular values and singular vectors. \n
 * @param[in] RANGE
          RANGE is CHARACTER*1 \n
          = 'A': all singular values will be found. \n
          = 'V': all singular values in the half-open interval [VL,VU)
                 will be found. \n
          = 'I': the IL-th through IU-th singular values will be found. \n
 * @param[in] N
          N is INTEGER \n
          The order of the bidiagonal matrix.  N >= 0. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The n diagonal elements of the bidiagonal matrix B. \n
 * @param[in] E
          E is REAL array, dimension (fla_max(1,N-1)) \n
          The (n-1) superdiagonal elements of the bidiagonal matrix
          B in elements 1 to N-1. \n
 * @param[in] VL
         VL is REAL \n
          If RANGE='V', the lower bound of the interval to
          be searched for singular values. VU > VL.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] VU
          VU is REAL \n
          If RANGE='V', the upper bound of the interval to
          be searched for singular values. VU > VL.
          Not referenced if RANGE = 'A' or 'I'. \n
 * @param[in] IL
          IL is INTEGER \n
          If RANGE='I', the index of the
          smallest singular value to be returned.
          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
          Not referenced if RANGE = 'A' or 'V'. \n
 * @param[in] IU
          IU is INTEGER
          If RANGE='I', the index of the
          largest singular value to be returned. \n
          1 <= IL <= IU <= min(M,N), if min(M,N) > 0. \n
          Not referenced if RANGE = 'A' or 'V'.
 * @param[out] NS
          NS is INTEGER \n
          The total number of singular values found.  0 <= NS <= N.
          If RANGE = 'A', NS = N, and if RANGE = 'I', NS = IU-IL+1. \n
 * @param[out] S
          S is REAL array, dimension (N) \n
          The first NS elements contain the selected singular values in
          ascending order. \n
 * @param[out] Z
          Z is REAL array, dimension (2*N,K) \n
          If JOBZ = 'V', then if INFO = 0 the first NS columns of Z
          contain the singular vectors of the matrix B corresponding to
          the selected singular values, with U in rows 1 to N and V
          in rows N+1 to N*2, i.e. \n
          Z = [ U ] \n
              [ V ] \n
          If JOBZ = 'N', then Z is not referenced. \n
          Note: The user must ensure that at least K = NS+1 columns are
          supplied in the array Z; if RANGE = 'V', the exact value of
          NS is not known in advance and an upper bound must be used. \n
 * @param[in] LDZ
          LDZ is INTEGER \n
          The leading dimension of the array Z. LDZ >= 1, and if
          JOBZ = 'V', LDZ >= fla_max(2,N*2). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (14*N) \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (12*N) \n
          If JOBZ = 'V', then if INFO = 0, the first NS elements of
          IWORK are zero. If INFO > 0, then IWORK contains the indices
          of the eigenvectors that failed to converge in DSTEVX. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n
          > 0:  if INFO = i, then i eigenvectors failed to converge
                   in SSTEVX. The indices of the eigenvectors
                   (as returned by SSTEVX) are stored in the
                   array IWORK. \n
                if INFO = N*2 + 1, an internal error occurred. \n

 *  * */
    template <typename T>
    void bdsvdx(char *uplo, char *jobz, char *range, integer *n, T *d, T *e, T vl, T vu,
                integer *il, integer *iu, integer *ns, T *s, T *z, integer *ldz, T *work,
                integer *iwork, integer *info)
    {
        bdsvdx(uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work, iwork, info);
    }
    /** @}*/ // end of bdsvdx

    /** @defgroup ggsvd3 ggsvd3
     * @ingroup SVD
     * @{
     */
    /*! @brief GGSVD3 computes the singular value decomposition (SVD) for OTHER matrices
 * @details
 * \b Purpose:
    \verbatim
      GGSVD3 computes the generalized singular value decomposition (GSVD)
      of an M-by-N real matrix A and P-by-N real matrix B:

            U**T*A*Q = D1*( 0 R),    V**T*B*Q = D2*( 0 R)

      where U, V and Q are orthogonal matrices.
      Let K+L = the effective numerical rank of the matrix (A**T,B**T)**T,
      then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and
      D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the
      following structures, respectively:

      If M-K-L >= 0,

                          K  L
             D1 =     K ( I  0)
                      L ( 0  C)
                  M-K-L ( 0  0)

                        K  L
             D2 =   L ( 0  S)
                  P-L ( 0  0)

                      N-K-L  K    L
        ( 0 R) = K (  0   R11  R12)
                  L (  0    0   R22)

      where

        C = diag( ALPHA(K+1), ... , ALPHA(K+L)),
        S = diag( BETA(K+1),  ... , BETA(K+L)),
        C**2 + S**2 = I.

        R is stored in A(1:K+L,N-K-L+1:N) on exit.

      If M-K-L < 0,

                        K M-K K+L-M
             D1 =   K ( I  0    0  )
                  M-K ( 0  C    0  )

                          K M-K K+L-M
             D2 =   M-K ( 0  S    0 )
                  K+L-M ( 0  0    I )
                    P-L ( 0  0    0 )

                         N-K-L  K   M-K  K+L-M
        ( 0 R) =     K ( 0    R11  R12  R13 )
                    M-K ( 0     0   R22  R23 )
                  K+L-M ( 0     0    0   R33 )

      where
        C = diag( ALPHA(K+1), ... , ALPHA(M)),
        S = diag( BETA(K+1),  ... , BETA(M)),
        C**2 + S**2 = I.
        (R11 R12 R13) is stored in A(1:M, N-K-L+1:N), and R33 is stored
        ( 0  R22 R23)
        in B(M-K+1:L,N+M-K-L+1:N) on exit.

      The routine computes C, S, R, and optionally the orthogonal
      transformation matrices U, V and Q.

      In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
      A and B implicitly gives the SVD of A*inv(B):
                           A*inv(B) = U*(D1*inv(D2))*V**T.
      If ( A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B is
      also equal to the CS decomposition of A and B. Furthermore, the GSVD
      can be used to derive the solution of the eigenvalue problem:
                           A**T*A x = lambda* B**T*B x.
      In some literature, the GSVD of A and B is presented in the form
                       U**T*A*X = ( 0 D1),   V**T*B*X = ( 0 D2)
      where U and V are orthogonal and X is nonsingular, D1 and D2 are
      ``diagonal''.  The former GSVD form can be converted to the latter
      form by taking the nonsingular matrix X as
                           X = Q*( I   0   )
                                 ( 0 inv(R)).
    \endverbatim

 * @param[in] JOBU
          JOBU is CHARACTER*1 \n
          = 'U':  Orthogonal matrix U is computed; \n
          = 'N':  U is not computed. \n
 * @param[in] JOBV
          JOBV is CHARACTER*1 \n
          = 'V':  Orthogonal matrix V is computed; \n
          = 'N':  V is not computed. \n
 * @param[in] JOBQ
          JOBQ is CHARACTER*1 \n
          = 'Q':  Orthogonal matrix Q is computed; \n
          = 'N':  Q is not computed. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrices A and B.  N >= 0. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows of the matrix B.  P >= 0. \n
 * @param[out] K
          K is INTEGER \n
 * @param[out] L
          L is INTEGER \n
          On exit, K and L specify the dimension of the subblocks
          described in Purpose.
          K + L = effective numerical rank of (A**T,B**T)**T. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, A contains the triangular matrix R, or part of R.
          See Purpose for details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the P-by-N matrix B. \n
          On exit, B contains the triangular matrix R if M-K-L < 0.
          See Purpose for details. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,P). \n
 * @param[out] ALPHA
          ALPHA is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, ALPHA and BETA contain the generalized singular
          value pairs of A and B; \n
            ALPHA(1:K) = 1, \n
            BETA(1:K)  = 0, \n
          and if M-K-L >= 0, \n
            ALPHA(K+1:K+L) = C, \n
            BETA(K+1:K+L)  = S, \n
          or if M-K-L < 0, \n
            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0 \n
            BETA(K+1:M) =S, BETA(M+1:K+L) =1 \n
          and \n
            ALPHA(K+L+1:N) = 0 \n
            BETA(K+L+1:N)  = 0 \n
 * @param[out] U
          U is REAL array, dimension (LDU,M) \n
          If JOBU = 'U', U contains the M-by-M orthogonal matrix U. \n
          If JOBU = 'N', U is not referenced. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U. LDU >= fla_max(1,M) if
          JOBU = 'U'; LDU >= 1 otherwise. \n
 * @param[out] V
          V is REAL array, dimension (LDV,P) \n
          If JOBV = 'V', V contains the P-by-P orthogonal matrix V. \n
          If JOBV = 'N', V is not referenced. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. LDV >= fla_max(1,P) if
          JOBV = 'V'; LDV >= 1 otherwise. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ,N) \n
          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q. \n
          If JOBQ = 'N', Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. LDQ >= fla_max(1,N) if
          JOBQ = 'Q'; LDQ >= 1 otherwise. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
 \n
          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
          On exit, IWORK stores the sorting information. More
          precisely, the following loop will sort ALPHA \n
             for I = K+1, min(M,K+L) \n
                 swap ALPHA(I) and ALPHA(IWORK(I)) \n
             endfor \n
          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, the Jacobi-type procedure failed to
                converge.  For further details, see subroutine STGSJA. \n

 *  * */
    template <typename T>
    void ggsvd3(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k,
                integer *l, T *a, integer *lda, T *b, integer *ldb, T *alpha, T *beta, T *u,
                integer *ldu, T *v, integer *ldv, T *q, integer *ldq, T *work, integer *lwork,
                integer *iwork, integer *info)
    {
        ggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq,
               work, lwork, iwork, info);
    }
    template <typename T, typename Ta>
    void ggsvd3(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k,
                integer *l, T *a, integer *lda, T *b, integer *ldb, Ta *alpha, Ta *beta, T *u,
                integer *ldu, T *v, integer *ldv, T *q, integer *ldq, T *work, integer *lwork,
                Ta *rwork, integer *iwork, integer *info)
    {
        ggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq,
               work, lwork, rwork, iwork, info);
    }
    /** @}*/ // end of ggsvd3

    /** @defgroup gebrd gebrd
     * @ingroup SVD
     * @{
     */
    /*! @brief Reduction to bidiagonal form
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a general real m-by-n matrix a to upper or lower bidiagonal form B by an
        orthogonal transformation: Q**T * A * P = B.

        If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows in the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns in the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the m-by-n general matrix to be reduced. \n
              On exit, \n
              if m >= n, the diagonal and the first superdiagonal are
              overwritten with the upper bidiagonal matrix b; the
              elements below the diagonal, with the array tauq, represent
              the orthogonal matrix Q as a product of elementary
              reflectors, and the elements above the first superdiagonal,
              with the array taup, represent the orthogonal matrix P as
              a product of elementary reflectors; \n
              if m < n, the diagonal and the first subdiagonal are
              overwritten with the lower bidiagonal matrix b; the
              elements below the first subdiagonal, with the array tauq,
              represent the orthogonal matrix Q as a product of
              elementary reflectors, and the elements above the diagonal,
              with the array taup, represent the orthogonal matrix P as
              a product of elementary reflectors. \n
              See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,m). \n
    * @param[out] d
              d is float/double array, dimension (min(m,n)) \n
              The diagonal elements of the bidiagonal matrix b:
              D(i) = A(i,i).
    * @param[out] e
              e is float/double array, dimension (min(m,n)-1) \n
              The off-diagonal elements of the bidiagonal matrix b: \n
              if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; \n
              if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1. \n
    * @param[out] tauq
              tauq is float/double array, dimension (min(m,n)) \n
              The scalar factors of the elementary reflectors which
              represent the orthogonal matrix Q. \n
              See Further Details.
    * @param[out] taup
              taup is float/double array, dimension (min(m,n)) \n
              The scalar factors of the elementary reflectors which
              represent the orthogonal matrix P. \n
              See Further Details.
              *
              * \n
              * **Further Details**
              * \verbatim
                      The matrices Q and P are represented as products of elementary reflectors:

                      If m >= n,

                      Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)

                      Each H(i) and G(i) has the form:

                      H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

                      where tauq and taup are real scalars, and V and u are real vectors;
                      V(1:i-1) = 0, V(i) = 1, and V(i+1:m) is stored on exit in A(i+1:m,i);
                      u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
                      tauq is stored in tauq(i) and taup in taup(i).

                      If m < n,

                      Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)

                      Each H(i) and G(i) has the form:

                      H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

                      where tauq and taup are real scalars, and V and u are real vectors;
                      V(1:i) = 0, V(i+1) = 1, and V(i+2:m) is stored on exit in A(i+2:m,i);
                      u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
                      tauq is stored in tauq(i) and taup in taup(i).

                      The contents of A on exit are illustrated by the following examples:

                      m = 6 and n = 5 (m > n): m = 5 and n = 6 (m < n):

                      (  d   e   u1  u1  u1)  (  d   u1  u1  u1  u1  u1)
                      (  v1  d   e   u2  u2)  (  e   d   u2  u2  u2  u2)
                      (  v1  v2  d   e   u3)  (  v1  e   d   u3  u3  u3)
                      (  v1  v2  v3  d   e )  (  v1  v2  e   d   u4  u4)
                      (  v1  v2  v3  v4  d )  (  v1  v2  v3  e   d   u5)
                      (  v1  v2  v3  v4  v5)

                      where d and e denote diagonal and off-diagonal elements of B, vi denotes an
    element of the vector defining H(i), and ui an element of the vector defining G(i). \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The length of the array WORK.  LWORK >= fla_max(1,M,N). \n
              For optimum performance LWORK >= (M+N)*NB, where NB
              is the optimal blocksize. \n

              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit. \n
              < 0:  if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void gebrd(integer *m, integer *n, T *a, integer *lda, T *d, T *e, T *tauq, T *taup, T *work,
               integer *lwork, integer *info)
    {
        gebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
    }
    template <typename T, typename Ta>
    void gebrd(integer *m, integer *n, T *a, integer *lda, Ta *d, Ta *e, T *tauq, T *taup, T *work,
               integer *lwork, integer *info)
    {
        gebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info);
    }
    /** @}*/ // end of gebrd

    /** @defgroup gebd2 gebd2
     * @ingroup SVD
     * @{
     */
    /*! @brief Reduction to bidiagonal form (unblocked algorithm)
    *
    * @details
    * \b Purpose:
    * \verbatim
        Reduction of a general real m-by-n matrix a to upper or lower bidiagonal form B by an
        orthogonal transformation: Q**T * A * P = B.

        If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
    \endverbatim

    * @param[in] m
              m is integer* \n
              The number of rows in the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns in the matrix a.  n >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the m-by-n general matrix to be reduced. \n
              On exit, \n
              if m >= n, the diagonal and the first superdiagonal are
              overwritten with the upper bidiagonal matrix b; the
              elements below the diagonal, with the array tauq, represent
              the orthogonal matrix Q as a product of elementary
              reflectors, and the elements above the first superdiagonal,
              with the array taup, represent the orthogonal matrix P as
              a product of elementary reflectors; \n
              if m < n, the diagonal and the first subdiagonal are
              overwritten with the lower bidiagonal matrix b; the
              elements below the first subdiagonal, with the array tauq,
              represent the orthogonal matrix Q as a product of
              elementary reflectors, and the elements above the diagonal,
              with the array taup, represent the orthogonal matrix P as
              a product of elementary reflectors. \n
              See Further Details. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,m). \n
    * @param[out] d
              d is float/double array, dimension (min(m,n)) \n
              The diagonal elements of the bidiagonal matrix b: \n
              D(i) = A(i,i).
    * @param[out] e
              e is float/double array, dimension (min(m,n)-1) \n
              The off-diagonal elements of the bidiagonal matrix b: \n
              if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1; \n
              if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
    * @param[out] tauq
              tauq is float/double array, dimension (min(m,n)) \n
              The scalar factors of the elementary reflectors which
              represent the orthogonal matrix Q.
              See Further Details. \n
    * @param[out] taup
              taup is float/double array, dimension (min(m,n)) \n
              The scalar factors of the elementary reflectors which
              represent the orthogonal matrix P. \n
              See Further Details.
              *
              * \n
              * **Further Details**
              * \verbatim
                      The matrices Q and P are represented as products of elementary reflectors:

                      If m >= n,

                      Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)

                      Each H(i) and G(i) has the form:

                      H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

                      where tauq and taup are real scalars, and V and u are real vectors;
                      V(1:i-1) = 0, V(i) = 1, and V(i+1:m) is stored on exit in A(i+1:m,i);
                      u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
                      tauq is stored in tauq(i) and taup in taup(i).

                      If m < n,

                      Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)

                      Each H(i) and G(i) has the form:

                      H(i) = I - tauq * V * V**T  and G(i) = I - taup * u * u**T

                      where tauq and taup are real scalars, and V and u are real vectors;
                      V(1:i) = 0, V(i+1) = 1, and V(i+2:m) is stored on exit in A(i+2:m,i);
                      u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
                      tauq is stored in tauq(i) and taup in taup(i).

                      The contents of A on exit are illustrated by the following examples:

                      m = 6 and n = 5 (m > n): m = 5 and n = 6 (m < n):

                      (  d   e   u1  u1  u1)  (  d   u1  u1  u1  u1  u1)
                      (  v1  d   e   u2  u2)  (  e   d   u2  u2  u2  u2)
                      (  v1  v2  d   e   u3)  (  v1  e   d   u3  u3  u3)
                      (  v1  v2  v3  d   e )  (  v1  v2  e   d   u4  u4)
                      (  v1  v2  v3  v4  d )  (  v1  v2  v3  e   d   u5)
                      (  v1  v2  v3  v4  v5)

                      where d and e denote diagonal and off-diagonal elements of B, vi denotes an
    element of the vector defining H(i), and ui an element of the vector defining G(i). \endverbatim
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (fla_max(M,N)) \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0: successful exit \n
              < 0: if INFO = -i, the i-th argument had an illegal value. \n

    *     *  */
    template <typename T>
    void gebd2(integer *m, integer *n, T *a, integer *lda, T *d, T *e, T *tauq, T *taup, T *work,
               integer *info)
    {
        gebd2(m, n, a, lda, d, e, tauq, taup, work, info);
    }
    template <typename T, typename Ta>
    void gebd2(integer *m, integer *n, T *a, integer *lda, Ta *d, Ta *e, T *tauq, T *taup, T *work,
               integer *info)
    {
        gebd2(m, n, a, lda, d, e, tauq, taup, work, info);
    }
    /** @}*/ // end of gebd2

    /** @defgroup labrd labrd
     * @ingroup SVD
     * @{
     */
    /*! @brief LABRD reduces the first nb rows and columns of a general matrix to a bidiagonal form

 * @details
 * \b Purpose:
    \verbatim
    LABRD reduces the first NB rows and columns of a real general
    m by n matrix A to upper or lower bidiagonal form by an orthogonal
    transformation Q**T * A * P, and   returns the matrices X and Y which
    are needed to apply the transformation to the unreduced part of A.

    If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
    bidiagonal form.

    This is an auxiliary routine called by GEBRD
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows in the matrix A. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns in the matrix A. \n
 * @param[in] NB
          NB is INTEGER \n
          The number of leading rows and columns of A to be reduced. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the m by n general matrix to be reduced.
          On exit, the first NB rows and columns of the matrix are
          overwritten; the rest of the array is unchanged. \n
          If m >= n, elements on and below the diagonal in the first NB
            columns, with the array TAUQ, represent the orthogonal
            matrix Q as a product of elementary reflectors; and
            elements above the diagonal in the first NB rows, with the
            array TAUP, represent the orthogonal matrix P as a product
            of elementary reflectors. \n
          If m < n, elements below the diagonal in the first NB
            columns, with the array TAUQ, represent the orthogonal
            matrix Q as a product of elementary reflectors, and
            elements on and above the diagonal in the first NB rows,
            with the array TAUP, represent the orthogonal matrix P as
            a product of elementary reflectors. \n
          See Further Details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] D
          D is REAL array, dimension (NB) \n
          The diagonal elements of the first NB rows and columns of
          the reduced matrix.  D(i) = A(i,i). \n
 * @param[out] E
          E is REAL array, dimension (NB) \n
          The off-diagonal elements of the first NB rows and columns of
          the reduced matrix. \n
 * @param[out] TAUQ
          TAUQ is REAL array, dimension (NB) \n
          The scalar factors of the elementary reflectors which
          represent the orthogonal matrix Q. See Further Details. \n
 * @param[out] TAUP
          TAUP is REAL array, dimension (NB) \n
          The scalar factors of the elementary reflectors which
          represent the orthogonal matrix P. See Further Details. \n
 * @param[out] X
          X is REAL array, dimension (LDX,NB) \n
          The m-by-nb matrix X required to update the unreduced part
          of A. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X. LDX >= fla_max(1,M). \n
 * @param[out] Y
          Y is REAL array, dimension (LDY,NB) \n
          The n-by-nb matrix Y required to update the unreduced part
          of A. \n
 * @param[in] LDY
          LDY is INTEGER \n
          The leading dimension of the array Y. LDY >= fla_max(1,N).  \n

 * */
    template <typename T, typename Ta>
    void labrd(integer *m, integer *n, integer *nb, T *a, integer *lda, Ta *d, Ta *e, T *tauq,
               T *taup, T *x, integer *ldx, T *y, integer *ldy)
    {
        labrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy);
    }
    /** @}*/ // end of labrd

    /** @defgroup gbbrd gbbrd
     * @ingroup SVD
     * @{
     */
    /*! @brief Reduces a general band matrix to bidiagonal form
 *
 * @details
 * \b Purpose :

   \verbatim
    GBBRD reduces a real general m-by-n band matrix A to upper
     bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.

     The routine computes B, and optionally forms Q or P**T, or computes
    Q**T*C for a given matrix C.
    \endverbatim

 * @param[in] VECT
          VECT is CHARACTER*1 \n
          Specifies whether or not the matrices Q and P**T are to be
          formed. \n
          = 'N': do not form Q or P**T; \n
          = 'Q': form Q only; \n
          = 'P': form P**T only; \n
          = 'B': form both. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] NCC
          NCC is INTEGER \n
          The number of columns of the matrix C.  NCC >= 0. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals of the matrix A. KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals of the matrix A. KU >= 0. \n
 * @param[in,out] AB
          AB is REAL array, dimension (LDAB,N) \n
          On entry, the m-by-n band matrix A, stored in rows 1 to
          KL+KU+1. The j-th column of A is stored in the j-th column of
          the array AB as follows: \n
          AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(m,j+kl).
          On exit, A is overwritten by values generated during the
          reduction. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array A. LDAB >= KL+KU+1. \n
 * @param[out] D
          D is REAL array, dimension (min(M,N)) \n
          The diagonal elements of the bidiagonal matrix B. \n
 * @param[out] E
          E is REAL array, dimension (min(M,N)-1) \n
          The superdiagonal elements of the bidiagonal matrix B. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ,M) \n
          If VECT = 'Q' or 'B', the m-by-m orthogonal matrix Q. \n
          If VECT = 'N' or 'P', the array Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. \n
          LDQ >= fla_max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise. \n
 * @param[out] PT
          PT is REAL array, dimension (LDPT,N) \n
          If VECT = 'P' or 'B', the n-by-n orthogonal matrix P'. \n
          If VECT = 'N' or 'Q', the array PT is not referenced. \n
 * @param[in] LDPT
          LDPT is INTEGER \n
          The leading dimension of the array PT. \n
          LDPT >= fla_max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise. \n
 * @param[in,out] C
          C is REAL array, dimension (LDC,NCC)
          On entry, an m-by-ncc matrix C. \n
          On exit, C is overwritten by Q**T*C. \n
          C is not referenced if NCC = 0. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. \n
          LDC >= fla_max(1,M) if NCC > 0; LDC >= 1 if NCC = 0. \n
 * @param[out]	WORK
          WORK is COMPLEX array, dimension (fla_max(M,N)) \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (fla_max(M,N)) \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n

 * * */
    template <typename T>
    void gbbrd(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, T *ab,
               integer *ldab, T *d, T *e, T *q, integer *ldq, T *pt, integer *ldpt, T *c,
               integer *ldc, T *work, integer *info)
    {
        gbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, info);
    }
    template <typename T, typename Ta>
    void gbbrd(char *vect, integer *m, integer *n, integer *ncc, integer *kl, integer *ku, T *ab,
               integer *ldab, Ta *d, Ta *e, T *q, integer *ldq, T *pt, integer *ldpt, T *c,
               integer *ldc, T *work, Ta *rwork, integer *info)
    {
        gbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, rwork, info);
    }
    /** @}*/ // end of gbbrd

    /** @defgroup ungbr ungbr
     * @ingroup SVD
     * @{
     */
    /*! @brief Form Q from bidiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
        Generate one of the complex unitary matrices Q or P**H determined by CGEBRD when reducing
        a complex matrix a to bidiagonal form: A = Q * B * P**H.  Q and P**H are defined as
        products of elementary reflectors H(i) or G(i) respectively.

        If vect = 'Q', A is assumed to have been an M-by-K matrix, and Q is of order M:
        if m >= k, Q = H(1) H(2) . . . H(k) and CUNGBR returns the first n columns of Q,
        where m >= n >= k;
        if m < k, Q = H(1) H(2) . . . H(m-1) and CUNGBR returns Q as an M-by-M matrix.

        If vect = 'P', A is assumed to have been a K-by-N matrix, and P**H is of order N:
        if k < n, P**H = G(k) . . . G(2) G(1) and CUNGBR returns the first m rows of P**H,
        where n >= m >= k;
        if k >= n, P**H = G(n-1) . . . G(2) G(1) and CUNGBR returns P**H as an n-by-n matrix.
    \endverbatim

    * @param[in] vect
              vect is char* \n
              Specifies whether the matrix Q or the matrix P**H is
              required, as defined in the transformation applied by CGEBRD: \n
              = 'Q':  generate Q; \n
              = 'P':  generate P**H. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix Q or P**H to be returned. \n
              m >= 0.
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix Q or P**H to be returned. \n
              n >= 0. \n
              If vect = 'Q', m >= n >= min(m,k); \n
              if vect = 'P', n >= m >= min(n,k). \n
    * @param[in] k
              k is integer* \n
              If vect = 'Q', the number of columns in the original M-by-K
              matrix reduced by CGEBRD. \n
              If vect = 'P', the number of rows in the original K-by-N
              matrix reduced by CGEBRD. \n
              k >= 0.
    * @param[in,out] a
              a is COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On entry, the vectors which define the elementary reflectors,
              as returned by CGEBRD. \n
              On exit, the m-by-n matrix Q or P**H. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. lda >= fla_max(1,m).
    * @param[in] tau
              tau is COMPLEX/COMPLEX*16 array, dimension \n
              (min(m,k)) if vect = 'Q' \n
              (min(n,k)) if vect = 'P' \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i) or G(i), which determines Q or P**H, as
              returned by CGEBRD in its array argument tauq or taup. \n
    * @param[out]	WORK
              WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= fla_max(1,min(M,N)).
              For optimum performance LWORK >= min(M,N)*NB, where NB
              is the optimal blocksize. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void ungbr(char *vect, integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        ungbr(vect, m, n, k, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of ungbr

    /** @defgroup orgbr orgbr
     * @ingroup SVD
     * @{
     */
    /*! @brief Form Q from bidiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
        Generate one of the real orthogonal matrices Q or P**T determined by SGEBRD when reducing
        a real matrix a to bidiagonal form: A = Q * B * P**T.  Q and P**T are defined as products
        of elementary reflectors H(i) or G(i) respectively.

        If vect = 'Q', A is assumed to have been an M-by-K matrix, and Q is of order M:
        if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n columns of Q,
        where m >= n >= k;
        if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an M-by-M matrix.

        If vect = 'P', A is assumed to have been a K-by-N matrix, and P**T is of order N:
        if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m rows of P**T,
        where n >= m >= k;
        if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as an n-by-n matrix.
    \endverbatim

    * @param[in] vect
              vect is char* \n
              Specifies whether the matrix Q or the matrix P**T is
              required, as defined in the transformation applied by SGEBRD:
              = 'Q':  generate Q; \n
              = 'P':  generate P**T. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix Q or P**T to be returned.
              m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix Q or P**T to be returned.
              n >= 0. \n
              If vect = 'Q', m >= n >= min(m,k); \n
              if vect = 'P', n >= m >= min(n,k). \n
    * @param[in] k
              k is integer* \n
              If vect = 'Q', the number of columns in the original M-by-K
              matrix reduced by SGEBRD. \n
              If vect = 'P', the number of rows in the original K-by-N
              matrix reduced by SGEBRD. \n
              k >= 0. \n
    * @param[in,out] a
              a is float/double array, dimension (lda,n) \n
              On entry, the vectors which define the elementary reflectors,
              as returned by SGEBRD. \n
              On exit, the m-by-n matrix Q or P**T. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. lda >= fla_max(1,m). \n
    * @param[in] tau
              tau is float/double array, dimension \n
              (min(m,k)) if vect = 'Q' \n
              (min(n,k)) if vect = 'P' \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i) or G(i), which determines Q or P**T, as
              returned by SGEBRD in its array argument tauq or taup. \n
    * @param[out]	WORK
              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. LWORK >= fla_max(1,min(M,N)).
              For optimum performance LWORK >= min(M,N)*NB, where NB
              is the optimal blocksize. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void orgbr(char *vect, integer *m, integer *n, integer *k, T *a, integer *lda, T *tau, T *work,
               integer *lwork, integer *info)
    {
        orgbr(vect, m, n, k, a, lda, tau, work, lwork, info);
    }
    /** @}*/ // end of orgbr

    /** @defgroup ormbr ormbr
     * @ingroup SVD
     * @{
     */
    /*! @brief Apply Q or Q' from bidiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
        If vect = 'Q', SORMBR overwrites the general real m-by-n matrix c with
            side = 'L'  side = 'R'
            trans = 'N':   Q * C C * Q
            trans = 'T':   Q**T * C C * Q**T

        If vect = 'P', SORMBR overwrites the general real m-by-n matrix c with
            side = 'L'  side = 'R'
            trans = 'N':   P * C C * P
            trans = 'T':   P**T * C C * P**T

        Here Q and P**T are the orthogonal matrices determined by SGEBRD when reducing a real
        matrix a to bidiagonal form: A = Q * B * P**T. Q and P**T are defined as products of
        elementary reflectors H(i) and G(i) respectively.

        Let nq = m if side = 'L' and nq = n if side = 'R'. Thus nq is the order of the orthogonal
        matrix Q or P**T that is applied.

        If vect = 'Q', A is assumed to have been an NQ-by-K matrix:
            if nq >= k, Q = H(1) H(2) . . . H(k);
            if nq < k, Q = H(1) H(2) . . . H(nq-1).

        If vect = 'P', A is assumed to have been a K-by-NQ matrix:
            if k < nq, P = G(1) G(2) . . . G(k);
            if k >= nq, P = G(1) G(2) . . . G(nq-1).
    \endverbatim

    * @param[in] vect
              vect is char* \n
              = 'Q': apply Q or Q**T; \n
              = 'P': apply P or P**T. \n
    * @param[in] side
              side is char* \n
              = 'L': apply Q, Q**T, P or P**T from the Left; \n
              = 'R': apply Q, Q**T, P or P**T from the Right. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q  or P; \n
              = 'T':  Transpose, apply Q**T or P**T. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] k
              k is integer* \n
              If vect = 'Q', the number of columns in the original
              matrix reduced by SGEBRD. \n
              If vect = 'P', the number of rows in the original
              matrix reduced by SGEBRD. \n
              k >= 0. \n
    * @param[in] a
              a is float/double array, dimension \n
              (lda,min(nq,k)) if vect = 'Q' \n
              (lda,nq)  if vect = 'P' \n
              The vectors which define the elementary reflectors H(i) and
              G(i), whose products determine the matrices Q and P, as
              returned by SGEBRD. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              If vect = 'Q', lda >= fla_max(1,nq); \n
              if vect = 'P', lda >= fla_max(1,min(nq,k)). \n
    * @param[in] tau
              tau is float/double array, dimension (min(nq,k)) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i) or G(i) which determines Q or P, as returned
              by SGEBRD in the array argument tauq or taup. \n
    * @param[in,out] c
              c is float/double array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
              or P*C or P**T*C or C*P or C*P**T. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. ldc >= fla_max(1,m). \n
    * @param[out]	WORK
              WORK is REAL array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              If SIDE = 'L', LWORK >= fla_max(1,N); \n
              if SIDE = 'R', LWORK >= fla_max(1,M). \n
              For optimum performance LWORK >= N*NB if SIDE = 'L', and
              LWORK >= M*NB if SIDE = 'R', where NB is the optimal
              blocksize. \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void ormbr(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        ormbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @}*/ // end of ormbr

    /** @defgroup unmbr unmbr
     * @ingroup SVD
     * @{
     */
    /*! @brief Apply Q or Q' from bidiagonal reduction
    *
    * @details
    * \b Purpose:
    * \verbatim
        If vect = 'Q', CUNMBR overwrites the general complex m-by-n matrix c with
        side = 'L'  side = 'R'
            trans = 'N':   Q * C C * Q
            trans = 'C':   Q**H * C C * Q**H

        If vect = 'P', CUNMBR overwrites the general complex m-by-n matrix c with
        side = 'L'  side = 'R'
            trans = 'N':   P * C C * P
            trans = 'C':   P**H * C C * P**H

        Here Q and P**H are the orthogonal matrices determined by CGEBRD when reducing a complex
        matrix a to bidiagonal form: A = Q * B * P**H. Q and P**H are defined as products of
        elementary reflectors H(i) and G(i) respectively.

        Let nq = m if side = 'L' and nq = n if side = 'R'. Thus nq is the order of the orthogonal
        matrix Q or P**H that is applied.

        If vect = 'Q', A is assumed to have been an NQ-by-K matrix:
            if nq >= k, Q = H(1) H(2) . . . H(k);
            if nq < k, Q = H(1) H(2) . . . H(nq-1).

        If vect = 'P', A is assumed to have been a K-by-NQ matrix:
            if k < nq, P = G(1) G(2) . . . G(k);
            if k >= nq, P = G(1) G(2) . . . G(nq-1).
    \endverbatim

    * @param[in] vect
              vect is char* \n
              = 'Q': apply Q or Q**H; \n
              = 'P': apply P or P**H. \n
    * @param[in] side
              side is char* \n
              = 'L': apply Q, Q**H, P or P**H from the Left; \n
              = 'R': apply Q, Q**H, P or P**H from the Right. \n
    * @param[in] trans
              trans is char* \n
              = 'N':  No transpose, apply Q  or P; \n
              = 'C':  Transpose, apply Q**H or P**H. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix c. m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix c. n >= 0. \n
    * @param[in] k
              k is integer* \n
              If vect = 'Q', the number of columns in the original
              matrix reduced by CGEBRD. \n
              If vect = 'P', the number of rows in the original
              matrix reduced by CGEBRD. \n
              k >= 0. \n
    * @param[in] a
              a is COMPLEX/COMPLEX*16 array, dimension \n
              (lda,min(nq,k)) if vect = 'Q' \n
              (lda,nq)  if vect = 'P' \n
              The vectors which define the elementary reflectors H(i) and
              G(i), whose products determine the matrices Q and P, as
              returned by CGEBRD. \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a. \n
              If vect = 'Q', lda >= fla_max(1,nq); \n
              if vect = 'P', lda >= fla_max(1,min(nq,k)). \n
    * @param[in] tau
              tau is COMPLEX/COMPLEX*16 array, dimension (min(nq,k)) \n
              tau(i) must contain the scalar factor of the elementary
              reflector H(i) or G(i) which determines Q or P, as returned
              by CGEBRD in the array argument tauq or taup. \n
    * @param[in,out] c
              c is COMPLEX/COMPLEX*16 array, dimension (ldc,n) \n
              On entry, the m-by-n matrix c. \n
              On exit, c is overwritten by Q*C or Q**H*C or C*Q**H or C*Q
              or P*C or P**H*C or C*P or C*P**H. \n
    * @param[in] ldc
              ldc is integer* \n
              The leading dimension of the array c. ldc >= fla_max(1,m). \n
    * @param[out]	WORK
              WORK is COMPLEX array, dimension (MAX(1,LWORK)) \n
              On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
    * @param[in]	LWORK
              LWORK is INTEGER \n
              The dimension of the array WORK. \n
              If SIDE = 'L', LWORK >= fla_max(1,N); \n
              if SIDE = 'R', LWORK >= fla_max(1,M); \n
              if N = 0 or M = 0, LWORK >= 1. \n
              For optimum performance LWORK >= fla_max(1,N*NB) if SIDE = 'L',
              and LWORK >= fla_max(1,M*NB) if SIDE = 'R', where NB is the
              optimal blocksize. (NB = 0 if M = 0 or N = 0.) \n
 \n
              If LWORK = -1, then a workspace query is assumed; the routine
              only calculates the optimal size of the WORK array, returns
              this value as the first entry of the WORK array, and no error
              message related to LWORK is issued by XERBLA. \n
    * @param[out]	INFO
              INFO is INTEGER \n
              = 0:  successful exit \n
              < 0:  if INFO = -i, the i-th argument had an illegal value \n

    *     *  */
    template <typename T>
    void unmbr(char *vect, char *side, char *trans, integer *m, integer *n, integer *k, T *a,
               integer *lda, T *tau, T *c, integer *ldc, T *work, integer *lwork, integer *info)
    {
        unmbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
    }
    /** @}*/ // end of unmbr

    /** @defgroup gesvj0 gesvj0
     * @ingroup SVD
     * @{
     */
    /*! @brief GSVJ0 pre-processor for the routine gesvj

 * @details
 * \b Purpose:
    \verbatim
    GSVJ0 is called from GESVJ as a pre-processor and that is its main
    purpose. It applies Jacobi rotations in the same way as GESVJ does, but
    it does not check convergence (stopping criterion). Few tuning
    parameters (marked by [TP]) are available for the implementer.
    \endverbatim

 * @param[in] JOBV
          JOBV is CHARACTER*1 \n
          Specifies whether the output from this procedure is used
          to compute the matrix V: \n
          = 'V': the product of the Jacobi rotations is accumulated
                 by postmulyiplying the N-by-N array V.
                (See the description of V.) \n
          = 'A': the product of the Jacobi rotations is accumulated
                 by postmulyiplying the MV-by-N array V.
                (See the descriptions of MV and V.) \n
          = 'N': the Jacobi rotations are not accumulated. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the input matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the input matrix A.
          M >= N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, M-by-N matrix A, such that A*diag(D) represents
          the input matrix. \n
          On exit,
          A_onexit * D_onexit represents the input matrix A*diag(D)
          post-multiplied by a sequence of Jacobi rotations, where the
          rotation threshold and the total number of sweeps are given in
          TOL and NSWEEP, respectively. \n
          (See the descriptions of D, TOL and NSWEEP.) \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          The array D accumulates the scaling factors from the fast scaled
          Jacobi rotations. \n
          On entry, A*diag(D) represents the input matrix. \n
          On exit, A_onexit*diag(D_onexit) represents the input matrix
          post-multiplied by a sequence of Jacobi rotations, where the
          rotation threshold and the total number of sweeps are given in
          TOL and NSWEEP, respectively. \n
          (See the descriptions of A, TOL and NSWEEP.) \n
 * @param[in,out] SVA
          SVA is REAL array, dimension (N) \n
          On entry, SVA contains the Euclidean norms of the columns of
          the matrix A*diag(D). \n
          On exit, SVA contains the Euclidean norms of the columns of
          the matrix onexit*diag(D_onexit). \n
 * @param[in] MV
          MV is INTEGER \n
          If JOBV = 'A', then MV rows of V are post-multipled by a
                           sequence of Jacobi rotations. \n
          If JOBV = 'N',  then MV is not referenced. \n
 * @param[in,out] V
          V is REAL array, dimension (LDV,N) \n
          If JOBV = 'V' then N rows of V are post-multipled by a
                           sequence of Jacobi rotations. \n
          If JOBV = 'A' then MV rows of V are post-multipled by a
                           sequence of Jacobi rotations. \n
          If JOBV = 'N',  then V is not referenced. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V, LDV >= 1. \n
          If JOBV = 'V', LDV >= N. \n
          If JOBV = 'A', LDV >= MV. \n
 * @param[in] EPS
          EPS is REAL \n
          EPS = SLAMCH('Epsilon') \n
 * @param[in] SFMIN
          SFMIN is REAL \n
          SFMIN = SLAMCH('Safe Minimum') \n
 * @param[in] TOL
          TOL is REAL \n
          TOL is the threshold for Jacobi rotations. For a pair
          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is
          applied only if ABS(COS(angle(A(:,p),A(:,q)))) > TOL. \n
 * @param[in] NSWEEP
          NSWEEP is INTEGER \n
          NSWEEP is the number of sweeps of Jacobi rotations to be
          performed. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          LWORK is the dimension of WORK. LWORK >= M. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, then the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gsvj0(char *jobv, integer *m, integer *n, T *a, integer *lda, T *d, T *sva, integer *mv,
               T *v, integer *ldv, T *eps, T *sfmin, T *tol, integer *nsweep, T *work,
               integer *lwork, integer *info)
    {
        gsvj0(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
    }
    template <typename T, typename Ta>
    void gsvj0(char *jobv, integer *m, integer *n, T *a, integer *lda, T *d, Ta *sva, integer *mv,
               T *v, integer *ldv, Ta *eps, Ta *sfmin, Ta *tol, integer *nsweep, T *work,
               integer *lwork, integer *info)
    {
        gsvj0(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork, info);
    }

    /** @}*/ // end of gesvj0

    /** @defgroup gesvj1 gesvj1
     * @ingroup SVD
     * @{
     */
    /*! @brief GSVJ1 pre-processor for the routine gesvj, applies Jacobi rotations targeting only
 particular pivots

 * @details
 * \b Purpose:
    \verbatim
    GSVJ1 is called from SGESVJ as a pre-processor and that is its main
    purpose. It applies Jacobi rotations in the same way as SGESVJ does, but
    it targets only particular pivots and it does not check convergence
    (stopping criterion). Few tunning parameters (marked by [TP]) are
    available for the implementer.

    Further Details
    ~~~~~~~~~~~~~~~
    SGSVJ1 applies few sweeps of Jacobi rotations in the column space of
    the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
    off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
    block-entries (tiles) of the (1,2) off-diagonal block are marked by the
    [x]'s in the following scheme:

       | *  *  * [x] [x] [x]|
       | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
       | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
       |[x] [x] [x] *  *  * |
       |[x] [x] [x] *  *  * |
       |[x] [x] [x] *  *  * |

    In terms of the columns of A, the first N1 columns are rotated 'against'
    the remaining N-N1 columns, trying to increase the angle between the
    corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
    tiled using quadratic tiles of side KBL. Here, KBL is a tunning parameter.
    The number of sweeps is given in NSWEEP and the orthogonality threshold
    is given in TOL.
    \endverbatim

 * @param[in] JOBV
          JOBV is CHARACTER*1 \n
          Specifies whether the output from this procedure is used
          to compute the matrix V: \n
          = 'V': the product of the Jacobi rotations is accumulated
                 by postmulyiplying the N-by-N array V.
                (See the description of V.) \n
          = 'A': the product of the Jacobi rotations is accumulated
                 by postmulyiplying the MV-by-N array V.
                (See the descriptions of MV and V.) \n
          = 'N': the Jacobi rotations are not accumulated. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the input matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the input matrix A.
          M >= N >= 0. \n
 * @param[in] N1
          N1 is INTEGER \n
          N1 specifies the 2 x 2 block partition, the first N1 columns are
          rotated 'against' the remaining N-N1 columns of A. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, M-by-N matrix A, such that A*diag(D) represents
          the input matrix. \n
          On exit,
          A_onexit * D_onexit represents the input matrix A*diag(D)
          post-multiplied by a sequence of Jacobi rotations, where the
          rotation threshold and the total number of sweeps are given in
          TOL and NSWEEP, respectively. \n
          (See the descriptions of N1, D, TOL and NSWEEP.) \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          The array D accumulates the scaling factors from the fast scaled
          Jacobi rotations. \n
          On entry, A*diag(D) represents the input matrix. \n
          On exit, A_onexit*diag(D_onexit) represents the input matrix
          post-multiplied by a sequence of Jacobi rotations, where the
          rotation threshold and the total number of sweeps are given in
          TOL and NSWEEP, respectively. \n
          (See the descriptions of N1, A, TOL and NSWEEP.) \n
 * @param[in,out] SVA
          SVA is REAL array, dimension (N) \n
          On entry, SVA contains the Euclidean norms of the columns of
          the matrix A*diag(D). \n
          On exit, SVA contains the Euclidean norms of the columns of
          the matrix onexit*diag(D_onexit). \n
 * @param[in] MV
          MV is INTEGER \n
          If JOBV = 'A', then MV rows of V are post-multipled by a
                           sequence of Jacobi rotations. \n
          If JOBV = 'N',  then MV is not referenced. \n
 * @param[in,out] V
          V is REAL array, dimension (LDV,N) \n
          If JOBV = 'V' then N rows of V are post-multipled by a
                           sequence of Jacobi rotations. \n
          If JOBV = 'A' then MV rows of V are post-multipled by a
                           sequence of Jacobi rotations. \n
          If JOBV = 'N',  then V is not referenced. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V, LDV >= 1. \n
          If JOBV = 'V', LDV >= N. \n
          If JOBV = 'A', LDV >= MV. \n
 * @param[in] EPS
          EPS is REAL \n
          EPS = SLAMCH('Epsilon') \n
 * @param[in] SFMIN
          SFMIN is REAL \n
          SFMIN = SLAMCH('Safe Minimum') \n
 * @param[in] TOL
          TOL is REAL \n
          TOL is the threshold for Jacobi rotations. For a pair
          A(:,p), A(:,q) of pivot columns, the Jacobi rotation is
          applied only if ABS(COS(angle(A(:,p),A(:,q)))) > TOL. \n
 * @param[in] NSWEEP
          NSWEEP is INTEGER \n
          NSWEEP is the number of sweeps of Jacobi rotations to be
          performed. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK) \n
 * @param[in] LWORK
          LWORK is INTEGER \n
          LWORK is the dimension of WORK. LWORK >= M. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, then the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void gsvj1(char *jobv, integer *m, integer *n, integer *n1, T *a, integer *lda, T *d, T *sva,
               integer *mv, T *v, integer *ldv, T *eps, T *sfmin, T *tol, integer *nsweep, T *work,
               integer *lwork, integer *info)
    {
        gsvj1(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork,
              info);
    }
    template <typename T, typename Ta>
    void gsvj1(char *jobv, integer *m, integer *n, integer *n1, T *a, integer *lda, T *d, Ta *sva,
               integer *mv, T *v, integer *ldv, Ta *eps, Ta *sfmin, Ta *tol, integer *nsweep,
               T *work, integer *lwork, integer *info)
    {
        gsvj1(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work, lwork,
              info);
    }
    /** @}*/ // end of gesvj1

    /** @defgroup las2 las2
     * @ingroup SVD
     * @{
     */
    /*! @brief LAS2 computes singular values of a 2-by-2 triangular matrix

 * @details
 * \b Purpose:
    \verbatim
    LAS2  computes the singular values of the 2-by-2 matrix
       [  F   G  ]
       [  0   H  ].
    On   return, SSMIN is the smaller singular value and SSMAX is the
    larger singular value.
    \endverbatim

 * @param[in] F
          F is REAL \n
          The (1,1) element of the 2-by-2 matrix. \n
 * @param[in] G
          G is REAL \n
          The (1,2) element of the 2-by-2 matrix. \n
 * @param[in] H
          H is REAL \n
          The (2,2) element of the 2-by-2 matrix. \n
 * @param[out] SSMIN
          SSMIN is REAL \n
          The smaller singular value. \n
 * @param[out] SSMAX
          SSMAX is REAL \n
          The larger singular value.  \n

 *  * */
    template <typename T>
    void las2(T *f, T *g, T *h, T *ssmin, T *ssmax)
    {
        las2(f, g, h, ssmin, ssmax);
    }
    /** @}*/ // end of las2

    /** @defgroup lasv2 lasv2
     * @ingroup SVD
     * @{
     */
    /*! @brief LASV2 computes the singular value decomposition of a 2-by-2 triangular matrix

 * @details
 * \b Purpose:
    \verbatim
    LASV2 computes the singular value decomposition of a 2-by-2
    triangular matrix
       [  F   G  ]
       [  0   H  ].
    On   return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
    smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
    right singular vectors for abs(SSMAX), giving the decomposition

       [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
       [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
    \endverbatim

 * @param[in] F
          F is REAL \n
          The (1,1) element of the 2-by-2 matrix. \n
 * @param[in] G
          G is REAL \n
          The (1,2) element of the 2-by-2 matrix. \n
 * @param[in] H
          H is REAL \n
          The (2,2) element of the 2-by-2 matrix. \n
 * @param[out] SSMIN
          SSMIN is REAL \n
          abs(SSMIN) is the smaller singular value. \n
 * @param[out] SSMAX
          SSMAX is REAL \n
          abs(SSMAX) is the larger singular value. \n
 * @param[out] SNL
          SNL is REAL \n
 * @param[out] CSL
          CSL is REAL \n
          The vector (CSL, SNL) is a unit left singular vector for the
          singular value abs(SSMAX). \n
 * @param[out] SNR
          SNR is REAL \n
 * @param[out] CSR
          CSR is REAL \n
          The vector (CSR, SNR) is a unit right singular vector for the
          singular value abs(SSMAX).  \n

 *  * */
    template <typename T>
    void lasv2(T *f, T *g, T *h, T *ssmin, T *ssmax, T *snr, T *csr, T *snl, T *csl)
    {
        lasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl);
    }
    /** @}*/ // end of lasv2

    /** @defgroup ggsvp3 ggsvp3
     * @ingroup SVD
     * @{
     */
    /*! @brief GGSVP3 computes orthogonal matrices U, V and Q
 * @details
 * \b Purpose:
    \verbatim
     GGSVP3 computes orthogonal matrices U, V and Q such that

                        N-K-L  K    L
      U**T*A*Q =     K ( 0    A12  A13)  if M-K-L >= 0;
                     L ( 0     0   A23)
                 M-K-L ( 0     0    0 )

                      N-K-L  K    L
             =     K ( 0    A12  A13)  if M-K-L < 0;
                 M-K ( 0     0   A23)

                      N-K-L  K    L
      V**T*B*Q =   L ( 0     0   B13)
                 P-L ( 0     0    0 )

     where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective
     numerical rank of the (M+P)-by-N matrix (A**T,B**T)**T.

     This decomposition is the preprocessing step for computing the
     Generalized Singular Value Decomposition (GSVD), see subroutine
     SGGSVD.
    \endverbatim

 * @param[in] JOBU
          JOBU is CHARACTER*1 \n
          = 'U':  Orthogonal matrix U is computed; \n
          = 'N':  U is not computed. \n
 * @param[in] JOBV
          JOBV is CHARACTER*1 \n
          = 'V':  Orthogonal matrix V is computed; \n
          = 'N':  V is not computed. \n
 * @param[in] JOBQ
          JOBQ is CHARACTER*1 \n
          = 'Q':  Orthogonal matrix Q is computed; \n
          = 'N':  Q is not computed. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows of the matrix B.  P >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrices A and B.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A. \n
          On exit, A contains the triangular (or trapezoidal) matrix
          described in the Purpose section. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the P-by-N matrix B. \n
          On exit, B contains the triangular matrix described in
          the Purpose section. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,P). \n
 * @param[in] TOLA
          TOLA is REAL \n
 * @param[in] TOLB
          TOLB is REAL \n
          TOLA and TOLB are the thresholds to determine the effective
          numerical rank of matrix B and a subblock of A. Generally,
          they are set to \n
             TOLA = MAX(M,N)*norm(A)*MACHEPS, \n
             TOLB = MAX(P,N)*norm(B)*MACHEPS. \n
          The size of TOLA and TOLB may affect the size of backward
          errors of the decomposition. \n
 * @param[out] K
          K is INTEGER \n
 * @param[out] L
          L is INTEGER \n
          On exit, K and L specify the dimension of the subblocks
          described in Purpose section. \n
          K + L = effective numerical rank of (A**T,B**T)**T. \n
 * @param[out] U
          U is REAL array, dimension (LDU,M) \n
          If JOBU = 'U', U contains the orthogonal matrix U. \n
          If JOBU = 'N', U is not referenced. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U. LDU >= fla_max(1,M) if
          JOBU = 'U'; LDU >= 1 otherwise. \n
 * @param[out] V
          V is REAL array, dimension (LDV,P) \n
          If JOBV = 'V', V contains the orthogonal matrix V. \n
          If JOBV = 'N', V is not referenced. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. LDV >= fla_max(1,P) if
          JOBV = 'V'; LDV >= 1 otherwise. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ,N) \n
          If JOBQ = 'Q', Q contains the orthogonal matrix Q. \n
          If JOBQ = 'N', Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. LDQ >= fla_max(1,N) if
          JOBQ = 'Q'; LDQ >= 1 otherwise. \n
 * @param[out]	IWORK
          IWORK is INTEGER array, dimension (N) \n
 * @param[out]	TAU
          TAU is REAL array, dimension (N) \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)) \n
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. \n
 * @param[in]	LWORK
          LWORK is INTEGER \n
          The dimension of the array WORK. \n
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
    void ggsvp3(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, T *a,
                integer *lda, T *b, integer *ldb, T *tola, T *tolb, integer *k, integer *l, T *u,
                integer *ldu, T *v, integer *ldv, T *q, integer *ldq, integer *iwork, T *tau,
                T *work, integer *lwork, integer *info)
    {
        ggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq,
               iwork, tau, work, lwork, info);
    }
    template <typename T, typename Ta>
    void ggsvp3(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, T *a,
                integer *lda, T *b, integer *ldb, Ta *tola, Ta *tolb, integer *k, integer *l, T *u,
                integer *ldu, T *v, integer *ldv, T *q, integer *ldq, integer *iwork, Ta *rwork,
                T *tau, T *work, integer *lwork, integer *info)
    {
        ggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq,
               iwork, rwork, tau, work, lwork, info);
    }
    /** @}*/ // end of ggsvp3

    /** @defgroup tgsja tgsja
     * @ingroup SVD
     * @{
     */
    /*! @brief TGSJA computes the generalized singular value decomposition (GSVD)  \n
     of two real upper triangular (or trapezoidal) matrices A and B.

 * @details
 * \b Purpose:
    \verbatim
     TGSJA computes the generalized singular value decomposition (GSVD)
     of two real upper triangular (or trapezoidal) matrices A and B.

     On entry, it is assumed that matrices A and B have the following
     forms, which may be obtained by the preprocessing subroutine SGGSVP
     from a general M-by-N matrix A and P-by-N matrix B:

                  N-K-L  K    L
        A =    K ( 0    A12  A13) if M-K-L >= 0;
               L ( 0     0   A23)
           M-K-L ( 0     0    0 )

                N-K-L  K    L
        A =  K ( 0    A12  A13) if M-K-L < 0;
           M-K ( 0     0   A23)

                N-K-L  K    L
        B =  L ( 0     0   B13)
           P-L ( 0     0    0 )

     where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     otherwise A23 is (M-K)-by-L upper trapezoidal.

     On exit,

            U**T *A*Q = D1*( 0 R),    V**T *B*Q = D2*( 0 R),

     where U, V and Q are orthogonal matrices.
     R is a nonsingular upper triangular matrix, and D1 and D2 are
     ``diagonal'' matrices, which are of the following structures:

     If M-K-L >= 0,

                         K  L
            D1 =     K ( I  0)
                     L ( 0  C)
                 M-K-L ( 0  0)

                       K  L
            D2 = L   ( 0  S)
                 P-L ( 0  0)

                    N-K-L  K    L
       ( 0 R) = K (  0   R11  R12) K
                 L (  0    0   R22) L

     where

       C = diag( ALPHA(K+1), ... , ALPHA(K+L)),
       S = diag( BETA(K+1),  ... , BETA(K+L)),
       C**2 + S**2 = I.

       R is stored in A(1:K+L,N-K-L+1:N) on exit.

     If M-K-L < 0,

                    K M-K K+L-M
         D1 =   K ( I  0    0  )
              M-K ( 0  C    0  )

                      K M-K K+L-M
         D2 =   M-K ( 0  S    0  )
              K+L-M ( 0  0    I  )
                P-L ( 0  0    0  )

                    N-K-L  K   M-K  K+L-M
     ( 0 R) =    K ( 0    R11  R12  R13 )
               M-K ( 0     0   R22  R23 )
             K+L-M ( 0     0    0   R33 )

     where
     C = diag( ALPHA(K+1), ... , ALPHA(M)),
     S = diag( BETA(K+1),  ... , BETA(M)),
     C**2 + S**2 = I.

     R = ( R11 R12 R13) is stored in A(1:M, N-K-L+1:N) and R33 is stored
         (  0  R22 R23)
     in B(M-K+1:L,N+M-K-L+1:N) on exit.

     The computation of the orthogonal transformation matrices U, V or Q
     is optional.  These matrices may either be formed explicitly, or they
     may be postmultiplied into input matrices U1, V1, or Q1.
    \endverbatim

 * @param[in] JOBU
          JOBU is CHARACTER*1 \n
          = 'U':  U must contain an orthogonal matrix U1 on entry, and
                  the product U1*U is returned; \n
          = 'I':  U is initialized to the unit matrix, and the
                  orthogonal matrix U is returned; \n
          = 'N':  U is not computed. \n
 * @param[in] JOBV
          JOBV is CHARACTER*1 \n
          = 'V':  V must contain an orthogonal matrix V1 on entry, and
                  the product V1*V is returned; \n
          = 'I':  V is initialized to the unit matrix, and the
                  orthogonal matrix V is returned; \n
          = 'N':  V is not computed. \n
 * @param[in] JOBQ
          JOBQ is CHARACTER*1 \n
          = 'Q':  Q must contain an orthogonal matrix Q1 on entry, and
                  the product Q1*Q is returned; \n
          = 'I':  Q is initialized to the unit matrix, and the
                  orthogonal matrix Q is returned; \n
          = 'N':  Q is not computed. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] P
          P is INTEGER \n
          The number of rows of the matrix B.  P >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrices A and B.  N >= 0. \n
 * @param[in] K
          K is INTEGER \n
 * @param[in] L
          L is INTEGER \n
          K and L specify the subblocks in the input matrices A and B:
          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N)
          of A and B, whose GSVD is going to be computed by STGSJA.
          See Further Details. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          On entry, the M-by-N matrix A.
          On exit, A(N-K+1:N,1:MIN(K+L,M)) contains the triangular
          matrix R or part of R.  See Purpose for details. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >= fla_max(1,M). \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          On entry, the P-by-N matrix B.
          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains
          a part of R.  See Purpose for details. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >= fla_max(1,P). \n
 * @param[in] TOLA
          TOLA is REAL \n
 * @param[in] TOLB
          TOLB is REAL \n
          TOLA and TOLB are the convergence criteria for the Jacobi-
          Kogbetliantz iteration procedure. Generally, they are the
          same as used in the preprocessing step, say \n
              TOLA = fla_max(M,N)*norm(A)*MACHEPS, \n
              TOLB = fla_max(P,N)*norm(B)*MACHEPS. \n
 * @param[out] ALPHA
          ALPHA is REAL array, dimension (N) \n
 * @param[out] BETA
          BETA is REAL array, dimension (N) \n
          On exit, ALPHA and BETA contain the generalized singular
          value pairs of A and B; \n
            ALPHA(1:K) = 1, \n
            BETA(1:K)  = 0, \n
          and if M-K-L >= 0, \n
            ALPHA(K+1:K+L) = diag(C), \n
            BETA(K+1:K+L)  = diag(S), \n
          or if M-K-L < 0, \n
            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0 \n
            BETA(K+1:M) = S, BETA(M+1:K+L) = 1. \n
          Furthermore, if K+L < N, \n
            ALPHA(K+L+1:N) = 0 and \n
            BETA(K+L+1:N)  = 0. \n
 * @param[in,out] U
          U is REAL array, dimension (LDU,M) \n
          On entry, if JOBU = 'U', U must contain a matrix U1 (usually
          the orthogonal matrix returned by SGGSVP). \n
          On exit, \n
          if JOBU = 'I', U contains the orthogonal matrix U; \n
          if JOBU = 'U', U contains the product U1*U. \n
          If JOBU = 'N', U is not referenced. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U. LDU >= fla_max(1,M) if
          JOBU = 'U'; LDU >= 1 otherwise. \n
 * @param[in,out] V
          V is REAL array, dimension (LDV,P) \n
          On entry, if JOBV = 'V', V must contain a matrix V1 (usually
          the orthogonal matrix returned by SGGSVP). \n
          On exit, \n
          if JOBV = 'I', V contains the orthogonal matrix V; \n
          if JOBV = 'V', V contains the product V1*V. \n
          If JOBV = 'N', V is not referenced. \n
 * @param[in] LDV
          LDV is INTEGER \n
          The leading dimension of the array V. LDV >= fla_max(1,P) if
          JOBV = 'V'; LDV >= 1 otherwise. \n
 * @param[in,out] Q
          Q is REAL array, dimension (LDQ,N) \n
          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually
          the orthogonal matrix returned by SGGSVP). \n
          On exit, \n
          if JOBQ = 'I', Q contains the orthogonal matrix Q; \n
          if JOBQ = 'Q', Q contains the product Q1*Q. \n
          If JOBQ = 'N', Q is not referenced. \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q. LDQ >= fla_max(1,N) if
          JOBQ = 'Q'; LDQ >= 1 otherwise. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (2*N) \n
 * @param[out]	NCYCLE
          NCYCLE is INTEGER \n
          The number of cycles required for convergence. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          = 1:  the procedure does not converge after MAXIT cycles. \n

 *  * */
    template <typename T>
    void tgsja(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k,
               integer *l, T *a, integer *lda, T *b, integer *ldb, T *tola, T *tolb, T *alpha,
               T *beta, T *u, integer *ldu, T *v, integer *ldv, T *q, integer *ldq, T *work,
               integer *ncycle, integer *info)
    {
        tgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v,
              ldv, q, ldq, work, ncycle, info);
    }
    template <typename T, typename Ta>
    void tgsja(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, integer *k,
               integer *l, T *a, integer *lda, T *b, integer *ldb, Ta *tola, Ta *tolb, Ta *alpha,
               Ta *beta, T *u, integer *ldu, T *v, integer *ldv, T *q, integer *ldq, T *work,
               integer *ncycle, integer *info)
    {
        tgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v,
              ldv, q, ldq, work, ncycle, info);
    }
    /** @}*/ // end of tgsja

    /** @defgroup lags2 lags2
     * @ingroup SVD
     * @{
     */
    /*! @brief LAGS2 computes 2-by-2 orthogonal matrices U, V, and Q, and applies them \n
     to matrices A and B such that the rows of the transformed A and B are parallel
 * @details
 * \b Purpose:
    \verbatim
    LAGS2 computes 2-by-2 orthogonal matrices U, V and Q, such
    that if (UPPER) then

              U**T *A*Q = U**T *(A1 A2)*Q = (x  0 )
                                (0  A3)     (x  x )
    and
              V**T*B*Q = V**T *(B1 B2)*Q = (x  0 )
                               (0  B3)     (x  x )

    or if (.NOT.UPPER) then

              U**T *A*Q = U**T *(A1 0 )*Q = (x  x )
                                (A2 A3)     (0  x )
    and
              V**T*B*Q = V**T*(B1 0 )*Q = (x  x )
                              (B2 B3)     (0  x )

    The rows of the transformed A and B are parallel, where

      U = ( CSU  SNU), V = ( CSV SNV), Q = ( CSQ   SNQ)
          (-SNU  CSU)      (-SNV CSV)      (-SNQ   CSQ)

    Z**T denotes the transpose of Z.
    \endverbatim


 * @param[in] UPPER
          UPPER is LOGICAL \n
          = .TRUE.: the input matrices A and B are upper triangular. \n
          = .FALSE.: the input matrices A and B are lower triangular. \n
 * @param[in] A1
          A1 is REAL \n
 * @param[in] A2
          A2 is REAL \n
 * @param[in] A3
          A3 is REAL \n
          On entry, A1, A2 and A3 are elements of the input 2-by-2
          upper (lower) triangular matrix A. \n
 * @param[in] B1
          B1 is REAL \n
 * @param[in] B2
          B2 is REAL \n
 * @param[in] B3
          B3 is REAL \n
          On entry, B1, B2 and B3 are elements of the input 2-by-2
          upper (lower) triangular matrix B. \n
 * @param[out] CSU
          CSU is REAL \n
 * @param[out] SNU
          SNU is REAL \n
          The desired orthogonal matrix U. \n
 * @param[out] CSV
          CSV is REAL \n
 * @param[out] SNV
          SNV is REAL \n
          The desired orthogonal matrix V. \n
 * @param[out] CSQ
          CSQ is REAL \n
 * @param[out] SNQ
          SNQ is REAL \n
          The desired orthogonal matrix Q.  \n

 * */
    template <typename T>
    void lags2(logical *upper, T *a1, T *a2, T *a3, T *b1, T *b2, T *b3, T *csu, T *snu, T *csv,
               T *snv, T *csq, T *snq)
    {
        lags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
    }
    template <typename T, typename Ta>
    void lags2(logical *upper, Ta *a1, T *a2, Ta *a3, Ta *b1, T *b2, Ta *b3, Ta *csu, T *snu,
               Ta *csv, T *snv, Ta *csq, T *snq)
    {
        lags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
    }
    /** @}*/ // end of lags2

    /** @defgroup lapll lapll
     * @ingroup SVD
     * @{
     */
    /*! @brief LAPLL measures the linear dependence of two vectors

 * @details
 * \b Purpose:
    \verbatim
    Given two column vectors X and Y, let
                         A = (X Y).
    The subroutine first computes the QR factorization of A = Q*R,
    and then computes the SVD of the 2-by-2 upper triangular matrix R.
    The smaller singular value of R is   returned in SSMIN, which is used
    as the measurement of the linear dependency of the vectors X and Y.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The length of the vectors X and Y. \n
 * @param[in,out] X
          X is REAL array,
                         dimension (1+(N-1)*INCX) \n
          On entry, X contains the N-vector X.
          On exit, X is overwritten. \n
 * @param[in] INCX
          INCX is INTEGER \n
          The increment between successive elements of X. INCX > 0. \n
 * @param[in,out] Y
          Y is REAL array,
                         dimension (1+(N-1)*INCY) \n
          On entry, Y contains the N-vector Y.
          On exit, Y is overwritten. \n
 * @param[in] INCY
          INCY is INTEGER \n
          The increment between successive elements of Y. INCY > 0. \n
 * @param[out] SSMIN
          SSMIN is REAL \n
          The smallest singular value of the N-by-2 matrix A = (X Y).  \n

 * */
    template <typename T>
    void lapll(integer *n, T *x, integer *incx, T *y, integer *incy, T *ssmin)
    {
        lapll(n, x, incx, y, incy, ssmin);
    }
    template <typename T, typename Ta>
    void lapll(integer *n, T *x, integer *incx, T *y, integer *incy, Ta *ssmin)
    {
        lapll(n, x, incx, y, incy, ssmin);
    }
    /** @}*/ // end of lapll

    /** @defgroup lasq1 lasq1
     * @ingroup SVD
     * @{
     */
    /*! @brief LASQ1 computes the singular values of a real square bidiagonal matrix. \n
     Used by sbdsqr
 * @details
 * \b Purpose:
    \verbatim
    LASQ1 computes the singular values of a real N-by-N bidiagonal
    matrix with diagonal D and off-diagonal E. The singular values
    are computed to high relative accuracy, in the absence of
    denormalization, underflow and overflow. The algorithm was first
    presented in

    "Accurate singular values and differential qd algorithms" by K. V.
    Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
    1994,

    and the present implementation is described in "An implementation of
    the dqds Algorithm (Positive Case)", LAPACK Working Note.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of rows and columns in the matrix. N >= 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, D contains the diagonal elements of the
          bidiagonal matrix whose SVD is desired. On normal exit,
          D contains the singular values in decreasing order. \n
 * @param[in,out] E
          E is REAL array, dimension (N) \n
          On entry, elements E(1:N-1) contain the off-diagonal elements
          of the bidiagonal matrix whose SVD is desired. \n
          On exit, E is overwritten. \n
 * @param[out] WORK
          WORK is REAL array, dimension (4*N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if INFO = -i, the i-th argument had an illegal value \n
          > 0: the algorithm failed \n
               = 1, a split was marked by a positive value in E \n
               = 2, current block of Z not diagonalized after 100*N
                    iterations (in inner while loop)  On exit D and E
                    represent a matrix with the same singular values
                    which the calling subroutine could use to finish the
                    computation, or even feed back into SLASQ1 \n
               = 3, termination criterion of outer while loop not met
                    (program created more than N unreduced blocks) \n

 *  * */
    template <typename T>
    void lasq1(integer *n, T *d, T *e, T *work, integer *info)
    {
        lasq1(n, d, e, work, info);
    }
    /** @}*/ // end of lasq1

    /** @defgroup lasq2 lasq2
     * @ingroup SVD
     * @{
     */
    /*! @brief LASQ2 computes all the eigenvalues of the symmetric positive definite   \n
     tridiagonal matrix associated with the qd Array Z to high relative accuracy.  \n
     Used by sbdsqr and sstegr.
 * @details
 * \b Purpose:
    \verbatim
    LASQ2 computes all the eigenvalues of the symmetric positive
    definite tridiagonal matrix associated with the qd array Z to high
    relative accuracy are computed to high relative accuracy, in the
    absence of denormalization, underflow and overflow.

    To see the relation of Z to the tridiagonal matrix, let L be a
    unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
    let U be an upper bidiagonal matrix with 1's above and diagonal
    Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
    symmetric tridiagonal to which it is similar.

    Note : SLASQ2 defines a logical variable, IEEE, which is true
    on machines which follow ieee-754 floating-point standard in their
    handling of infinities and NaNs, and false otherwise. This variable
    is passed to SLASQ3.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of rows and columns in the matrix. N >= 0. \n
 * @param[in,out] Z
          Z is REAL array, dimension (4*N) \n
          On entry Z holds the qd array. On exit, entries 1 to N hold
          the eigenvalues in decreasing order, Z(2*N+1) holds the
          trace, and Z(2*N+2) holds the sum of the eigenvalues. If
          N > 2, then Z(2*N+3) holds the iteration count, Z(2*N+4)
          holds NDIVS/NIN^2, and Z(2*N+5) holds the percentage of
          shifts that failed. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0: successful exit \n
          < 0: if the i-th argument is a scalar and had an illegal
               value, then INFO = -i, if the i-th argument is an
               array and the j-entry had an illegal value, then
               INFO = -(i*100+j) \n
          > 0: the algorithm failed \n
                = 1, a split was marked by a positive value in E \n
                = 2, current block of Z not diagonalized after 100*N \n
                     iterations (in inner while loop).  On exit Z holds
                     a qd array with the same eigenvalues as the given Z. \n
                = 3, termination criterion of outer while loop not met
                     (program created more than N unreduced blocks)  \n

 *  * */
    template <typename T>
    void lasq2(integer *n, T *z, integer *info)
    {
        lasq2(n, z, info);
    }
    /** @}*/ // end of lasq2

    /** @defgroup lasq3 lasq3
     * @ingroup SVD
     * @{
     */
    /*! @brief LASQ3 checks for deflation, computes a shift and calls dqds. Used by sbdsqr.

 * @details
 * \b Purpose:
    \verbatim
    LASQ3 checks for deflation, computes a shift (TAU) and calls dqds.
    In case of failure it changes shifts, and tries again until output
    is positive.
    \endverbatim

 * @param[in] I0
          I0 is INTEGER \n
          First index. \n
 * @param[in,out] N0
          N0 is INTEGER \n
          Last index. \n
 * @param[in,out] Z
          Z is REAL array, dimension (4*N0) \n
          Z holds the qd array. \n
 * @param[in,out] PP
          PP is INTEGER \n
          PP=0 for ping, PP=1 for pong. \n
          PP=2 indicates that flipping was applied to the Z array
          and that the initial tests for deflation should not be
          performed. \n
 * @param[out] DMIN
          DMIN is REAL \n
          Minimum value of d. \n
 * @param[out] SIGMA
          SIGMA is REAL \n
          Sum of shifts used in current segment. \n
 * @param[in,out] DESIG
          DESIG is REAL \n
          Lower order part of SIGMA \n
 * @param[in] QMAX
          QMAX is REAL \n
          Maximum value of q. \n
 * @param[in,out] NFAIL
          NFAIL is INTEGER \n
          Increment NFAIL by 1 each time the shift was too big. \n
 * @param[in,out] ITER
          ITER is INTEGER \n
          Increment ITER by 1 for each iteration. \n
 * @param[in,out] NDIV
          NDIV is INTEGER \n
          Increment NDIV by 1 for each division. \n
 * @param[in] IEEE
          IEEE is LOGICAL \n
          Flag for IEEE or non IEEE arithmetic (passed to SLASQ5). \n
 * @param[in,out] TTYPE
          TTYPE is INTEGER \n
          Shift type. \n
 * @param[in,out] DMIN1
          DMIN1 is REAL \n
 * @param[in,out] DMIN2
          DMIN2 is REAL \n
 * @param[in,out] DN
          DN is REAL \n
 * @param[in,out] DN1
          DN1 is REAL \n
 * @param[in,out] DN2
          DN2 is REAL \n
 * @param[in,out] G
          G is REAL \n
 * @param[in,out] TAU
          TAU is REAL \n
          These are passed as arguments in order to save their values
          between calls to LASQ3. \n

 *  * */
    template <typename T>
    void lasq3(integer *i0, integer *n0, T *z, integer *pp, T *dmin, T *sigma, T *desig, T *qmax,
               integer *nfail, integer *iter, integer *ndiv, logical *ieee, integer *ttype,
               T *dmin1, T *dmin2, T *dn, T *dn1, T *dn2, T *g, T *tau)
    {
        lasq3(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2,
              dn, dn1, dn2, g, tau);
    }
    /** @}*/ // end of lasq3

    /** @defgroup lasq4 lasq4
     * @ingroup SVD
     * @{
     */
    /*! @brief LASQ4 computes an approximation to the smallest eigenvalue using   \n
     values of d from the previous transform. Used by sbdsqr
 * @details
 * \b Purpose:
    \verbatim
    LASQ4 computes an approximation TAU to the smallest eigenvalue
    using values of d from the previous transform.
    \endverbatim

 * @param[in] I0
          I0 is INTEGER \n
          First index. \n
 * @param[in] N0
          N0 is INTEGER \n
          Last index. \n
 * @param[in] Z
          Z is REAL array, dimension (4*N0) \n
          Z holds the qd array. \n
 * @param[in] PP
          PP is INTEGER \n
          PP=0 for ping, PP=1 for pong. \n
 * @param[in] N0IN
          N0IN is INTEGER \n
          The value of N0 at start of EIGTEST. \n
 * @param[in] DMIN
          DMIN is REAL \n
          Minimum value of d. \n
 * @param[in] DMIN1
          DMIN1 is REAL \n
          Minimum value of d, excluding D(N0). \n
 * @param[in] DMIN2
          DMIN2 is REAL \n
          Minimum value of d, excluding D(N0) and D(N0-1). \n
 * @param[in] DN
          DN is REAL \n
          d(N) \n
 * @param[in] DN1
          DN1 is REAL \n
          d(N-1) \n
 * @param[in] DN2
          DN2 is REAL \n
          d(N-2) \n
 * @param[out] TAU
          TAU is REAL \n
          This is the shift. \n
 * @param[out] TTYPE
          TTYPE is INTEGER \n
          Shift type. \n
 * @param[in,out] G
          G is REAL \n
          G is passed as an argument in order to save its value between
          calls to SLASQ4.  \n

 *  * */
    template <typename T>
    void lasq4(integer *i0, integer *n0, T *z, integer *pp, integer *n0in, T *dmin, T *dmin1,
               T *dmin2, T *dn, T *dn1, T *dn2, T *tau, integer *ttype, T *g)
    {
        lasq4(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
    }
    /** @}*/ // end of lasq4

    /** @defgroup lasq5 lasq5
     * @ingroup SVD
     * @{
     */
    /*! @brief LASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr

 * @details
 * \b Purpose:
    \verbatim
    LASQ5 computes one dqds transform in ping-pong form, one
    version for IEEE machines another for non IEEE machines.
    \endverbatim

 * @param[in] I0
          I0 is INTEGER \n
          First index. \n
 * @param[in] N0
          N0 is INTEGER \n
          Last index. \n
 * @param[in] Z
          Z is REAL array, dimension (4*N) \n
          Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
          an extra argument. \n
 * @param[in] PP
          PP is INTEGER \n
          PP=0 for ping, PP=1 for pong. \n
 * @param[in] TAU
          TAU is REAL \n
          This is the shift. \n
 * @param[in] SIGMA
          SIGMA is REAL \n
          This is the accumulated shift up to this step. \n
 * @param[out] DMIN
          DMIN is REAL \n
          Minimum value of d. \n
 * @param[out] DMIN1
          DMIN1 is REAL \n
          Minimum value of d, excluding D(N0). \n
 * @param[out] DMIN2
          DMIN2 is REAL \n
          Minimum value of d, excluding D(N0) and D(N0-1). \n
 * @param[out] DN
          DN is REAL \n
          d(N0), the last value of d. \n
 * @param[out] DNM1
          DNM1 is REAL \n
          d(N0-1). \n
 * @param[out] DNM2
          DNM2 is REAL \n
          d(N0-2). \n
 * @param[in] IEEE
          IEEE is LOGICAL \n
          Flag for IEEE or non IEEE arithmetic. \n
 * @param[in] EPS
          EPS is REAL \n
          This is the value of epsilon used. \n

 *  * */
    template <typename T>
    void lasq5(integer *i0, integer *n0, T *z, integer *pp, T *tau, T *sigma, T *dmin, T *dmin1,
               T *dmin2, T *dn, T *dnm1, T *dnm2, logical *ieee, T *eps)
    {
        lasq5(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps);
    }

    /** @}*/ // end of lasq5

    /** @defgroup lasq6 lasq6
     * @ingroup SVD
     * @{
     */
    /*! @brief LASQ6 computes one dqd transform in ping-pong form. Used by sbdsqr and sstegr.

 * @details
 * \b Purpose:
    \verbatim
    LASQ6 computes one dqd (shift equal to zero) transform in
    ping-pong form, with protection against underflow and overflow.
    \endverbatim

 * @param[in] I0
          I0 is INTEGER \n
          First index. \n
 * @param[in] N0
          N0 is INTEGER \n
          Last index. \n
 * @param[in] Z
          Z is REAL array, dimension (4*N) \n
          Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
          an extra argument. \n
 * @param[in] PP
          PP is INTEGER \n
          PP=0 for ping, PP=1 for pong.  \n
 * @param[out] DMIN
          DMIN is REAL \n
          Minimum value of d. \n
 * @param[out] DMIN1
          DMIN1 is REAL \n
          Minimum value of d, excluding D(N0). \n
 * @param[out] DMIN2
          DMIN2 is REAL \n
          Minimum value of d, excluding D(N0) and D(N0-1). \n
 * @param[out] DN
          DN is REAL \n
          d(N0), the last value of d. \n
 * @param[out] DNM1
          DNM1 is REAL \n
          d(N0-1). \n
 * @param[out] DNM2
          DNM2 is REAL \n
          d(N0-2).  \n

 *  * */
    template <typename T>
    void lasq6(integer *i0, integer *n0, T *z, integer *pp, T *dmin, T *dmin1, T *dmin2, T *dn,
               T *dnm1, T *dnm2)
    {
        lasq6(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2);
    }
    /** @}*/ // end of lasq6

    /** @defgroup lasd0 lasd0
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD0 computes the singular values of a real upper bidiagonal n-by-m  \n
     matrix B with diagonal d and off-diagonal e. Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    Using a divide and conquer approach, SLASD0 computes the singular
    value decomposition (SVD) of a real upper bidiagonal N-by-M
    matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
    The algorithm computes orthogonal matrices U and VT such that
    B = U * S * VT. The singular values S are overwritten on D.

    A related subroutine, SLASDA, computes only the singular values,
    and optionally, the singular vectors in compact form.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          On entry, the row dimension of the upper bidiagonal matrix.
          This is also the dimension of the main diagonal array D. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          Specifies the column dimension of the bidiagonal matrix. \n
          = 0: The bidiagonal matrix has column dimension M = N; \n
          = 1: The bidiagonal matrix has column dimension M = N+1; \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry D contains the main diagonal of the bidiagonal
          matrix. \n
          On exit D, if INFO = 0, contains its singular values. \n
 * @param[in,out] E
          E is REAL array, dimension (M-1) \n
          Contains the subdiagonal entries of the bidiagonal matrix.
          On exit, E has been destroyed. \n
 * @param[out] U
          U is REAL array, dimension (LDU, N) \n
          On exit, U contains the left singular vectors. \n
 * @param[in] LDU
          LDU is INTEGER \n
          On entry, leading dimension of U. \n
 * @param[out] VT
          VT is REAL array, dimension (LDVT, M) \n
          On exit, VT**T contains the right singular vectors. \n
 * @param[in] LDVT
          LDVT is INTEGER \n
          On entry, leading dimension of VT. \n
 * @param[in] SMLSIZ
          SMLSIZ is INTEGER \n
          On entry, maximum size of the subproblems at the
          bottom of the computation tree. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (8*N) \n
 * @param[out] WORK
          WORK is REAL array, dimension (3*M**2+2*M) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, a singular value did not converge  \n

 *  * */
    template <typename T>
    void lasd0(integer *n, integer *sqre, T *d, T *e, T *u, integer *ldu, T *vt, integer *ldvt,
               integer *smlsiz, integer *iwork, T *work, integer *info)
    {
        lasd0(n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info);
    }
    /** @}*/ // end of lasd0

    /** @defgroup lasdt lasdt
     * @ingroup SVD
     * @{
     */
    /*! @brief LASDT creates a tree of subproblems for bidiagonal divide and conquer. \n
     Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    LASDT creates a tree of subproblems for bidiagonal divide and
    conquer.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          On entry, the number of diagonal elements of the
          bidiagonal matrix. \n
 * @param[out] LVL
          LVL is INTEGER \n
          On exit, the number of levels on the computation tree. \n
 * @param[out] ND
          ND is INTEGER \n
          On exit, the number of nodes on the tree. \n
 * @param[out] INODE
          INODE is INTEGER array, dimension (N) \n
          On exit, centers of subproblems. \n
 * @param[out] NDIML
          NDIML is INTEGER array, dimension (N) \n
          On exit, row dimensions of left children. \n
 * @param[out] NDIMR
          NDIMR is INTEGER array, dimension (N) \n
          On exit, row dimensions of right children. \n
 * @param[in] MSUB
          MSUB is INTEGER \n
          On entry, the maximum row dimension each subproblem at the
          bottom of the tree can be of.  \n

 *  * */
    inline void slasdt(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml,
                       integer *ndimr, integer *msub)
    {
        slasdt_(n, lvl, nd, inode, ndiml, ndimr, msub);
    }
    inline void dlasdt(integer *n, integer *lvl, integer *nd, integer *inode, integer *ndiml,
                       integer *ndimr, integer *msub)
    {
        dlasdt_(n, lvl, nd, inode, ndiml, ndimr, msub);
    }
    /** @}*/ // end of lasdt

    /** @defgroup lasd1 lasd1
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD1 computes the SVD of an upper bidiagonal matrix B of the specified size. Used by
 sbdsdc

 * @details
 * \b Purpose:
    \verbatim
    LASD1 computes the SVD of an upper bidiagonal N-by-M matrix B,
    where N = NL + NR + 1 and M = N + SQRE. SLASD1 is called from SLASD0.

    A related subroutine SLASD7 handles the case in which the singular
    values (and the singular vectors in factored form) are desired.

    SLASD1 computes the SVD as follows:

                  (D1(in)    0    0       0)
      B = U(in) * (  Z1**T   a   Z2**T    b) * VT(in)
                  (  0       0   D2(in)   0)

        = U(out) * (D(out) 0) * VT(out)

    where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
    with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
    elsewhere; and the entry b is empty if SQRE = 0.

    The left singular vectors of the original matrix are stored in U, and
    the transpose of the right singular vectors are stored in VT, and the
    singular values are in D.  The algorithm consists of three stages:

       The first stage consists of deflating the size of the problem
       when there are multiple singular values or when there are zeros in
       the Z vector.  For each such occurrence the dimension of the
       secular equation problem is reduced by one.  This stage is
       performed by the routine SLASD2.

       The second stage consists of calculating the updated
       singular values. This is done by finding the square roots of the
       roots of the secular equation via the routine SLASD4 (as called
       by SLASD3). This routine also calculates the singular vectors of
       the current problem.

       The final stage consists of computing the updated singular vectors
       directly using the updated singular values.  The singular vectors
       for the current problem are multiplied with the singular vectors
       from the overall problem.
    \endverbatim

 * @param[in] NL
          NL is INTEGER \n
          The row dimension of the upper block.  NL >= 1. \n
 * @param[in] NR
          NR is INTEGER \n
          The row dimension of the lower block.  NR >= 1. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          = 0: the lower block is an NR-by-NR square matrix. \n
          = 1: the lower block is an NR-by-(NR+1) rectangular matrix. \n
 \n
          The bidiagonal matrix has row dimension N = NL + NR + 1,
          and column dimension M = N + SQRE. \n
 * @param[in,out] D
          D is REAL array, dimension (NL+NR+1). \n
          N = NL+NR+1 \n
          On entry D(1:NL,1:NL) contains the singular values of the
          upper block; and D(NL+2:N) contains the singular values of
          the lower block. On exit D(1:N) contains the singular values
          of the modified matrix. \n
 * @param[in,out] ALPHA
          ALPHA is REAL \n
          Contains the diagonal element associated with the added row. \n
 * @param[in,out] BETA
          BETA is REAL \n
          Contains the off-diagonal element associated with the added
          row. \n
 * @param[in,out] U
          U is REAL array, dimension (LDU,N) \n
          On entry U(1:NL, 1:NL) contains the left singular vectors of
          the upper block; U(NL+2:N, NL+2:N) contains the left singular
          vectors of the lower block. On exit U contains the left
          singular vectors of the bidiagonal matrix. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U.  LDU >= fla_max(1, N). \n
 * @param[in,out] VT
          VT is REAL array, dimension (LDVT,M) \n
          where M = N + SQRE. \n
          On entry VT(1:NL+1, 1:NL+1)**T contains the right singular
          vectors of the upper block; VT(NL+2:M, NL+2:M)**T contains
          the right singular vectors of the lower block. On exit
          VT**T contains the right singular vectors of the
          bidiagonal matrix. \n
 * @param[in] LDVT
          LDVT is INTEGER \n
          The leading dimension of the array VT.  LDVT >= fla_max(1, M). \n
 * @param[in,out] IDXQ
          IDXQ is INTEGER array, dimension (N) \n
          This contains the permutation which will reintegrate the
          subproblem just solved back into sorted order, i.e.
          D(IDXQ(I = 1, N)) will be in ascending order. \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (4*N) \n
 * @param[out] WORK
          WORK is REAL array, dimension (3*M**2+2*M) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, a singular value did not converge \n

 *  * */
    template <typename T>
    void lasd1(integer *nl, integer *nr, integer *sqre, T *d, T *alpha, T *beta, T *u, integer *ldu,
               T *vt, integer *ldvt, integer *idxq, integer *iwork, T *work, integer *info)
    {
        lasd1(nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt, idxq, iwork, work, info);
    }
    /** @}*/ // end of lasd1

    /** @defgroup lasd2 lasd2
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD2 merges the two sets of singular values together into a single sorted set. Used
 by sbdsdc

 * @details
 * \b Purpose:
    \verbatim
    LASD2 merges the two sets of singular values together into a single
    sorted set.  Then it tries to deflate the size of the problem.
    There are two ways in which deflation can occur:  when two or more
    singular values are close together or if there is a tiny entry in the
    Z vector.  For each such occurrence the order of the related secular
    equation problem is reduced by one.

    SLASD2 is called from SLASD1.
    \endverbatim

 * @param[in] NL
          NL is INTEGER \n
          The row dimension of the upper block.  NL >= 1. \n
 * @param[in] NR
          NR is INTEGER \n
          The row dimension of the lower block.  NR >= 1. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          = 0: the lower block is an NR-by-NR square matrix. \n
          = 1: the lower block is an NR-by-(NR+1) rectangular matrix. \n
 \n
          The bidiagonal matrix has N = NL + NR + 1 rows and
          M = N + SQRE >= N columns. \n
 * @param[out] K
          K is INTEGER \n
          Contains the dimension of the non-deflated matrix,
          This is the order of the related secular equation. 1 <= K <=N. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry D contains the singular values of the two submatrices
          to be combined.  On exit D contains the trailing (N-K) updated
          singular values (those which were deflated) sorted into
          increasing order. \n
 * @param[out] Z
          Z is REAL array, dimension (N) \n
          On exit Z contains the updating row vector in the secular
          equation. \n
 * @param[in] ALPHA
          ALPHA is REAL \n
          Contains the diagonal element associated with the added row. \n
 * @param[in] BETA
          BETA is REAL \n
          Contains the off-diagonal element associated with the added
          row. \n
 * @param[in,out] U
          U is REAL array, dimension (LDU,N) \n
          On entry U contains the left singular vectors of two
          submatrices in the two square blocks with corners at (1,1),
          (NL, NL), and (NL+2, NL+2), (N,N). \n
          On exit U contains the trailing (N-K) updated left singular
          vectors (those which were deflated) in its last N-K columns. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U.  LDU >= N. \n
 * @param[in,out] VT
          VT is REAL array, dimension (LDVT,M) \n
          On entry VT**T contains the right singular vectors of two
          submatrices in the two square blocks with corners at (1,1),
          (NL+1, NL+1), and (NL+2, NL+2), (M,M). \n
          On exit VT**T contains the trailing (N-K) updated right singular
          vectors (those which were deflated) in its last N-K columns.
          In case SQRE =1, the last row of VT spans the right null
          space. \n
 * @param[in] LDVT
          LDVT is INTEGER \n
          The leading dimension of the array VT.  LDVT >= M. \n
 * @param[out] DSIGMA
          DSIGMA is REAL array, dimension (N) \n
          Contains a copy of the diagonal elements (K-1 singular values
          and one zero) in the secular equation. \n
 * @param[out] U2
          U2 is REAL array, dimension (LDU2,N) \n
          Contains a copy of the first K-1 left singular vectors which
          will be used by SLASD3 in a matrix multiply (SGEMM) to solve
          for the new left singular vectors. U2 is arranged into four
          blocks. The first block contains a column with 1 at NL+1 and
          zero everywhere else; the second block contains non-zero
          entries only at and above NL; the third contains non-zero
          entries only below NL+1; and the fourth is dense. \n
 * @param[in] LDU2
          LDU2 is INTEGER \n
          The leading dimension of the array U2.  LDU2 >= N. \n
 * @param[out] VT2
          VT2 is REAL array, dimension (LDVT2,N) \n
          VT2**T contains a copy of the first K right singular vectors
          which will be used by SLASD3 in a matrix multiply (SGEMM) to
          solve for the new right singular vectors. VT2 is arranged into
          three blocks. The first block contains a row that corresponds
          to the special 0 diagonal element in SIGMA; the second block
          contains non-zeros only at and before NL +1; the third block
          contains non-zeros only at and after  NL +2. \n
 * @param[in] LDVT2
          LDVT2 is INTEGER \n
          The leading dimension of the array VT2.  LDVT2 >= M. \n
 * @param[out] IDXP
          IDXP is INTEGER array, dimension (N) \n
          This will contain the permutation used to place deflated
          values of D at the end of the array. On output IDXP(2:K)
          points to the nondeflated D-values and IDXP(K+1:N)
          points to the deflated singular values. \n
 * @param[out] IDX
          IDX is INTEGER array, dimension (N) \n
          This will contain the permutation used to sort the contents of
          D into ascending order. \n
 * @param[out] IDXC
          IDXC is INTEGER array, dimension (N) \n
          This will contain the permutation used to arrange the columns
          of the deflated U matrix into three groups:  the first group
          contains non-zero entries only at and above NL, the second
          contains non-zero entries only below NL+2, and the third is
          dense. \n
 * @param[in,out] IDXQ
          IDXQ is INTEGER array, dimension (N) \n
          This contains the permutation which separately sorts the two
          sub-problems in D into ascending order.  Note that entries in
          the first hlaf of this permutation must first be moved one
          position backward; and entries in the second half
          must first have NL+1 added to their values. \n
 * @param[out] COLTYP
          COLTYP is INTEGER array, dimension (N) \n
          As workspace, this will contain a label which will indicate
          which of the following types a column in the U2 matrix or a
          row in the VT2 matrix is: \n
          1 : non-zero in the upper half only \n
          2 : non-zero in the lower half only \n
          3 : dense \n
          4 : deflated \n
 \n
          On exit, it is an array of dimension 4, with COLTYP(I) being
          the dimension of the I-th type columns. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value.  \n

 *  * */
    template <typename T>
    void lasd2(integer *nl, integer *nr, integer *sqre, integer *k, T *d, T *z, T *alpha, T *beta,
               T *u, integer *ldu, T *vt, integer *ldvt, T *dsigma, T *u2, integer *ldu2, T *vt2,
               integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *idxq,
               integer *coltyp, integer *info)
    {
        lasd2(nl, nr, sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2, ldvt2,
              idxp, idx, idxc, idxq, coltyp, info);
    }
    /** @}*/ // end of lasd2

    /** @defgroup lasd3 lasd3
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD3 finds all square roots of the roots of the secular equation, \n
     as defined by the values in D and Z, and then updates the singular       \n
     vectors by matrix multiplication. Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    LASD3 finds all the square roots of the roots of the secular
    equation, as defined by the values in D and Z.  It makes the
    appropriate calls to SLASD4 and then updates the singular
    vectors by matrix multiplication.

    This code makes very mild assumptions about floating point
    arithmetic. It will work on machines with a guard digit in
    add/subtract, or on those binary machines without guard digits
    which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
    It could conceivably fail on hexadecimal or decimal machines
    without guard digits, but we know of none.

    SLASD3 is called from SLASD1.
    \endverbatim

 * @param[in] NL
          NL is INTEGER \n
          The row dimension of the upper block.  NL >= 1. \n
 * @param[in] NR
          NR is INTEGER \n
          The row dimension of the lower block.  NR >= 1. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          = 0: the lower block is an NR-by-NR square matrix. \n
          = 1: the lower block is an NR-by-(NR+1) rectangular matrix. \n
 \n
          The bidiagonal matrix has N = NL + NR + 1 rows and
          M = N + SQRE >= N columns. \n
 * @param[in] K
          K is INTEGER \n
          The size of the secular equation, 1 =< K = < N. \n
 * @param[out] D
          D is REAL array, dimension(K) \n
          On exit the square roots of the roots of the secular equation,
          in ascending order. \n
 * @param[out] Q
          Q is REAL array, dimension (LDQ,K) \n
 * @param[in] LDQ
          LDQ is INTEGER \n
          The leading dimension of the array Q.  LDQ >= K. \n
 * @param[in,out] DSIGMA
          DSIGMA is REAL array, dimension(K) \n
          The first K elements of this array contain the old roots
          of the deflated updating problem.  These are the poles
          of the secular equation. \n
 * @param[out] U
          U is REAL array, dimension (LDU, N) \n
          The last N - K columns of this matrix contain the deflated
          left singular vectors. \n
 * @param[in] LDU
          LDU is INTEGER \n
          The leading dimension of the array U.  LDU >= N. \n
 * @param[in] U2
          U2 is REAL array, dimension (LDU2, N) \n
          The first K columns of this matrix contain the non-deflated
          left singular vectors for the split problem. \n
 * @param[in] LDU2
          LDU2 is INTEGER \n
          The leading dimension of the array U2.  LDU2 >= N. \n
 * @param[out] VT
          VT is REAL array, dimension (LDVT, M) \n
          The last M - K columns of VT**T contain the deflated
          right singular vectors. \n
 * @param[in] LDVT
          LDVT is INTEGER \n
          The leading dimension of the array VT.  LDVT >= N. \n
 * @param[in,out] VT2
          VT2 is REAL array, dimension (LDVT2, N) \n
          The first K columns of VT2**T contain the non-deflated
          right singular vectors for the split problem. \n
 * @param[in] LDVT2
          LDVT2 is INTEGER \n
          The leading dimension of the array VT2.  LDVT2 >= N. \n
 * @param[in] IDXC
          IDXC is INTEGER array, dimension (N) \n
          The permutation used to arrange the columns of U (and rows of
          VT) into three groups:  the first group contains non-zero
          entries only at and above (or before) NL +1; the second
          contains non-zero entries only at and below (or after) NL+2;
          and the third is dense. The first column of U and the row of
          VT are treated separately, however. \n
 \n
          The rows of the singular vectors found by SLASD4
          must be likewise permuted before the matrix multiplies can
          take place. \n
 * @param[in] CTOT
          CTOT is INTEGER array, dimension (4) \n
          A count of the total number of the various types of columns
          in U (or rows in VT), as described in IDXC. The fourth column
          type is any column which has been deflated. \n
 * @param[in,out] Z
          Z is REAL array, dimension (K) \n
          The first K elements of this array contain the components
          of the deflation-adjusted updating row vector. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, a singular value did not converge  \n

 *  * */
    template <typename T>
    void lasd3(integer *nl, integer *nr, integer *sqre, integer *k, T *d, T *q, integer *ldq,
               T *dsigma, T *u, integer *ldu, T *u2, integer *ldu2, T *vt, integer *ldvt, T *vt2,
               integer *ldvt2, integer *idxc, integer *ctot, T *z, integer *info)
    {
        lasd3(nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2, ldu2, vt, ldvt, vt2, ldvt2, idxc,
              ctot, z, info);
    }
    /** @}*/ // end of lasd3

    /** @defgroup lasd4 lasd4
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD4 computes the square root of the i-th updated eigenvalue  \n
     of a positive symmetric rank-one modification to a positive diagonal \n
     matrix. Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    This subroutine computes the square root of the I-th updated
    eigenvalue of a positive symmetric rank-one modification to
    a positive diagonal matrix whose entries are given as the squares
    of the corresponding entries in the array d, and that

           0 <= D(i) < D(j)  for  i < j

    and that RHO > 0. This is arranged by the calling routine, and is
    no loss in generality.  The rank-one modified system is thus

           diag(D) * diag(D) +  RHO * Z * Z_transpose.

    where we assume the Euclidean norm of Z is 1.

    The method consists of approximating the rational functions in the
    secular equation by simpler interpolating rational functions.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The length of all arrays. \n
 * @param[in] I
          I is INTEGER \n
          The index of the eigenvalue to be computed.  1 <= I <= N. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The original eigenvalues.  It is assumed that they are in
          order, 0 <= D(I) < D(J)  for I < J. \n
 * @param[in] Z
          Z is REAL array, dimension (N) \n
          The components of the updating vector. \n
 * @param[out] DELTA
          DELTA is REAL array, dimension (N) \n
          If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th
          component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
          contains the information necessary to construct the
          (singular) eigenvectors. \n
 * @param[in] RHO
          RHO is REAL \n
          The scalar in the symmetric updating formula. \n
 * @param[out] SIGMA
          SIGMA is REAL \n
          The computed sigma_I, the I-th updated eigenvalue. \n
 * @param[out] WORK
          WORK is REAL array, dimension (N) \n
          If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th
          component.  If N = 1, then WORK(1) = 1. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          > 0:  if INFO = 1, the updating process failed.  \n

 *  * */
    template <typename T>
    void lasd4(integer *n, integer *i, T *d, T *z, T *delta, T *rho, T *sigma, T *work,
               integer *info)
    {
        lasd4(n, i, d, z, delta, rho, sigma, work, info);
    }
    /** @}*/ // end of lasd4

    /** @defgroup lasdq lasdq
     * @ingroup SVD
     * @{
     */
    /*! @brief LASDQ computes the SVD of a real bidiagonal matrix with diagonal \n
     d and off-diagonal e. Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    LASDQ computes the singular value decomposition (SVD) of a real
    (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
    E, accumulating the transformations if desired. Letting B denote
    the input bidiagonal matrix, the algorithm computes orthogonal
    matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose
    of P). The singular values S are overwritten on D.

    The input matrix U  is changed to U  * Q  if desired.
    The input matrix VT is changed to P**T * VT if desired.
    The input matrix C  is changed to Q**T * C  if desired.

    See "Computing  Small Singular Values of Bidiagonal Matrices With
    Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
    LAPACK Working Note #3, for a detailed description of the algorithm.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          On entry, UPLO specifies whether the input bidiagonal matrix
          is upper or lower bidiagonal, and whether it is square are
          not. \n
             UPLO = 'U' or 'u'   B is upper bidiagonal. \n
             UPLO = 'L' or 'l'   B is lower bidiagonal. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          = 0: then the input matrix is N-by-N. \n
          = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and
               (N+1)-by-N if UPLU = 'L'. \n
 \n
          The bidiagonal matrix has \n
          N = NL + NR + 1 rows and \n
          M = N + SQRE >= N columns. \n
 * @param[in] N
          N is INTEGER \n
          On entry, N specifies the number of rows and columns
          in the matrix. N must be at least 0. \n
 * @param[in] NCVT
          NCVT is INTEGER \n
          On entry, NCVT specifies the number of columns of
          the matrix VT. NCVT must be at least 0. \n
 * @param[in] NRU
          NRU is INTEGER \n
          On entry, NRU specifies the number of rows of
          the matrix U. NRU must be at least 0. \n
 * @param[in] NCC
          NCC is INTEGER \n
          On entry, NCC specifies the number of columns of
          the matrix C. NCC must be at least 0. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, D contains the diagonal entries of the
          bidiagonal matrix whose SVD is desired. On normal exit,
          D contains the singular values in ascending order. \n
 * @param[in,out] E
          E is REAL array. \n
          dimension is (N-1) if SQRE = 0 and N if SQRE = 1. \n
          On entry, the entries of E contain the offdiagonal entries
          of the bidiagonal matrix whose SVD is desired. On normal
          exit, E will contain 0. If the algorithm does not converge,
          D and E will contain the diagonal and superdiagonal entries
          of a bidiagonal matrix orthogonally equivalent to the one
          given as input. \n
 * @param[in,out] VT
          VT is REAL array, dimension (LDVT, NCVT) \n
          On entry, contains a matrix which on exit has been
          premultiplied by P**T, dimension N-by-NCVT if SQRE = 0
          and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0). \n
 * @param[in] LDVT
          LDVT is INTEGER \n
          On entry, LDVT specifies the leading dimension of VT as
          declared in the calling (sub) program. LDVT must be at
          least 1. If NCVT is nonzero LDVT must also be at least N. \n
 * @param[in,out] U
          U is REAL array, dimension (LDU, N) \n
          On entry, contains a  matrix which on exit has been
          postmultiplied by Q, dimension NRU-by-N if SQRE = 0
          and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0). \n
 * @param[in] LDU
          LDU is INTEGER \n
          On entry, LDU  specifies the leading dimension of U as
          declared in the calling (sub) program. LDU must be at
          least fla_max(1, NRU). \n
 * @param[in,out] C
          C is REAL array, dimension (LDC, NCC) \n
          On entry, contains an N-by-NCC matrix which on exit
          has been premultiplied by Q**T  dimension N-by-NCC if SQRE = 0
          and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0). \n
 * @param[in] LDC
          LDC is INTEGER \n
          On entry, LDC  specifies the leading dimension of C as
          declared in the calling (sub) program. LDC must be at
          least 1. If NCC is nonzero, LDC must also be at least N. \n
 * @param[out] WORK
          WORK is REAL array, dimension (4*N) \n
          Workspace. Only referenced if one of NCVT, NRU, or NCC is
          nonzero, and if N is at least 2. \n
 * @param[out] INFO
          INFO is INTEGER \n
          On exit, a value of 0 indicates a successful exit. \n
          If INFO < 0, argument number -INFO is illegal. \n
          If INFO > 0, the algorithm did not converge, and INFO
          specifies how many superdiagonals did not converge.  \n

 *  * */
    template <typename T>
    void lasdq(char *uplo, integer *sqre, integer *n, integer *ncvt, integer *nru, integer *ncc,
               T *d, T *e, T *vt, integer *ldvt, T *u, integer *ldu, T *c, integer *ldc, T *work,
               integer *info)
    {
        lasdq(uplo, sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info);
    }
    /** @}*/ // end of lasdq

    /** @defgroup lasda lasda
     * @ingroup SVD
     * @{
     */
    /*! @brief LASDA computes the singular value decomposition (SVD) of a \n
     real upper bidiagonal matrix with diagonal d and off-diagonal e. \n
     Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    Using a divide and conquer approach, LASDA computes the singular
    value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
    B with diagonal D and offdiagonal E, where M = N + SQRE. The
    algorithm computes the singular values in the SVD B = U * S * VT.
    The orthogonal matrices U and VT are optionally computed in
    compact form.

    A related subroutine, SLASD0, computes the singular values and
    the singular vectors in explicit form.
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          Specifies whether singular vectors are to be computed
          in compact form, as follows \n
          = 0: Compute singular values only. \n
          = 1: Compute singular vectors of upper bidiagonal
               matrix in compact form. \n
 * @param[in] SMLSIZ
          SMLSIZ is INTEGER \n
          The maximum size of the subproblems at the bottom of the
          computation tree. \n
 * @param[in] N
          N is INTEGER \n
          The row dimension of the upper bidiagonal matrix. This is
          also the dimension of the main diagonal array D. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          Specifies the column dimension of the bidiagonal matrix. \n
          = 0: The bidiagonal matrix has column dimension M = N; \n
          = 1: The bidiagonal matrix has column dimension M = N + 1. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry D contains the main diagonal of the bidiagonal
          matrix. On exit D, if INFO = 0, contains its singular values. \n
 * @param[in] E
          E is REAL array, dimension (M-1) \n
          Contains the subdiagonal entries of the bidiagonal matrix.
          On exit, E has been destroyed. \n
 * @param[out] U
          U is REAL array, \n
          dimension (LDU, SMLSIZ) if ICOMPQ = 1, and not referenced
          if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left
          singular vector matrices of all subproblems at the bottom
          level. \n
 * @param[in] LDU
          LDU is INTEGER, LDU = > N. \n
          The leading dimension of arrays U, VT, DIFL, DIFR, POLES,
          GIVNUM, and Z. \n
 * @param[out] VT
          VT is REAL array, \n
          dimension (LDU, SMLSIZ+1) if ICOMPQ = 1, and not referenced
          if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT**T contains the right
          singular vector matrices of all subproblems at the bottom
          level. \n
 * @param[out] K
          K is INTEGER array, dimension (N) \n
          if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0. \n
          If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th
          secular equation on the computation tree. \n
 * @param[out] DIFL
          DIFL is REAL array, dimension (LDU, NLVL), \n
          where NLVL = floor(log_2 (N/SMLSIZ))). \n
 * @param[out] DIFR
          DIFR is REAL array, \n
                  dimension (LDU, 2 * NLVL) if ICOMPQ = 1 and \n
                  dimension (N) if ICOMPQ = 0. \n
          If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1)
          record distances between singular values on the I-th
          level and singular values on the (I -1)-th level, and
          DIFR(1:N, 2 * I) contains the normalizing factors for
          the right singular vector matrix. See SLASD8 for details. \n
 * @param[out] Z
          Z is REAL array, \n
                  dimension (LDU, NLVL) if ICOMPQ = 1 and \n
                  dimension (N) if ICOMPQ = 0. \n
          The first K elements of Z(1, I) contain the components of
          the deflation-adjusted updating row vector for subproblems
          on the I-th level. \n
 * @param[out] POLES
          POLES is REAL array, \n
          dimension (LDU, 2 * NLVL) if ICOMPQ = 1, and not referenced
          if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and
          POLES(1, 2*I) contain  the new and old singular values
          involved in the secular equations on the I-th level. \n
 * @param[out] GIVPTR
          GIVPTR is INTEGER array, \n
          dimension (N) if ICOMPQ = 1, and not referenced if
          ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR(I) records
          the number of Givens rotations performed on the I-th
          problem on the computation tree. \n
 * @param[out] GIVCOL
          GIVCOL is INTEGER array, \n
          dimension (LDGCOL, 2 * NLVL) if ICOMPQ = 1, and not
          referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
          GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations
          of Givens rotations performed on the I-th level on the
          computation tree. \n
 * @param[in] LDGCOL
          LDGCOL is INTEGER, LDGCOL = > N. \n
          The leading dimension of arrays GIVCOL and PERM. \n
 * @param[out] PERM
          PERM is INTEGER array, dimension (LDGCOL, NLVL) \n
          if ICOMPQ = 1, and not referenced \n
          if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records
          permutations done on the I-th level of the computation tree. \n
 * @param[out] GIVNUM
          GIVNUM is REAL array, \n
          dimension (LDU, 2 * NLVL) if ICOMPQ = 1, and not \n
          referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
          GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S-
          values of Givens rotations performed on the I-th level on
          the computation tree. \n
 * @param[out] C
          C is REAL array, \n
          dimension (N) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0.
          If ICOMPQ = 1 and the I-th subproblem is not square, on exit,
          C(I) contains the C-value of a Givens rotation related to
          the right null space of the I-th subproblem. \n
 * @param[out] S
          S is REAL array, dimension (N) if \n
          ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1
          and the I-th subproblem is not square, on exit, S(I)
          contains the S-value of a Givens rotation related to
          the right null space of the I-th subproblem. \n
 * @param[out] WORK
          WORK is REAL array, dimension
          (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)). \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (7*N). \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, a singular value did not converge  \n

 *  * */
    template <typename T>
    void lasda(integer *icompq, integer *smlsiz, integer *n, integer *sqre, T *d, T *e, T *u,
               integer *ldu, T *vt, integer *k, T *difl, T *difr, T *z, T *poles, integer *givptr,
               integer *givcol, integer *ldgcol, integer *perm, T *givnum, T *c, T *s, T *work,
               integer *iwork, integer *info)
    {
        lasda(icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol,
              ldgcol, perm, givnum, c, s, work, iwork, info);
    }
    /** @}*/ // end of lasda

    /** @defgroup lasd6 lasd6
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD6 computes the SVD of an updated upper bidiagonal matrix \n
     obtained by merging two smaller ones by appending a row.           \n
     Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    LASD6 computes the SVD of an updated upper bidiagonal matrix B
    obtained by merging two smaller ones by appending a row. This
    routine is used only for the problem which requires all singular
    values and optionally singular vector matrices in factored form.
    B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.
    A related subroutine, SLASD1, handles the case in which all singular
    values and singular vectors of the bidiagonal matrix are desired.

    SLASD6 computes the SVD as follows:

                  (D1(in)    0    0       0)
      B = U(in) * (  Z1**T   a   Z2**T    b) * VT(in)
                  (  0       0   D2(in)   0)

        = U(out) * (D(out) 0) * VT(out)

    where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
    with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
    elsewhere; and the entry b is empty if SQRE = 0.

    The singular values of B can be computed using D1, D2, the first
    components of all the right singular vectors of the lower block, and
    the last components of all the right singular vectors of the upper
    block. These components are stored and updated in VF and VL,
    respectively, in SLASD6. Hence U and VT are not explicitly
    referenced.

    The singular values are stored in D. The algorithm consists of two
    stages:

          The first stage consists of deflating the size of the problem
          when there are multiple singular values or if there is a zero
          in the Z vector. For each such occurrence the dimension of the
          secular equation problem is reduced by one. This stage is
          performed by the routine SLASD7.

          The second stage consists of calculating the updated
          singular values. This is done by finding the roots of the
          secular equation via the routine SLASD4 (as called by SLASD8).
          This routine also updates VF and VL and computes the distances
          between the updated singular values and the old singular
          values.

    SLASD6 is called from SLASDA.
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          Specifies whether singular vectors are to be computed in
          factored form: \n
          = 0: Compute singular values only. \n
          = 1: Compute singular vectors in factored form as well. \n
 * @param[in] NL
          NL is INTEGER \n
          The row dimension of the upper block.  NL >= 1. \n
 * @param[in] NR
          NR is INTEGER \n
          The row dimension of the lower block.  NR >= 1. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          = 0: the lower block is an NR-by-NR square matrix. \n
          = 1: the lower block is an NR-by-(NR+1) rectangular matrix. \n
 \n
          The bidiagonal matrix has row dimension N = NL + NR + 1,
          and column dimension M = N + SQRE. \n
 * @param[in,out] D
          D is REAL array, dimension (NL+NR+1). \n
          On entry D(1:NL,1:NL) contains the singular values of the
          upper block, and D(NL+2:N) contains the singular values
          of the lower block. On exit D(1:N) contains the singular
          values of the modified matrix. \n
 * @param[in,out] VF
          VF is REAL array, dimension (M) \n
          On entry, VF(1:NL+1) contains the first components of all
          right singular vectors of the upper block; and VF(NL+2:M)
          contains the first components of all right singular vectors
          of the lower block. On exit, VF contains the first components
          of all right singular vectors of the bidiagonal matrix. \n
 * @param[in,out] VL
          VL is REAL array, dimension (M) \n
          On entry, VL(1:NL+1) contains the  last components of all
          right singular vectors of the upper block; and VL(NL+2:M)
          contains the last components of all right singular vectors of
          the lower block. On exit, VL contains the last components of
          all right singular vectors of the bidiagonal matrix. \n
 * @param[in,out] ALPHA
          ALPHA is REAL \n
          Contains the diagonal element associated with the added row. \n
 * @param[in,out] BETA
          BETA is REAL \n
          Contains the off-diagonal element associated with the added
          row. \n
 * @param[in,out] IDXQ
          IDXQ is INTEGER array, dimension (N) \n
          This contains the permutation which will reintegrate the
          subproblem just solved back into sorted order, i.e.
          D(IDXQ(I = 1, N)) will be in ascending order. \n
 * @param[out] PERM
          PERM is INTEGER array, dimension (N) \n
          The permutations (from deflation and sorting) to be applied
          to each block. Not referenced if ICOMPQ = 0. \n
 * @param[out] GIVPTR
          GIVPTR is INTEGER \n
          The number of Givens rotations which took place in this
          subproblem. Not referenced if ICOMPQ = 0. \n
 * @param[out] GIVCOL
          GIVCOL is INTEGER array, dimension (LDGCOL, 2) \n
          Each pair of numbers indicates a pair of columns to take place
          in a Givens rotation. Not referenced if ICOMPQ = 0. \n
 * @param[in] LDGCOL
          LDGCOL is INTEGER \n
          leading dimension of GIVCOL, must be at least N. \n
 * @param[out] GIVNUM
          GIVNUM is REAL array, dimension (LDGNUM, 2) \n
          Each number indicates the C or S value to be used in the
          corresponding Givens rotation. Not referenced if ICOMPQ = 0. \n
 * @param[in] LDGNUM
          LDGNUM is INTEGER \n
          The leading dimension of GIVNUM and POLES, must be at least N. \n
 * @param[out] POLES
          POLES is REAL array, dimension (LDGNUM, 2) \n
          On exit, POLES(1,*) is an array containing the new singular
          values obtained from solving the secular equation, and
          POLES(2,*) is an array containing the poles in the secular
          equation. Not referenced if ICOMPQ = 0. \n
 * @param[out] DIFL
          DIFL is REAL array, dimension (N) \n
          On exit, DIFL(I) is the distance between I-th updated
          (undeflated) singular value and the I-th (undeflated) old
          singular value. \n
 * @param[out] DIFR
          DIFR is REAL array, \n
                   dimension (LDDIFR, 2) if ICOMPQ = 1 and
                   dimension (K) if ICOMPQ = 0. \n
          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not
          defined and will not be referenced. \n
 \n
          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the
          normalizing factors for the right singular vector matrix.
 \n
          See SLASD8 for details on DIFL and DIFR. \n
 * @param[out] Z
          Z is REAL array, dimension (M) \n
          The first elements of this array contain the components
          of the deflation-adjusted updating row vector. \n
 * @param[out] K
          K is INTEGER \n
          Contains the dimension of the non-deflated matrix,
          This is the order of the related secular equation. 1 <= K <=N. \n
 * @param[out] C
          C is REAL \n
          C contains garbage if SQRE =0 and the C-value of a Givens
          rotation related to the right null space if SQRE = 1. \n
 * @param[out] S
          S is REAL \n
          S contains garbage if SQRE =0 and the S-value of a Givens
          rotation related to the right null space if SQRE = 1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (4 * M) \n
 * @param[out] IWORK
          IWORK is INTEGER array, dimension (3 * N) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, a singular value did not converge  \n

 *  * */
    template <typename T>
    void lasd6(integer *icompq, integer *nl, integer *nr, integer *sqre, T *d, T *vf, T *vl,
               T *alpha, T *beta, integer *idxq, integer *perm, integer *givptr, integer *givcol,
               integer *ldgcol, T *givnum, integer *ldgnum, T *poles, T *difl, T *difr, T *z,
               integer *k, T *c, T *s, T *work, integer *iwork, integer *info)
    {
        lasd6(icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq, perm, givptr, givcol, ldgcol,
              givnum, ldgnum, poles, difl, difr, z, k, c, s, work, iwork, info);
    }
    /** @}*/ // end of lasd6

    /** @defgroup lasd7 lasd7
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD7 merges the two sets of singular values together into a   \n
     single sorted set. Then it tries to deflate the size of the problem. \n
     Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    LASD7 merges the two sets of singular values together into a single
    sorted set. Then it tries to deflate the size of the problem. There
    are two ways in which deflation can occur:  when two or more singular
    values are close together or if there is a tiny entry in the Z
    vector. For each such occurrence the order of the related
    secular equation problem is reduced by one.

    SLASD7 is called from SLASD6.
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          Specifies whether singular vectors are to be computed
          in compact form, as follows: \n
          = 0: Compute singular values only. \n
          = 1: Compute singular vectors of upper
               bidiagonal matrix in compact form. \n
 * @param[in] NL
          NL is INTEGER \n
          The row dimension of the upper block. NL >= 1. \n
 * @param[in] NR
          NR is INTEGER \n
          The row dimension of the lower block. NR >= 1. \n
 * @param[in] SQRE
          SQRE is INTEGER \n
          = 0: the lower block is an NR-by-NR square matrix. \n
          = 1: the lower block is an NR-by-(NR+1) rectangular matrix. \n
 \n
          The bidiagonal matrix has \n
          N = NL + NR + 1 rows and \n
          M = N + SQRE >= N columns. \n
 * @param[out] K
          K is INTEGER \n
          Contains the dimension of the non-deflated matrix, this is
          the order of the related secular equation. 1 <= K <=N. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry D contains the singular values of the two submatrices
          to be combined. On exit D contains the trailing (N-K) updated
          singular values (those which were deflated) sorted into
          increasing order. \n
 * @param[out] Z
          Z is REAL array, dimension (M) \n
          On exit Z contains the updating row vector in the secular
          equation. \n
 * @param[out] ZW
          ZW is REAL array, dimension (M) \n
          Workspace for Z. \n
 * @param[in,out] VF
          VF is REAL array, dimension (M) \n
          On entry, VF(1:NL+1) contains the first components of all
          right singular vectors of the upper block; and VF(NL+2:M)
          contains the first components of all right singular vectors
          of the lower block. On exit, VF contains the first components
          of all right singular vectors of the bidiagonal matrix. \n
 * @param[out] VFW
          VFW is REAL array, dimension (M) \n
          Workspace for VF. \n
 * @param[in,out] VL
          VL is REAL array, dimension (M) \n
          On entry, VL(1:NL+1) contains the  last components of all
          right singular vectors of the upper block; and VL(NL+2:M)
          contains the last components of all right singular vectors
          of the lower block. On exit, VL contains the last components
          of all right singular vectors of the bidiagonal matrix. \n
 * @param[out] VLW
          VLW is REAL array, dimension (M) \n
          Workspace for VL. \n
 * @param[in] ALPHA
          ALPHA is REAL \n
          Contains the diagonal element associated with the added row. \n
 * @param[in] BETA
          BETA is REAL \n
          Contains the off-diagonal element associated with the added
          row. \n
 * @param[out] DSIGMA
          DSIGMA is REAL array, dimension (N) \n
          Contains a copy of the diagonal elements (K-1 singular values
          and one zero) in the secular equation. \n
 * @param[out] IDX
          IDX is INTEGER array, dimension (N) \n
          This will contain the permutation used to sort the contents of
          D into ascending order. \n
 * @param[out] IDXP
          IDXP is INTEGER array, dimension (N) \n
          This will contain the permutation used to place deflated
          values of D at the end of the array. On output IDXP(2:K)
          points to the nondeflated D-values and IDXP(K+1:N)
          points to the deflated singular values. \n
 * @param[in] IDXQ
          IDXQ is INTEGER array, dimension (N) \n
          This contains the permutation which separately sorts the two
          sub-problems in D into ascending order.  Note that entries in
          the first half of this permutation must first be moved one
          position backward; and entries in the second half
          must first have NL+1 added to their values. \n
 * @param[out] PERM
          PERM is INTEGER array, dimension (N) \n
          The permutations (from deflation and sorting) to be applied
          to each singular block. Not referenced if ICOMPQ = 0. \n
 * @param[out] GIVPTR
          GIVPTR is INTEGER \n
          The number of Givens rotations which took place in this
          subproblem. Not referenced if ICOMPQ = 0. \n
 * @param[out] GIVCOL
          GIVCOL is INTEGER array, dimension (LDGCOL, 2) \n
          Each pair of numbers indicates a pair of columns to take place
          in a Givens rotation. Not referenced if ICOMPQ = 0. \n
 * @param[in] LDGCOL
          LDGCOL is INTEGER \n
          The leading dimension of GIVCOL, must be at least N. \n
 * @param[out] GIVNUM
          GIVNUM is REAL array, dimension (LDGNUM, 2) \n
          Each number indicates the C or S value to be used in the
          corresponding Givens rotation. Not referenced if ICOMPQ = 0. \n
 * @param[in] LDGNUM
          LDGNUM is INTEGER \n
          The leading dimension of GIVNUM, must be at least N. \n
 * @param[out] C
          C is REAL \n
          C contains garbage if SQRE =0 and the C-value of a Givens
          rotation related to the right null space if SQRE = 1. \n
 * @param[out] S
          S is REAL \n
          S contains garbage if SQRE =0 and the S-value of a Givens
          rotation related to the right null space if SQRE = 1. \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value.  \n

 *  * */
    template <typename T>
    void lasd7(integer *icompq, integer *nl, integer *nr, integer *sqre, integer *k, T *d, T *z,
               T *zw, T *vf, T *vfw, T *vl, T *vlw, T *alpha, T *beta, T *dsigma, integer *idx,
               integer *idxp, integer *idxq, integer *perm, integer *givptr, integer *givcol,
               integer *ldgcol, T *givnum, integer *ldgnum, T *c, T *s, integer *info)
    {
        lasd7(icompq, nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma, idx, idxp,
              idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info);
    }
    /** @}*/ // end of lasd7

    /** @defgroup lasd8 lasd8
     * @ingroup SVD
     * @{
     */
    /*! @brief LASD8 finds the square roots of the roots of the secular equation, \n
     and stores, for each element in D, the distance to its two nearest poles. \n
     Used by sbdsdc
 * @details
 * \b Purpose:
    \verbatim
    LASD8 finds the square roots of the roots of the secular equation,
    as defined by the values in DSIGMA and Z. It makes the appropriate
    calls to SLASD4, and stores, for each  element in D, the distance
    to its two nearest poles (elements in DSIGMA). It also updates
    the arrays VF and VL, the first and last components of all the
    right singular vectors of the original bidiagonal matrix.

    LASD8 is called from LASD6.
    \endverbatim

 * @param[in] ICOMPQ
          ICOMPQ is INTEGER \n
          Specifies whether singular vectors are to be computed in
          factored form in the calling routine: \n
          = 0: Compute singular values only. \n
          = 1: Compute singular vectors in factored form as well. \n
 * @param[in] K
          K is INTEGER \n
          The number of terms in the rational function to be solved
          by SLASD4.  K >= 1. \n
 * @param[out] D
          D is REAL array, dimension (K) \n
          On output, D contains the updated singular values. \n
 * @param[in,out] Z
          Z is REAL array, dimension (K) \n
          On entry, the first K elements of this array contain the
          components of the deflation-adjusted updating row vector. \n
          On exit, Z is updated. \n
 * @param[in,out] VF
          VF is REAL array, dimension (K) \n
          On entry, VF contains  information passed through DBEDE8. \n
          On exit, VF contains the first K components of the first
          components of all right singular vectors of the bidiagonal
          matrix. \n
 * @param[in,out] VL
          VL is REAL array, dimension (K) \n
          On entry, VL contains  information passed through DBEDE8. \n
          On exit, VL contains the first K components of the last
          components of all right singular vectors of the bidiagonal
          matrix. \n
 * @param[out] DIFL
          DIFL is REAL array, dimension (K) \n
          On exit, DIFL(I) = D(I) - DSIGMA(I). \n
 * @param[out] DIFR
          DIFR is REAL array, \n
                   dimension (LDDIFR, 2) if ICOMPQ = 1 and
                   dimension (K) if ICOMPQ = 0. \n
          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not
          defined and will not be referenced. \n
 \n
          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the
          normalizing factors for the right singular vector matrix. \n
 * @param[in] LDDIFR
          LDDIFR is INTEGER \n
          The leading dimension of DIFR, must be at least K. \n
 * @param[in,out] DSIGMA
          DSIGMA is REAL array, dimension (K) \n
          On entry, the first K elements of this array contain the old
          roots of the deflated updating problem.  These are the poles
          of the secular equation. \n
          On exit, the elements of DSIGMA may be very slightly altered
          in value. \n
 * @param[out] WORK
          WORK is REAL array, dimension (3*K) \n
 * @param[out] INFO
          INFO is INTEGER \n
          = 0:  successful exit. \n
          < 0:  if INFO = -i, the i-th argument had an illegal value. \n
          > 0:  if INFO = 1, a singular value did not converge  \n

 *  * */
    template <typename T>
    void lasd8(integer *icompq, integer *k, T *d, T *z, T *vf, T *vl, T *difl, T *difr,
               integer *lddifr, T *dsigma, T *work, integer *info)
    {
        lasd8(icompq, k, d, z, vf, vl, difl, difr, lddifr, dsigma, work, info);
    }
    /** @}*/ // end of lasd8
    /** @}*/ // end of SVD Computational
    /** @}*/ // end of SVD
}
#endif