/******************************************************************************
 * Copyright (C) 2021-2025, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file libflame_interface_blas.hh
 *  libflame_interface.hh defines all the BLAS-like routines for Libflame CPP templated public
 *  interfaces.
 *  */
#ifndef LIBFLAME_INTERFACE_BLAS_HH
#define LIBFLAME_INTERFACE_BLAS_HH


namespace libflame
{
     /** @defgroup BLAS BLAS-like
     * @ingroup AOCL_LAPACK
     * @{
     */

    /** @defgroup general Initialize, copy, convert
     * @ingroup BLAS BLAS-like
     * @{
     */
    /** @defgroup laset laset
     * @ingroup general
     * @{
     */
    /*! @brief Initialize the off-diagonal elements and the diagonal elements of a matrix to given
    values
    *
    * @details
    * \b Purpose:
    * \verbatim
      Initialize an m-by-n matrix a to beta on the diagonal and alpha on the offdiagonals.
      \endverbatim

    * @param[in] uplo
              uplo is char* \n
              Specifies the part of the matrix a to be set. \n
              = 'U':   Upper triangular part is set; the strictly lower
              triangular part of a is not changed. \n
              = 'L':   Lower triangular part is set; the strictly upper
              triangular part of a is not changed. \n
              Otherwise:  All of the matrix a is set. \n
    * @param[in] m
              m is integer* \n
              The number of rows of the matrix a.  m >= 0. \n
    * @param[in] n
              n is integer* \n
              The number of columns of the matrix a.  n >= 0. \n
    * @param[in] alpha
              alpha is float/double/COMPLEX/COMPLEX*16* \n
              The constant to which the offdiagonal elements are to be set. \n
    * @param[in] beta
              beta is float/double/COMPLEX/COMPLEX*16* \n
              The constant to which the diagonal elements are to be set. \n
    * @param[out] a
              a is float/double/COMPLEX/COMPLEX*16 array, dimension (lda,n) \n
              On exit, the leading m-by-n submatrix of A is set as follows: \n \n
              if uplo = 'U', A(i,j) = alpha, 1<=i<=j-1, 1<=j<=n, \n
              if uplo = 'L', A(i,j) = alpha, j+1<=i<=m, 1<=j<=n, \n
              otherwise,  A(i,j) = alpha, 1<=i<=m, 1<=j<=n, i.ne.j, \n \n
              and, for all uplo, A(i,i) = beta, 1<=i<=min(m,n). \n
    * @param[in] lda
              lda is integer* \n
              The leading dimension of the array a.  lda >= fla_max(1,m). \n

    *  */
    template <typename T>
    void laset(char *uplo, integer *m, integer *n, T *alpha, T *beta, T *a, integer *lda)
    {
        laset(uplo, m, n, alpha, beta, a, lda);
    }
    /** @}*/ // end of laset

    /** @defgroup larnv larnv
     * @ingroup general
     * @{
     */
    /*! @brief LARNV returns a vector of random numbers from a uniform or normal distribution
 * @details
 * \b Purpose:
    \verbatim
     LARNV returns a vector of n random real numbers from a uniform or
     normal distribution.
    \endverbatim

 * @param[in] IDIST
          IDIST is INTEGER \n
          Specifies the distribution of the random numbers: \n
          = 1:  uniform (0,1) \n
          = 2:  uniform (-1,1) \n
          = 3:  normal (0,1) \n
 * @param[in,out] ISEED
          ISEED is INTEGER array, dimension (4) \n
          On entry, the seed of the random number generator; the array
          elements must be between 0 and 4095, and ISEED(4) must be
          odd. \n
          On exit, the seed is updated. \n
 * @param[in] N
          N is INTEGER \n
          The number of random numbers to be generated. \n
 * @param[out] X
          X is REAL array, dimension (N) \n
          The generated random numbers. \n

 * */
    template <typename T>
    void larnv(integer *idist, integer *iseed, integer *n, T *x)
    {
        larnv(idist, iseed, n, x);
    }
    /** @}*/ // end of larnv

    /** @defgroup laruv laruv
     * @ingroup general
     * @{
     */
    /*! @brief LARUV returns a vector of n random real numbers from a uniform distribution

 * @details
 * \b Purpose:
    \verbatim
    LARUV returns a vector of n random real numbers from a uniform (0,1)
    distribution (n <= 128).

    This is an auxiliary routine called by SLARNV and CLARNV.
    \endverbatim

 * @param[in,out] ISEED
          ISEED is INTEGER array, dimension (4) \n
          On entry, the seed of the random number generator; the array
          elements must be between 0 and 4095, and ISEED(4) must be
          odd. \n
          On exit, the seed is updated. \n
 * @param[in] N
          N is INTEGER \n
          The number of random numbers to be generated. N <= 128. \n
 * @param[out] X
          X is REAL array, dimension (N) \n
          The generated random numbers. \n

 *  * */
    template <typename T>
    void laruv(integer *iseed, integer *n, T *x)
    {
        laruv(iseed, n, x);
    }
    /** @}*/ // end of laruv

    /** @defgroup lacpy lacpy
     * @ingroup general
     * @{
     */
    /*! @brief LACPY copies all or part of one two-dimensional array to another.

 * @details
 * \b Purpose:
    \verbatim
     LACPY copies all or part of a two-dimensional matrix A to another
     matrix B.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies the part of the matrix A to be copied to B. \n
          = 'U':      Upper triangular part \n
          = 'L':      Lower triangular part \n
          Otherwise:  All of the matrix A \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The m by n matrix A.  If UPLO = 'U', only the upper triangle
          or trapezoid is accessed; if UPLO = 'L', only the lower
          triangle or trapezoid is accessed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] B
          B is REAL array, dimension (LDB,N) \n
          On exit, B = A in the locations specified by UPLO. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,M). \n

 * */
    template <typename T>
    void lacpy(char *uplo, integer *m, integer *n, T *a, integer *lda, T *b, integer *ldb)
    {
        lacpy(uplo, m, n, a, lda, b, ldb);
    }
    /** @}*/ // end of lacpy

    /** @defgroup lacp2 lacp2
     * @ingroup general
     * @{
     */
    /*! @brief LACP2 copies all or part of a real two-dimensional array to a complex array

 * @details
 * \b Purpose:
    \verbatim
     LACP2 copies all or part of a real two-dimensional matrix A to a
     complex matrix B.
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies the part of the matrix A to be copied to B. \n
          = 'U':      Upper triangular part \n
          = 'L':      Lower triangular part \n
          Otherwise:  All of the matrix A \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
          is accessed; if UPLO = 'L', only the lower trapezium is
          accessed. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,M). \n
 * @param[out] B
          B is COMPLEX array, dimension (LDB,N) \n
          On exit, B = A in the locations specified by UPLO. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(1,M). \n

 * */
    template <typename T, typename Ta>
    void lacp2(char *uplo, integer *m, integer *n, T *a, integer *lda, Ta *b, integer *ldb)
    {
        lacp2(uplo, m, n, a, lda, b, ldb);
    }
    /** @}*/ // end of lacp2

    /** @defgroup tfttp tfttp
     * @ingroup general
     * @{
     */
    /*! @brief TFTTP copies a triangular matrix from the rectangular full \n
     packed format (TF) to the standard packed format (TP).

 * @details
 * \b Purpose:
    \verbatim
     TFTTP copies a triangular matrix A from rectangular full packed
     format (TF) to standard packed format (TP).
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  ARF is in Normal format; \n
          = 'T':  ARF is in Transpose format; \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in] ARF
          ARF is REAL array, dimension ( N*(N+1)/2), \n
          On entry, the upper or lower triangular matrix A stored in
          RFP format. For a further discussion see Notes below. \n
 * @param[out] AP
          AP is REAL array, dimension ( N*(N+1)/2), \n
          On exit, the upper or lower triangular matrix A, packed
          columnwise in a linear array. The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tfttp(char *transr, char *uplo, integer *n, T *arf, T *ap, integer *info)
    {
        tfttp(transr, uplo, n, arf, ap, info);
    }
    /** @}*/ // end of tfttp

    /** @defgroup tfttr tfttr
     * @ingroup general
     * @{
     */
    /*! @brief TFTTR copies a triangular matrix from the rectangular full \n
     packed format (TF) to the standard full format (TR).
 * @details
 * \b Purpose:
    \verbatim
     TFTTR copies a triangular matrix A from rectangular full packed
     format (TF) to standard full format (TR).
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  ARF is in Normal format; \n
          = 'T':  ARF is in Transpose format. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices ARF and A. N >= 0. \n
 * @param[in] ARF
          ARF is REAL array, dimension (N*(N+1)/2). \n
          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')
          matrix A in RFP format. See the "Notes" below for more
          details. \n
 * @param[out] A
          A is REAL array, dimension (LDA,N) \n
          On exit, the triangular matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of the array A contains
          the upper triangular matrix, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of the array A contains
          the lower triangular matrix, and the strictly upper
          triangular part of A is not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tfttr(char *transr, char *uplo, integer *n, T *arf, T *a, integer *lda, integer *info)
    {
        tfttr(transr, uplo, n, arf, a, lda, info);
    }
    /** @}*/ // end of tfttr

    /** @defgroup tpttf tpttf
     * @ingroup general
     * @{
     */
    /*! @brief TPTTF copies a triangular matrix from the standard packed \n
     format (TP) to the rectangular full packed format (TF).
 * @details
 * \b Purpose:
    \verbatim
     TPTTF copies a triangular matrix A from standard packed format (TP)
     to rectangular full packed format (TF).
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  ARF in Normal format is wanted; \n
          = 'T':  ARF in Conjugate-transpose format is wanted. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular; \n
          = 'L':  A is lower triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension ( N*(N+1)/2), \n
          On entry, the upper or lower triangular matrix A, packed
          columnwise in a linear array. The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[out] ARF
          ARF is REAL array, dimension ( N*(N+1)/2), \n
          On exit, the upper or lower triangular matrix A stored in
          RFP format. For a further discussion see Notes below. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tpttf(char *transr, char *uplo, integer *n, T *ap, T *arf, integer *info)
    {
        tpttf(transr, uplo, n, ap, arf, info);
    }
    /** @}*/ // end of tpttf

    /** @defgroup tpttr tpttr
     * @ingroup general
     * @{
     */
    /*! @brief TPTTR copies a triangular matrix from the standard    \n
     packed format (TP) to the standard full format (TR).
 * @details
 * \b Purpose:
    \verbatim
     TPTTR copies a triangular matrix A from standard packed format (TP)
     to standard full format (TR).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular. \n
          = 'L':  A is lower triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in] AP
          AP is REAL array, dimension ( N*(N+1)/2), \n
          On entry, the upper or lower triangular matrix A, packed
          columnwise in a linear array. The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[out] A
          A is REAL array, dimension ( LDA, N) \n
          On exit, the triangular matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void tpttr(char *uplo, integer *n, T *ap, T *a, integer *lda, integer *info)
    {
        tpttr(uplo, n, ap, a, lda, info);
    }
    /** @}*/ // end of tpttr

    /** @defgroup trttf trttf
     * @ingroup general
     * @{
     */
    /*! @brief TRTTF copies a triangular matrix from the standard full \n
     format (TR) to the rectangular full packed format (TF).
 * @details
 * \b Purpose:
    \verbatim
     TRTTF copies a triangular matrix A from standard full format (TR)
     to rectangular full packed format (TF) .
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  ARF in Normal form is wanted; \n
          = 'T':  ARF in Transpose form is wanted. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  Upper triangle of A is stored; \n
          = 'L':  Lower triangle of A is stored. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N). \n
          On entry, the triangular matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of the array A contains
          the upper triangular matrix, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of the array A contains
          the lower triangular matrix, and the strictly upper
          triangular part of A is not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the matrix A. LDA >= fla_max(1,N). \n
 * @param[out] ARF
          ARF is REAL array, dimension (NT). \n
          NT=N*(N+1)/2. On exit, the triangular matrix A in RFP format. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void trttf(char *transr, char *uplo, integer *n, T *a, integer *lda, T *arf, integer *info)
    {
        trttf(transr, uplo, n, a, lda, arf, info);
    }
    /** @}*/ // end of trttf

    /** @defgroup trttp trttp
     * @ingroup general
     * @{
     */
    /*! @brief TRTTP copies a triangular matrix from the standard full \n
     format (TR) to the standard packed format (TP).
 * @details
 * \b Purpose:
    \verbatim
    TRTTP copies a triangular matrix A from full format (TR) to standard
    packed format (TP).
    \endverbatim

 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          = 'U':  A is upper triangular. \n
          = 'L':  A is lower triangular. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrices AP and A.  N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          On exit, the triangular matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(1,N). \n
 * @param[out] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          On exit, the upper or lower triangular matrix A, packed
          columnwise in a linear array. The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 *  * */
    template <typename T>
    void trttp(char *uplo, integer *n, T *a, integer *lda, T *ap, integer *info)
    {
        trttp(uplo, n, a, lda, ap, info);
    }
    /** @}*/ // end of trttp
    /** @}*/ // end of general

    /** @defgroup MN Matrix Norm
     * @ingroup BLAS
     * @{
     */
    /** @defgroup lange lange
     * @ingroup MN
     * @{
     */
    /*! @brief LANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, \n
     or the largest absolute value of any element of a general rectangular matrix.
 * @details
 * \b Purpose:
    \verbatim
     LANGE  returns the value of the one norm,  or the Frobenius norm, or
     the  infinity norm,  or the  element of  largest absolute value  of a
     real matrix A.
     \return LANGE
         LANGE = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                 (
                 ( norm1(A),         NORM = '1', 'O' or 'o'
                 (
                 ( normI(A),         NORM = 'I' or 'i'
                 (
                 ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

     where  norm1  denotes the  one norm of a matrix (maximum column sum),
     normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
     normF  denotes the  Frobenius norm of a matrix (square root of sum of
     squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be returned in SLANGE as described
          above. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0.  When M = 0,
          SLANGE is set to zero. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0.  When N = 0,
          SLANGE is set to zero. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The m by n matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(M,1). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
          referenced. \n
 * @return Returns value of the norm.
 * */
    template <typename T>
    T lange(char *norm, integer *m, integer *n, T *a, integer *lda, T *work)
    {
        lange(norm, m, n, a, lda, work);
    }
    template <typename T, typename Ta>
    Ta lange(char *norm, integer *m, integer *n, T *a, integer *lda, Ta *work)
    {
        lange(norm, m, n, a, lda, work);
    }
    /** @}*/ // end of lange

    /** @defgroup langb langb
     * @ingroup  MN
     * @{
     */
    /*! @brief LANGB   returns the value of the 1-norm, Frobenius norm, infinity-norm, \n
     or the largest absolute value of any element of general band matrix
 * @details
 * \b Purpose:
    \verbatim
    LANGB    returns the value of the one norm, or the Frobenius norm, or
    the  infinity norm, or the element of  largest absolute value  of an
    n by n band matrix  A, with kl sub-diagonals and ku super-diagonals.


       SLANGB = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in SLANGB as described
          above. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANGB is
          set to zero. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of sub-diagonals of the matrix A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of super-diagonals of the matrix A.  KU >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th
          column of A is stored in the j-th column of the array AB as
          follows: \n
          AB(ku+1+i-j,j) = A(i,j) for fla_max(1,j-ku)<=i<=min(n,j+kl). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= KL+KU+1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
          referenced.   \n

 * @return Returns the value of the norm.
 * */
    template <typename T>
    T langb(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *work)
    {
        return langb(norm, n, kl, ku, ab, ldab, work);
    }
    template <typename T, typename Ta>
    Ta langb(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, Ta *work)
    {
        return langb(norm, n, kl, ku, ab, ldab, work);
    }
    /** @}*/ // end of langb

    /** @defgroup langt langt
     * @ingroup  MN
     * @{
     */
    /*! @brief LANGT   returns the value of the 1-norm, Frobenius norm, infinity-norm, \n
     or the largest absolute value of any element of a general tridiagonal matrix
 * @details
 * \b Purpose:
    \verbatim
    LANGT    returns the value of the one norm, or the Frobenius norm, or
    the  infinity norm, or the  element of  largest absolute value  of a
    real tridiagonal matrix A.

       SLANGT = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in SLANGT as described
          above. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANGT is
          set to zero. \n
 * @param[in] DL
          DL is REAL array, dimension (N-1) \n
          The (n-1) sub-diagonal elements of A. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The diagonal elements of A. \n
 * @param[in] DU
          DU is REAL array, dimension (N-1) \n
          The (n-1) super-diagonal elements of A. \n

 * @return Returns the value of norm.
 * */
    template <typename T>
    T langt(char *norm, integer *n, T *dl, T *d, T *du)
    {
        return langt(norm, n, dl, d, du);
    }
    template <typename T, typename Ta>
    Ta langt(char *norm, integer *n, T *dl, T *d, T *du)
    {
        return langt(norm, n, dl, d, du);
    }
    /** @}*/ // end of langt

    /** @defgroup lanhs lanhs
     * @ingroup  MN
     * @{
     */
    /*! @brief LANHS returns the value of the 1-norm, Frobenius norm, infinity-norm, or the \n
     largest absolute value of any element of an upper Hessenberg matrix
 * @details
 * \b Purpose:
    \verbatim
    LANHS returns the value of the one norm, or the Frobenius norm, or
    the infinity norm, or the  element of  largest absolute value  of a
    Hessenberg matrix A.

       SLANHS = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in SLANHS as described
          above. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANHS is
          set to zero. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The n by n upper Hessenberg matrix A; the part of A below the
          first sub-diagonal is not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(N,1). \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
          referenced. \n

 * @return Returns the value of norm.
 * */
    template <typename T>
    T lanhs(char *norm, integer *n, T *a, integer *lda, T *work)
    {
        return lanhs(norm, n, a, lda, work);
    }
    template <typename T, typename Ta>
    Ta lanhs(char *norm, integer *n, T *a, integer *lda, Ta *work)
    {
        return lanhs(norm, n, a, lda, work);
    }
    /** @}*/ // end of lanhs

    /** @defgroup lanhe lanhe
     * @ingroup  MN
     * @{
     */
    /*! @brief LANHE returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a complex Hermitian matrix.
 * @details
 * \b Purpose:
    \verbatim
    LANHE  returns the value of the one norm,  or the Frobenius norm, or
    the  infinity norm,  or the  element of  largest absolute value  of a
    complex hermitian matrix A.
 * \b Returns LANHE

       CLANHE = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                ( norm1(A),         NORM = '1', 'O' or 'o'
                (
                ( normI(A),         NORM = 'I' or 'i'
                (
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be returned in CLANHE as described
          above. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          hermitian matrix A is to be referenced. \n
          = 'U':  Upper triangular part of A is referenced \n
          = 'L':  Lower triangular part of A is referenced \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, CLANHE is
          set to zero. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA,N) \n
          The hermitian matrix A.  If UPLO = 'U', the leading n by n
          upper triangular part of A contains the upper triangular part
          of the matrix A, and the strictly lower triangular part of A
          is not referenced.  If UPLO = 'L', the leading n by n lower
          triangular part of A contains the lower triangular part of
          the matrix A, and the strictly upper triangular part of A is
          not referenced. Note that the imaginary parts of the diagonal
          elements need not be set and are assumed to be zero. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(N,1). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced. \n

 * @return Returns value of the norm.
 * */
    template <typename T, typename Ta>
    Ta lanhe(char *norm, char *uplo, integer *n, T *a, integer *lda, Ta *work)
    {
        return lanhe(norm, uplo, n, a, lda, work);
    }

    /** @}*/ // end of lanhe

    /** @defgroup lansy lansy
     * @ingroup  MN
     * @{
     */
    /*! @brief LANSY returns the value of the 1-norm, or the Frobenius norm, or the infinity
      norm, or the element of largest absolute value of a real symmetric matrix.
 * @details
 * \b Purpose:
    \verbatim
     LANSY  returns the value of the one norm,  or the Frobenius norm, or
     the  infinity norm,  or the  element of  largest absolute value  of a
     real symmetric matrix A.

     \return LANSY
        LANSY = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                 (
                 ( norm1(A),         NORM = '1', 'O' or 'o'
                 (
                 ( normI(A),         NORM = 'I' or 'i'
                 (
                 ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

     where  norm1  denotes the  one norm of a matrix (maximum column sum),
     normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
     normF  denotes the  Frobenius norm of a matrix (square root of sum of
     squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be returned in SLANSY as described
          above. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is to be referenced. \n
          = 'U':  Upper triangular part of A is referenced \n
          = 'L':  Lower triangular part of A is referenced \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANSY is
          set to zero. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The symmetric matrix A.  If UPLO = 'U', the leading n by n
          upper triangular part of A contains the upper triangular part
          of the matrix A, and the strictly lower triangular part of A
          is not referenced.  If UPLO = 'L', the leading n by n lower
          triangular part of A contains the lower triangular part of
          the matrix A, and the strictly upper triangular part of A is
          not referenced. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(N,1). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced. \n

 * @return Returns value of the norm.
 * */
    template <typename T>
    T lansy(char *norm, char *uplo, integer *n, T *a, integer *lda, T *work)
    {
        return lansy(norm, uplo, n, a, lda, work);
    }
    template <typename T, typename Ta>
    Ta lansy(char *norm, char *uplo, integer *n, T *a, integer *lda, Ta *work)
    {
        return lansy(norm, uplo, n, a, lda, work);
    }
    /** @}*/ // end of lansy

    /** @defgroup lanhf lanhf
     * @ingroup  MN
     * @{
     */
    /*! @brief LANHF   returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a Hermitian matrix in RFP format
 * @details
 * \b Purpose:
    \verbatim
    LANHF    returns the value of the one norm, or the Frobenius norm, or
    the  infinity norm, or the  element of  largest absolute value  of a
    complex Hermitian matrix A in RFP format.

       CLANHF = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a  matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER \n
          Specifies the value to be   returned in CLANHF as described
          above. \n
 * @param[in] TRANSR
          TRANSR is CHARACTER \n
          Specifies whether the RFP format of A is normal or
          conjugate-transposed format. \n
          = 'N':  RFP format is Normal \n
          = 'C':  RFP format is Conjugate-transposed \n
 * @param[in] UPLO
          UPLO is CHARACTER \n
          On entry, UPLO specifies whether the RFP matrix A came from
          an upper or lower triangular matrix as follows: \n
 \n
          UPLO = 'U' or 'u' RFP A came from an upper triangular
          matrix
 \n
          UPLO = 'L' or 'l' RFP A came from a  lower triangular
          matrix \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, CLANHF is
          set to zero. \n
 * @param[in] A
          A is COMPLEX array, dimension (N*(N+1)/2); \n
          On entry, the matrix A in RFP Format.
          RFP Format is described by TRANSR, UPLO and N as follows:
          If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;
          K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If
          TRANSR = 'C' then RFP is the Conjugate-transpose of RFP A
          as defined when TRANSR = 'N'. The contents of RFP A are
          defined by UPLO as follows: If UPLO = 'U' the RFP A
          contains the (N*(N+1)/2) elements of upper packed A
          either in normal or conjugate-transpose Format. If
          UPLO = 'L' the RFP A contains the (N*(N+1) /2) elements
          of lower packed A either in normal or conjugate-transpose
          Format. The LDA of RFP A is (N+1)/2 when TRANSR = 'C'. When
          TRANSR is 'N' the LDA is N+1 when N is even and is N when
          is odd. See the Note below for more details. \n
          Unchanged on exit. \n
 * @param[out] WORK
          WORK is REAL array, dimension (LWORK), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced. \n

 * @return Returns value of norm.
 * */
    template <typename T, typename Ta>
    Ta lanhf(char *norm, char *transr, char *uplo, integer *n, T *a, Ta *work)
    {
        return lanhf(norm, transr, uplo, n, a, work);
    }
    /** @}*/ // end of lanhf

    /** @defgroup lansf lansf
     * @ingroup  MN
     * @{
     */
    /*! @brief LANSF returns the value of the one norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a real symmetric matrix A in RFP format
 * @details
 * \b Purpose:
    \verbatim
    LANSF returns the value of the one norm, or the Frobenius norm, or
    the infinity norm, or the element of largest absolute value of a
    real symmetric matrix A in RFP format.

    \  return SLANSF

       SLANSF = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a  matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in SLANSF as described
          above. \n
 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          Specifies whether the RFP format of A is normal or
          transposed format. \n
          = 'N':  RFP format is Normal; \n
          = 'T':  RFP format is Transpose. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          On entry, UPLO specifies whether the RFP matrix A came from
          an upper or lower triangular matrix as follows: \n
          = 'U': RFP A came from an upper triangular matrix; \n
          = 'L': RFP A came from a lower triangular matrix. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A. N >= 0. When N = 0, SLANSF is
          set to zero. \n
 * @param[in] A
          A is REAL array, dimension (N*(N+1)/2); \n
          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')
          part of the symmetric matrix A stored in RFP format. See the
          "Notes" below for more details. \n
          Unchanged on exit. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced.  \n

 * @return Returns the value of norm.
 * */
    template <typename T>
    T lansf(char *norm, char *transr, char *uplo, integer *n, T *a, T *work)
    {
        return lansf(norm, transr, uplo, n, a, work);
    }
    /** @}*/ // end of lansf

    /** @defgroup lanhp lanhp
     * @ingroup  MN
     * @{
     */
    /*! @brief LANHP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a complex Hermitian matrix supplied in packed form
 * @details
 * \b Purpose:
    \verbatim
    LANHP returns the value of the one norm, or the Frobenius norm, or
    the  infinity norm, or the  element of  largest absolute value  of a
    complex hermitian matrix A, supplied in packed form.

       CLANHP = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in CLANHP as described
          above. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          hermitian matrix A is supplied. \n
          = 'U':  Upper triangular part of A is supplied \n
          = 'L':  Lower triangular part of A is supplied \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, CLANHP is
          set to zero. \n
 * @param[in] AP
          AP is COMPLEX array, dimension (N*(N+1)/2) \n
          The upper or lower triangle of the hermitian matrix A, packed
          columnwise in a linear array.  The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
          Note that the  imaginary parts of the diagonal elements need
          not be set and are assumed to be zero. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced. \n

 * @return Returns the value of norm.
 * */
    template <typename T, typename Ta>
    Ta lanhp(char *norm, char *uplo, integer *n, T *ap, Ta *work)
    {
        return lanhp(norm, uplo, n, ap, work);
    }
    /** @}*/ // end of lanhp

    /** @defgroup lansp lansp
     * @ingroup  MN
     * @{
     */
    /*! @brief LANSP returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a symmetric matrix supplied in packed form
 * @details
 * \b Purpose:
    \verbatim
    LANSP returns the value of the one norm, or the Frobenius norm, or
    the infinity norm, or the element of largest absolute value of a
    real symmetric matrix A, supplied in packed form.

    \  return SLANSP

       SLANSP = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in SLANSP as described
          above. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          symmetric matrix A is supplied. \n
          = 'U':  Upper triangular part of A is supplied \n
          = 'L':  Lower triangular part of A is supplied \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANSP is
          set to zero. \n
 * @param[in] AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangle of the symmetric matrix A, packed
          columnwise in a linear array.  The j-th column of A is stored
          in the array AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced. \n

 * @return Returns the value of norm.
 * */
    template <typename T>
    T lansp(char *norm, char *uplo, integer *n, T *ap, T *work)
    {
        return lansp(norm, uplo, n, ap, work);
    }
    template <typename T, typename Ta>
    Ta lansp(char *norm, char *uplo, integer *n, T *ap, Ta *work)
    {
        return lansp(norm, uplo, n, ap, work);
    }
    /** @}*/ // end of lansp

    /** @defgroup lanhb lanhb
     * @ingroup  MN
     * @{
     */
    /*! @brief LANHB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a Hermitian band matrix
 * @details
 * \b Purpose:
    \verbatim
    LANHB returns the value of the one norm, or the Frobenius norm, or
    the  infinity norm, or the element of  largest absolute value  of an
    n by n hermitian band matrix A, with k super-diagonals.

       LANHB = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in CLANHB as described
          above. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          band matrix A is supplied. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, CLANHB is
          set to zero. \n
 * @param[in] K
          K is INTEGER \n
          The number of super-diagonals or sub-diagonals of the
          band matrix A.  K >= 0. \n
 * @param[in] AB
          AB is COMPLEX array, dimension (LDAB,N) \n
          The upper or lower triangle of the hermitian band matrix A,
          stored in the first K+1 rows of AB.  The j-th column of A is
          stored in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for fla_max(1,j-k)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k). \n
          Note that the imaginary parts of the diagonal elements need
          not be set and are assumed to be zero. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= K+1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced. \n

 * @return Returns the value of norm.
 * */
    template <typename T, typename Ta>
    Ta lanhb(char *norm, char *uplo, integer *n, integer *k, T *ab, integer *ldab, Ta *work)
    {
        return lanhb(norm, uplo, n, k, ab, ldab, work);
    }
    /** @}*/ // end of lanhb

    /** @defgroup lansb lansb
     * @ingroup  MN
     * @{
     */
    /*! @brief LANSB returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a symmetric band matrix
 * @details
 * \b Purpose:
    \verbatim
    LANSB returns the value of the one norm, or the Frobenius norm, or
    the infinity norm, or the element of largest absolute value of an
    n by n symmetric band matrix A, with k super-diagonals.

    \  return SLANSB

       SLANSB = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in SLANSB as described
          above. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the upper or lower triangular part of the
          band matrix A is supplied. \n
          = 'U':  Upper triangular part is supplied \n
          = 'L':  Lower triangular part is supplied \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANSB is
          set to zero. \n
 * @param[in] K
          K is INTEGER \n
          The number of super-diagonals or sub-diagonals of the
          band matrix A.  K >= 0. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB,N) \n
          The upper or lower triangle of the symmetric band matrix A,
          stored in the first K+1 rows of AB.  The j-th column of A is
          stored in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for fla_max(1,j-k)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k). \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= K+1. \n
 * @param[out] WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
          WORK is not referenced.  \n

 * @return Returns value of norm.
 * */
    template <typename T, typename Ta>
    Ta lansb(char *norm, char *uplo, integer *n, integer *k, T *ab, integer *ldab, Ta *work)
    {
        return lansb(norm, uplo, n, k, ab, ldab, work);
    }
    /** @}*/ // end of lansb

    /** @defgroup lanht lanht
     * @ingroup  MN
     * @{
     */
    /*! @brief LANHT returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a complex Hermitian tridiagonal matrix
 * @details
 * \b Purpose:
    \verbatim
    LANHT returns the value of the one norm, or the Frobenius norm, or
    the infinity norm, or the element of largest absolute value of a
    complex Hermitian tridiagonal matrix A.

    \  return CLANHT


       CLANHT = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in CLANHT as described
          above. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, CLANHT is
          set to zero. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The diagonal elements of A. \n
 * @param[in] E
          E is COMPLEX array, dimension (N-1) \n
          The (n-1) sub-diagonal or super-diagonal elements of A.  \n

 * @return Returns the value of norm.
 * */
    template <typename T, typename Ta>
    Ta lanht(char *norm, integer *n, Ta *d, T *e)
    {
        return lanht(norm, n, d, e);
    }
    /** @}*/ // end of lanht

    /** @defgroup lanst lanst
     * @ingroup  MN
     * @{
     */
    /*! @brief LANST returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
     or the element of largest absolute value of a real symmetric tridiagonal matrix
 * @details
 * \b Purpose:
    \verbatim
    LANST returns the value of the one norm, or the Frobenius norm, or
    the infinity norm, or the element of largest absolute value  of a
    real symmetric tridiagonal matrix A.

    \  return SLANST

       SLANST = (fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                (
                (norm1(A),        NORM = '1', 'O' or 'o'
                (
                (normI(A),        NORM = 'I' or 'i'
                (
                (normF(A),        NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be   returned in SLANST as described
          above. \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANST is
          set to zero. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The diagonal elements of A. \n
 * @param[in] E
          E is REAL array, dimension (N-1) \n
          The (n-1) sub-diagonal or super-diagonal elements of A. \n

 * @return Returns the value of norm.
 * */
    template <typename T>
    T lanst(char *norm, integer *n, T *d, T *e)
    {
        return lanst(norm, n, d, e);
    }
    /** @}*/ // end of lanst

    /** @defgroup lantr lantr
     * @ingroup  MN
     * @{
     */
    /*! @brief LANTR returns the value of the 1-norm, or the Frobenius norm, or the infinity norm,
 \n or the element of largest absolute value of a trapezoidal or triangular matrix.
 * @details
 * \b Purpose:
    \verbatim
     LANTR  returns the value of the one norm,  or the Frobenius norm, or
     the  infinity norm,  or the  element of  largest absolute value  of a
     trapezoidal or triangular matrix A.
     \return LANTR
        LANTR = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
                 (
                 ( norm1(A),         NORM = '1', 'O' or 'o'
                 (
                 ( normI(A),         NORM = 'I' or 'i'
                 (
                 ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
     where  norm1  denotes the  one norm of a matrix (maximum column sum),
     normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
     normF  denotes the  Frobenius norm of a matrix (square root of sum of
     squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in] NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be returned in SLANTR as described
          above. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the matrix A is upper or lower trapezoidal. \n
          = 'U':  Upper trapezoidal \n
          = 'L':  Lower trapezoidal \n
          Note that A is triangular instead of trapezoidal if M = N. \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          Specifies whether or not the matrix A has unit diagonal. \n
          = 'N':  Non-unit diagonal \n
          = 'U':  Unit diagonal \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0, and if
          UPLO = 'U', M <= N.  When M = 0, SLANTR is set to zero. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0, and if
          UPLO = 'L', N <= M.  When N = 0, SLANTR is set to zero. \n
 * @param[in] A
          A is REAL array, dimension (LDA,N) \n
          The trapezoidal matrix A (A is triangular if M = N). \n
          If UPLO = 'U', the leading m by n upper trapezoidal part of
          the array A contains the upper trapezoidal matrix, and the
          strictly lower triangular part of A is not referenced. \n
          If UPLO = 'L', the leading m by n lower trapezoidal part of
          the array A contains the lower trapezoidal matrix, and the
          strictly upper triangular part of A is not referenced.  Note
          that when DIAG = 'U', the diagonal elements of A are not
          referenced and are assumed to be one. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A.  LDA >= fla_max(M,1). \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
          referenced. \n

 * @return Returns value of the norm.
 * */
    template <typename T>
    T lantr(char *norm, char *uplo, char *diag, integer *m, integer *n, T *a, integer *lda, T *work)
    {
        return lantr(norm, uplo, diag, m, n, a, lda, work);
    }
    template <typename T, typename Ta>
    Ta lantr(char *norm, char *uplo, char *diag, integer *m, integer *n, T *a, integer *lda,
             Ta *work)
    {
        return lantr(norm, uplo, diag, m, n, a, lda, work);
    }
    /** @}*/ // end of lantr

    /** @defgroup lantp lantp
     * @ingroup  MN
     * @{
     */
    /*! @brief LANTP returns the value of the one norm, or the Frobenius norm, or
      the infinity norm, or the element of largest absolute value  of a
      triangular matrix A, supplied in packed form.

 * @details
 * \b Purpose:
    \verbatim
    SLANTP = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
             (
             ( norm1(A),         NORM = '1', 'O' or 'o'
             (
             ( normI(A),         NORM = 'I' or 'i'
             (
             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in]	NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be returned in SLANTP as described
          above. \n
 * @param[in]	UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the matrix A is upper or lower triangular. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in]	DIAG
          DIAG is CHARACTER*1 \n
          Specifies whether or not the matrix A is unit triangular. \n
          = 'N':  Non-unit triangular \n
          = 'U':  Unit triangular \n
 * @param[in]	N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANTP is
          set to zero. \n
 * @param[in]	AP
          AP is REAL array, dimension (N*(N+1)/2) \n
          The upper or lower triangular matrix A, packed columnwise in
          a linear array.  The j-th column of A is stored in the array
          AP as follows: \n
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; \n
          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. \n
          Note that when DIAG = 'U', the elements of the array AP
          corresponding to the diagonal elements of the matrix A are
          not referenced, but are assumed to be one. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
          referenced. \n

 * @return Returns value of the norm.
 * */
    template <typename T>
    T lantp(char *norm, char *uplo, char *diag, integer *n, T *ap, T *work)
    {
        return lantp(norm, uplo, diag, n, ap, work);
    }
    template <typename T, typename Ta>
    Ta lantp(char *norm, char *uplo, char *diag, integer *n, T *ap, Ta *work)
    {
        return lantp(norm, uplo, diag, n, ap, work);
    }
    /** @}*/ // end of lantp

    /** @defgroup lantb lantb
     * @ingroup  MN
     * @{
     */
    /*! @brief LANTB returns the value of the one norm, or the Frobenius norm, or
     the infinity norm, or the element of largest absolute value of an
     n by n triangular band matrix A,  with ( k + 1 ) diagonals.

 * @details
 * \b Purpose:
    \verbatim
    SLANTB = ( fla_max(abs(A(i,j))), NORM = 'M' or 'm'
             (
             ( norm1(A),         NORM = '1', 'O' or 'o'
             (
             ( normI(A),         NORM = 'I' or 'i'
             (
             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

    where  norm1  denotes the  one norm of a matrix (maximum column sum),
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    normF  denotes the  Frobenius norm of a matrix (square root of sum of
    squares).  Note that  fla_max(abs(A(i,j)))  is not a consistent matrix norm.
    \endverbatim

 * @param[in]	NORM
          NORM is CHARACTER*1 \n
          Specifies the value to be returned in SLANTB as described
          above. \n
 * @param[in]	UPLO
          UPLO is CHARACTER*1 \n
          Specifies whether the matrix A is upper or lower triangular. \n
          = 'U':  Upper triangular \n
          = 'L':  Lower triangular \n
 * @param[in]	DIAG
          DIAG is CHARACTER*1 \n
          Specifies whether or not the matrix A is unit triangular. \n
          = 'N':  Non-unit triangular \n
          = 'U':  Unit triangular \n
 * @param[in]	N
          N is INTEGER \n
          The order of the matrix A.  N >= 0.  When N = 0, SLANTB is
          set to zero. \n
 * @param[in]	K
          K is INTEGER \n
          The number of super-diagonals of the matrix A if UPLO = 'U',
          or the number of sub-diagonals of the matrix A if UPLO = 'L'.
          K >= 0. \n
 * @param[in]	AB
          AB is REAL array, dimension (LDAB,N) \n
          The upper or lower triangular band matrix A, stored in the
          first k+1 rows of AB.  The j-th column of A is stored
          in the j-th column of the array AB as follows: \n
          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for fla_max(1,j-k)<=i<=j; \n
          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k). \n
          Note that when DIAG = 'U', the elements of the array AB
          corresponding to the diagonal elements of the matrix A are
          not referenced, but are assumed to be one. \n
 * @param[in]	LDAB
          LDAB is INTEGER \n
          The leading dimension of the array AB.  LDAB >= K+1. \n
 * @param[out]	WORK
          WORK is REAL array, dimension (MAX(1,LWORK)), \n
          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
          referenced. \n

 * @return Returns value of the norm.
 * */
    template <typename T>
    T lantb(char *norm, char *uplo, char *diag, integer *n, integer *k, T *ab, integer *ldab,
            T *work)
    {
        return lantb(norm, uplo, diag, n, k, ab, ldab, work);
    }
    template <typename T, typename Ta>
    Ta lantb(char *norm, char *uplo, char *diag, integer *n, integer *k, T *ab, integer *ldab,
             Ta *work)
    {
        return lantb(norm, uplo, diag, n, k, ab, ldab, work);
    }
    /** @}*/ // end of lantb
    /** @}*/ // end of Matrix norm

    /** @defgroup Scalar Scalar operations
     * @ingroup BLAS
     * @{
     */

    /** @defgroup isnan isnan
     * @ingroup Scalar
     * @{
     */
    /*! @brief ISNAN tests input for NaN

 * @details
 * \b Purpose:
    \verbatim
    ISNAN   returns .TRUE. if its argument is NaN, and .FALSE.
    otherwise.  To be replaced by the Fortran 2003 intrinsic in the
    future.
    \endverbatim

 * @param[in] SIN
          SIN is REAL \n
          Input to test for NaN. \n
 * @return LOGICAL Boolean. Return TRUE if argument is NAN.
 * */
    template <typename T>
    logical isnan(T *sin)
    {
        return isnan(sin);
    }
    /** @}*/ // end of isNAN

    /** @defgroup laisnan laisnan
     * @ingroup Scalar
     * @{
     */
    /*! @brief LAISNAN tests input for NaN by comparing two arguments for inequality

 * @details
 * \b Purpose:
    \verbatim
    This routine is not for general use.  It exists solely to avoid
    over-optimization in ISNAN.

    LAISNAN checks for NaNs by comparing its two arguments for
    inequality.  NaN is the only floating-point value where NaN != NaN
      returns .TRUE.  To check for NaNs, pass the same variable as both
    arguments.

    A compiler must assume that the two arguments are
    not the same variable, and the test will not be optimized away.
    Interprocedural or whole-program optimization may delete this
    test.  The ISNAN functions will be replaced by the correct
    Fortran 03 intrinsic once the intrinsic is widely available.
    \endverbatim

 * @param[in] SIN1
          SIN1 is REAL \n
 * @param[in] SIN2
          SIN2 is REAL \n
          Two numbers to compare for inequality. \n

 * @return LOGICAL Boolean. Return TRUE if argument is NAN.
 * */
    template <typename T>
    logical laisnan(real *sin1, real *sin2)
    {
        return laisnan(sin1, sin2);
    }

    /** @}*/ // end of laisnan

    /** @defgroup ladiv ladiv
     * @ingroup Scalar
     * @{
     */
    /*! @brief LADIV performs complex division in real arithmetic, avoiding unnecessary overflow

 * @details
 * \b Purpose:
    \verbatim
    LADIV performs complex division in  real arithmetic

                          a + i*b
               p + i*q = ---------
                          c + i*d

    The algorithm is due to Michael Baudin and Robert L. Smith
    and can be found in the paper
    "A Robust Complex Division in Scilab"
    \endverbatim

 * @param[in] A
          A is REAL \n
 * @param[in] B
          B is REAL \n
 * @param[in] C
          C is REAL \n
 * @param[in] D
          D is REAL \n
          The scalars a, b, c, and d in the above expression. \n
 * @param[out] P
          P is REAL \n
 * @param[out] Q
          Q is REAL \n
          The scalars p and q in the above expression.  \n

 * */
    template <typename T>
    void ladiv(T *a, T *b, T *c, T *d, T *p, T *q)
    {
        ladiv(a, b, c, d, p, q);
    }
    /** @}*/ // end of ladiv

    /** @defgroup lapy2 lapy2
     * @ingroup Scalar
     * @{
     */
    /*! @brief LAPY2 returns sqrt(x2+y2)

 * @details
 * \b Purpose:
    \verbatim
    LAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
    overflow.
    \endverbatim

 * @param[in] X
          X is REAL \n
 * @param[in] Y
          Y is REAL \n
          X and Y specify the values x and y. \n

* @return T Template based return value of the function.
 * */
    template <typename T>
    T lapy2(T *x, T *y)
    {
        return lapy2(x, y);
    }
    /** @}*/ // end of lapy2

    /** @defgroup lapy3 lapy3
     * @ingroup Scalar
     * @{
     */
    /*! @brief LAPY3 returns sqrt(x2+y2+z2)

 * @details
 * \b Purpose:
    \verbatim
     LAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
     unnecessary overflow.
    \endverbatim

 * @param[in] X
          X is REAL \n
 * @param[in] Y
          Y is REAL \n
 * @param[in] Z
          Z is REAL \n
          X, Y and Z specify the values x, y and z. \n

 * @return T Template based return value of the function.
 * */
    template <typename T>
    T lapy3(T *x, T *y, T *z)
    {
        return lapy3(x, y, z);
    }
    /** @}*/ // end of lapy3

    /** @}*/ // end of scalar operations

    /** @defgroup L1 Level 1 BLAS-like vector ops
     * @ingroup BLAS
     * @{
     */

    /** @defgroup lacgv lacgv
     * @ingroup L1
     * @{
     */
    /*! @brief LACGV conjugates a complex vector

 * @details
 * \b Purpose:
    \verbatim
    CLACGV conjugates a complex vector of length N.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The length of the vector X.  N >= 0. \n
 * @param[in,out] X
          X is COMPLEX array, dimension
                         (1+(N-1)*abs(INCX)) \n
          On entry, the vector of length N to be conjugated.
          On exit, X is overwritten with conjg(X). \n
 * @param[in] INCX
          INCX is INTEGER \n
          The spacing between successive elements of X. \n

 * */
    template <typename T>
    void lacgv(integer *n, T *x, integer *incx)
    {
        lacgv(n, x, incx);
    }
    /** @}*/ // end of lacgv

    /** @defgroup lasrt lasrt
     * @ingroup L1
     * @{
     */
    /*! @brief LASRT sorts numbers in increasing or decreasing order

 * @details
 * \b Purpose:
    \verbatim
     Sort the numbers in D in increasing order (if ID = 'I') or
     in decreasing order (if ID = 'D').

     Use Quick Sort, reverting to Insertion sort on arrays of
     size <= 20. Dimension of STACK limits N to about 2**32.
    \endverbatim

 * @param[in] ID
          ID is CHARACTER*1 \n
          = 'I': sort D in increasing order; \n
          = 'D': sort D in decreasing order. \n
 * @param[in] N
          N is INTEGER \n
          The length of the array D. \n
 * @param[in,out] D
          D is REAL array, dimension (N) \n
          On entry, the array to be sorted. \n
          On exit, D has been sorted into increasing order \n
          (D(1) <= ... <= D(N)) or into decreasing order \n
          (D(1) >= ... >= D(N)), depending on ID. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          = 0:  successful exit \n
          < 0:  if INFO = -i, the i-th argument had an illegal value \n

 * */
    template <typename T>
    void lasrt(char *id, integer *n, T *d, integer *info)
    {
        lasrt(id, n, d, info);
    }
    /** @}*/ // end of lasrt

    /** @defgroup lassq lassq
     * @ingroup L1
     * @{
     */
    /*! @brief LASSQ updates a sum of squares represented in scaled form

 * @details
 * \b Purpose:
    \verbatim
     LASSQ  returns the values  scl  and  smsq  such that

        ( scl**2)*smsq = x( 1)**2 +...+ x( n)**2 + ( scale**2)*sumsq,

     where  x( i) = X( 1 + ( i - 1)*INCX). The value of  sumsq  is
     assumed to be non-negative and  scl  returns the value

        scl = fla_max( scale, abs( x( i))).

     scale and sumsq must be supplied in SCALE and SUMSQ and
     scl and smsq are overwritten on SCALE and SUMSQ respectively.

     The routine makes only one pass through the vector x.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of elements to be used from the vector X. \n
 * @param[in] X
          X is REAL array, dimension (1+(N-1)*INCX) \n
          The vector for which a scaled sum of squares is computed. \n
             x( i)  = X( 1 + ( i - 1)*INCX), 1 <= i <= n. \n
 * @param[in] INCX
          INCX is INTEGER \n
          The increment between successive values of the vector X.
          INCX > 0. \n
 * @param[in,out] SCALE
          SCALE is REAL \n
          On entry, the value  scale  in the equation above. \n
          On exit, SCALE is overwritten with  scl , the scaling factor
          for the sum of squares. \n
 * @param[in,out] SUMSQ
          SUMSQ is REAL \n
          On entry, the value  sumsq  in the equation above. \n
          On exit, SUMSQ is overwritten with  smsq , the basic sum of
          squares from which  scl  has been factored out. \n

 * */
    template <typename T>
    void lassq(integer *n, T *x, integer *incx, T *scale, T *sumsq)
    {
        lassq(n, x, incx, scale, sumsq);
    }
    template <typename T, typename Ta>
    void lassq(integer *n, T *x, integer *incx, Ta *scale, Ta *sumsq)
    {
        lassq(n, x, incx, scale, sumsq);
    }
    /** @}*/ // end of lassq

    /** @defgroup rot rot
     * @ingroup L1
     * @{
     */
     /*! @brief ROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors

      * @details
      * \b Purpose:
          \verbatim 
           ROT   applies a plane rotation, where the cos (C) is real and the
           sin (S) is complex, and the vectors CX and CY are complex.
          \endverbatim  

      * @param[in] N
               N is INTEGER \n
               The number of elements in the vectors CX and CY. \n
      * @param[in,out] CX
               CX is COMPLEX array, dimension (N) \n
               On input, the vector X.
               On output, CX is overwritten with C*X + S*Y. \n
      * @param[in] INCX
               INCX is INTEGER \n
               The increment between successive values of CY.  INCX <> 0. \n
      * @param[in,out] CY
               CY is COMPLEX array, dimension (N) \n
               On input, the vector Y.
               On output, CY is overwritten with -CONJG(S)*X + C*Y. \n
      * @param[in] INCY
               INCY is INTEGER \n
               The increment between successive values of CY.  INCX <> 0. \n
      * @param[in] C
               C is REAL \n
      * @param[in] S
               S is COMPLEX \n
               C and S define a rotation \n
                  [  C          S  ] \n
                  [ -conjg(S)   C  ] \n
               where C*C + S*CONJG(S) = 1.0. \n

      *  * */
    template< typename T, typename Ta >
    void rot(integer *n, T *cx, integer *incx, T * cy, integer *incy, Ta *c, T *s)
    {
        rot(n, cx, incx, cy, incy, c, s);
    }
    /** @}*/ // end of rot

    /** @defgroup rscl rscl
     * @ingroup L1
     * @{
     */
    /*! @brief RSCL multiplies a vector by the reciprocal of a real scalar

 * @details
 * \b Purpose:
    \verbatim
    RSCL multiplies an n-element real vector x by the real scalar 1/a.
    This is done without overflow or underflow as long as
    the final result x/a does not overflow or underflow.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The number of components of the vector x. \n
 * @param[in] SA
          SA is REAL \n
          The scalar a which is used to divide each component of x.
          SA must be >= 0, or the subroutine will divide by zero. \n
 * @param[in,out] SX
          SX is REAL array, dimension
                         (1+(N-1)*abs(INCX)) \n
          The n-element vector x. \n
 * @param[in] INCX
          INCX is INTEGER \n
          The increment between successive values of the vector SX. \n
          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),    1< i<= n \n

 *  * */
    template <typename T>
    void rscl(integer *n, T *sa, T *sx, integer *incx)
    {
        rscl(n, sa, sx, incx);
    }
    template <typename T, typename Ta>
    void rscl(integer *n, Ta *sa, T *sx, integer *incx)
    {
        rscl(n, sa, sx, incx);
    }
    /** @}*/ // end of rscl

    /** @}*/ // end of l1

    /** @defgroup L2 Level 2 BLAS-like matrix-vector ops
     * @ingroup BLAS
     * @{
     */

    /** @defgroup lascl lascl
     * @ingroup L2
     * @{
     */
    /*! @brief LASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom

 * @details
 * \b Purpose:
    \verbatim
     LASCL multiplies the M by N real matrix A by the real scalar
     CTO/CFROM.  This is done without over/underflow as long as the final
     result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     A may be full, upper triangular, lower triangular, upper Hessenberg,
     or banded.
    \endverbatim

 * @param[in] TYPE
          TYPE is CHARACTER*1 \n
          TYPE indices the storage type of the input matrix. \n
          = 'G':  A is a full matrix. \n
          = 'L':  A is a lower triangular matrix. \n
          = 'U':  A is an upper triangular matrix. \n
          = 'H':  A is an upper Hessenberg matrix. \n
          = 'B':  A is a symmetric band matrix with lower bandwidth KL
                  and upper bandwidth KU and with the only the lower
                  half stored. \n
          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
                  and upper bandwidth KU and with the only the upper
                  half stored. \n
          = 'Z':  A is a band matrix with lower bandwidth KL and upper
                  bandwidth KU. See SGBTRF for storage details. \n
 * @param[in] KL
          KL is INTEGER \n
          The lower bandwidth of A.  Referenced only if TYPE = 'B',
          'Q' or 'Z'. \n
 * @param[in] KU
          KU is INTEGER \n
          The upper bandwidth of A.  Referenced only if TYPE = 'B',
          'Q' or 'Z'. \n
 * @param[in] CFROM
          CFROM is REAL \n
 * @param[in] CTO
          CTO is REAL \n
          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
          without over/underflow if the final result CTO*A(I,J)/CFROM
          can be represented without over/underflow.  CFROM must be
          nonzero. \n
 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A.  M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of the matrix A.  N >= 0. \n
 * @param[in,out] A
          A is REAL array, dimension (LDA,N) \n
          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
          storage type. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. \n
          If TYPE = 'G', 'L', 'U', 'H', LDA >= fla_max(1,M); \n
             TYPE = 'B', LDA >= KL+1; \n
             TYPE = 'Q', LDA >= KU+1; \n
             TYPE = 'Z', LDA >= 2*KL+KU+1. \n
 * @param[out]	INFO
          INFO is INTEGER \n
          0  - successful exit \n
          <0 - if INFO = -i, the i-th argument had an illegal value. \n

 * */
    template <typename T>
    void lascl(char *type, integer *kl, integer *ku, T *cfrom, T *cto, integer *m, integer *n, T *a,
               integer *lda, integer *info)
    {
        lascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }
    template <typename T, typename Ta>
    void lascl(char *type, integer *kl, integer *ku, Ta *cfrom, Ta *cto, integer *m, integer *n,
               T *a, integer *lda, integer *info)
    {
        lascl(type, kl, ku, cfrom, cto, m, n, a, lda, info);
    }

    /** @}*/ // end of lascl

    /** @defgroup lageamv la_geamv: matrix-vector multiply |A| * |x|, general
     * @ingroup L2
     * @{
     */
    /*! @brief LA_GEAMV computes a matrix-vector product using a general matrix to calculate error
 bounds

 * @details
 * \b Purpose:
    \verbatim
    LA_GEAMV  performs one of the matrix-vector operations

            y := alpha*abs(A)*abs(x) + beta*abs(y),
       or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),

    where alpha and beta are scalars, x and y are vectors and A is an
    m by n matrix.

    This function is primarily used in calculating error bounds.
    To protect against underflow during evaluation, components in
    the resulting vector are perturbed away from zero by (N+1)
    times the underflow threshold.  To prevent unnecessarily large
    errors for block-structure embedded in general matrices,
    "symbolically" zero components are not perturbed.  A zero
    entry is considered "symbolic" if all multiplications involved
    in computing that entry have at least one zero multiplicand.

    Level 2 Blas routine.
    \endverbatim

 * @param[in] TRANS
          TRANS is INTEGER \n
          On entry, TRANS specifies the operation to be performed as
          follows: \n
 \n
            BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y) \n
            BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y) \n
            BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y) \n
 \n
          Unchanged on exit. \n
 * @param[in] M
          M is INTEGER \n
          On entry, M specifies the number of rows of the matrix A.
          M must be at least zero. \n
          Unchanged on exit. \n
 * @param[in] N
          N is INTEGER \n
          On entry, N specifies the number of columns of the matrix A.
          N must be at least zero. \n
          Unchanged on exit. \n
 * @param[in] ALPHA
          ALPHA is REAL \n
          On entry, ALPHA specifies the scalar alpha.
          Unchanged on exit. \n
 * @param[in] A
          A is REAL array, dimension (LDA, n) \n
          Before entry, the leading m by n part of the array A must
          contain the matrix of coefficients. \n
          Unchanged on exit. \n
 * @param[in] LDA
          LDA is INTEGER \n
          On entry, LDA specifies the first dimension of A as declared
          in the calling (sub) program. LDA must be at least
          fla_max(1, m). \n
          Unchanged on exit. \n
 * @param[in] X
          X is REAL array, dimension
          (1 + (n - 1)*abs(INCX)) when TRANS = 'N' or 'n' \n
          and at least
          (1 + (m - 1)*abs(INCX)) otherwise. \n
          Before entry, the incremented array X must contain the
          vector x. \n
          Unchanged on exit. \n
 * @param[in] INCX
          INCX is INTEGER \n
          On entry, INCX specifies the increment for the elements of
          X. INCX must not be zero. \n
          Unchanged on exit. \n
 * @param[in] BETA
          BETA is REAL \n
          On entry, BETA specifies the scalar beta. When BETA is
          supplied as zero then Y need not be set on input. \n
          Unchanged on exit. \n
 * @param[in,out] Y
          Y is REAL array,
          dimension at least
          (1 + (m - 1)*abs(INCY)) when TRANS = 'N' or 'n' \n
          and at least
          (1 + (n - 1)*abs(INCY)) otherwise. \n
          Before entry with BETA non-zero, the incremented array Y
          must contain the vector y. On exit, Y is overwritten by the
          updated vector y. \n
 * @param[in] INCY
          INCY is INTEGER \n
          On entry, INCY specifies the increment for the elements of
          Y. INCY must not be zero. \n
          Unchanged on exit. \n

 * */
    template <typename T>
    void la_geamv(integer *trans, integer *m, integer *n, T *alpha, T *a, integer *lda, T *x,
                  integer *incx, T *beta, T *y, integer *incy)
    {
        la_geamv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
    template <typename T, typename Ta>
    void la_geamv(integer *trans, integer *m, integer *n, Ta *alpha, T *a, integer *lda, T *x,
                  integer *incx, Ta *beta, Ta *y, integer *incy)
    {
        la_geamv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
    }
    /** @}*/ // end of lageamc

    /** @defgroup lageamb la_gbamv: matrix-vector multiply |A| * |x|, general banded
     * @ingroup L2
     * @{
     */
    /*! @brief LA_GBAMV performs a matrix-vector operation to calculate error bounds

 * @details
 * \b Purpose:
    \verbatim
    LA_GBAMV  performs one of the matrix-vector operations

            y := alpha*abs(A)*abs(x) + beta*abs(y),
       or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),

    where alpha and beta are scalars, x and y are vectors and A is an
    m by n matrix.

    This function is primarily used in calculating error bounds.
    To protect against underflow during evaluation, components in
    the resulting vector are perturbed away from zero by (N+1)
    times the underflow threshold.  To prevent unnecessarily large
    errors for block-structure embedded in general matrices,
    "symbolically" zero components are not perturbed.  A zero
    entry is considered "symbolic" if all multiplications involved
    in computing that entry have at least one zero multiplicand.
    Level 2 Blas routine.
    \endverbatim

 * @param[in] TRANS
          TRANS is INTEGER \n
          On entry, TRANS specifies the operation to be performed as
          follows: \n
 \n
            BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y) \n
            BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y) \n
            BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y) \n
 \n
          Unchanged on exit. \n
 * @param[in] M
          M is INTEGER \n
          On entry, M specifies the number of rows of the matrix A.
          M must be at least zero. \n
          Unchanged on exit. \n
 * @param[in] N
          N is INTEGER \n
          On entry, N specifies the number of columns of the matrix A.
          N must be at least zero. \n
          Unchanged on exit. \n
 * @param[in] KL
          KL is INTEGER \n
          The number of subdiagonals within the band of A.  KL >= 0. \n
 * @param[in] KU
          KU is INTEGER \n
          The number of superdiagonals within the band of A.  KU >= 0. \n
 * @param[in] ALPHA
          ALPHA is REAL \n
          On entry, ALPHA specifies the scalar alpha.
          Unchanged on exit. \n
 * @param[in] AB
          AB is REAL array, dimension (LDAB, n) \n
          Before entry, the leading m by n part of the array AB must
          contain the matrix of coefficients. \n
          Unchanged on exit. \n
 * @param[in] LDAB
          LDAB is INTEGER \n
          On entry, LDA specifies the first dimension of AB as declared
          in the calling (sub) program. LDAB must be at least
          fla_max(1, m). \n
          Unchanged on exit. \n
 * @param[in] X
          X is REAL array, dimension
          (1 + (n - 1)*abs(INCX)) when TRANS = 'N' or 'n' \n
          and at least
          (1 + (m - 1)*abs(INCX)) otherwise.
          Before entry, the incremented array X must contain the
          vector x. \n
          Unchanged on exit. \n
 * @param[in] INCX
          INCX is INTEGER \n
          On entry, INCX specifies the increment for the elements of
          X. INCX must not be zero. \n
          Unchanged on exit. \n
 * @param[in] BETA
          BETA is REAL \n
          On entry, BETA specifies the scalar beta. When BETA is
          supplied as zero then Y need not be set on input. \n
          Unchanged on exit. \n
 * @param[in,out] Y
          Y is REAL array, dimension
          (1 + (m - 1)*abs(INCY)) when TRANS = 'N' or 'n' \n
          and at least
          (1 + (n - 1)*abs(INCY)) otherwise. \n
          Before entry with BETA non-zero, the incremented array Y
          must contain the vector y. On exit, Y is overwritten by the
          updated vector y. \n
 * @param[in] INCY
          INCY is INTEGER \n
          On entry, INCY specifies the increment for the elements of
          Y. INCY must not be zero. \n
          Unchanged on exit. \n

 * */
    template <typename T>
    void la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, T *alpha, T *ab,
                  integer *ldab, T *x, integer *incx, T *beta, T *y, integer *incy)
    {
        la_gbamv(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
    }
    template <typename T, typename Ta>
    void la_gbamv(integer *trans, integer *m, integer *n, integer *kl, integer *ku, Ta *alpha,
                  T *ab, integer *ldab, T *x, integer *incx, Ta *beta, Ta *y, integer *incy)
    {
        la_gbamv(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy);
    }
    /** @}*/ // end of lageamb

    /** @defgroup la_heamv la_heamv: matrix-vector multiply |A| * |x|, Hermitian/symmetric
     * @ingroup L2
     * @{
     */
    /*! @brief LA_HEAMV computes a matrix-vector product using a Hermitian indefinite matrix to
 calculate error bounds

 * @details
 * \b Purpose:
    \verbatim
    LA_HEAMV  performs the matrix-vector operation

            y := alpha*abs(A)*abs(x) + beta*abs(y),

    where alpha and beta are scalars, x and y are vectors and A is an
    n by n symmetric matrix.

    This function is primarily used in calculating error bounds.
    To protect against underflow during evaluation, components in
    the resulting vector are perturbed away from zero by (N+1)
    times the underflow threshold.  To prevent unnecessarily large
    errors for block-structure embedded in general matrices,
    "symbolically" zero components are not perturbed.  A zero
    entry is considered "symbolic" if all multiplications involved
    in computing that entry have at least one zero multiplicand.
    \endverbatim

 * @param[in] UPLO
          UPLO is INTEGER \n
          On entry, UPLO specifies whether the upper or lower
          triangular part of the array A is to be referenced as
          follows: \n
 \n
              UPLO = BLAS_UPPER   Only the upper triangular part of A
                                  is to be referenced.
 \n
              UPLO = BLAS_LOWER   Only the lower triangular part of A
                                  is to be referenced.
 \n
          Unchanged on exit. \n
 * @param[in] N
          N is INTEGER \n
          On entry, N specifies the number of columns of the matrix A.
          N must be at least zero. \n
          Unchanged on exit. \n
 * @param[in] ALPHA
          ALPHA is REAL . \n
          On entry, ALPHA specifies the scalar alpha.
          Unchanged on exit. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA, n). \n
          Before entry, the leading m by n part of the array A must
          contain the matrix of coefficients. \n
          Unchanged on exit. \n
 * @param[in] LDA
          LDA is INTEGER \n
          On entry, LDA specifies the first dimension of A as declared
          in the calling (sub) program. LDA must be at least
          fla_max(1, n). \n
          Unchanged on exit. \n
 * @param[in] X
          X is COMPLEX array, dimension
          (1 + (n - 1)*abs(INCX)) \n
          Before entry, the incremented array X must contain the
          vector x. \n
          Unchanged on exit. \n
 * @param[in] INCX
          INCX is INTEGER \n
          On entry, INCX specifies the increment for the elements of
          X. INCX must not be zero. \n
          Unchanged on exit. \n
 * @param[in] BETA
          BETA is REAL. \n
          On entry, BETA specifies the scalar beta. When BETA is
          supplied as zero then Y need not be set on input.
          Unchanged on exit. \n
 * @param[in,out] Y
          Y is REAL array, dimension
          (1 + (n - 1)*abs(INCY)) \n
          Before entry with BETA non-zero, the incremented array Y
          must contain the vector y. On exit, Y is overwritten by the
          updated vector y. \n
 * @param[in] INCY
          INCY is INTEGER \n
          On entry, INCY specifies the increment for the elements of
          Y. INCY must not be zero. \n
          Unchanged on exit. \n

 * */
    template <typename T, typename Ta>
    void la_heamv(integer *uplo, integer *n, Ta *alpha, T *a, integer *lda, T *x, integer *incx,
                  Ta *beta, Ta *y, integer *incy)
    {
        la_heamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
    }
    /** @}*/ // end of la_heamv

    /** @defgroup lascl2 lascl2
     * @ingroup L2 l
     * @{
     */
    /*! @brief LASCL2 performs diagonal scaling on a vector

 * @details
 * \b Purpose:
    \verbatim
    LASCL2 performs a diagonal scaling on a vector:
      x <-- D * x
    where the diagonal matrix D is stored as a vector.

    Eventually to be replaced by BLAS_sge_diag_scale in the new BLAS
    standard.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of D and X. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of X. N >= 0. \n
 * @param[in] D
          D is REAL array, length M \n
          Diagonal matrix D, stored as a vector of length M. \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,N) \n
          On entry, the vector X to be scaled by D.
          On exit, the scaled vector. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the vector X. LDX >= M.  \n

 *  * */
    template <typename T>
    void lascl2(integer *m, integer *n, T *d, T *x, integer *ldx)
    {
        lascl2(m, n, d, x, ldx);
    }
    template <typename T, typename Ta>
    void lascl2(integer *m, integer *n, Ta *d, T *x, integer *ldx)
    {
        lascl2(m, n, d, x, ldx);
    }
    /** @}*/ // end of lascl2

    /** @defgroup larscl2 larscl2
     * @ingroup L2
     * @{
     */
    /*! @brief LARSCL2 performs reciprocal diagonal scaling on a vector

 * @details
 * \b Purpose:
    \verbatim
    LARSCL2 performs a reciprocal diagonal scaling on an vector:
      x <-- inv(D) * x
    where the diagonal matrix D is stored as a vector.

    Eventually to be replaced by BLAS_sge_diag_scale in the new BLAS
    standard.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of D and X. M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns of X. N >= 0. \n
 * @param[in] D
          D is REAL array, length M \n
          Diagonal matrix D, stored as a vector of length M. \n
 * @param[in,out] X
          X is REAL array, dimension (LDX,N) \n
          On entry, the vector X to be scaled by D.
          On exit, the scaled vector. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the vector X. LDX >= M.  \n

 *  * */
    template <typename T>
    void larscl2(integer *m, integer *n, T *d, T *x, integer *ldx)
    {
        larscl2(m, n, d, x, ldx);
    }
    template <typename T, typename Ta>
    void larscl2(integer *m, integer *n, Ta *d, T *x, integer *ldx)
    {
        larscl2(m, n, d, x, ldx);
    }
    /** @}*/ // end of larscl2

    /** @defgroup la_wwaddw la_wwaddw
     * @ingroup L2
     * @{
     */
    /*! @brief LA_WWADDW adds a vector into a doubled-single vector

 * @details
 * \b Purpose:
    \verbatim
    LA_WWADDW adds a vector W into a doubled-single vector (X, Y).

    This works for all extant IBM's hex and binary floating point
    arithmetic, but not for decimal.
    \endverbatim

 * @param[in] N
          N is INTEGER \n
          The length of vectors X, Y, and W. \n
 * @param[in,out] X
          X is REAL array, dimension (N) \n
          The first part of the doubled-single accumulation vector. \n
 * @param[in,out] Y
          Y is REAL array, dimension (N) \n
          The second part of the doubled-single accumulation vector. \n
 * @param[in] W
          W is REAL array, dimension (N) \n
          The vector to be added.  \n

 * */
    template <typename T>
    void la_wwaddw(integer *n, T *x, T *y, T *w)
    {
        la_wwaddw(n, x, y, w);
    }
    /** @}*/ // end of la_wwaddw
    /** @}*/ // end of L2

    /** @defgroup L3 Level 3 BLAS-like matrix-matrix ops
     * @ingroup BLAS
     * @{
     */
    /** @defgroup lagtm lagtm
     * @ingroup L3
     * @{
     */
    /*! @brief LAGTM performs a matrix-matrix product of the form C = Î±AB+Î²C, where A is a \n
     tridiagonal matrix, B and C are rectangular matrices, and Î± and Î² are scalars,   \n
     which may be 0, 1, or -1
 * @details
 * \b Purpose:
    \verbatim
    LAGTM performs a matrix-vector product of the form

       B := alpha * A * X + beta * B

    where A is a tridiagonal matrix of order N, B and X are N by NRHS
    matrices, and alpha and beta are real scalars, each of which may be
    0., 1., or -1.
    \endverbatim

 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          Specifies the operation applied to A. \n
          = 'N':  No transpose, B := alpha * A * X + beta * B \n
          = 'T':  Transpose,   B := alpha * A'* X + beta * B \n
          = 'C':  Conjugate transpose = Transpose \n
 * @param[in] N
          N is INTEGER \n
          The order of the matrix A.  N >= 0. \n
 * @param[in] NRHS
          NRHS is INTEGER \n
          The number of right hand sides, i.e., the number of columns
          of the matrices X and B. \n
 * @param[in] ALPHA
          ALPHA is REAL \n
          The scalar alpha. ALPHA must be 0., 1., or -1.; otherwise,
          it is assumed to be 0. \n
 * @param[in] DL
          DL is REAL array, dimension (N-1) \n
          The (n-1) sub-diagonal elements of T. \n
 * @param[in] D
          D is REAL array, dimension (N) \n
          The diagonal elements of T. \n
 * @param[in] DU
          DU is REAL array, dimension (N-1) \n
          The (n-1) super-diagonal elements of T. \n
 * @param[in] X
          X is REAL array, dimension (LDX,NRHS) \n
          The N by NRHS matrix X. \n
 * @param[in] LDX
          LDX is INTEGER \n
          The leading dimension of the array X.  LDX >= fla_max(N,1). \n
 * @param[in] BETA
          BETA is REAL \n
          The scalar beta. BETA must be 0., 1., or -1.; otherwise,
          it is assumed to be 1. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,NRHS) \n
          On entry, the N by NRHS matrix B.
          On exit, B is overwritten by the matrix expression \n
          B := alpha * A * X + beta * B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B.  LDB >= fla_max(N,1).  \n

 * */
    template <typename T>
    void lagtm(char *trans, integer *n, integer *nrhs, T *alpha, T *dl, T *d, T *du, T *x,
               integer *ldx, T *beta, T *b, integer *ldb)
    {
        lagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
    }
    template <typename T, typename Ta>
    void lagtm(char *trans, integer *n, integer *nrhs, Ta *alpha, T *dl, T *d, T *du, T *x,
               integer *ldx, Ta *beta, T *b, integer *ldb)
    {
        lagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb);
    }
    /** @}*/ // end of lagtm

    /** @defgroup lacrm lacrm
     * @ingroup L3
     * @{
     */
    /*! @brief LACRM multiplies a complex matrix by a square real matrix

 * @details
 * \b Purpose:
    \verbatim
    LACRM performs a very simple matrix-matrix multiplication:
             C := A * B,
    where A is M by N and complex; B is N by N and real;
    C is M by N and complex.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A and of the matrix C.
          M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns and rows of the matrix B and
          the number of columns of the matrix C.
          N >= 0. \n
 * @param[in] A
          A is COMPLEX array, dimension (LDA, N) \n
          On entry, A contains the M by N matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >=max(1,M). \n
 * @param[in] B
          B is REAL array, dimension (LDB, N) \n
          On entry, B contains the N by N matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >=max(1,N). \n
 * @param[out] C
          C is COMPLEX array, dimension (LDC, N) \n
          On exit, C contains the M by N matrix C. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >=max(1,N). \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*M*N) \n

 * */
    template <typename T, typename Ta>
    void lacrm(integer *m, integer *n, T *a, integer *lda, Ta *b, integer *ldb, T *c, integer *ldc,
               T *rwork)
    {
        lacrm(m, n, a, lda, b, ldb, c, ldc, rwork);
    }
    /** @}*/ // end of lacrm

    /** @defgroup larcm larcm
     * @ingroup L3
     * @{
     */
    /*! @brief LARCM copies all or part of a real two-dimensional array to a complex array.

 * @details
 * \b Purpose:
    \verbatim
    LARCM performs a very simple matrix-matrix multiplication:
             C := A * B,
    where A is M by M and real; B is M by N and complex;
    C is M by N and complex.
    \endverbatim

 * @param[in] M
          M is INTEGER \n
          The number of rows of the matrix A and of the matrix C.
          M >= 0. \n
 * @param[in] N
          N is INTEGER \n
          The number of columns and rows of the matrix B and
          the number of columns of the matrix C.
          N >= 0. \n
 * @param[in] A
          A is REAL array, dimension (LDA, M) \n
          On entry, A contains the M by M matrix A. \n
 * @param[in] LDA
          LDA is INTEGER \n
          The leading dimension of the array A. LDA >=max(1,M). \n
 * @param[in] B
          B is COMPLEX array, dimension (LDB, N) \n
          On entry, B contains the M by N matrix B. \n
 * @param[in] LDB
          LDB is INTEGER \n
          The leading dimension of the array B. LDB >=max(1,M). \n
 * @param[out] C
          C is COMPLEX array, dimension (LDC, N) \n
          On exit, C contains the M by N matrix C. \n
 * @param[in] LDC
          LDC is INTEGER \n
          The leading dimension of the array C. LDC >=max(1,M). \n
 * @param[out]	RWORK
          RWORK is REAL array, dimension (2*M*N) \n

 * */
    template <typename T, typename Ta>
    void larcm(integer *m, integer *n, Ta *a, integer *lda, T *b, integer *ldb, T *c, integer *ldc,
               Ta *rwork)
    {
        larcm(m, n, a, lda, b, ldb, c, ldc, rwork);
    }
    /** @}*/ // end of larcm

    /** @defgroup hfrk hfrk
     * @ingroup L3
     * @{
     */
    /*! @brief SFRK performs a symmetric rank-k operation for matrix in RFP format

 * @details
 * \b Purpose:
    \verbatim
     Level 3 BLAS like routine for C in RFP Format.

     SSFRK performs one of the symmetric rank--k operations

        C := alpha*A*A**T + beta*C,

     or

        C := alpha*A**T*A + beta*C,

     where alpha and beta are real scalars, C is an n--by--n symmetric
     matrix and A is an n--by--k matrix in the first case and a k--by--n
     matrix in the second case.
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  The Normal Form of RFP A is stored; \n
          = 'T':  The Transpose Form of RFP A is stored. \n
 * @param[in] UPLO
           UPLO is CHARACTER*1 \n
           On  entry, UPLO specifies whether the upper or lower
           triangular part of the array C is to be referenced as
           follows: \n
              UPLO = 'U' or 'u'   Only the upper triangular part of C
                                  is to be referenced. \n
              UPLO = 'L' or 'l'   Only the lower triangular part of C
                                  is to be referenced. \n
           Unchanged on exit. \n
 * @param[in] TRANS
           TRANS is CHARACTER*1 \n
           On entry, TRANS specifies the operation to be performed as
           follows: \n
              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C. \n
              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C. \n
           Unchanged on exit. \n
 * @param[in] N
           N is INTEGER \n
           On entry, N specifies the order of the matrix C. N must be
           at least zero. \n
           Unchanged on exit. \n
 * @param[in] K
           K is INTEGER \n
           On entry with TRANS = 'N' or 'n', K specifies the number
           of  columns of the matrix A, and on entry with TRANS = 'T'
           or 't', K specifies the number of rows of the matrix A. K
           must be at least zero. \n
           Unchanged on exit. \n
 * @param[in] ALPHA
           ALPHA is REAL \n
           On entry, ALPHA specifies the scalar alpha.
           Unchanged on exit. \n
 * @param[in] A
           A is REAL array, dimension (LDA,ka) \n
           where KA
           is K  when TRANS = 'N' or 'n', and is N otherwise. Before
           entry with TRANS = 'N' or 'n', the leading N--by--K part of
           the array A must contain the matrix A, otherwise the leading
           K--by--N part of the array A must contain the matrix A.
           Unchanged on exit. \n
 * @param[in] LDA
           LDA is INTEGER \n
           On entry, LDA specifies the first dimension of A as declared
           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
           then  LDA must be at least  fla_max( 1, n), otherwise  LDA must
           be at least  fla_max( 1, k).
           Unchanged on exit. \n
 * @param[in] BETA
           BETA is REAL \n
           On entry, BETA specifies the scalar beta.
           Unchanged on exit. \n
 * @param[in,out] C
           C is REAL array, dimension (NT) \n
           NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP
           Format. RFP Format is described by TRANSR, UPLO and N.  \n

 *  * */
    template <typename T>
    void sfrk(char *transr, char *uplo, char *trans, integer *n, integer *k, T *alpha, T *a,
              integer *lda, T *beta, T *c)
    {
        sfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
    }
    template <typename T, typename Ta>
    void hfrk(char *transr, char *uplo, char *trans, integer *n, integer *k, Ta *alpha, T *a,
              integer *lda, Ta *beta, T *c)
    {
        hfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c);
    }
    /** @}*/ // end of

    /** @defgroup tfsm tfsm
     * @ingroup L3
     * @{
     */
    /*! @brief TFSM solves a matrix equation (one operand is a triangular matrix in RFP format).

 * @details
 * \b Purpose:
    \verbatim
     TFSM  solves the matrix equation

        op( A)*X = alpha*B  or  X*op( A) = alpha*B

     where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     non-unit,  upper or lower triangular matrix  and  op( A)  is one  of

        op( A) = A   or   op( A) = A**T.

     A is in Rectangular Full Packed (RFP) Format.

     The matrix X is overwritten on B.
    \endverbatim

 * @param[in] TRANSR
          TRANSR is CHARACTER*1 \n
          = 'N':  The Normal Form of RFP A is stored; \n
          = 'T':  The Transpose Form of RFP A is stored. \n
 * @param[in] SIDE
          SIDE is CHARACTER*1 \n
          On entry, SIDE specifies whether op( A) appears on the left
          or right of X as follows: \n
            SIDE = 'L' or 'l'   op( A)*X = alpha*B. \n
            SIDE = 'R' or 'r'   X*op( A) = alpha*B. \n
          Unchanged on exit. \n
 * @param[in] UPLO
          UPLO is CHARACTER*1 \n
          On entry, UPLO specifies whether the RFP matrix A came from
          an upper or lower triangular matrix as follows: \n
          UPLO = 'U' or 'u' RFP A came from an upper triangular matrix \n
          UPLO = 'L' or 'l' RFP A came from a  lower triangular matrix \n
          Unchanged on exit. \n
 * @param[in] TRANS
          TRANS is CHARACTER*1 \n
          On entry, TRANS  specifies the form of op( A) to be used
          in the matrix multiplication as follows: \n
            TRANS  = 'N' or 'n'   op( A) = A. \n
            TRANS  = 'T' or 't'   op( A) = A'. \n
          Unchanged on exit. \n
 * @param[in] DIAG
          DIAG is CHARACTER*1 \n
          On entry, DIAG specifies whether or not RFP A is unit
          triangular as follows: \n
            DIAG = 'U' or 'u'   A is assumed to be unit triangular. \n
            DIAG = 'N' or 'n'   A is not assumed to be unit
                                triangular. \n
          Unchanged on exit. \n
 * @param[in] M
          M is INTEGER \n
          On entry, M specifies the number of rows of B. M must be at
          least zero.
          Unchanged on exit. \n
 * @param[in] N
          N is INTEGER \n
          On entry, N specifies the number of columns of B.  N must be
          at least zero. \n
          Unchanged on exit. \n
 * @param[in] ALPHA
          ALPHA is REAL \n
          On entry,  ALPHA specifies the scalar  alpha. When  alpha is
          zero then  A is not referenced and  B need not be set before
          entry. \n
          Unchanged on exit. \n
 * @param[in] A
          A is REAL array, dimension (NT) \n
          NT = N*(N+1)/2. On entry, the matrix A in RFP Format.
          RFP Format is described by TRANSR, UPLO and N as follows:
          If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;
          K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If
          TRANSR = 'T' then RFP is the transpose of RFP A as
          defined when TRANSR = 'N'. The contents of RFP A are defined
          by UPLO as follows: If UPLO = 'U' the RFP A contains the NT
          elements of upper packed A either in normal or
          transpose Format. If UPLO = 'L' the RFP A contains
          the NT elements of lower packed A either in normal or
          transpose Format. The LDA of RFP A is (N+1)/2 when
          TRANSR = 'T'. When TRANSR is 'N' the LDA is N+1 when N is
          even and is N when is odd. \n
          See the Note below for more details. Unchanged on exit. \n
 * @param[in,out] B
          B is REAL array, dimension (LDB,N) \n
          Before entry,  the leading  m by n part of the array  B must
          contain  the  right-hand  side  matrix  B,  and  on exit  is
          overwritten by the solution matrix  X. \n
 * @param[in] LDB
          LDB is INTEGER \n
          On entry, LDB specifies the first dimension of B as declared
          in  the  calling  (sub)  program.   LDB  must  be  at  least
          fla_max( 1, m). \n
          Unchanged on exit. \n

 *  * */
    template <typename T>
    void tfsm(char *transr, char *side, char *uplo, char *trans, char *diag, integer *m, integer *n,
              T *alpha, T *a, T *b, integer *ldb)
    {
        tfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb);
    }
    /** @}*/ // end of trfsm

    /** @}*/ // end of L3
    /** @}*/ // end of BLAS
}
#endif